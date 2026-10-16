# SPDX-License-Identifier: GPL-3.0-or-later
# Top-level driver mirroring the upstream `HctCco -d domain.vol -i base.xml ...` invocation.

"""
    generate(domain_vol, base_xml; nterm, aperf, min_dist = 5.0, root = nothing,
             gamma = 3.0, root_radius_max = 1.588, seed = 0, T = Float32,
             volume = :upstream, metric = :upstream, parallel = true, progress = 0,
             export_xml = nothing) -> Engine

Grow a coronary tree from the segmented base tree `base_xml` inside the myocardial mask
`domain_vol`, as upstream `HctCco -d domain_vol -i base_xml -n nterm -a aperf -m min_dist
-p root...` does. `seed` replaces upstream's wall-clock seed (see `CRand`).

Two upstream quirks are reproduced by default and can be corrected:
- `volume = :upstream` counts tree volume the way upstream's computeTotalVolume does, which
  skips everything below a single-child node; `volume = :full` counts the whole tree.
- `metric = :upstream` uses upstream's L3 wall-distance field; `metric = :euclidean` uses the
  Euclidean distance the HCT-CCO paper describes.
- `sampling = :voxel` proposes terminals at integer voxel points as upstream does, which caps
  the terminal count near the number of admissible voxels; `sampling = :continuous` samples
  anywhere inside admissible voxels.
- Upstream moves the end of a split segment to the optimized bifurcation point, which can bend
  imaged base-tree vessels sharply; `fix_base = true` keeps junctions on the imaged
  centerlines.
- Collision and spacing tests use each vessel's radius at insertion time, but early in
  growth all radii are a fraction of their final size, so vessels that passed the test
  later overlap; `collision = :final` tests at radii scaled to the final root radius.
- Upstream places no limit on branching angles; `max_branch_angle = 90` rejects bifurcations
  in which a daughter (or the moved parent segment) turns back by more than 90°, the limit
  Murray's minimum-work analysis of branching angles implies.

**Forests.** With a vector of base trees (`base_xml = ["lca.xml", "rca.xml"]`), the trees grow
together in one domain (a vector of `.vol` masks is merged into their union) and compete for
every terminal: a candidate may join any tree, and the one that adds the least total volume
wins, so the territories emerge from the optimization. All trees share the perfusion pressure
drop and the terminal flow, so each root radius follows from the same law; `root_radius_max`
may be a vector, one cap per tree (in file order). Forests need `mode = :fast`.
With `root_fixed = true` the caps are the imaged root radii rather than upper bounds: after
growth every radius is multiplied by the forest's shared scale min_k cap_k / r_k, so the binding
tree's root sits exactly at its imaged radius and the absolute radii follow from anatomy instead
of the growth hydraulics (upstream's perfusion flow and pressure drop, which are coronary values).
Growth itself is unchanged; use `collision = :final`, which already tests at that scale.
Unconstrained competition lets the volume objective starve a tree with a long imaged trunk (an
XCAT RCA won 165 of 2,000 terminals). To keep anatomical territories without a hard seam, a
domain voxel value ≥ `0xC0` lists in its low six bits which trees may claim a terminal there
(bit k−1 for tree k): give each territory its own tree and a band around the watershed both
trees, and the trees interleave across the band. Any other inside value (e.g. 255) allows all.

**Hierarchical territories.** For finer control (up to 15 groups), pass `territories`, an
`Array{UInt16,3}` on the domain grid (see [`set_territories!`](@ref)), and `base_groups`, a
function mapping an imaged segment's midpoint to its group (1–15). Grown segments inherit their
parent's group, and a terminal joins only groups its voxel allows, so each stretch of an imaged
artery keeps supplying the myocardium nearest to it, as the epicardial segment above a wall
region does in real hearts. Bit 15 marks voxels vessels may cross but where no terminal goes.

Growth is batched by default (`growth = :batched`, fast mode): after the first `n_seq`
terminals it adds terminals in parallel rounds on CPU threads or any GPU (`device =`
`CuArray`, `ROCArray`, `MtlArray`, ...), about 15× faster than sequential growth at 800k
terminals and statistically equivalent to it, not bitwise identical (see `docs/deviations.md`
§12 and [`grow!`](@ref)). `growth = :sequential` is the reference path; `mode = :exact`, the
bit-for-bit C++ replica, is always sequential.
"""
function generate(
        domain_vol, base_xml; nterm::Integer,
        aperf::Real, min_dist::Real = 5.0, root = nothing, gamma::Real = 3.0,
        root_radius_max = 1.588, seed::Integer = 0, T::Type = Float32,
        parallel::Bool = true, progress::Integer = 0, export_xml = nothing,
        trace = nothing, volume::Symbol = :upstream, metric::Symbol = :upstream,
        mode::Symbol = :fast, sampling::Symbol = :voxel, fix_base::Bool = false,
        collision::Symbol = :current, max_branch_angle = nothing,
        growth::Symbol = mode === :exact ? :sequential : :batched, n_seq::Integer = 2_000, batch_frac::Real = 0.01, batch_max::Integer = 16_384, device = identity,
        territories = nothing, base_groups = nothing, forbidden = nothing, root_fixed::Bool = false, kwargs...
    )
    p = CCOParams(;
        nterm = Int(nterm), aperf = Float64(aperf), gamma = Float64(gamma),
        root_radius_max = Float64(first(root_radius_max)),
        min_distance_to_border = Float64(min_dist), kwargs...
    )
    dom = MaskDomain{T}(
        domain_vol; root = root === nothing ? nothing : Tuple(Int32.(root)),
        min_dist = min_dist, metric = metric
    )
    territories === nothing || set_territories!(dom, territories)
    forbidden === nothing || set_forbidden!(dom, forbidden...)       # (mask, origin mm, cell mm)
    m = import_gxl(T, base_xml, p)
    E = Engine(
        m, dom; seed = seed, mode = mode, volume = volume, sampling = sampling, fix_base = fix_base,
        collision = collision, max_branch_angle = max_branch_angle, root_caps = root_radius_max
    )  # mode = :exact is test-only
    if base_groups !== nothing
        # territory group of each imaged segment, from its midpoint; grown segments inherit it
        for j in 2:nsegments(m)
            a = m.pos[m.parent[j] == 1 ? j : m.parent[j]]; b = m.pos[j]
            E.treeid[j] = UInt8(base_groups(((a[1] + b[1]) / 2, (a[2] + b[2]) / 2, (a[3] + b[3]) / 2)))
        end
    end
    grow!(
        E; parallel = parallel, progress = Int(progress), trace = trace, growth = growth,
        n_seq = n_seq, batch_frac = batch_frac, batch_max = batch_max, device = device
    )
    if root_fixed
        # one shared scale for the whole forest (all trees at the same pressure), set by the tree
        # whose cap binds; collision = :final already tested at this scale during growth
        s = Inf
        for (k, c) in enumerate(E.roots)
            E.rcaps[k] > 0 && (s = min(s, Float64(E.rcaps[k]) / Float64(m.radius[c])))
        end
        isfinite(s) && (m.radius .*= eltype(m.radius)(s))
    end
    export_xml === nothing || export_gxl(E.m, export_xml)
    return E
end
