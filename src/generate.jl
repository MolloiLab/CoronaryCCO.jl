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
"""
function generate(
        domain_vol::AbstractString, base_xml::AbstractString; nterm::Integer,
        aperf::Real, min_dist::Real = 5.0, root = nothing, gamma::Real = 3.0,
        root_radius_max::Real = 1.588, seed::Integer = 0, T::Type = Float32,
        parallel::Bool = true, progress::Integer = 0, export_xml = nothing,
        trace = nothing, volume::Symbol = :upstream, metric::Symbol = :upstream,
        mode::Symbol = :fast, sampling::Symbol = :voxel, fix_base::Bool = false,
        collision::Symbol = :current, max_branch_angle = nothing, kwargs...
    )
    p = CCOParams(;
        nterm = Int(nterm), aperf = Float64(aperf), gamma = Float64(gamma),
        root_radius_max = Float64(root_radius_max),
        min_distance_to_border = Float64(min_dist), kwargs...
    )
    dom = MaskDomain{T}(
        domain_vol; root = root === nothing ? nothing : Tuple(Int32.(root)),
        min_dist = min_dist, metric = metric
    )
    m = import_gxl(T, base_xml, p)
    E = Engine(
        m, dom; seed = seed, mode = mode, volume = volume, sampling = sampling, fix_base = fix_base,
        collision = collision, max_branch_angle = max_branch_angle
    )  # mode = :exact is test-only
    grow!(E; parallel = parallel, progress = Int(progress), trace = trace)
    export_xml === nothing || export_gxl(E.m, export_xml)
    return E
end
