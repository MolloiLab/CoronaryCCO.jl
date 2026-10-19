# SPDX-License-Identifier: GPL-3.0-or-later
# The fourth stage, after grow, flow and transport: a shareable dynamic CT phantom of the transported bed, in the
# DynamicPhantom format (one compartment volume + one CSV per scenario), the same for every organ. Only the voxel
# classifier, which voxel of the subject's label volume is what, comes from the caller: it holds the subject's (licensed)
# anatomy.

using .DynamicPhantom: compartment_classes, VolumeWriter, write_slice!, write_phantom_csv
using .OrganCT: tissue_iodine, vessel_volumes

"""
    Pool(name, material, segment, what; scale = 1, ml = 0)

A compartment filled from the whole-body model: `segment` names a `Body` segment, `what` which of its curves: `:out`
(outflow blood, shifted by the voxel's transit delay), `:mean` (the segment's mean blood), `:organ` (an organ bed: its
blood and extracellular iodine over its tissue volume `ml`), or `:tissue` (the organ's own perfused tissue, its mean).
`scale` multiplies the curve (e.g. a blood fraction). `material` is the `XrayAttenuation.Materials` base.
"""
struct Pool
    name::String
    material::String
    segment::Symbol
    what::Symbol
    scale::Float32
    ml::Float32
end
Pool(name, material, segment, what; scale = 1, ml = 0) = Pool(name, material, segment, what, scale, ml)

"""
    PhantomSpec(; times, tissue_material, blood_material, tissue_classes, vessel_classes, vessel_rmin_mm, delay_bin_s,
                energy_keV, materials, pools)
    PhantomSpec(spec::PhantomSpec; kwargs...)

Everything that sets an organ's dynamic CT phantom. Every field is a keyword and every one is required: an organ is an
instance (`HEART_PHANTOM`, `LUNG_PHANTOM`, organs.jl); a variant is the second form.

- `times`: the frame times (s) of the iodine curves.
- `tissue_material`: the `XrayAttenuation.Materials` name of the perfused tissue (its base; iodine is added on top). A
  `:tissue` pool takes it too.
- `blood_material`: the base of the grown and imaged vessels' classes.
- `tissue_classes`: tissue voxels are grouped into this many kinetic classes (k-means on their joint curves over every
  state), so one compartment volume serves every state; at most 999 (ids 1001-1999).
- `vessel_classes`, `vessel_rmin_mm`: grown vessels of radius ≥ `vessel_rmin_mm` are drawn into the tissue as blood, and
  the vessels a classifier names (imaged ones) take their grown segment's curve; all grouped into this many classes.
- `delay_bin_s`: `:out` pools are binned by the voxel's transit delay (s).
- `energy_keV`: energy of the base HU written for each material.
- `materials`: static label → material name (labels 0-99; no iodine).
- `pools`: the `Pool`s a classifier can name, by index.
"""
struct PhantomSpec
    times::Vector{Float32}
    tissue_material::String
    blood_material::String
    tissue_classes::Int
    vessel_classes::Int
    vessel_rmin_mm::Float32
    delay_bin_s::Float32
    energy_keV::Float32
    materials::Dict{Int, String}
    pools::Vector{Pool}
end
function PhantomSpec(; kwargs...)
    bad = setdiff(keys(kwargs), fieldnames(PhantomSpec))
    isempty(bad) || throw(ArgumentError("unknown PhantomSpec fields $(collect(bad))"))
    miss = setdiff(fieldnames(PhantomSpec), keys(kwargs))
    isempty(miss) || throw(ArgumentError("PhantomSpec: no value for $(collect(miss))"))
    s = PhantomSpec((kwargs[f] for f in fieldnames(PhantomSpec))...)
    all(k -> 0 <= k <= 99, keys(s.materials)) || throw(ArgumentError("PhantomSpec: static labels must be 0-99"))
    all(p -> p.what in (:out, :mean, :organ, :tissue), s.pools) || throw(ArgumentError("PhantomSpec: pool `what` must be :out, :mean, :organ or :tissue"))
    all(p -> p.what !== :organ || p.ml > 0, s.pools) || throw(ArgumentError("PhantomSpec: an :organ pool needs its tissue volume `ml` > 0"))
    1 <= s.tissue_classes <= 999 || throw(ArgumentError("PhantomSpec: tissue_classes must be 1-999 (ids 1001-1999)"))
    1 <= s.vessel_classes <= 63535 || throw(ArgumentError("PhantomSpec: vessel_classes must be 1-63535 (ids from 2001)"))
    return s
end
PhantomSpec(s::PhantomSpec; kwargs...) = PhantomSpec(; merge((; (f => getfield(s, f) for f in fieldnames(PhantomSpec))...), values(kwargs))...)

# what a classifier returns for a voxel centre: (kind, a, d)
const VOX_STATIC = 0x00 # a = static label (PhantomSpec.materials)
const VOX_POOL = 0x01 # a = pool index (PhantomSpec.pools); d = transit delay (s) for `:out` pools
const VOX_TISSUE = 0x02 # a = the territories' tissue voxel number (0: tissue the territories do not reach)
const VOX_VESSEL = 0x03 # a = grown segment: > 0 an arterial segment, < 0 a venous one (-a)

"""
    phantom(spec, scans, classify, grid; pairs, terr, nvox, out, name = "phantom", facts = []) -> NamedTuple

Writes the phantom of a transported bed into folder `out`: `<name>_compartments.u16.zst` (UInt16 ids, column-major,
zstd) and `label_healthy.csv`, in the DynamicPhantom format (`DynamicPhantom.load`, `DynamicPhantom.bs_phantom`).

- `scans`: one `(state, C, T)` per state (`"rest"`, `"stress"`): the solved circuit and its `transport` result on
  `spec.times`.
- `classify(p) -> (kind, a, d)` for a voxel centre `p` (mm, the trees' frame): one of `VOX_STATIC`, `VOX_POOL`,
  `VOX_TISSUE`, `VOX_VESSEL` (see above). It is the subject-specific step; everything else here is the same for every
  organ.
- `grid`: the output `VoxelGrid` (voxel centres, isotropic `h` mm).
- `pairs`, `terr`, `nvox`: the territories (`flow`'s `pairs` and `terr`) over the `nvox` tissue voxels the classifier
  numbers.

Iodine is in mg per mL of voxel. Ids: 0-99 static labels, 100-999 pools (and delay bins), 1000 tissue without a territory
(static), 1001… tissue classes (1 = least enhancing), 2001… vessel classes. The result holds the ids, their voxel counts
and each state's tissue flow (mL/min per mL of tissue), the ground truth.
"""
function phantom(
        spec::PhantomSpec, scans, classify, grid::VoxelGrid; pairs, terr, nvox::Integer, out::AbstractString,
        name::AbstractString = "phantom", facts = Pair{String, Any}[]
    )
    mkpath(out); nT = length(spec.times); nS = length(scans); states = [String(s[1]) for s in scans]
    # ---- tissue: every tissue voxel's curve in every state, k-means into kinetic classes
    Xs = [Array(tissue_iodine(Array(T.tissue), pairs, terr, collect(Int32(1):Int32(nvox)))) for (_, _, T) in scans] # nvox × nT each
    X = Matrix{Float32}(reduce(vcat, permutedims.(Xs))) # (nS nT) × nvox
    a, cen = compartment_classes(X, spec.tissue_classes)
    ord = sortperm(vec(sum(cen; dims = 1))); rank = invperm(ord)
    tid = Int32[1000 + rank[a[v]] for v in 1:nvox]
    tmean = [vec(sum(Xs[s]; dims = 1)) ./ nvox for s in 1:nS]
    # ---- vessels: grown segments of radius ≥ vessel_rmin_mm and the imaged (base) ones, k-means on their curves
    C1 = scans[1][2]; trees = (C1.tree, C1.veins); vol = vessel_volumes(trees...); na = length(C1.tree.r)
    vcurve_of(t, j) = Float32[(T = scans[s][3]; (t == 1 ? T.artery[f, j] : T.vein[f, j]) / vol[t == 1 ? j : na + j]) for s in 1:nS for f in 1:nT]
    seg = Tuple{Int, Int}[(t, j) for t in 1:2 for j in eachindex(trees[t].r) if trees[t].r[j] >= spec.vessel_rmin_mm || trees[t].isbase[j]]
    V = isempty(seg) ? zeros(Float32, nS * nT, 0) : reduce(hcat, (vcurve_of(t, j) for (t, j) in seg))
    nvc = clamp(spec.vessel_classes, 1, max(1, length(seg)))
    av, vcen = isempty(seg) ? (Int[], zeros(Float32, nS * nT, 1)) : compartment_classes(V, nvc)
    vord = sortperm(vec(sum(vcen; dims = 1))); vrank = invperm(vord)
    vclass = Dict{Tuple{Int, Int}, Int32}(seg[k] => Int32(2000 + vrank[av[k]]) for k in eachindex(seg))
    # a segment a classifier names outside that set: its nearest class centre
    function vessel_id(t, j)
        return get!(vclass, (t, j)) do
            c = vcurve_of(t, j); Int32(2000 + vrank[argmin([sum(abs2, c .- vcen[:, k]) for k in axes(vcen, 2)])])
        end
    end
    # the grown lumens drawn into the tissue, by the z range they cover (slice-by-slice stamp)
    stamped = [(t, j) for (t, j) in seg if trees[t].r[j] >= spec.vessel_rmin_mm]
    segab(t, j) = (trees[t].a[:, j], trees[t].b[:, j], trees[t].r[j])
    zr = [(sa = segab(t, j); (min(sa[1][3], sa[2][3]) - sa[3], max(sa[1][3], sa[2][3]) + sa[3])) for (t, j) in stamped]
    zo = sortperm(first.(zr))
    # ---- pools, ids from 100 in order of first use (an :out pool once per delay bin)
    pid = Dict{Tuple{Int, Int}, Int}()
    function poolid(i, d)
        id = get!(pid, (i, spec.pools[i].what === :out ? round(Int, d / spec.delay_bin_s) : 0), 100 + length(pid))
        id <= 999 || error("phantom: more than 900 pool compartments (pools x delay bins); widen delay_bin_s")
        return id
    end
    function pcurve(i, b, s, f)
        P = spec.pools[i]; t = spec.times[f]; r = scans[s][3].body
        c = P.what === :tissue ? tmean[s][f] :
            P.what === :out ? _interp(t - b * spec.delay_bin_s, r.t, r.out[P.segment]) :
            P.what === :mean ? _interp(t, r.t, r.mean[P.segment]) :
            (
                g = r.circ.segs[findfirst(x -> x.name === P.segment, r.circ.segs)];
                (g.v * _interp(t, r.t, r.mean[P.segment]) + g.ve * _interp(t, r.t, r.ecv[P.segment])) / P.ml
            )
        return Float32(P.scale * c)
    end
    # ---- stream the volume, one z slice at a time
    nx, ny, nz = grid.dims; h = grid.h; o = grid.origin
    file = "$(name)_compartments.u16.zst"; w = VolumeWriter(joinpath(out, file), (nx, ny, nz))
    kind = zeros(UInt8, nx, ny); A = zeros(Int32, nx, ny); D = zeros(Float32, nx, ny); sl = zeros(UInt16, nx, ny)
    act = 1
    for kz in 1:nz
        z = o[3] + (kz - 1) * h
        AK.foreachindex(1:(nx * ny); min_elems = 4096) do q
            ix = 1 + (q - 1) % nx; jy = 1 + (q - 1) ÷ nx
            k, aa, d = classify((o[1] + (ix - 1) * h, o[2] + (jy - 1) * h, z))
            kind[q] = k; A[q] = aa; D[q] = d
        end
        for q in 1:(nx * ny)
            k = kind[q]; aa = A[q]
            sl[q] = UInt16(
                k == VOX_STATIC ? aa : k == VOX_POOL ? poolid(aa, D[q]) : k == VOX_TISSUE ? (aa > 0 ? tid[aa] : 1000) :
                    vessel_id(aa > 0 ? 1 : 2, abs(aa))
            )
        end
        # grown lumens replace tissue where they cover a tissue voxel's centre
        while act <= length(zo) && zr[zo[act]][2] < z
            act += 1
        end
        for kk in act:length(zo)
            k = zo[kk]; zr[k][1] > z && break
            zr[k][2] < z && continue
            (t, j) = stamped[k]; pa, pb, r = segab(t, j); ab = pb .- pa; L2 = sum(abs2, ab); id = UInt16(vessel_id(t, j))
            lo1 = max(1, floor(Int, (min(pa[1], pb[1]) - r - o[1]) / h) + 1); hi1 = min(nx, ceil(Int, (max(pa[1], pb[1]) + r - o[1]) / h) + 1)
            lo2 = max(1, floor(Int, (min(pa[2], pb[2]) - r - o[2]) / h) + 1); hi2 = min(ny, ceil(Int, (max(pa[2], pb[2]) + r - o[2]) / h) + 1)
            for jy in lo2:hi2, ix in lo1:hi1
                kind[ix, jy] == VOX_TISSUE || continue
                p1 = o[1] + (ix - 1) * h - pa[1]; p2 = o[2] + (jy - 1) * h - pa[2]; p3 = z - pa[3]
                u = L2 > 0 ? clamp((p1 * ab[1] + p2 * ab[2] + p3 * ab[3]) / L2, 0, 1) : 0.0
                (p1 - u * ab[1])^2 + (p2 - u * ab[2])^2 + (p3 - u * ab[3])^2 <= r^2 && (sl[ix, jy] = id)
            end
        end
        write_slice!(w, sl)
    end
    close(w)
    cnt = Dict(i - 1 => n for (i, n) in enumerate(w.counts) if n > 0)
    # ---- the scenario table
    pkey = Dict(v => k for (k, v) in pid)
    vinv = Dict{Int32, Int}(Int32(2000 + vrank[k]) => k for k in axes(vcen, 2))
    cols = NamedTuple[]
    for id in sort!(collect(keys(cnt)))
        c = if id < 100
            haskey(spec.materials, id) || error("phantom: the classifier gave static label $id, which spec.materials does not name")
            (key = id, name = "static $id", base_material = spec.materials[id], kinetics = "static", voxels = cnt[id])
        elseif id < 1000
            (i, b) = pkey[id]; P = spec.pools[i]
            (
                key = id, name = P.what === :out ? "$(P.name) +$(round(b * spec.delay_bin_s; digits = 2)) s" : P.name,
                base_material = P.what === :tissue ? spec.tissue_material : P.material, kinetics = String(P.what), voxels = cnt[id],
            )
        elseif id == 1000
            (key = id, name = "tissue (no territory)", base_material = spec.tissue_material, kinetics = "static", voxels = cnt[id])
        elseif id < 2000
            (key = id, name = "tissue class $(id - 1000)", base_material = spec.tissue_material, kinetics = "tissue", voxels = cnt[id])
        else
            (key = id, name = "vessel class $(id - 2000)", base_material = spec.blood_material, kinetics = "vessel", voxels = cnt[id])
        end
        push!(cols, c)
    end
    flow = [s => sum(C.q[j] for j in eachindex(C.tree.parent) if C.tree.parent[j] == 0) * 60 / sum(terr) for (s, C, _) in scans]
    allfacts = Pair{String, Any}[
        "compartment_volume" => file, "dims" => "$nx $ny $nz", "voxel_mm" => h, "origin_mm" => join(string.(o), " "),
        "states" => join(states, " "), ("flow_ml_min_per_ml_" * s => round(f; digits = 4) for (s, f) in flow)..., facts...,
    ]
    function curve(j, st, f)
        id = cols[j].key; s = findfirst(==(st), states)
        (id < 100 || id == 1000) && return nothing
        id < 1000 && return (ib = pkey[id]; pcurve(ib[1], ib[2], s, f))
        id < 2000 && return cen[(s - 1) * nT + f, ord[id - 1000]]
        return vcen[(s - 1) * nT + f, vinv[Int32(id)]]
    end
    csv = write_phantom_csv(joinpath(out, "label_healthy.csv"); facts = allfacts, columns = cols, times = spec.times, states, curve, energy_keV = spec.energy_keV)
    return (; csv, volume = joinpath(out, file), ids = sort!(collect(keys(cnt))), counts = cnt, flow)
end
