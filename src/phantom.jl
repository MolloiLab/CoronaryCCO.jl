# SPDX-License-Identifier: GPL-3.0-or-later
# The fourth stage, after grow, flow and transport: a shareable dynamic CT phantom of the transported bed, one folder per
# scenario in the DynamicPhantom format (labels, ROIs, ground-truth MBF maps and their tables), the same for every organ.
# Only the voxel classifier and the ROI rule, which voxel of the subject is what, come from the caller: they hold the
# subject's (licensed) anatomy.

import XrayAttenuation as XA
using Unitful: ustrip, @u_str
using .DynamicPhantom: compartment_classes, VolumeWriter, write_slice!, write_bundle
using .OrganCT: tissue_iodine, blood_volumes, unresolved_masses
using .Flow: sinus_outflow

"""
    Pool(name, material, segment, what; scale = 1, ml = 0, through = segment, lead_ml = 0, lead_through = segment,
         upstream = false)

A compartment filled from the whole-body model: `segment` names a `Body` segment, `what` which of its curves: `:out`
(its outflow blood, shifted by the voxel's transit, `pool_delay`), `:mean` (the segment's mean blood), `:organ` (an organ
bed: its blood and extracellular iodine over its tissue volume `ml`), or `:tissue` (the organ's own perfused tissue, its
mean). `scale` multiplies the curve (e.g. a blood fraction). `material` is the `XrayAttenuation.Materials` base.

An `:out` pool is positional: each voxel carries the blood volume (mL) between it and the segment's outlet along its vessel.
That volume moves at the flow of segment `through`; `lead_ml` more lie before the pool's vessel begins, at the flow of
`lead_through` (the aortic arch before the descending aorta). An `upstream` pool lies before the outlet (a vein draining
into a chamber), so its voxels see the outlet's curve earlier.
"""
struct Pool
    name::String
    material::String
    segment::Symbol
    what::Symbol
    scale::Float32
    ml::Float32
    through::Symbol
    lead_ml::Float32
    lead_through::Symbol
    upstream::Bool
end
function Pool(
        name, material, segment, what; scale = 1, ml = 0, through = segment, lead_ml = 0, lead_through = segment,
        upstream = false
    )
    return Pool(name, material, segment, what, scale, ml, through, lead_ml, lead_through, upstream)
end

"""
    pool_delay(P::Pool, r, ml) -> Float32

The transit (s) of a voxel of `:out` pool `P` lying `ml` of blood from the outlet of `P.segment` along its vessel, in the
state of `Body.simulate` result `r`: each stretch of blood over the flow through it (central volume principle, plug flow),
`P.lead_ml` at the flow of `P.lead_through` and `ml` at the flow of `P.through`. Negative for an `upstream` pool, whose
voxels see the outlet's curve earlier. One anatomy thus gives each state its own transits, shorter where the flow is higher.
"""
function pool_delay(P::Pool, r, ml)
    flow(name) = (i = findfirst(x -> x.name === name, r.circ.segs); r.circ.q[i] * r.circ.co)
    delay = P.lead_ml / flow(P.lead_through) + ml / flow(P.through)
    return Float32(P.upstream ? -delay : delay)
end

"""
    pool_curve(P::Pool, r, t; delay = 0, tissue = 0) -> Float32

Pool `P`'s iodine (mg per mL) at time `t` (s) in a `Body.simulate` result `r` (a `transport` result's `body`): `:out` is
the segment's outflow `delay` s earlier (`pool_delay`), `:mean` its mean blood, `:organ` its blood and extracellular iodine over `P.ml`,
`:tissue` the organ's own tissue value `tissue`. The one definition of the body's compartments, for phantoms and pictures.
"""
function pool_curve(P::Pool, r, t; delay = 0.0f0, tissue = 0.0f0)
    c = P.what === :tissue ? tissue :
        P.what === :out ? _interp(t - delay, r.t, r.out[P.segment]) :
        P.what === :mean ? _interp(t, r.t, r.mean[P.segment]) :
        (
            g = r.circ.segs[findfirst(x -> x.name === P.segment, r.circ.segs)];
            (g.v * _interp(t, r.t, r.mean[P.segment]) + g.ve * _interp(t, r.t, r.ecv[P.segment])) / P.ml
        )
    return Float32(P.scale * c)
end

"""
    PhantomSpec(; times, tissue_material, blood_material, tissue_classes, vessel_classes, vessel_rmin_mm, delay_bin_ml,
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
- `delay_bin_ml`: `:out` pools are binned by the blood volume (mL) between the voxel and the segment's outlet.
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
    delay_bin_ml::Float32
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
# the largest blood volume (mL) an `:out` pool voxel can lie downstream of its segment's outlet
const MAX_UPSTREAM_ML = 600

"""
A voxel kind `classify(p)` returns: a static material; `a` is its label in `PhantomSpec.materials`.
"""
const VOX_STATIC = 0x00 # a = static label (PhantomSpec.materials)
"""
A voxel kind `classify(p)` returns: a blood pool of the body model; `a` is its index in `PhantomSpec.pools` and `d` the
blood (mL) upstream of it along its vessel.
"""
const VOX_POOL = 0x01 # a = pool index (PhantomSpec.pools); d = blood volume (mL) upstream along the vessel, for `:out` pools
"""
A voxel kind `classify(p)` returns: perfused tissue; `a` is its tissue voxel's number (0 where no territory reaches).
"""
const VOX_TISSUE = 0x02 # a = the territories' tissue voxel number (0: tissue the territories do not reach)
"""
A voxel kind `classify(p)` returns: a grown vessel; `a > 0` is an arterial segment, `a < 0` the venous segment `-a`.
"""
const VOX_VESSEL = 0x03 # a = grown segment: > 0 an arterial segment, < 0 a venous one (-a)


"""
    phantom(spec, scans, classify, grid; pairs, terr, nvox, out, roi = p -> 0, roi_names = [], roi_groups = [],
            tracking = [], facts = [], check = true) -> NamedTuple

Writes the phantom of a transported bed into the scenario folder `out` (the `DynamicPhantom` bundle: `labels.u16.zst`,
`rois.u8.zst`, `mbf_<state>.f32.zst`, `label_<state>.csv`, `ledger_<state>.csv`, `phantom.toml`).

- `scans`: one `(state, C, T)` per state (`"rest"`, `"stress"`): the solved circuit and its `transport` result on
  `spec.times`.
- `classify(p) -> (kind, a, d)` for a voxel centre `p` (mm, the trees' frame): one of `VOX_STATIC`, `VOX_POOL`,
  `VOX_TISSUE`, `VOX_VESSEL` (see above).
- `roi(p) -> 0:length(roi_names)`: the ground-truth compartment of the voxel at `p` (0: none). A
  compartment is a region of the organ with everything in it: tissue and every artery, vein and lumen that lies in it.
  Only tissue and vessel voxels take a compartment, so every compartment is a set of whole labels.
  `roi_groups`: `name => ROI ids` (e.g. the LV as the union of its territories).
- `grid`: the output `VoxelGrid` (voxel centres, isotropic `h` mm).
- `pairs`, `terr`, `nvox`: the territories (`flow`'s `pairs` and `terr`) over the `nvox` tissue voxels the classifier
  numbers.
- `tracking`: `name => (x, y, z)` (mm): bolus-tracking points; `DynamicPhantom.scan_times` reads their curves. A point in
  a vessel segment gives that segment a class of its own (one of `vessel_classes`), painted at the segment's blood
  concentration, so the point reads that segment's curve alone. (Every other vessel's mass is spread over its painted
  voxels; a tracked segment's painted iodine differs from its mass by its lumen's volume beyond the segment's.)

Labels: 0-99 static, 100-999 pools (and their delay bins), 1000 tissue without a territory (static), 1001… tissue
(`tissue_classes` kinetic classes in each ROI), 2001… drawn vessels (`vessel_classes` classes per ROI, arteries and veins
apart). A tissue or vessel label lies in one ROI. Its iodine is the mass of its members (beds spread over the tissue
sub-voxels that remain after lumens are drawn; each vessel segment over its painted voxels) over their volume, so
painting conserves mass, and its ground-truth MBF is its members' bed flow over their mass (a vessel's is 0). An artery or
vein with no painted voxel gives its iodine to the beds it feeds and drains. A compartment's iodine and ground truth are
therefore sums over its labels (`DynamicPhantom.load` gives them as `roi_table`): its MBF is its beds' flow over the mass
of everything in it, vessels included.

The ledger (`ledger_<state>.csv`) is the organ's iodine at each time: delivered through the inlets, held in its arteries,
beds and veins, out through its venous outlet, and painted (a tracked segment counted at its mass). With `check`, the
build fails when delivered ≠ held + outflow or painted ≠ held (beyond 1% of the peak delivered mass). `verbose` logs each
stage's seconds.
"""
function phantom(
        spec::PhantomSpec, scans, classify, grid::VoxelGrid; pairs, terr, nvox::Integer, out::AbstractString,
        roi = (p -> 0), roi_names::AbstractVector{<:AbstractString} = String[], roi_groups = Pair{String, Vector{Int}}[],
        tracking = Pair{String, Tuple}[], facts = Pair{String, Any}[], check::Bool = true, verbose::Bool = false
    )
    started = Ref(time()) # `verbose`: each stage's seconds
    lap(stage) = (verbose && @info "phantom: $stage $(round(time() - started[]; digits = 1)) s"; started[] = time())
    mkpath(out); nT = length(spec.times); nS = length(scans); states = [String(s[1]) for s in scans]
    R = length(roi_names); K = spec.tissue_classes; KV = spec.vessel_classes
    R <= 255 || throw(ArgumentError("phantom: at most 255 ROIs"))
    (R + 1) * K <= 999 || throw(ArgumentError("phantom: (ROIs + 1) × tissue_classes = $((R + 1) * K) > 999 (ids 1001-1999)"))
    2000 + 2 * (R + 1) * KV <= 65535 || throw(ArgumentError("phantom: (ROIs + 1) × 2 × vessel_classes does not fit UInt16 ids"))
    vml = Float64(grid.h)^3 / 1000; voxml = Float64(sum(terr)) / nvox
    vox = collect(Int32(1):Int32(nvox))
    # ---- tissue: every territory voxel's curves (bed + the iodine of vessels too thin to draw), k-means into kinetic classes
    beds(U) = [Array(tissue_iodine(Array(scans[s][3].tissue) .+ U[s], pairs, terr, vox)) for s in 1:nS] # nvox × nT, mg/mL
    Xs = beds([unresolved_masses(C, T.artery, T.vein; rmin = spec.vessel_rmin_mm) for (_, C, T) in scans])
    a, cen = compartment_classes(Matrix{Float32}(reduce(vcat, permutedims.(Xs))), K)
    trank = invperm(sortperm(vec(sum(cen; dims = 1)))); tclass = Int[trank[a[v]] for v in 1:nvox] # 1 = least enhancing
    tlabel(r, c) = 1000 + r * K + c
    lap("tissue classes")
    # ---- vessels: grown segments of radius ≥ vessel_rmin_mm and the imaged (base) ones, k-means on their blood curves; a
    # segment a tracking point lies in keeps a class of its own, so the point reads that segment's curve
    C1 = scans[1][2]; trees = (C1.tree, C1.veins); na = length(C1.tree.r); nv = length(C1.veins.r)
    vols = [blood_volumes(C) for (_, C, _) in scans]
    vcurve(t, j) = Float32[(T = scans[s][3]; (t == 1 ? T.artery[f, j] : T.vein[f, j]) / vols[s][t == 1 ? j : na + j]) for s in 1:nS for f in 1:nT]
    tracked = Dict{Tuple{Int, Int}, String}()
    for (tn, p) in tracking
        k, aa, _ = classify(Tuple(1.0 .* p))
        k == VOX_VESSEL && (tracked[(aa > 0 ? 1 : 2, abs(aa))] = tn)
    end
    KS = KV - length(tracked) # classes left for the other segments
    KS >= 1 || throw(ArgumentError("phantom: vessel_classes must exceed the $(length(tracked)) tracked segments"))
    seg = Tuple{Int, Int}[(t, j) for t in 1:2 for j in eachindex(trees[t].r) if (trees[t].r[j] >= spec.vessel_rmin_mm || trees[t].isbase[j]) && !haskey(tracked, (t, j))]
    av, vcen = isempty(seg) ? (Int[], zeros(Float32, nS * nT, 1)) : compartment_classes(reduce(hcat, (vcurve(t, j) for (t, j) in seg)), clamp(KS, 1, length(seg)))
    vrank = invperm(sortperm(vec(sum(vcen; dims = 1))))
    vcls = Dict{Tuple{Int, Int}, Int}(seg[k] => vrank[av[k]] for k in eachindex(seg))
    tracked_class = Dict{Tuple{Int, Int}, String}() # (artery 1 / vein 2, class) => the tracking point's name
    for (i, (s, tn)) in enumerate(sort!(collect(tracked)))
        vcls[s] = KS + i; tracked_class[(s[1], KS + i)] = tn
    end
    seg = [seg; sort!(collect(keys(tracked)))]
    vclass(t, j) = get!(vcls, (t, j)) do # a segment the classifier names outside that set: its nearest class
        c = vcurve(t, j); vrank[argmin([sum(abs2, c .- vcen[:, k]) for k in axes(vcen, 2)])]
    end
    vlabel(r, t, c) = 2000 + (2r + t - 1) * KV + c
    # the grown lumens drawn into the tissue, by the z range they cover (slice-by-slice stamp)
    stamped = [(t, j) for (t, j) in seg if trees[t].r[j] >= spec.vessel_rmin_mm]
    segab(t, j) = (trees[t].a[:, j], trees[t].b[:, j], trees[t].r[j])
    zr = [(sa = segab(t, j); (min(sa[1][3], sa[2][3]) - sa[3], max(sa[1][3], sa[2][3]) + sa[3])) for (t, j) in stamped]
    zo = sortperm(first.(zr))
    lap("vessel classes")
    # ---- pools, ids from 100 in order of first use (an :out pool once per delay bin)
    pid = Dict{Tuple{Int, Int}, Int}()
    nbin = round(Int, MAX_UPSTREAM_ML / spec.delay_bin_ml)
    pbin(i, d) = spec.pools[i].what === :out ? clamp(round(Int, d / spec.delay_bin_ml), 0, nbin) : 0
    function poolid(i, d)
        id = get!(pid, (i, pbin(i, d)), 100 + length(pid))
        id <= 999 || error("phantom: more than 900 pool compartments (pools x delay bins); widen delay_bin_ml")
        return id
    end
    tmean = [vec(sum(Xs[s]; dims = 1)) ./ nvox for s in 1:nS]
    pcurve(i, b, s, f) = (
        P = spec.pools[i]; r = scans[s][3].body;
        pool_curve(P, r, spec.times[f]; delay = P.what === :out ? pool_delay(P, r, b * spec.delay_bin_ml) : 0.0f0, tissue = tmean[s][f])
    )
    function label_at(k, aa, d, r) # a voxel's label from the classifier's answer and its ROI
        k == VOX_STATIC && return Int(aa)
        k == VOX_POOL && return poolid(aa, d)
        k == VOX_TISSUE && return aa > 0 ? tlabel(r, tclass[aa]) : 1000
        t = aa > 0 ? 1 : 2; return vlabel(r, t, vclass(t, abs(aa)))
    end
    lap("pools")
    # ---- stream the volumes one z slice at a time, counting every voxel by (ROI, label), tissue by (ROI, territory voxel)
    # and vessel voxels by (ROI, segment)
    nx, ny, nz = grid.dims; h = grid.h; o = grid.origin
    cnt_l = zeros(Int, R + 1, 65536); cnt_t = zeros(Int32, R + 1, nvox); cnt_s = (zeros(Int32, R + 1, na), zeros(Int32, R + 1, nv))
    wl = VolumeWriter(joinpath(out, "labels.u16.zst"), (nx, ny, nz); level = 3); wr = VolumeWriter(joinpath(out, "rois.u8.zst"), (nx, ny, nz); T = UInt8, level = 3)
    kind = zeros(UInt8, nx, ny); A = zeros(Int32, nx, ny); D = zeros(Float32, nx, ny); RO = zeros(UInt8, nx, ny)
    sl = zeros(UInt16, nx, ny); vt = zeros(Int8, nx, ny); vj = zeros(Int32, nx, ny)
    over = Dict{Tuple{Int32, UInt16}, Int}() # lumens drawn over a territory voxel's tissue, by (voxel, label)
    act = 1
    # every voxel is labelled in the parallel kernel from array tables: a tissue voxel's class, each segment's class (0: not
    # yet known), each (pool, delay bin)'s id (0: not yet used, set in the serial pass on its first voxel)
    vcl = (zeros(Int32, na), zeros(Int32, nv)); for ((t, j), c) in vcls
        vcl[t][j] = c
    end
    plut = zeros(Int32, length(spec.pools), nbin + 1)
    for kz in 1:nz
        z = o[3] + (kz - 1) * h
        AK.foreachindex(1:(nx * ny); min_elems = 4096) do q
            ix = 1 + (q - 1) % nx; jy = 1 + (q - 1) ÷ nx; p = (o[1] + (ix - 1) * h, o[2] + (jy - 1) * h, z)
            k, aa, d = classify(p)
            r = (k == VOX_TISSUE && aa > 0) || k == VOX_VESSEL ? Int(roi(p)) : 0 # a compartment holds tissue and vessels only
            kind[q] = k; A[q] = aa; D[q] = d; RO[q] = UInt8(r)
            vt[q] = k == VOX_VESSEL ? (aa > 0 ? 1 : 2) : 0; vj[q] = abs(aa)
            sl[q] = if k == VOX_STATIC
                UInt16(aa)
            elseif k == VOX_TISSUE
                UInt16(aa > 0 ? tlabel(r, tclass[aa]) : 1000)
            elseif k == VOX_POOL
                UInt16(plut[aa, pbin(aa, d) + 1]) # 0: first use, set below
            else
                c = vcl[vt[q]][abs(aa)]; c > 0 ? UInt16(vlabel(r, vt[q], c)) : 0xffff
            end
        end
        for q in 1:(nx * ny) # first uses of a pool bin, and segments outside the drawn set (their nearest class)
            if kind[q] == VOX_POOL && sl[q] == 0
                aa = A[q]; b = pbin(aa, D[q]); plut[aa, b + 1] == 0 && (plut[aa, b + 1] = poolid(aa, D[q])); sl[q] = UInt16(plut[aa, b + 1])
            elseif sl[q] == 0xffff
                c = vclass(vt[q], vj[q]); vcl[vt[q]][vj[q]] = c; sl[q] = UInt16(vlabel(Int(RO[q]), vt[q], c))
            end
        end
        # grown lumens replace tissue where they cover a tissue voxel's centre
        while act <= length(zo) && zr[zo[act]][2] < z
            act += 1
        end
        for kk in act:length(zo)
            k = zo[kk]; zr[k][1] > z && break
            zr[k][2] < z && continue
            (t, j) = stamped[k]; pa, pb, r = segab(t, j); ab = pb .- pa; L2 = sum(abs2, ab)
            lo1 = max(1, floor(Int, (min(pa[1], pb[1]) - r - o[1]) / h) + 1); hi1 = min(nx, ceil(Int, (max(pa[1], pb[1]) + r - o[1]) / h) + 1)
            lo2 = max(1, floor(Int, (min(pa[2], pb[2]) - r - o[2]) / h) + 1); hi2 = min(ny, ceil(Int, (max(pa[2], pb[2]) + r - o[2]) / h) + 1)
            for jy in lo2:hi2, ix in lo1:hi1
                kind[ix, jy] == VOX_TISSUE || continue
                p1 = o[1] + (ix - 1) * h - pa[1]; p2 = o[2] + (jy - 1) * h - pa[2]; p3 = z - pa[3]
                u = L2 > 0 ? clamp((p1 * ab[1] + p2 * ab[2] + p3 * ab[3]) / L2, 0, 1) : 0.0
                if (p1 - u * ab[1])^2 + (p2 - u * ab[2])^2 + (p3 - u * ab[3])^2 <= r^2
                    sl[ix, jy] = UInt16(vlabel(Int(RO[ix, jy]), t, vclass(t, j))); vt[ix, jy] = t; vj[ix, jy] = j
                end
            end
        end
        for q in 1:(nx * ny)
            r1 = Int(RO[q]) + 1; cnt_l[r1, Int(sl[q]) + 1] += 1
            if vt[q] > 0
                cnt_s[vt[q]][r1, vj[q]] += 1
                kind[q] == VOX_TISSUE && A[q] > 0 && (over[(A[q], sl[q])] = get(over, (A[q], sl[q]), 0) + 1)
            elseif kind[q] == VOX_TISSUE && A[q] > 0
                cnt_t[r1, A[q]] += 1
            end
        end
        write_slice!(wl, sl); write_slice!(wr, RO)
    end
    close(wl); close(wr)
    lap("slices (classify, label, compress)")
    # ---- the tracking points' labels (a pool bin a point names is kept even if no voxel holds it)
    track = Dict{String, Any}()
    for (tn, p) in tracking
        pt = Tuple(1.0 .* p); k, aa, d = classify(pt)
        track[tn] = Dict("point_mm" => collect(Float64.(pt)), "label" => label_at(k, aa, d, Int(roi(pt))))
    end
    lap("tracking points")
    # ---- masses: every artery and vein without a painted voxel gives its iodine to its beds; each bed is spread over the
    # territory voxels that kept tissue sub-voxels (its own territory, so its iodine and flow stay where it is) and each of
    # those over its tissue sub-voxels; each segment over its painted voxels
    painted = (vec(sum(cnt_s[1]; dims = 1)) .> 0, vec(sum(cnt_s[2]; dims = 1)) .> 0)
    nta = vec(sum(cnt_t; dims = 1))
    nterm = length(terr); wsum = zeros(nterm); wkeep = zeros(nterm)
    for (v, i, w) in pairs
        wsum[i] += w; nta[v] > 0 && (wkeep[i] += w)
    end
    pk = [p for p in pairs if nta[p[1]] > 0]
    tk = Float32[wkeep[i] > 0 ? terr[i] * wkeep[i] / wsum[i] : terr[i] for i in 1:nterm]
    keptbeds(U) = [Array(tissue_iodine(Array(scans[s][3].tissue) .+ U[s], pk, tk, vox)) for s in 1:nS]
    U = [unresolved_masses(C, T.artery, T.vein; undrawn = (.!painted[1], .!painted[2])) for (_, C, T) in scans]
    M = [x .* voxml for x in keptbeds(U)]
    kflow = [vec(Array(tissue_iodine(reshape(Float32[C.q[j] for j in C.terminals], 1, :), pk, tk, vox))) .* 60 for (_, C, _) in scans]
    nl = 65536; lm = zeros(Float64, nl, nS * nT); lq = zeros(Float64, nl, nS) # label iodine (mg) and bed flow (mL/min)
    # a terminal none of whose territory voxels kept tissue: its bed goes to the lumens drawn over its territory (same place,
    # same ROI); with none there either, it is lost and the ledger shows it
    lost = zeros(Float64, nS * nT); ov = Dict{Int32, Vector{Tuple{UInt16, Int}}}()
    for ((v, l), n) in over
        push!(get!(ov, v, Tuple{UInt16, Int}[]), (l, n))
    end
    bedpairs = Dict{Int32, Vector{Tuple{Int32, Float32}}}()
    for (v, i, w) in pairs
        wkeep[i] == 0 && push!(get!(bedpairs, i, Tuple{Int32, Float32}[]), (v, w))
    end
    for (i, vw) in bedpairs
        tot = sum(w * sum(last, get(ov, v, [(0x0000, 0)])) for (v, w) in vw)
        for s in 1:nS, f in 1:nT
            mi = Float64(scans[s][3].tissue[f, i]) + Float64(U[s][f, i])
            if tot == 0
                lost[(s - 1) * nT + f] += mi
            else
                for (v, w) in vw, (l, n) in get(ov, v, Tuple{UInt16, Int}[])
                    lm[Int(l) + 1, (s - 1) * nT + f] += mi * w * n / tot
                end
            end
        end
    end
    for v in 1:nvox
        nta[v] == 0 && continue
        for r in 0:R
            c = cnt_t[r + 1, v]; c == 0 && continue
            l = tlabel(r, tclass[v]) + 1; w = c / nta[v]
            for s in 1:nS
                lq[l, s] += w * kflow[s][v] * voxml
                for f in 1:nT
                    lm[l, (s - 1) * nT + f] += w * M[s][v, f]
                end
            end
        end
    end
    # a tracked segment is painted at its blood concentration (a tracking point measures a concentration, as a pool's
    # does), not with its mass spread over its painted voxels: an imaged lumen can hold more or less than the segment
    extra = zeros(Float64, nS * nT) # the iodine a tracked segment's lumen holds beyond its mass (the ledger counts its mass)
    for t in 1:2, j in axes(cnt_s[t], 2)
        painted[t][j] || continue
        n = sum(@view cnt_s[t][:, j]); c = vclass(t, j); own = haskey(tracked, (t, j))
        for r in 0:R
            cr = cnt_s[t][r + 1, j]; cr == 0 && continue
            l = vlabel(r, t, c) + 1
            for s in 1:nS
                T = scans[s][3]; w = own ? cr * vml / vols[s][t == 1 ? j : na + j] : cr / n
                for f in 1:nT
                    m = t == 1 ? T.artery[f, j] : T.vein[f, j]
                    lm[l, (s - 1) * nT + f] += w * m; own && (extra[(s - 1) * nT + f] += (w - cr / n) * m)
                end
            end
        end
    end
    lap("masses")
    # ---- the label table
    nvl = vec(sum(cnt_l; dims = 1)); pkey = Dict(v => k for (k, v) in pid)
    ids = sort!(unique!([[l - 1 for l in 1:nl if nvl[l] > 0]; [tr["label"] for tr in values(track)]]))
    dens(m) = hasproperty(XA.Materials, Symbol(m)) ? ustrip(u"g/cm^3", getproperty(XA.Materials, Symbol(m)).density) : NaN
    lroi(l) = (r = findall(>(0), @view cnt_l[:, l + 1]); length(r) == 1 ? r[1] - 1 : -1) # -1: a pool or static label spanning ROIs
    labels = NamedTuple[]; conc = Dict{Int, Vector{Vector{Float64}}}()
    for id in ids
        n = nvl[id + 1]
        name, base, kd, curve, mbf = if id < 100
            haskey(spec.materials, id) || error("phantom: the classifier gave static label $id, which spec.materials does not name")
            ("static $id", spec.materials[id], "static", nothing, nothing)
        elseif id < 1000
            (i, b) = pkey[id]; P = spec.pools[i]
            (
                P.what === :out ? "$(P.name) +$(round(b * spec.delay_bin_ml; digits = 2)) mL" : P.name, P.what === :tissue ? spec.tissue_material : P.material,
                "pool", [[Float64(pcurve(i, b, s, f)) for f in 1:nT] for s in 1:nS], nothing,
            )
        elseif id == 1000
            ("tissue (no territory)", spec.tissue_material, "static", nothing, nothing)
        elseif id < 2000
            r, c = divrem(id - 1001, K)
            (
                "$(r == 0 ? "outside the ROIs" : roi_names[r]) · tissue class $(c + 1)", spec.tissue_material, "tissue",
                [[lm[id + 1, (s - 1) * nT + f] / (n * vml) for f in 1:nT] for s in 1:nS], [lq[id + 1, s] / (n * vml * dens(spec.tissue_material)) for s in 1:nS],
            )
        else
            g, c = divrem(id - 2001, KV); r, t = divrem(g, 2)
            (
                "$(r == 0 ? "outside the ROIs" : roi_names[r]) · $(t == 0 ? "artery" : "vein") " *
                    (haskey(tracked_class, (t + 1, c + 1)) ? "at the $(tracked_class[(t + 1, c + 1)]) tracking point" : "class $(c + 1)"),
                spec.blood_material, t == 0 ? "artery" : "vein",
                [[lm[id + 1, (s - 1) * nT + f] / (n * vml) for f in 1:nT] for s in 1:nS], zeros(nS),
            )
        end
        curve === nothing || (conc[id] = curve)
        push!(labels, (; id, name, roi = lroi(id), kind = kd, base_material = base, voxels = n, mass_g = n * vml * dens(base), mbf, curve))
    end
    lap("label table")
    # ---- the iodine ledger
    organ = [l for l in labels if l.kind in ("tissue", "artery", "vein")]
    ledger = map(1:nS) do s
        _, C, T = scans[s]; inp = T.input; ts = Float32.(spec.times)
        Ai(t) = (x = t / inp.dt + 1; k = clamp(floor(Int, x), 1, length(inp.A) - 1); w = clamp(x - k, 0, 1); (1 - w) * inp.A[k] + w * inp.A[k + 1])
        delivered = [sum(Float64, C.fp) * Ai(t) for t in ts]
        art = vec(sum(Float64.(T.artery); dims = 2)); bed = vec(sum(Float64.(T.tissue); dims = 2)); ven = vec(sum(Float64.(T.vein); dims = 2))
        held = art .+ bed .+ ven; outflow = Float64.(sinus_outflow(C, inp, ts))
        pnt = [sum(l.voxels * vml * l.curve[s][f] for l in organ; init = 0.0) - extra[(s - 1) * nT + f] for f in 1:nT]
        (;
            t_s = Float64.(spec.times), delivered_mg = delivered, artery_mg = art, bed_mg = bed, vein_mg = ven, held_mg = held,
            outflow_mg = outflow, balance_mg = delivered .- held .- outflow, painted_mg = pnt, unpainted_mg = held .- pnt,
        )
    end
    worst = [(b = maximum(abs, L.balance_mg) / max(maximum(L.delivered_mg), eps()), p = maximum(abs, L.unpainted_mg) / max(maximum(L.delivered_mg), eps())) for L in ledger]
    if check
        for (s, w) in enumerate(worst)
            w.b <= 0.01 || error("phantom: $(states[s]) iodine does not balance: delivered − held − outflow reaches $(round(100w.b; digits = 2))% of the peak delivered mass")
            w.p <= 0.01 || error("phantom: $(states[s]) painted iodine differs from the held iodine by $(round(100w.p; digits = 2))% of the peak delivered mass")
        end
    end
    lap("ledger")
    # ---- tables, then each state's ground-truth MBF map from the labels
    flow = [s => sum(C.q[j] for j in eachindex(C.tree.parent) if C.tree.parent[j] == 0) * 60 / sum(terr) for (s, C, _) in scans]
    files = Dict{String, Any}("labels" => "labels.u16.zst", "rois" => "rois.u8.zst", ("mbf_" * st => "mbf_$st.f32.zst" for st in states)...)
    toml = Dict{String, Any}(
        "rois" => Dict("names" => collect(String.(roi_names)), "groups" => Dict(String(g) => collect(Int, r) for (g, r) in roi_groups)),
        "tracking" => track, "facts" => Dict{String, Any}(string(k) => (v isa Number || v isa AbstractString ? v : string(v)) for (k, v) in facts),
        "flow_ml_min_per_ml" => Dict(String(s) => Float64(f) for (s, f) in flow),
    )
    write_bundle(out; grid, files, states, times = spec.times, labels, ledger, toml, energy_keV = spec.energy_keV)
    luts = [
        (
            lut = zeros(Float32, nl); for l in labels
                l.kind == "tissue" && (lut[l.id + 1] = Float32(l.mbf[s]))
            end; lut
        ) for s in 1:nS
    ]
    lab = DynamicPhantom.ZstdDecompressorStream(open(joinpath(out, "labels.u16.zst")))
    ws = [VolumeWriter(joinpath(out, "mbf_$st.f32.zst"), (nx, ny, nz); T = Float32, level = 3) for st in states]
    buf = Matrix{UInt16}(undef, nx, ny); m = Matrix{Float32}(undef, nx, ny)
    for _ in 1:nz
        read!(lab, buf)
        for s in 1:nS
            lut = luts[s]; AK.foreachindex(m; min_elems = 1 << 16) do i
                m[i] = lut[Int(buf[i]) + 1]
            end
            write_slice!(ws[s], m)
        end
    end
    close(lab); foreach(close, ws)
    lap("tables and MBF maps")
    return (; dir = out, ids, counts = Dict(l.id => l.voxels for l in labels), flow, ledger, worst, lost)
end
