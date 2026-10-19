# SPDX-License-Identifier: GPL-3.0-or-later
# CT frames of any perfused organ from a solved circuit (heart, lungs, brain): iodine in the tissue
# voxels each terminal perfuses, grown vessels drawn into display planes as tubes, and a plane
# renderer in Hounsfield units. Nothing here knows which organ it is: the trees, the voxels and the
# base attenuation of each pixel come from the caller.
#
# Typical use (xcat-data/work/lung/lung_frames.jl does this for the lungs):
#   am, vm, tm = frame_masses(C, inp, times)                 # mg per artery, vein and terminal bed
#   P = Plane(centre, (1, 0, 0), (0, 1, 0), 520, 400; px = 0.5)
#   G = sample_grid(P, vox, lower)                           # organ voxels under each pixel
#   S = stamp_vessels(P, C.tree, C.veins)                    # lumens and sub-pixel vessels
#   hu = render_plane(P, base, tissue_iodine(tm, pairs, terr, G), S, hcat(am, vm),
#                     blood_volumes(C); hu_per_mg, blood_hu)
#
# Geometry (planes, samples, stamps) is built once on the host; the per-frame kernels take `to`.

module OrganCT

import AcceleratedKernels as AK
using Atomix: @atomic
using ..Flow: Circuit, InputCurve, VTree, frame_masses, _lengths, _bfs

export Plane, pixel_point, GridSample, sample_grid, tissue_iodine, region_iodine, VesselStamp,
    stamp_vessels, blood_volumes, unresolved_masses, render_plane

const F32 = Float32

@inline _dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _unit(a) = (s = sqrt(_dot(a, a)); (a[1] / s, a[2] / s, a[3] / s))
@inline _cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

# ---------------------------------------------------------------------------------------------
# Planes
# ---------------------------------------------------------------------------------------------

"""
    Plane(centre, e1, e2, n1, n2; px = 0.5, slab = px)

A display plane of `n1 × n2` pixels of `px` mm, centred on `centre` (mm). `e1` runs along the
columns (left → right in the image) and `e2` along the rows (top → bottom); `e2` is made
orthogonal to `e1`. Pixel `(i, j)` (linear index `i + n1·(j − 1)`) has its centre at
`centre + (i − ½ − n1/2)·px·e1 + (j − ½ − n2/2)·px·e2`. Vessels too thin to draw as tubes are
averaged over a slab `slab` mm thick.
"""
struct Plane
    c::NTuple{3, F32}
    e1::NTuple{3, F32}
    e2::NTuple{3, F32}
    nrm::NTuple{3, F32}        # e1 × e2
    n1::Int
    n2::Int
    px::F32
    slab::F32
end
function Plane(centre, e1, e2, n1::Integer, n2::Integer; px::Real = 0.5, slab::Real = px)
    u = _unit(F32.(Tuple(e1))); v = F32.(Tuple(e2)); d = _dot(v, u)
    v = _unit((v[1] - d * u[1], v[2] - d * u[2], v[3] - d * u[3]))
    return Plane(F32.(Tuple(centre)), u, v, _cross(u, v), Int(n1), Int(n2), F32(px), F32(slab))
end

"""
    pixel_point(P, x, y) -> (x, y, z) mm

The point of plane `P` at `x` columns and `y` rows from its top-left corner (pixel `(i, j)`
has its centre at `x = i − ½`, `y = j − ½`).
"""
@inline function pixel_point(P::Plane, x, y)
    a = (F32(x) - P.n1 / 2.0f0) * P.px; b = (F32(y) - P.n2 / 2.0f0) * P.px
    return (P.c[1] + a * P.e1[1] + b * P.e2[1], P.c[2] + a * P.e1[2] + b * P.e2[2], P.c[3] + a * P.e1[3] + b * P.e2[3])
end

# in-plane pixel coordinates (x, y) and signed distance from the plane of a point
@inline function _project(P::Plane, p)
    d = (p[1] - P.c[1], p[2] - P.c[2], p[3] - P.c[3])
    return _dot(d, P.e1) / P.px + P.n1 / 2.0f0, _dot(d, P.e2) / P.px + P.n2 / 2.0f0, _dot(d, P.nrm)
end

# ---------------------------------------------------------------------------------------------
# Tissue: voxels → pixels, terminal beds → voxels
# ---------------------------------------------------------------------------------------------

"""
    GridSample

Trilinear weights from organ voxels to the pixels of a plane: pixel `k` takes
`Σₛ w[s, k] · value(v[s, k])` over its 8 neighbouring voxel centres; `v` is the organ voxel
number (0 = none).
"""
struct GridSample
    v::Matrix{Int32}           # (8, npixels)
    w::Matrix{F32}
end

"""
    sample_grid(P, vox, lower; h = 1) -> GridSample

Organ voxels under every pixel of `P`. `vox` is a 3-D array over a grid of `h` mm voxels whose
voxel `I` is centred at `lower + (I − 1)·h` mm, and `vox[I]` is that voxel's organ number (its
column in the voxel centres given to `territory_shares`), 0 outside the organ. A pixel whose
nearest voxel is in the organ takes trilinear weights over its neighbours in the organ, so tissue
values vary smoothly instead of in voxel blocks. Other pixels get no weights.
"""
function sample_grid(P::Plane, vox::AbstractArray{<:Integer, 3}, lower; h::Real = 1)
    npx = P.n1 * P.n2; V = zeros(Int32, 8, npx); W = zeros(F32, 8, npx)
    lo = F32.(Tuple(lower)); hh = F32(h); dims = size(vox)
    AK.foreachindex(1:npx; min_elems = 4096) do k
        i = (k - 1) % P.n1 + 1; j = (k - 1) ÷ P.n1 + 1
        p = pixel_point(P, i - 0.5f0, j - 0.5f0)
        u = ((p[1] - lo[1]) / hh + 1.0f0, (p[2] - lo[2]) / hh + 1.0f0, (p[3] - lo[3]) / hh + 1.0f0)
        c = (round(Int, u[1]), round(Int, u[2]), round(Int, u[3]))
        inside = all(d -> 1 <= c[d] <= dims[d], 1:3) && vox[c[1], c[2], c[3]] > 0
        if inside
            i0 = (floor(Int, u[1]), floor(Int, u[2]), floor(Int, u[3])); f = (u[1] - i0[1], u[2] - i0[2], u[3] - i0[3])
            tot = 0.0f0
            for s in 0:7
                o = (s & 1, (s >> 1) & 1, (s >> 2) & 1); I = (i0[1] + o[1], i0[2] + o[2], i0[3] + o[3])
                all(d -> 1 <= I[d] <= dims[d], 1:3) || continue
                v = vox[I[1], I[2], I[3]]; v > 0 || continue
                wt = (o[1] == 1 ? f[1] : 1.0f0 - f[1]) * (o[2] == 1 ? f[2] : 1.0f0 - f[2]) * (o[3] == 1 ? f[3] : 1.0f0 - f[3])
                V[s + 1, k] = v; W[s + 1, k] = wt; tot += wt
            end
            if tot > 0
                for s in 1:8
                    W[s, k] /= tot
                end
            end
        end
    end
    return GridSample(V, W)
end

# the (voxel row, terminal, weight / territory) entries of `pairs` whose voxel has a row
function _bed_tables(pairs, terr, row)
    keep = [k for (k, (v, i, w)) in enumerate(pairs) if v <= length(row) && row[v] > 0]
    pv = Int32[row[pairs[k][1]] for k in keep]; pi = Int32[pairs[k][2] for k in keep]
    pw = F32[pairs[k][3] / terr[pairs[k][2]] for k in keep]
    return pv, pi, pw
end

"""
    tissue_iodine(tm, pairs, terr, voxels; to = identity) -> Matrix (length(voxels), ntimes)
    tissue_iodine(tm, pairs, terr, G::GridSample; to = identity) -> Matrix (npixels, ntimes)
    tissue_iodine(c::Circuit, inp::InputCurve, times, pairs, terr, voxels; to = identity)

Iodine in tissue voxels, mg per mL of voxel, at every frame. Each terminal's lumped bed (`tm`,
the third result of `frame_masses`: mg, ntimes × nterminals) is spread over its territory, so
voxel `v` holds `Σ (w / terr[i]) · tm[f, i]` over its entries `(v, i, w)` in `pairs`, and `terr`
are the territory volumes (mL), both from `territory_shares`. `voxels` selects the voxels
(columns of the voxel centres given to `territory_shares`), and rows follow its order. With a
[`GridSample`](@ref) in place of `voxels` the result is per pixel of that plane. With a circuit and
an input curve, `tm` comes from `frame_masses(c, inp, times)`.
"""
function tissue_iodine(tm::AbstractMatrix, pairs, terr, voxels::AbstractVector{<:Integer}; to = identity)
    nf = size(tm, 1); row = zeros(Int32, maximum(voxels; init = 0))
    for (k, v) in enumerate(voxels)
        row[v] = k
    end
    pv, pi, pw = to.(_bed_tables(pairs, terr, row)); np = length(pv)
    out = to(zeros(F32, length(voxels), nf))
    AK.foreachindex(1:(np * nf); backend = AK.get_backend(out), min_elems = 4096) do e
        p = (e - 1) % np + 1; f = (e - 1) ÷ np + 1
        @atomic out[pv[p], f] += pw[p] * tm[f, pi[p]]
    end
    return out
end
function tissue_iodine(tm::AbstractMatrix, pairs, terr, G::GridSample; to = identity)
    voxels = sort!(unique(filter(>(0), vec(G.v)))); nf = size(tm, 1)
    X = tissue_iodine(tm, pairs, terr, voxels; to)
    row = zeros(Int32, maximum(voxels; init = 0))
    for (k, v) in enumerate(voxels)
        row[v] = k
    end
    gv = to(Int32[v > 0 ? row[v] : 0 for v in G.v]); gw = to(G.w); npx = size(G.v, 2)
    out = to(zeros(F32, npx, nf))
    AK.foreachindex(1:(npx * nf); backend = AK.get_backend(out), min_elems = 4096) do e
        k = (e - 1) % npx + 1; f = (e - 1) ÷ npx + 1; s_ = 0.0f0
        for s in 1:8
            r = gv[s, k]; r > 0 && (s_ += gw[s, k] * X[r, f])
        end
        out[k, f] = s_
    end
    return out
end
tissue_iodine(c::Circuit, inp::InputCurve, times, pairs, terr, voxels; to = identity) =
    tissue_iodine(frame_masses(c, inp, times; to)[3], pairs, terr, voxels; to)

"""
    region_iodine(tm, pairs, terr, voxels; to = identity) -> Vector (ntimes)

Mean tissue iodine (mg/mL) over the voxels `voxels` at every frame: the mean of the rows of
[`tissue_iodine`](@ref), without building them (one weight per terminal, then `tm` times it).
"""
function region_iodine(tm::AbstractMatrix, pairs, terr, voxels::AbstractVector{<:Integer}; to = identity)
    row = zeros(Int32, maximum(voxels; init = 0))
    for v in voxels
        row[v] = 1
    end
    _, pi, pw = to.(_bed_tables(pairs, terr, row)); s = to(zeros(F32, size(tm, 2)))
    AK.foreachindex(1:length(pi); backend = AK.get_backend(s), min_elems = 4096) do p
        @atomic s[pi[p]] += pw[p]
    end
    return vec(Array(tm * s)) ./ F32(length(voxels))
end

# ---------------------------------------------------------------------------------------------
# Vessels
# ---------------------------------------------------------------------------------------------

"""
    VesselStamp

Grown vessels in one plane. `seg[k]` is the vessel whose lumen covers pixel `k` (0 = none) and
`cover[k]` the covered fraction of the pixel (2 × 2 sub-samples). Vessels thinner than the tube
threshold are deposits instead: pixel `dpix[e]` gains `dw[e] · mass` mg/mL from vessel `dseg[e]`.
"""
struct VesselStamp
    seg::Vector{Int32}
    cover::Vector{F32}
    dpix::Vector{Int32}
    dseg::Vector{Int32}
    dw::Vector{F32}            # 1/mL
end

# distance from point p to the segment a–b
@inline function _segdist(p, a, b)
    ab = (b[1] - a[1], b[2] - a[2], b[3] - a[3]); ap = (p[1] - a[1], p[2] - a[2], p[3] - a[3])
    L2 = _dot(ab, ab); t = L2 > 0 ? clamp(_dot(ap, ab) / L2, 0.0f0, 1.0f0) : 0.0f0
    d = (ap[1] - t * ab[1], ap[2] - t * ab[2], ap[3] - t * ab[3])
    return sqrt(_dot(d, d))
end

# Count (emit = nothing) or emit the tube candidates of segment j: (pixel, 4-bit sub-sample mask, d/r)
function _tube(P::Plane, a, b, r, emit)
    xa, ya, da = _project(P, a); xb, yb, db = _project(P, b)
    (min(da, db) > r || max(da, db) < -r) && return 0
    rp = r / P.px
    i0 = max(1, floor(Int, min(xa, xb) - rp)); i1 = min(P.n1, ceil(Int, max(xa, xb) + rp) + 1)
    j0 = max(1, floor(Int, min(ya, yb) - rp)); j1 = min(P.n2, ceil(Int, max(ya, yb) + rp) + 1)
    k = 0
    for j in j0:j1, i in i0:i1
        m = 0x00
        for s in 0:3
            p = pixel_point(P, i - 0.75f0 + 0.5f0 * (s & 1), j - 0.75f0 + 0.5f0 * (s >> 1))
            _segdist(p, a, b) <= r && (m |= UInt8(1) << s)
        end
        m == 0x00 && continue
        k += 1
        emit === nothing || emit(k, i + P.n1 * (j - 1), m, _segdist(pixel_point(P, i - 0.5f0, j - 0.5f0), a, b) / r)
    end
    return k
end

# Count or emit the sub-pixel samples of segment j inside the slab: (pixel, weight 1/mL)
function _deposit(P::Plane, a, b, step, emit)
    _, _, da = _project(P, a); _, _, db = _project(P, b); h = P.slab / 2
    (min(da, db) > h || max(da, db) < -h) && return 0
    L = sqrt((b[1] - a[1])^2 + (b[2] - a[2])^2 + (b[3] - a[3])^2); ns = max(1, ceil(Int, L / step))
    w = (1.0f0 / ns) / (P.px * P.px * P.slab * 1.0f-3)      # mg per mL of the pixel's voxel
    k = 0
    for s in 1:ns
        u = (s - 0.5f0) / ns
        x, y, d = _project(P, (a[1] + u * (b[1] - a[1]), a[2] + u * (b[2] - a[2]), a[3] + u * (b[3] - a[3])))
        abs(d) <= h || continue
        i = floor(Int, x) + 1; j = floor(Int, y) + 1
        (1 <= i <= P.n1 && 1 <= j <= P.n2) || continue
        k += 1
        emit === nothing || emit(k, i + P.n1 * (j - 1), w)
    end
    return k
end

"""
    stamp_vessels(P, trees...; rtube = P.px / 2, step = 0.1, deposit = true, weight = nothing) -> VesselStamp

Draws the segments of `trees` (`VTree`s, numbered one tree after the other, the order of `hcat`
of their mass matrices) into plane `P`. A segment of radius ≥ `rtube` mm is a tube: its lumen
holds the pixel sub-samples within its radius of its axis (3-D distance), and each pixel shows the
segment that covers most of it. A thinner segment is sampled every `step` mm along its axis, and
each sample inside the slab deposits its share of the segment's mass into its pixel, per mL of
the pixel's voxel (`px² · slab`), so its iodine is conserved. With `deposit = false` thinner segments draw nothing (their
iodine then belongs in the tissue: [`unresolved_masses`](@ref)). `weight` (per pixel, 0-1) is the share of a pixel a lumen
may replace, e.g. the share of its sub-samples that are tissue, as a phantom replaces only tissue voxels.
"""
function stamp_vessels(P::Plane, trees::VTree...; rtube::Real = P.px / 2, step::Real = 0.1, deposit::Bool = true, weight = nothing)
    rt = F32(rtube); st = F32(step)
    tp = Int32[]; ts = Int32[]; tm = UInt8[]; tq = F32[]; dp = Int32[]; ds = Int32[]; dw = F32[]
    off = 0
    for t in trees
        n = length(t.r); A = t.a; B = t.b; R = t.r
        ends(j) = ((A[1, j], A[2, j], A[3, j]), (B[1, j], B[2, j], B[3, j]))
        visit(j, et, ed) = (e = ends(j); R[j] >= rt ? _tube(P, e[1], e[2], R[j], et) : _deposit(P, e[1], e[2], st, ed))
        cnt = zeros(Int32, n)
        AK.foreachindex(1:n; min_elems = 1024) do j
            cnt[j] = visit(j, nothing, nothing)
        end
        tube = R .>= rt; ct = ifelse.(tube, cnt, Int32(0)); cd = ifelse.(tube, Int32(0), cnt)
        ot = AK.accumulate(+, ct; init = Int32(0)) .- ct .+ Int32(length(tp))
        od = AK.accumulate(+, cd; init = Int32(0)) .- cd .+ Int32(length(dp))
        append!(tp, zeros(Int32, sum(ct))); append!(ts, zeros(Int32, sum(ct))); append!(tm, zeros(UInt8, sum(ct))); append!(tq, zeros(F32, sum(ct)))
        append!(dp, zeros(Int32, sum(cd))); append!(ds, zeros(Int32, sum(cd))); append!(dw, zeros(F32, sum(cd)))
        AK.foreachindex(1:n; min_elems = 1024) do j
            cnt[j] > 0 || return
            visit(
                j, (k, q, m, dr) -> (e = ot[j] + k; tp[e] = q; ts[e] = off + j; tm[e] = m; tq[e] = dr),
                (k, q, w) -> (e = od[j] + k; dp[e] = q; ds[e] = off + j; dw[e] = w)
            )
        end
        off += n
    end
    # each pixel: the segment covering most sub-samples (ties: nearest axis), cover = union of all
    npx = P.n1 * P.n2; seg = zeros(Int32, npx); bits = zeros(UInt8, npx); best = zeros(Int, npx); bq = fill(Inf32, npx)
    for e in eachindex(tp)
        q = tp[e]; c = count_ones(tm[e]); bits[q] |= tm[e]
        if c > best[q] || (c == best[q] && tq[e] < bq[q])
            best[q] = c; bq[q] = tq[e]; seg[q] = ts[e]
        end
    end
    cover = F32.(count_ones.(bits)) ./ 4.0f0
    weight === nothing || (cover .*= F32.(weight); seg[cover .== 0] .= 0)
    deposit || (dp = Int32[]; ds = Int32[]; dw = F32[])
    return VesselStamp(seg, cover, dp, ds, dw)
end

# a tree's segments grouped by depth (roots first), so each level can be processed in parallel after its parents
function _levels(parent::AbstractVector{<:Integer})
    order = _bfs(Vector{Int32}(parent)); d = zeros(Int32, length(parent))
    for j in order
        p = parent[j]; d[j] = p == 0 ? 1 : d[p] + 1
    end
    L = [Int32[] for _ in 1:maximum(d; init = 0)]
    for j in order
        push!(L[d[j]], Int32(j))
    end
    return L
end

"""
    unresolved_masses(C::Circuit, artery, vein; rmin) -> Matrix (ntimes × nterminals), mg

The iodine in the grown arteries and veins of `C` narrower than `rmin` mm (`frame_masses`' `artery` and `vein`), given to
the terminal beds each one feeds or drains in proportion to their flow (an artery: a bed's forward flow over the forward
flow of all beds below it, so a reversed bed gets none; a vein: terminal plasma over the plasma routed through it). Added to the beds' masses, it keeps in the tissue the iodine of the
vessels a picture or phantom does not draw as lumens. Imaged segments (`isbase`) are left out: the scan draws them.
"""
function unresolved_masses(C::Circuit, am::AbstractMatrix, vm::AbstractMatrix; rmin::Real)
    nf = size(am, 1); t = C.tree; v = C.veins; nt = length(C.terminals); rm = F32(rmin)
    H = zeros(F32, nf, length(t.r))                 # Σ over thin ancestors (and itself) of mass / forward flow below
    par = t.parent; q = C.q; tr = t.r; tb = t.isbase; term = C.terminals
    order = _bfs(Vector{Int32}(par)); qf = zeros(F32, length(t.r))   # forward bed flow below each artery
    for i in 1:nt
        qf[term[i]] = max(q[term[i]], 0.0f0)
    end
    for h in length(order):-1:1
        j = order[h]; p = par[j]; p > 0 && (qf[p] += qf[j])
    end
    for lv in _levels(t.parent)                     # a level's segments are independent; each holds its frames contiguous
        AK.foreachindex(lv; min_elems = 2048) do h
            j = lv[h]; p = par[j]; thin = tr[j] < rm && !tb[j] && q[j] > 0 && qf[j] > 0
            for f in 1:nf
                H[f, j] = (p == 0 ? 0.0f0 : H[f, p]) + (thin ? am[f, j] / qf[j] : 0.0f0)
            end
        end
    end
    vorder = _bfs(v.parent); qr = zeros(F32, length(v.r))  # plasma routed through each vein
    for i in 1:nt
        qr[C.vleaf[i]] += C.fp[i]
    end
    for h in length(vorder):-1:1
        k = vorder[h]; p = v.parent[k]; p > 0 && (qr[p] += qr[k])
    end
    Hv = zeros(F32, nf, length(v.r)); vpar = v.parent; vr = v.r; vb = v.isbase
    for lv in _levels(v.parent)
        AK.foreachindex(lv; min_elems = 2048) do h
            k = lv[h]; p = vpar[k]; thin = vr[k] < rm && !vb[k] && qr[k] > 0
            for f in 1:nf
                Hv[f, k] = (p == 0 ? 0.0f0 : Hv[f, p]) + (thin ? vm[f, k] / qr[k] : 0.0f0)
            end
        end
    end
    X = zeros(F32, nf, nt); vl = C.vleaf; fpv = C.fp
    AK.foreachindex(1:nt; min_elems = 4096) do i
        j = term[i]; k = vl[i]; qa = qf[j]; fp = fpv[i]
        for f in 1:nf
            X[f, i] = qa * H[f, j] + fp * Hv[f, k]
        end
    end
    return X
end

"""
    blood_volumes(C::Circuit) -> Vector (mL)

The volume to divide every arterial and venous segment's `frame_masses` mass by for the whole-blood iodine concentration
in its lumen, in the numbering of [`stamp_vessels`](@ref). Iodine is in the plasma. A grown segment, drawn at its own
radius: its blood volume as the circuit holds it (geometric lumen times its tone's `xd²`), so the concentration carries
its own tube hematocrit. An imaged segment (`isbase`), drawn at the lumen the scan shows: its plasma volume over the
plasma fraction of discharged blood, `vp / (1 - hd)`, the concentration of blood at the size the scan shows, whatever
radius growth left the segment.
"""
function blood_volumes(C::Circuit)
    pf = 1.0f0 - C.phys.hd
    va = _lengths(C.tree)[2] .* C.xd .^ 2; vv = _lengths(C.veins)[2]
    for j in eachindex(va)
        C.tree.isbase[j] && (va[j] = C.vp_a[j] / pf)
    end
    for k in eachindex(vv)
        C.veins.isbase[k] && (vv[k] = C.vp_v[k] / pf)
    end
    return vcat(va, vv)
end

# ---------------------------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------------------------

"""
    render_plane(P, base, tissue, S, masses, vol; hu_per_mg, blood_hu, to = identity) -> Array (n1, n2, ntimes)

HU of every pixel of `P` at every frame. A pixel shows `base[k] + hu_per_mg · tissue[k, f]`
(`base` in HU; `tissue` in mg/mL, npixels × ntimes, or `nothing`) plus its sub-pixel vessel
deposits. Where a tube covers the fraction `c` of a pixel, that fraction shows its lumen,
`blood_hu + hu_per_mg · masses[f, g] / vol[g]`, instead. `masses` (ntimes × nsegments) and `vol`
(mL) follow the numbering of `S` ([`stamp_vessels`](@ref)); for a circuit `C` they are
`hcat(am, vm)` of `frame_masses` and [`blood_volumes`](@ref)`(C)`.
"""
function render_plane(
        P::Plane, base::AbstractVector, tissue, S::VesselStamp, masses::AbstractMatrix, vol::AbstractVector;
        hu_per_mg::Real, blood_hu::Real, to = identity
    )
    npx = P.n1 * P.n2; nf = size(masses, 1); k_ = F32(hu_per_mg); bh = F32(blood_hu)
    length(base) == npx || throw(DimensionMismatch("base has $(length(base)) pixels, the plane $npx"))
    m = to(masses); v = to(F32.(vol)); b = to(F32.(base)); seg = to(S.seg); cov = to(S.cover)
    tis = tissue === nothing ? to(zeros(F32, npx, nf)) : tissue
    dep = to(zeros(F32, npx, nf)); dpix = to(S.dpix); dseg = to(S.dseg); dw = to(S.dw); nd = length(S.dpix)
    AK.foreachindex(1:(nd * nf); backend = AK.get_backend(dep), min_elems = 4096) do e
        p = (e - 1) % nd + 1; f = (e - 1) ÷ nd + 1; x = m[f, dseg[p]]
        x > 0 && (@atomic dep[dpix[p], f] += dw[p] * x)
    end
    img = to(zeros(F32, npx, nf))
    AK.foreachindex(img; min_elems = 4096) do e
        k = (e - 1) % npx + 1; f = (e - 1) ÷ npx + 1; g = seg[k]
        h = b[k] + k_ * (tis[k, f] + dep[k, f])
        g > 0 && (h = (1.0f0 - cov[k]) * h + cov[k] * (bh + k_ * m[f, g] / v[g]))
        img[e] = h
    end
    return reshape(Array(img), P.n1, P.n2, nf)
end

end # module
