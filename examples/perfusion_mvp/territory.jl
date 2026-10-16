# SPDX-License-Identifier: GPL-3.0-or-later
# Nearest-terminal (Voronoi) territories on a voxel grid, with a uniform bucket grid.

module Territory

import AcceleratedKernels as AK
using Random

export nearest_terminal, territory_shares, demand_field

const F32 = Float32

struct Buckets
    lo::NTuple{3, F32}
    h::F32
    dims::NTuple{3, Int}
    start::Vector{Int32}   # CSR offsets, length prod(dims) + 1
    items::Vector{Int32}
end

function Buckets(pts::Matrix{F32}, h::F32)
    lo = ntuple(d -> minimum(@view pts[d, :]) - h, 3)
    hi = ntuple(d -> maximum(@view pts[d, :]) + h, 3)
    dims = ntuple(d -> max(1, ceil(Int, (hi[d] - lo[d]) / h)), 3)
    np = size(pts, 2)
    cell = Vector{Int32}(undef, np)
    for i in 1:np
        c = ntuple(d -> clamp(floor(Int, (pts[d, i] - lo[d]) / h), 0, dims[d] - 1), 3)
        cell[i] = 1 + c[1] + dims[1] * (c[2] + dims[2] * c[3])
    end
    nc = prod(dims)
    start = zeros(Int32, nc + 1)
    for c in cell
        start[c + 1] += 1
    end
    start[1] = 1
    for c in 1:nc
        start[c + 1] += start[c]
    end
    fill_ = copy(start); items = Vector{Int32}(undef, np)
    for i in 1:np
        items[fill_[cell[i]]] = i; fill_[cell[i]] += 1
    end
    return Buckets(lo, h, dims, start, items)
end

@inline function _nearest(B::Buckets, pts, x, y, z)
    c0 = (floor(Int, (x - B.lo[1]) / B.h), floor(Int, (y - B.lo[2]) / B.h), floor(Int, (z - B.lo[3]) / B.h))
    best = typemax(F32); bi = Int32(0)
    for ring in 0:maximum(B.dims)
        # every point in a cell at Chebyshev ring ≥ ring is at least (ring − 1)·h away
        ring > 0 && bi != 0 && (F32(ring - 1) * B.h)^2 > best && break
        for k in (c0[3] - ring):(c0[3] + ring), j in (c0[2] - ring):(c0[2] + ring), i in (c0[1] - ring):(c0[1] + ring)
            max(abs(i - c0[1]), abs(j - c0[2]), abs(k - c0[3])) == ring || continue
            (0 <= i < B.dims[1] && 0 <= j < B.dims[2] && 0 <= k < B.dims[3]) || continue
            c = 1 + i + B.dims[1] * (j + B.dims[2] * k)
            for s in B.start[c]:(B.start[c + 1] - 1)
                p = B.items[s]
                d = (pts[1, p] - x)^2 + (pts[2, p] - y)^2 + (pts[3, p] - z)^2
                if d < best
                    best = d; bi = p
                end
            end
        end
    end
    return bi
end

"""
    nearest_terminal(pts, q; h = 2) -> Vector{Int32}

For every query point `q[:, v]`, the index of the nearest column of `pts` (both 3×n, mm).
"""
function nearest_terminal(pts::Matrix{F32}, q::Matrix{F32}; h::Real = 2.0f0)
    B = Buckets(pts, F32(h))
    out = Vector{Int32}(undef, size(q, 2))
    AK.foreachindex(out; min_elems = 4096) do v
        out[v] = _nearest(B, pts, q[1, v], q[2, v], q[3, v])
    end
    return out
end

"""
    territory_shares(tpts, q; voxel_ml = 1e-3) -> (terr_ml, pairs)

Tissue volume served by each terminal, robust to any terminal density. Each voxel (centers `q`,
3×V, on a unit grid) is shared equally among the terminals inside it; a voxel with no terminal
goes to its nearest terminal. `pairs` lists (voxel, terminal, fraction of the voxel), so a
per-voxel quantity is Σ fraction × (terminal value / terminal volume).
"""
function territory_shares(tpts::Matrix{F32}, q::Matrix{F32}; voxel_ml::Real = 1.0e-3)
    nt = size(tpts, 2); nv = size(q, 2)
    lo = ntuple(d -> minimum(@view q[d, :]), 3)
    key(x, y, z) = (round(Int, x - lo[1]), round(Int, y - lo[2]), round(Int, z - lo[3]))
    vid = Dict{NTuple{3, Int}, Int32}(key(q[1, v], q[2, v], q[3, v]) => Int32(v) for v in 1:nv)
    tv = zeros(Int32, nt); cnt = zeros(Int32, nv)
    for i in 1:nt
        v = get(vid, key(tpts[1, i], tpts[2, i], tpts[3, i]), Int32(0))
        tv[i] = v; v > 0 && (cnt[v] += 1)
    end
    # terminals that round to a voxel outside the mask join their nearest mask voxel
    miss = findall(==(0), tv)
    if !isempty(miss)
        nv_ = nearest_terminal(q, tpts[:, miss])
        for (k, i) in enumerate(miss)
            tv[i] = nv_[k]; cnt[nv_[k]] += 1
        end
    end
    empty_ = findall(==(0), cnt)
    near = isempty(empty_) ? Int32[] : nearest_terminal(tpts, q[:, empty_])
    pairs = Tuple{Int32, Int32, F32}[]
    sizehint!(pairs, nt + length(empty_))
    for i in 1:nt
        tv[i] > 0 && push!(pairs, (tv[i], Int32(i), 1.0f0 / cnt[tv[i]]))
    end
    for (k, v) in enumerate(empty_)
        push!(pairs, (Int32(v), near[k], 1.0f0))
    end
    terr = zeros(F32, nt)
    for (v, i, f) in pairs
        terr[i] += f * F32(voxel_ml)
    end
    return terr, pairs
end

"""
    demand_field(pos, vol; rd1g = 0.15, D = 1.2, seed = 1) -> Vector{Float32}

Relative O₂ demand (mean 1, volume-weighted) at points `pos` (3 × n, mm) holding tissue volumes
`vol`. Local O₂ use in the heart is uneven and stable, and local flow follows it
(Bassingthwaighte, Beard & Li 2001; Deussen 2001); its spread falls with piece mass as the flow
spread does, RD(m) = RD(1 g)·m^(1−D) (Bassingthwaighte et al. 1989: D = 1.2, RD = 15% at 1 g).

The field is lognormal, exp(z), with z a sum of white noise smoothed at octave scales
`octaves` (default 0.5–32 mm, so the fractal law of Bassingthwaighte et al. 1989, measured on
0.1–10 g pieces, is followed down to the 1 mm voxel).
Each octave's variance scales as σ^(6(1−D)), so the variance of a cube average falls as
L^(6(1−D)) with cube side L, as RD ∝ m^(1−D) requires. The amplitude is set so cubes of 1 cm³
(≈ 1 g) have relative dispersion `rd1g`.
"""
function demand_field(pos::AbstractMatrix, vol::AbstractVector; rd1g::Real = 0.15, D::Real = 1.2, seed::Integer = 1, octaves = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0))
    n = size(pos, 2)
    lo = floor.(Int, vec(minimum(pos; dims = 2))) .- 2; hi = ceil.(Int, vec(maximum(pos; dims = 2))) .+ 2
    dims = Tuple(hi .- lo .+ 1)
    rng = Xoshiro(seed)
    # RD ∝ m^(1−D) means Var(cube average) ∝ L^(6(1−D)) for side L (m ∝ L³); octave components of
    # correlation length σ with variance ∝ σ^(6(1−D)) give that.
    α = 6 * (1 - D)
    z = zeros(Float32, dims)
    for σ in octaves
        f = randn(rng, Float32, dims)
        _blur!(f, Float32(σ))
        f ./= sqrt(sum(abs2, f) / length(f))
        z .+= Float32(σ^(α / 2)) .* f
    end
    zp = Vector{Float32}(undef, n)
    for k in 1:n
        I = CartesianIndex(Tuple(clamp.(round.(Int, pos[:, k]) .- lo .+ 1, 1, collect(dims))))
        zp[k] = z[I]
    end
    zp .-= sum(zp .* vol) / sum(vol)
    # amplitude: relative dispersion of 1 cm³ cubes equals rd1g
    cube = [Tuple(floor.(Int, pos[:, k] ./ 10)) for k in 1:n]
    function rd(a)
        m = exp.(a .* zp); m ./= sum(m .* vol) / sum(vol)
        acc = Dict{NTuple{3, Int}, Tuple{F32, F32}}()
        for k in 1:n
            s, w = get(acc, cube[k], (0.0f0, 0.0f0)); acc[cube[k]] = (s + F32(m[k] * vol[k]), w + F32(vol[k]))
        end
        wt = [w for (s, w) in values(acc)]; mv = [s / w for (s, w) in values(acc)]
        μ = sum(mv .* wt) / sum(wt)
        return sqrt(sum(wt .* (mv .- μ) .^ 2) / sum(wt)) / μ, m
    end
    a0, a1 = 0.0f0, 10.0f0
    for _ in 1:40
        am = (a0 + a1) / 2
        rd(am)[1] < rd1g ? (a0 = am) : (a1 = am)
    end
    return Float32.(rd((a0 + a1) / 2)[2])
end

# separable Gaussian blur, σ in grid cells
function _blur!(f::Array{Float32, 3}, σ::Float32)
    r = max(1, ceil(Int, 3σ)); w = Float32[exp(-(i / σ)^2 / 2) for i in -r:r]; w ./= sum(w)
    for d in 1:3
        g = similar(f); n = size(f, d)
        CI = CartesianIndices(f)
        AK.foreachindex(g; min_elems = 4096) do li
            I = CI[li]; i = I[d]; s = 0.0f0
            for o in -r:r
                ii = clamp(i + o, 1, n)
                J = Base.setindex(Tuple(I), ii, d)
                s += w[o + r + 1] * f[CartesianIndex(J)]
            end
            g[I] = s
        end
        f .= g
    end
    return f
end

end # module
