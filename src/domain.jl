# SPDX-License-Identifier: GPL-3.0-or-later
# Ports ImageMaskDomainCtrl (src/MyocardialDomain.h) and CoronaryGeometry::getImageDistance
# from HCT-CCO commit 9bc141f. The distance map reproduces DGtal's exact Euclidean distance
# transformation (ExactPredicateLpSeparableMetric<.,2>): pow(squared distance, 0.5) to the
# nearest background voxel inside the image domain.

const MASK_THRESHOLD = 0x80   # myMaskThreshold = 128

# One line of the separable exact Lp distance transform: out[x] = min_i g[i] + |x - i|^p over
# integer raw distances (DGtal's RawValue). Translates of |t|^p cross at most once, so the
# minimum is exact; lines are short (<= a few hundred voxels), so a direct scan is used.
function _lp_line!(out, g, n, p, INF)
    @inbounds for x in 1:n
        best = INF
        for i in 1:n
            gi = g[i]
            gi >= INF && continue
            v = gi + abs(x - i)^p
            v < best && (best = v)
        end
        out[x] = best
    end
    return
end

"""
    distance_map(T, img; p = 3) -> Array{T,3}

DGtal DistanceTransformation over the image domain with foreground `128 < v <= 255`, using
`ExactPredicateLpSeparableMetric<Space, p>`: each voxel gets `pow(raw, 1/p)` where `raw` is
the minimum over background voxels of `sum(|dx_k|^p)` (0 on background).

Upstream (CoronaryGeometry::getImageDistance) instantiates the metric with p equal to the
space dimension, so in 3D its "L2Metric" is the L3 metric; `p = 3` reproduces that, and
`p = 2` gives the true Euclidean distance the HCT-CCO paper describes.
"""
function distance_map(::Type{T}, img::AbstractArray{UInt8, 3}; p::Integer = 3) where {T}
    INF = Int64(1) << 50
    D = Int64[v > MASK_THRESHOLD ? INF : 0 for v in img]
    nmax = maximum(size(img))
    for dim in 1:3
        n = size(D, dim)
        others = [size(D, d) for d in 1:3 if d != dim]
        Threads.@threads for j in 1:prod(others)
            g = Vector{Int64}(undef, nmax); out = similar(g)
            a, b = (j - 1) % others[1] + 1, (j - 1) ÷ others[1] + 1
            idx(u) = dim == 1 ? (u, a, b) : dim == 2 ? (a, u, b) : (a, b, u)
            for u in 1:n
                g[u] = D[idx(u)...]
            end
            _lp_line!(out, g, n, p, INF)
            for u in 1:n
                D[idx(u)...] = out[u]
            end
        end
    end
    dist = Array{T, 3}(undef, size(img))
    for i in eachindex(D)
        dist[i] = upow(T(D[i]), one(T) / T(p))
    end
    return dist
end

"""
    MaskDomain{T}

Image-mask domain controller (upstream ImageMaskDomainCtrl<3>): the mask, its exact distance
map, and the clearance `min_dist` (`-m/--minDistanceToBorder`).
"""
mutable struct MaskDomain{T}
    img::Array{UInt8, 3}
    dist::Array{T, 3}
    lower::NTuple{3, Int32}
    upper::NTuple{3, Int32}
    center::NTuple{3, Int32}
    min_dist::T
    nbtry::Int
    cache::Vector{NTuple{3, Int32}}
    cache_dist::T
end

"""
    read_vol_union(paths) -> (img, lower)

The union of several masks (voxelwise maximum) on the box that covers all of them, for a
forest grown in one domain (for example the left and right coronary territories).
"""
function read_vol_union(paths::AbstractVector{<:AbstractString})
    vols = [read_vol(p) for p in paths]
    lo = ntuple(d -> minimum(v[2][d] for v in vols), 3)
    hi = ntuple(d -> maximum(v[2][d] + Int32(size(v[1], d)) - Int32(1) for v in vols), 3)
    img = zeros(UInt8, (hi .- lo .+ Int32(1))...)
    for (v, l) in vols
        o = l .- lo
        for k in axes(v, 3), j in axes(v, 2), i in axes(v, 1)
            I = (i + o[1], j + o[2], k + o[3])
            img[I...] = max(img[I...], v[i, j, k])
        end
    end
    return img, Int32.(lo)
end

function MaskDomain{T}(
        path::Union{AbstractString, AbstractVector{<:AbstractString}};
        root::Union{Nothing, NTuple{3, <:Integer}} = nothing,
        min_dist::Real = 5.0, nbtry::Integer = 100, metric::Symbol = :upstream
    ) where {T}
    metric in (:upstream, :euclidean) || throw(ArgumentError("metric must be :upstream (L3, as HCT-CCO) or :euclidean"))
    img, lower = path isa AbstractString ? read_vol(path) : read_vol_union(path)
    upper = lower .+ Int32.(size(img)) .- Int32(1)
    dist = distance_map(T, img; p = metric === :euclidean ? 2 : 3)
    d = MaskDomain{T}(
        img, dist, lower, upper, (Int32(0), Int32(0), Int32(0)), T(min_dist),
        nbtry, NTuple{3, Int32}[], T(-1)
    )
    if root !== nothing && is_inside(d, Int32.(root))
        d.center = Int32.(root)
    else
        d.center = max_distant_point(d)
    end
    return d
end

@inline in_bounds(d, p::NTuple{3, Int32}) =
    d.lower[1] <= p[1] <= d.upper[1] && d.lower[2] <= p[2] <= d.upper[2] && d.lower[3] <= p[3] <= d.upper[3]
@inline _ix(d, p) = (p[1] - d.lower[1] + 1, p[2] - d.lower[2] + 1, p[3] - d.lower[3] + 1)
@inline img_at(d, p::NTuple{3, Int32}) = @inbounds d.img[_ix(d, p)...]
@inline dist_at(d, p::NTuple{3, Int32}) = @inbounds d.dist[_ix(d, p)...]

# ImageMaskDomainCtrl::isInside: strictly above the threshold.
@inline is_inside(d, p::NTuple{3, Int32}) = in_bounds(d, p) && img_at(d, p) > MASK_THRESHOLD

# std::min / std::max on doubles.
@inline cmin(a, b) = b < a ? b : a
@inline cmax(a, b) = a < b ? b : a

# ImageMaskDomainCtrl::interpolateDistance (trilinear, corners in bit order x, y, z).
@inline function _interp_axis(x::T, lo::Int32, up::Int32) where {T}
    l = T(lo); h = T(up)
    coord = cmin(cmax(x, l), h)
    base = floor(coord)
    i0 = itrunc(base)
    if i0 < lo
        i0 = lo; base = T(i0)
    end
    if i0 >= up
        i0 = up == lo ? up : up - Int32(1)
        base = T(i0)
    end
    i1 = min(i0 + Int32(1), up)
    t = coord - base
    i1 == i0 && (t = zero(T))
    return i0, i1, t
end

@inline function interpolate_distance(d, pt::P3{T}) where {T}
    lo, up = d.lower, d.upper
    @inbounds for i in 1:3
        (pt[i] < T(lo[i]) || pt[i] > T(up[i])) && return T(-Inf)
    end
    a0, a1, fa = _interp_axis(pt[1], lo[1], up[1])
    b0, b1, fb = _interp_axis(pt[2], lo[2], up[2])
    c0, c1, fc = _interp_axis(pt[3], lo[3], up[3])
    value = zero(T)
    for mask in 0:7
        w = one(T)
        ci = (mask & 1) != 0 ? (w *= fa; a1) : (w *= (one(T) - fa); a0)
        cj = (mask & 2) != 0 ? (w *= fb; b1) : (w *= (one(T) - fb); b0)
        ck = (mask & 4) != 0 ? (w *= fc; c1) : (w *= (one(T) - fc); c0)
        value += w * dist_at(d, (ci, cj, ck))
    end
    return value
end

sample_distance(d::MaskDomain{T}, p::NTuple{3, Int32}) where {T} = interpolate_distance(d, (T(p[1]), T(p[2]), T(p[3])))

# ImageMaskDomainCtrl::randomPoint
# Admissible voxels, listed once; the fallback when rejection sampling fails.
function point_cache!(d::MaskDomain)
    if isempty(d.cache) || d.cache_dist != d.min_dist
        empty!(d.cache)
        for k in d.lower[3]:d.upper[3], j in d.lower[2]:d.upper[2], i in d.lower[1]:d.upper[1]
            p = (i, j, k)
            dist = sample_distance(d, p)
            img_at(d, p) >= MASK_THRESHOLD && dist >= d.min_dist && push!(d.cache, p)
        end
        d.cache_dist = d.min_dist
    end
    return d
end

function random_point!(d::MaskDomain{T}, rng::CRand) where {T}
    dp = d.upper .- d.lower
    n = 0
    while n < d.nbtry
        p = (d.lower[1] + crand!(rng) % dp[1], d.lower[2] + crand!(rng) % dp[2], d.lower[3] + crand!(rng) % dp[3])
        dist = sample_distance(d, p)
        img_at(d, p) >= MASK_THRESHOLD && dist >= d.min_dist && return p
        n += 1
    end
    point_cache!(d)
    isempty(d.cache) && return (Int32(0), Int32(0), Int32(0))
    idx = reinterpret(UInt32, crand!(rng)) % UInt32(length(d.cache))
    return d.cache[idx + 1]
end

# ImageMaskDomainCtrl::checkNoIntersectDomain: sample the segment every unit of length.
function check_no_intersect(d, p1::NTuple{3, Int32}, p2::NTuple{3, Int32})
    T = typeof(d.min_dist)
    (in_bounds(d, p1) && in_bounds(d, p2)) || return false
    st = (T(p1[1]), T(p1[2]), T(p1[3]))
    en = (T(p2[1]), T(p2[2]), T(p2[3]))
    dir = sub3(en, st)
    len = norm3(dir)
    len == zero(T) && return interpolate_distance(d, st) >= d.min_dist
    dir = div3(dir, len)
    s = zero(T)
    while s <= len
        interpolate_distance(d, add3(st, mul3(dir, s))) < d.min_dist && return false
        s += one(T)
    end
    return true
end

# ImageMaskDomainCtrl::maxDistantPointFromBorder: first voxel (domain order) with the maximum.
function max_distant_point(d::MaskDomain{T}) where {T}
    m = zero(T); best = (Int32(0), Int32(0), Int32(0))
    for k in d.lower[3]:d.upper[3], j in d.lower[2]:d.upper[2], i in d.lower[1]:d.upper[1]
        v = dist_at(d, (i, j, k))
        v > m && (m = v; best = (i, j, k))
    end
    return best
end
