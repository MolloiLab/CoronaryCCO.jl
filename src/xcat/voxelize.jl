# SPDX-License-Identifier: GPL-3.0-or-later
# Inside/outside voxelization of closed triangle meshes by parity ray casting.
#
# For each (x, y) voxel column a ray along +z crosses the surface at sorted heights; voxel
# centers between crossing 2k−1 and 2k are inside. Triangles are binned by the columns their
# xy bounding box covers, and columns are processed in parallel with AcceleratedKernels.
# Rays through a vertex or edge are nudged by a fixed sub-voxel offset, so no column is
# counted twice.

"""
    VoxelGrid(origin, h, dims)

Voxel `(i, j, k)` has its center at `origin .+ h .* (i - 1, j - 1, k - 1)` (mm).
"""
struct VoxelGrid
    origin::NTuple{3, Float64}
    h::Float64
    dims::NTuple{3, Int}
end

center(g::VoxelGrid, i, j, k) = (g.origin[1] + g.h * (i - 1), g.origin[2] + g.h * (j - 1), g.origin[3] + g.h * (k - 1))

"""
    grid_around(points, h; margin = 0) -> VoxelGrid

Smallest grid with spacing `h` covering the 3×K `points` plus `margin` mm.
"""
function grid_around(X::AbstractMatrix, h::Real; margin::Real = 0)
    lo = ntuple(d -> minimum(@view X[d, :]) - margin, 3)
    hi = ntuple(d -> maximum(@view X[d, :]) + margin, 3)
    dims = ntuple(d -> max(1, ceil(Int, (hi[d] - lo[d]) / h) + 1), 3)
    return VoxelGrid(lo, Float64(h), dims)
end

"""
    voxelize(X, T, g::VoxelGrid) -> Array{Bool, 3}

Voxels of `g` whose centers lie inside the closed triangle mesh with vertices `X` (3×K, mm)
and triangles `T` (3×M, 1-based).
"""
function voxelize(X::AbstractMatrix, T::AbstractMatrix{<:Integer}, g::VoxelGrid)
    nx, ny, nz = g.dims
    # bin triangles by column
    ncol = nx * ny
    cnt = zeros(Int32, ncol + 1)
    rng = Vector{NTuple{4, Int}}(undef, size(T, 2))
    for t in 1:size(T, 2)
        a, b, c = T[1, t], T[2, t], T[3, t]
        x0 = min(X[1, a], X[1, b], X[1, c]); x1 = max(X[1, a], X[1, b], X[1, c])
        y0 = min(X[2, a], X[2, b], X[2, c]); y1 = max(X[2, a], X[2, b], X[2, c])
        i0 = max(1, ceil(Int, (x0 - g.origin[1]) / g.h) + 1); i1 = min(nx, floor(Int, (x1 - g.origin[1]) / g.h) + 1)
        j0 = max(1, ceil(Int, (y0 - g.origin[2]) / g.h) + 1); j1 = min(ny, floor(Int, (y1 - g.origin[2]) / g.h) + 1)
        rng[t] = (i0, i1, j0, j1)
        for j in j0:j1, i in i0:i1
            cnt[i + nx * (j - 1) + 1] += 1
        end
    end
    cnt[1] = 1
    for c in 1:ncol
        cnt[c + 1] += cnt[c]
    end
    fill_ = copy(cnt); items = Vector{Int32}(undef, cnt[end] - 1)
    for t in 1:size(T, 2)
        i0, i1, j0, j1 = rng[t]
        for j in j0:j1, i in i0:i1
            c = i + nx * (j - 1)
            items[fill_[c]] = t; fill_[c] += 1
        end
    end
    inside = zeros(Bool, g.dims)   # Bool, not BitArray: columns share BitArray chunks across threads
    cols = collect(1:ncol)
    AK.foreachindex(cols; min_elems = 256) do c
        lo = cnt[c]; hi = cnt[c + 1] - 1
        lo > hi && return
        i = mod1(c, nx); j = cld(c, nx)
        # fixed irrational nudge avoids rays through vertices and edges
        px = g.origin[1] + g.h * (i - 1) + 1.0e-7 * g.h * 0.7548776662466927
        py = g.origin[2] + g.h * (j - 1) + 1.0e-7 * g.h * 0.5698402909980532
        zs = Float64[]
        for s in lo:hi
            t = items[s]
            a, b, cc = T[1, t], T[2, t], T[3, t]
            z = _ray_z(X, a, b, cc, px, py)
            z === nothing || push!(zs, z)
        end
        sort!(zs)
        for q in 1:2:(length(zs) - 1)
            k0 = max(1, ceil(Int, (zs[q] - g.origin[3]) / g.h) + 1)
            k1 = min(nz, floor(Int, (zs[q + 1] - g.origin[3]) / g.h) + 1)
            for k in k0:k1
                inside[i, j, k] = true
            end
        end
    end
    return inside
end

# Height where the vertical line through (px, py) meets triangle abc, or nothing.
@inline function _ray_z(X, a, b, c, px, py)
    x1, y1 = X[1, a] - px, X[2, a] - py
    x2, y2 = X[1, b] - px, X[2, b] - py
    x3, y3 = X[1, c] - px, X[2, c] - py
    d1 = x1 * y2 - x2 * y1; d2 = x2 * y3 - x3 * y2; d3 = x3 * y1 - x1 * y3
    ((d1 >= 0 && d2 >= 0 && d3 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0)) || return nothing
    s = d1 + d2 + d3
    s == 0 && return nothing
    # barycentric weights of a, b, c are d2/s, d3/s, d1/s
    return (d2 * X[3, a] + d3 * X[3, b] + d1 * X[3, c]) / s
end
