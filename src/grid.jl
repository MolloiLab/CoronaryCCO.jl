# SPDX-License-Identifier: GPL-3.0-or-later
# Uniform grid over segment bounding boxes. It only narrows which segments are tested; every
# answer is decided by upstream's own predicate (is_intersecting, proj_distance, norms), and
# query boxes are padded by the largest radius plus a rounding margin, so each query returns
# exactly what upstream's scan over all segments returns.

mutable struct SegGrid{T}
    origin::P3{T}
    h::T
    dims::NTuple{3, Int}
    cells::Vector{Vector{Int32}}
    box::Vector{NTuple{6, Int32}}   # per segment: cell index range it is registered in
    margin::T
end

function SegGrid{T}(lower::P3{T}, upper::P3{T}, h::Real) where {T}
    pad = T(4h)
    origin = (lower[1] - pad, lower[2] - pad, lower[3] - pad)
    dims = ntuple(i -> Int(cld(upper[i] - lower[i] + 2pad, T(h))) + 1, 3)
    cells = [Int32[] for _ in 1:prod(dims)]
    return SegGrid{T}(origin, T(h), dims, cells, NTuple{6, Int32}[], T(1.0e-3))
end

@inline _cell(g::SegGrid, x, i) = clamp(Int(fld(x - g.origin[i], g.h)) + 1, 1, g.dims[i])
@inline _lin(g::SegGrid, i, j, k) = i + (j - 1) * g.dims[1] + (k - 1) * g.dims[1] * g.dims[2]

function _cellbox(g::SegGrid, a::P3, b::P3, pad)
    return (
        Int32(_cell(g, min(a[1], b[1]) - pad, 1)), Int32(_cell(g, max(a[1], b[1]) + pad, 1)),
        Int32(_cell(g, min(a[2], b[2]) - pad, 2)), Int32(_cell(g, max(a[2], b[2]) + pad, 2)),
        Int32(_cell(g, min(a[3], b[3]) - pad, 3)), Int32(_cell(g, max(a[3], b[3]) + pad, 3)),
    )
end

function grid_insert!(g::SegGrid, m::CCOModel, j)
    length(g.box) < j && resize!(g.box, nsegments(m))
    bx = _cellbox(g, m.pos[j], m.pos[m.parent[j]], zero(g.h))
    g.box[j] = bx
    for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        push!(g.cells[_lin(g, i, jj, k)], Int32(j))
    end
    return
end

function grid_remove!(g::SegGrid, j)
    bx = g.box[j]
    for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        c = g.cells[_lin(g, i, jj, k)]
        idx = findfirst(==(Int32(j)), c)
        c[idx] = c[end]; pop!(c)
    end
    return
end

function grid_build!(g::SegGrid, m::CCOModel)
    foreach(empty!, g.cells)
    resize!(g.box, nsegments(m))
    for j in 2:nsegments(m)
        grid_insert!(g, m, j)
    end
    return g
end

# Call f(j) for every segment registered in a cell overlapping box [a, b] padded by `pad`
# (a segment can be visited more than once). Stops early, returning true, if f returns true.
@inline function grid_any(f::F, g::SegGrid, a::P3, b::P3, pad) where {F}
    bx = _cellbox(g, a, b, pad + g.margin)
    @inbounds for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        for j in g.cells[_lin(g, i, jj, k)]
            f(j) && return true
        end
    end
    return false
end
