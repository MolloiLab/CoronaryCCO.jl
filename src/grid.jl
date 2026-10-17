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
    box::Vector{NTuple{6, Int16}}   # per segment: cell index range it is registered in (Int16: 12 bytes per segment)
    margin::T
    cellr::Vector{T}                 # per cell: largest reference radius among its segments
    rref::Vector{T}                  # per segment: radius when its cells' bounds were last set
    growth::T                        # max over segments of current radius / rref (≥ 1)
    prune::Bool                      # skip cells whose radius bound cannot reach a query
end

function SegGrid{T}(lower::P3{T}, upper::P3{T}, h::Real) where {T}
    pad = T(4h)
    origin = (lower[1] - pad, lower[2] - pad, lower[3] - pad)
    dims = ntuple(i -> Int(cld(upper[i] - lower[i] + 2pad, T(h))) + 1, 3)
    maximum(dims) <= typemax(Int16) || throw(ArgumentError("grid of $(dims) cells: at most $(typemax(Int16)) per axis"))
    cells = [Int32[] for _ in 1:prod(dims)]
    return SegGrid{T}(origin, T(h), dims, cells, NTuple{6, Int16}[], T(1.0e-3), zeros(T, prod(dims)), T[], one(T), false)
end

@inline _cell(g, x, i) = clamp(unsafe_trunc(Int, floor((x - g.origin[i]) / g.h)) + 1, 1, g.dims[i])
@inline _lin(g, i, j, k) = i + (j - 1) * g.dims[1] + (k - 1) * g.dims[1] * g.dims[2]

@inline function _cellbox(g, a::P3, b::P3, pad)
    return (
        Int32(_cell(g, min(a[1], b[1]) - pad, 1)), Int32(_cell(g, max(a[1], b[1]) + pad, 1)),
        Int32(_cell(g, min(a[2], b[2]) - pad, 2)), Int32(_cell(g, max(a[2], b[2]) + pad, 2)),
        Int32(_cell(g, min(a[3], b[3]) - pad, 3)), Int32(_cell(g, max(a[3], b[3]) + pad, 3)),
    )
end

function grid_insert!(g::SegGrid, m::CCOModel, j)
    length(g.box) < j && resize!(g.box, nsegments(m))
    bx = _cellbox(g, m.pos[j], m.pos[m.parent[j]], zero(g.h))
    g.box[j] = Int16.(bx)
    r = m.radius[j]
    length(g.rref) < j && resize!(g.rref, nsegments(m))
    g.rref[j] = r
    for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        c = _lin(g, i, jj, k)
        push!(g.cells[c], Int32(j))
        r > g.cellr[c] && (g.cellr[c] = r)
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

# After a radius refresh: every stored radius is at most growth x its reference radius, so
# cell bounds scaled by `growth` stay valid. Recompute the bounds (and reset the references)
# only when `growth` exceeds `tol`, since that costs a pass over every cell.
function grid_refresh_radii!(g::SegGrid{T}, m::CCOModel{T}; tol = 1.25) where {T}
    rad = m.radius; rref = g.rref
    resize!(rref, nsegments(m))
    n = nsegments(m)
    nch = max(1, min(4 * Threads.nthreads(), n ÷ 32768))
    part = ones(T, nch)
    AK.foreachindex(1:nch; max_tasks = nch, min_elems = 1) do c
        q = one(T)
        len = cld(n - 1, nch)
        @inbounds for j in (2 + (c - 1) * len):min(n, 1 + c * len)
            rref[j] > 0 && (x = rad[j] / rref[j]; x > q && (q = x))
        end
        part[c] = q
    end
    gr = maximum(part)
    g.growth = gr
    gr <= T(tol) && return g
    cr = g.cellr; cells = g.cells
    AK.foreachindex(cells; max_tasks = Threads.nthreads(), min_elems = 4096) do c
        r = zero(T)
        @inbounds for j in cells[c]
            rad[j] > r && (r = rad[j])
        end
        @inbounds cr[c] = r
    end
    copyto!(rref, rad)
    g.growth = one(T)
    return g
end

# Segments registered in cell c (a SegGrid cell list or a CSR slice of a GridView).
@inline cell_items(g::SegGrid, c) = @inbounds g.cells[c]

# Like grid_any, for tests whose reach depends on each segment's radius: `need(r)` is the
# padding a segment of radius r requires. Cells are searched within need(rmax); with pruning
# on, a cell is skipped when its distance to box [a, b] exceeds need(its radius bound).
@inline function grid_any_r(f::F, g, a::P3{T}, b::P3, need::N, rmax) where {F, N, T}
    pad = need(rmax)
    bx = _cellbox(g, a, b, pad + g.margin)
    lo = (min(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))
    hi = (max(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3]))
    h = g.h; o = g.origin
    @inbounds for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        c = _lin(g, i, jj, k)
        if g.prune
            cx = o[1] + (i - 1) * h; cy = o[2] + (jj - 1) * h; cz = o[3] + (k - 1) * h
            gx = max(zero(T), lo[1] - (cx + h), cx - hi[1])
            gy = max(zero(T), lo[2] - (cy + h), cy - hi[2])
            gz = max(zero(T), lo[3] - (cz + h), cz - hi[3])
            gx * gx + gy * gy + gz * gz > (need(g.cellr[c] * g.growth) + g.margin)^2 && continue
        end
        for j in cell_items(g, c)
            f(j) && return true
        end
    end
    return false
end

# Call f(j) for every segment registered in a cell overlapping box [a, b] padded by `pad`
# (a segment can be visited more than once). Stops early, returning true, if f returns true.
@inline function grid_any(f::F, g, a::P3, b::P3, pad) where {F}
    bx = _cellbox(g, a, b, pad + g.margin)
    @inbounds for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
        for j in cell_items(g, _lin(g, i, jj, k))
            f(j) && return true
        end
    end
    return false
end
