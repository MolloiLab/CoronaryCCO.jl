# SPDX-License-Identifier: GPL-3.0-or-later
# Flat, device-capable snapshot of the frozen tree used by one batched round.
#
# A batched round evaluates every location against the same tree, so the tree can be packed
# once per round into plain arrays and a compressed (CSR) spatial grid. The candidate code
# (nearest segments, isAddable, Kamiya solve, tentative volume) reads the snapshot through the
# same field names as the Engine, so one implementation serves the sequential engine on the
# CPU and batched rounds on multithreaded CPUs or GPUs (AcceleratedKernels). Symbols and
# Float64 parameters are replaced by flags and `T` scalars so the snapshot is isbits once its
# arrays live on a device.

using Adapt: Adapt

# Model arrays and scalars the candidate code reads (subset of CCOModel).
struct ModelView{T, VP, VI, VT}
    pos::VP
    parent::VI
    left::VI
    right::VI
    flow::VT
    glen::VT
    beta::VT
    radius::VT
    qterm::T
    lf::T
    pdrop::T
    root_radius_max::T
    gamma::T
    kterm_total::Int
end
Adapt.@adapt_structure ModelView

struct FastView{VT}
    Rt::VT
    W::VT
end
Adapt.@adapt_structure FastView

# CSR grid: segments of cell c are items[start[c]:start[c + 1] - 1].
struct GridView{T, VI, VT, VB}
    origin::P3{T}
    h::T
    dims::NTuple{3, Int}
    margin::T
    cellr::VT
    growth::T
    prune::Bool
    start::VI
    items::VI
    box::VB               # per segment: registered cell range (for de-duplication)
end
Adapt.@adapt_structure GridView

struct DomainView{T, A8, AT, A16, AF}
    img::A8
    dist::AT
    allow::A16
    lower::NTuple{3, Int32}
    upper::NTuple{3, Int32}
    min_dist::T
    fine::AF
    fine_origin::NTuple{3, T}
    fine_h::T
end
Adapt.@adapt_structure DomainView

# Run parameters the candidate code reads, in T.
struct ParamView{T}
    num_neighbors::Int
    upstream_attach::Bool
    max_connection::T          # size-scaled (CCOParams.size_scale)
    reach::T                   # TERRITORY_REACH, size-scaled
    kamiya_iters::Int
    kamiya_tol::T
    side_branch_bias::T
    term_flow_cv::T
    sym_n::Int                 # symmetry floor knots (see CCOParams.symmetry_floor)
    sym_ld::NTuple{12, T}
    sym_s::NTuple{12, T}
end
symmetry_min(p::ParamView{T}, d::T) where {T} = _symmetry_min(p.sym_ld, p.sym_s, p.sym_n, d)

struct TreeView{T, M, F, G, D, VI8, VB, VU8, VR, VT}
    m::M
    F::F
    g::G
    dom::D
    p::ParamView{T}
    tcount::VI8
    middle::VB
    isbase::VB
    treeid::VU8
    cscale::T
    lscale::T
    max_angle_cos::T
    rmax::T
    fix_base::Bool
    full::Bool
    level::T              # connectivity restriction level of this round
    roots::VR             # root segment of each tree
    rcaps::VT             # root radius cap of each tree
    rest::VT              # per tree: unscaled volume of the other trees
    rmin::VT              # per tree: tightest cap ratio of the other trees
end
Adapt.@adapt_structure TreeView

@inline cell_items(g::GridView, c) = @inbounds view(g.items, g.start[c]:(g.start[c + 1] - 1))

@inline _isfast(::TreeView) = true
@inline _fullvol(V::TreeView) = V.full
@inline _level(V::TreeView) = V.level
@inline _maxconn(::Type{T}, p::ParamView) where {T} = T(p.max_connection)
@inline _reach(::Type{T}, p::ParamView) where {T} = T(p.reach)
@inline _ktol(::Type{T}, p::ParamView) where {T} = T(p.kamiya_tol)

# Host-side buffers reused across rounds (CSR grid and, for devices, the static domain).
mutable struct ViewCache
    start::Vector{Int32}
    items::Vector{Int32}
    dom::Any               # (img, dist) already moved by `to`, keyed by the `to` function
    to::Any
    bufs::Dict{Symbol, Any}   # device buffers reused across rounds (grown geometrically)
end
ViewCache() = ViewCache(Int32[], Int32[], nothing, nothing, Dict{Symbol, Any}())

# x on the device: aliases x on the CPU; otherwise refills a reusable device buffer that is at
# least as long as x (kernels index explicitly, so trailing capacity is never read).
function _up!(vc::ViewCache, key::Symbol, x::AbstractVector, to)
    to === identity && return x
    buf = get(vc.bufs, key, nothing)
    if buf === nothing || length(buf) < length(x)
        buf = to(similar(x, max(length(x), 1) + length(x) ÷ 2))
        vc.bufs[key] = buf
    end
    copyto!(buf, 1, x, 1, length(x))
    return buf
end

function _csr!(vc::ViewCache, g::SegGrid)
    nc = length(g.cells)
    resize!(vc.start, nc + 1)
    vc.start[1] = 1
    @inbounds for c in 1:nc
        vc.start[c + 1] = vc.start[c] + Int32(length(g.cells[c]))
    end
    resize!(vc.items, vc.start[end] - 1)
    AK.foreachindex(1:nc; max_tasks = Threads.nthreads(), min_elems = 4096) do c
        o = vc.start[c] - 1
        @inbounds for (i, j) in enumerate(g.cells[c])
            vc.items[o + i] = j
        end
    end
    return vc
end

"""
    tree_view(E, vc; to = identity) -> TreeView

Snapshot of the frozen tree for one batched round. With `to = identity` the views alias the
engine's arrays (no copy); with a device array constructor (for example `MtlArray`) they are
uploaded. The domain is uploaded once per `ViewCache`.
"""
function tree_view(E::Engine{T}, vc::ViewCache; to = identity) where {T}
    m = E.m; g = E.g; p = E.p
    up(k, x) = _up!(vc, k, x, to)
    M = ModelView(
        up(:pos, m.pos), up(:parent, m.parent), up(:left, m.left), up(:right, m.right), up(:flow, m.flow),
        up(:glen, m.glen), up(:beta, m.beta), up(:radius, m.radius), m.qterm, m.lf, m.pdrop,
        m.root_radius_max, m.gamma, m.kterm_total
    )
    F = FastView(up(:Rt, E.F.Rt), up(:W, E.F.W))
    _csr!(vc, g)
    G = GridView(
        g.origin, g.h, g.dims, g.margin, up(:cellr, g.cellr), g.growth, g.prune, up(:start, vc.start),
        up(:items, vc.items), up(:box, g.box)
    )
    if vc.dom === nothing || vc.to !== to
        vc.dom = (to(E.dom.img), to(E.dom.dist), to(E.dom.allow), to(E.dom.fine)); vc.to = to
    end
    D = DomainView(vc.dom[1], vc.dom[2], vc.dom[3], E.dom.lower, E.dom.upper, E.dom.min_dist, vc.dom[4], E.dom.fine_origin, E.dom.fine_h)
    P = ParamView{T}(
        p.num_neighbors, p.upstream_attach, T(p.max_connection * p.size_scale), T(TERRITORY_REACH * p.size_scale), p.kamiya_iters,
        T(p.kamiya_tol * p.size_scale^3), T(p.side_branch_bias),
        T(p.term_flow_cv), _symmetry_knots(T, p.symmetry_floor)...
    )
    return TreeView(
        M, F, G, D, P, up(:tcount, E.tcount), up(:middle, E.middle), up(:isbase, E.isbase), up(:treeid, E.treeid), E.cscale,
        E.lscale, E.max_angle_cos, E.rmax,
        E.fix_base, E.volume === :full, T(_level(E)), up(:roots, E.roots), up(:rcaps, E.rcaps),
        up(:rest, E.rest), up(:rmin, E.rmin)
    )
end
