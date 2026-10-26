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

# The grid's cell lists as on the host (see SegGrid): head per cell, next and item per entry.
struct GridView{T, VI, VT, VB}
    origin::P3{T}
    h::T
    dims::NTuple{3, Int}
    margin::T
    cellr::VT
    growth::T
    prune::Bool
    head::VI
    next::VI
    item::VI
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
    big::VR                   # segments tested directly by every grid query (Engine.big; buffer may be longer)
    nbig::Int
    fix_base::Bool
    full::Bool
    keepout::T            # vessel keep-out factor for new terminals (see distance_valid)
    level::T              # connectivity restriction level of this round
    roots::VR             # root segment of each tree
    rcaps::VT             # root radius cap of each tree
    rest::VT              # per tree: unscaled volume of the other trees
    rmin::VT              # per tree: tightest cap ratio of the other trees
    open::UInt16          # trees at or below their share of the flow now (Engine.open)
end
Adapt.@adapt_structure TreeView

@inline cell_items(g::GridView, c) = @inbounds CellItems(g.next, g.item, g.head[c])

@inline _isfast(::TreeView) = true
@inline _fullvol(V::TreeView) = V.full
@inline _level(V::TreeView) = V.level
@inline function _any_big(f::F, V::TreeView) where {F}
    @inbounds for i in 1:V.nbig
        f(V.big[i]) && return true
    end
    return false
end
@inline _keepout(::Type{T}, V::TreeView) where {T} = T(V.keepout)
@inline _maxconn(::Type{T}, p::ParamView) where {T} = T(p.max_connection)
@inline _reach(::Type{T}, p::ParamView) where {T} = T(p.reach)
@inline _ktol(::Type{T}, p::ParamView) where {T} = T(p.kamiya_tol)

"""
    Unified(wrap, copy)

A device that shares memory with the host (Apple silicon, integrated GPUs). `wrap(x)` returns a device array over the same
memory as the host vector `x` without copying; `copy(x)` moves any other array there. Each batched round then hands the
growing tree to the device in place instead of copying it, which halves the memory of a large tree. With `GPUArray` the
GPU package's array type and its no-copy wrapper of host memory:

    Unified(x -> unsafe_wrap(GPUArray{eltype(x), ndims(x)}, x), x -> GPUArray(x))

`wrap` and `copy` must return the same array type (the round's views mix wrapped and copied arrays).
"""
struct Unified{W, C}
    wrap::W
    copy::C
end
(u::Unified)(x) = u.copy(x)

# Host-side buffers reused across rounds (for devices: the static domain and the tree's device arrays).
mutable struct ViewCache
    dom::Any               # (img, dist) already moved by `to`, keyed by the `to` function
    to::Any
    bufs::Dict{Symbol, Any}   # device buffers reused across rounds (grown geometrically)
end
ViewCache() = ViewCache(nothing, nothing, Dict{Symbol, Any}())

# x on the device: aliases x on the CPU; otherwise refills a reusable device buffer that is at
# least as long as x (kernels index explicitly, so trailing capacity is never read).
function _up!(vc::ViewCache, key::Symbol, x::AbstractVector, to)
    to === identity && return x
    to isa Unified && return _wrap!(vc, key, x, to)
    return _copyup!(vc, key, x, to)
end
function _copyup!(vc::ViewCache, key::Symbol, x::AbstractVector, to)
    buf = get(vc.bufs, key, nothing)
    if buf === nothing || length(buf) < length(x)
        buf = to(similar(x, max(length(x), 1) + length(x) ÷ 2))
        vc.bufs[key] = buf
    end
    copyto!(buf, 1, x, 1, length(x))
    return buf
end

# A unified-memory device takes x in place: a device array over the first m elements of x's memory, the fewest that fill
# whole pages (a no-copy device buffer must start and end on a page boundary; large host arrays start on one). Elements
# past length(x) are never read. Re-wrapped only when x was reallocated or outgrew the wrapped length.
const _PAGE = Mmap.PAGESIZE
function _pagelen(T, n)
    m = max(n, 1)
    while (m * sizeof(T)) % _PAGE != 0
        m += 1
    end
    return m
end
function _wrap!(vc::ViewCache, key::Symbol, x::Vector, to::Unified)
    c = get(vc.bufs, key, nothing)
    UInt(pointer(x)) % _PAGE == 0 || return _copyup!(vc, Symbol(key, :_copy), x, to.copy) # not page-aligned: copy
    if c === nothing || c[1] !== x.ref.mem || c[3] < length(x)
        # wrap all whole pages of the capacity (not just the length), so a growing tree is re-wrapped only when its host
        # array is reallocated (geometric growth); reserve half again first when even that would not hold x
        m = _pagefloor(eltype(x), length(x.ref.mem))
        if m < max(length(x), 1) # also: at least one page, even for an empty array
            sizehint!(x, _pagelen(eltype(x), length(x) + length(x) ÷ 2))
            m = _pagefloor(eltype(x), length(x.ref.mem))
        end
        UInt(pointer(x)) % _PAGE == 0 || return _copyup!(vc, Symbol(key, :_copy), x, to.copy) # moved by sizehint!: check again
        d = to.wrap(unsafe_wrap(Array, pointer(x), m))
        c = (x.ref.mem, d, m); vc.bufs[key] = c
    end
    return c[2]
end
# the most elements of type T, at most n, that fill whole pages
function _pagefloor(T, n)
    m = n
    while m > 0 && (m * sizeof(T)) % _PAGE != 0
        m -= 1
    end
    return m
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
    G = GridView(
        g.origin, g.h, g.dims, g.margin, up(:cellr, g.cellr), g.growth, g.prune, up(:head, g.head),
        up(:next, g.next), up(:item, g.item), up(:box, g.box)
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
        E.lscale, E.max_angle_cos, E.rmax, up(:big, E.big), length(E.big),
        E.fix_base, E.volume === :full, _keepout(T, E), T(_level(E)), up(:roots, E.roots), up(:rcaps, E.rcaps),
        up(:rest, E.rest), up(:rmin, E.rmin), E.open
    )
end
