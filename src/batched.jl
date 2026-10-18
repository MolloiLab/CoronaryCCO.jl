# SPDX-License-Identifier: GPL-3.0-or-later
# Batched growth: many terminals per round instead of one per step.
#
# A round draws `nb` terminal locations from the same random stream as sequential growth,
# evaluates every location's candidate bifurcations in parallel against the same (frozen)
# tree, then accepts, in draw order, every location whose best bifurcation does not conflict
# with one already accepted in the round. Conflicts are: a host segment adjacent to another
# accepted host (same segment, parent, child or sibling), a new terminal closer than the
# current spacing threshold to another new terminal, and new segments that intersect, or pass
# too close to, the new segments of another insertion (tested at the same scaled radii as
# every other collision test). Accepted insertions are committed structurally, then the
# baseline resistances, subtree volume factors and betas are recomputed once over the union
# of their root paths, and all radii are refreshed once.
#
# Insertions within a round do not see each other's flow when their bifurcation is
# optimized. Once the tree is large, each insertion changes the flow of its path by about
# 1/N, so this is a small perturbation; Cury et al. (2021) grew CCO networks concurrently on
# this premise and compared them against sequential growth. The selection is deterministic
# and independent of the thread count. See docs/deviations.md.

# Set VASCULATURE_ROUND_LOG=1 to print one timing line per batched round to stderr.


# Exact radii of this round's changed segments `D` (new ones and their paths to the root), in parallel; new segments
# (index > n0) get their grid cells' radius bounds set, older ones raise the grid's drift factor if they grew.
function _path_radii!(E::Engine{T}, D, n0) where {T}
    m = E.m; g = E.g
    AK.foreachindex(D; max_tasks = Threads.nthreads(), min_elems = 2048) do i
        j = D[i]; m.radius[j] = fast_radius(E, j)
    end
    gr = g.growth
    for j in D
        r = m.radius[j]; note_radius!(E, j, r)
        if j > n0
            g.rref[j] = r; bx = g.box[j]
            for k in bx[5]:bx[6], jj in bx[3]:bx[4], i in bx[1]:bx[2]
                c = _lin(g, i, jj, k); r > g.cellr[c] && (g.cellr[c] = r)
            end
        elseif g.rref[j] > 0
            gr = max(gr, r / g.rref[j])
        end
    end
    g.growth = gr
    return
end
_round_log() = get(ENV, "VASCULATURE_ROUND_LOG", "") == "1"

mutable struct BatchPool{T}
    pts::Vector{P3{T}}
    seeds::Vector{UInt64}      # per-location random streams
    res::Vector{Tuple{Int32, T, P3{T}, T, Int, Int, Int}}   # (s, score, pOpt, qNew, iters, ncand, nvalid)
    stamp::Vector{Int32}       # per segment: round in which it was claimed by an accepted host
    dirty::Vector{Int32}       # per segment: round in which its hydraulics were invalidated
    cnt::Vector{Int32}         # per segment: dirty children still to process
    round::Int32
    acc::Vector{Int}           # accepted location indices of the round
    cells::Dict{NTuple{3, Int32}, Vector{Int32}}   # accepted-insertion hash (index into acc)
    vc::ViewCache
    to::Any                    # array constructor for the evaluation device (identity: CPU)
end
# `device` is an array constructor (`MtlArray`, `CuArray`, ...; `identity` or `Array` for the
# CPU), a KernelAbstractions backend, or a `Unified` device; all become a function that moves an array there.
_mover(device) = device === identity || device === Array || device isa KA.CPU ? identity :
    device isa KA.Backend ? (x -> copyto!(KA.allocate(device, eltype(x), size(x)), x)) : device

BatchPool{T}(; to = identity) where {T} = BatchPool{T}(
    P3{T}[], UInt64[], Tuple{Int32, T, P3{T}, T, Int, Int, Int}[], Int32[], Int32[],
    Int32[], Int32(0), Int[], Dict{NTuple{3, Int32}, Vector{Int32}}(), ViewCache(), to
)

# Candidate evaluation of a round against the frozen tree V, in two kernels that run wherever
# V's arrays live (CPU threads or a GPU): the nearest segments of every location, then one
# work item per (candidate, location) pair. Short work items keep GPU threads in step.
function nearest_locations!(nc, bd, bj, V, pts)
    nn = V.p.num_neighbors
    AK.foreachindex(nc; max_tasks = Threads.nthreads(), min_elems = 64) do k
        nc[k] = Int32(nearest_k!(bd, bj, k, V, pts[k], nn))
    end
    return
end

# In two passes, so that the costly part (Kamiya iterations and the path walk) runs only on candidates that survive the
# bifurcation collision test (31-47% of a 1 M lung's candidates end before any Kamiya iteration, most of them here): first the test for every pair, then the survivors packed into one
# list, so a GPU's threads stay busy instead of idling beside rejected neighbours. Each candidate's result is the same as
# in one pass (eval_candidate).
function evaluate_candidates!(cok, cscore, copt, cq, cit, bj, nc, V, pts)
    nn = size(bj, 1)
    AK.foreachindex(cok; max_tasks = Threads.nthreads(), min_elems = 16) do i
        T = typeof(V.m.lf)
        c = (i - 1) % nn + 1; k = (i - 1) ÷ nn + 1
        live = false
        c <= nc[k] && (live = bifurcation_clear(V, pts[k], bj[c, k]))
        cok[i] = live; cit[i] = Int32(0)
    end
    hlive = Array(cok); idx = Int32[i for i in eachindex(hlive) if hlive[i]]
    isempty(idx) && return
    didx = similar(cok, Int32, length(idx)); copyto!(didx, idx)
    AK.foreachindex(didx; max_tasks = Threads.nthreads(), min_elems = 16) do j
        T = typeof(V.m.lf)
        i = didx[j]; c = (i - 1) % nn + 1; k = (i - 1) ÷ nn + 1
        p = pts[k]; s = bj[c, k]
        ok, vol, pOpt, qNew, it = try_candidate!(nothing, V, p, s)
        cok[i] = ok; cscore[i] = ok ? candidate_score(V, s, vol) : zero(T)
        copt[i] = pOpt; cq[i] = qNew; cit[i] = Int32(max(it, 0))
    end
    return
end

# The three segments an insertion at s creates: parent end → bifurcation, bifurcation → old
# distal end of s, bifurcation → new terminal.
@inline function _new_segments(m::CCOModel{T}, s, pOpt::P3{T}, p::P3{T}) where {T}
    return ((m.pos[m.parent[s]], pOpt), (pOpt, m.pos[s]), (pOpt, p))
end

@inline _hkey(x::P3{T}, h::T) where {T} = (Int32(fld(x[1], h)), Int32(fld(x[2], h)), Int32(fld(x[3], h)))

# Cells (of size h) touched by the box around segment (a, b) padded by `pad`.
@inline function _hrange(a::P3{T}, b::P3{T}, pad::T, h::T) where {T}
    lo = _hkey((min(a[1], b[1]) - pad, min(a[2], b[2]) - pad, min(a[3], b[3]) - pad), h)
    hi = _hkey((max(a[1], b[1]) + pad, max(a[2], b[2]) + pad, max(a[3], b[3]) + pad), h)
    return lo, hi
end

# Does location k conflict with an insertion already accepted this round?
function _conflicts(E::Engine{T}, pool::BatchPool{T}, k, thr::T, h::T) where {T}
    m = E.m; c = E.cscale
    s, _, pOpt, _, _, _, _ = pool.res[k]
    p = pool.pts[k]
    rnd = pool.round
    st = pool.stamp
    # hosts must not be adjacent to another accepted host
    for j in (s, m.parent[s], m.left[s], m.right[s])
        j > 1 && st[j] == rnd && return true
    end
    rs = fast_radius(E, s) * c
    segs = _new_segments(m, s, pOpt, p)
    seen = 0
    for (a, b) in segs
        lo, hi = _hrange(a, b, max(thr, T(2) * rs) + rs, h)
        for kz in lo[3]:hi[3], ky in lo[2]:hi[2], kx in lo[1]:hi[1]
            lst = get(pool.cells, (kx, ky, kz), nothing)
            lst === nothing && continue
            for ai in lst
                k2 = pool.acc[ai]
                s2, _, pOpt2, _, _, _, _ = pool.res[k2]
                p2 = pool.pts[k2]
                # terminal spacing (generate_a_location's threshold)
                norm3(sub3(p, p2)) <= thr && return true
                r2 = fast_radius(E, s2) * c
                for (a2, b2) in _new_segments(m, s2, pOpt2, p2)
                    is_intersecting(a, b, rs, a2, b2, r2) && return true
                    # a new terminal must keep clear of the other insertion's new vessels
                    proj_distance(a2, b2, p) < max(T(0.05), T(2) * r2) && return true
                end
                for (a1, b1) in segs
                    proj_distance(a1, b1, p2) < max(T(0.05), T(2) * rs) && return true
                end
            end
        end
    end
    return false
end

function _accept!(E::Engine{T}, pool::BatchPool{T}, k, thr::T, h::T) where {T}
    m = E.m; c = E.cscale
    s, _, pOpt, _, _, _, _ = pool.res[k]
    rnd = pool.round
    for j in (s, m.parent[s], m.left[s], m.right[s])
        j > 1 && (pool.stamp[j] = rnd)
    end
    push!(pool.acc, k)
    ai = Int32(length(pool.acc))
    rs = fast_radius(E, s) * c
    for (a, b) in _new_segments(m, s, pOpt, pool.pts[k])
        lo, hi = _hrange(a, b, rs, h)
        for kz in lo[3]:hi[3], ky in lo[2]:hi[2], kx in lo[1]:hi[1]
            lst = get!(() -> Int32[], pool.cells, (kx, ky, kz))
            (isempty(lst) || lst[end] != ai) && push!(lst, ai)
        end
    end
    return
end

# Recompute Rt, W and betas over the union of the root paths of this round's insertions,
# children before parents.
function _refresh_dirty!(E::Engine{T}, pool::BatchPool{T}, starts) where {T}
    m = E.m; rnd = pool.round
    n = nsegments(m)
    if length(pool.dirty) < n
        k = length(pool.dirty)
        resize!(pool.dirty, n + n ÷ 2); resize!(pool.cnt, n + n ÷ 2)
        fill!(@view(pool.dirty[(k + 1):end]), Int32(0))
    end
    D = Int32[]
    for j0 in starts
        j = j0
        while j != 1 && pool.dirty[j] != rnd
            pool.dirty[j] = rnd; pool.cnt[j] = 0
            push!(D, j)
            j = m.parent[j]
        end
    end
    for j in D
        p = m.parent[j]
        p != 1 && (pool.cnt[p] += Int32(1))
    end
    # wavefronts: every node whose dirty children are done, in parallel (a node writes only
    # its own Rt, W and its children's betas)
    front = Int32[j for j in D if pool.cnt[j] == 0]
    nxt = Int32[]
    while !isempty(front)
        AK.foreachindex(front; max_tasks = Threads.nthreads(), min_elems = 2048) do i
            fast_node!(E, front[i])
        end
        empty!(nxt)
        for j in front
            p = m.parent[j]
            if p != 1
                pool.cnt[p] -= Int32(1)
                pool.cnt[p] == 0 && push!(nxt, p)
            end
        end
        front, nxt = nxt, front
    end
    return D
end

"""
    batched_round!(E, pool, nb; parallel = true) -> accepted

One round of batched growth with `nb` locations. Returns the number of terminals added. The
result does not depend on `parallel` or on the thread count.
"""
function batched_round!(E::Engine{T}, pool::BatchPool{T}, nb::Int; parallel::Bool = true) where {T}
    m = E.m; st = E.stats
    update_cscale!(E)
    pool.round += Int32(1)
    n0 = nsegments(m)
    length(pool.stamp) < n0 && append!(pool.stamp, zeros(Int32, n0 + n0 ÷ 2 - length(pool.stamp)))
    # 1. locations: each draws from its own stream, seeded in order from the main stream, so
    #    they can be generated in parallel and still do not depend on the thread count
    ta = time_ns()
    resize!(pool.pts, nb); resize!(pool.seeds, nb)
    for k in 1:nb
        pool.seeds[k] = (UInt64(reinterpret(UInt32, crand!(E.rng))) << 31) ⊻ UInt64(reinterpret(UInt32, crand!(E.rng)))
    end
    point_cache!(E.dom)
    pts = pool.pts; seeds = pool.seeds
    if parallel
        AK.foreachindex(pts; max_tasks = Threads.nthreads(), min_elems = 16) do k
            pts[k] = generate_new_location!(E, CRand(seeds[k]))
        end
    else
        for k in 1:nb
            pts[k] = generate_new_location!(E, CRand(seeds[k]))
        end
    end
    st.locations += nb
    tb = time_ns()
    # 2. every location's best candidate against the frozen tree (on pool.to's device)
    to = parallel ? pool.to : identity
    tv0 = time_ns()
    V = tree_view(E, pool.vc; to = to)
    st.t_view += (time_ns() - tv0) / 1.0e9
    nn = E.p.num_neighbors
    dpts = to(pts)
    nc = to(zeros(Int32, nb)); bd = to(zeros(T, nn, nb)); bj = to(zeros(Int32, nn, nb))
    cok = to(zeros(Bool, nn * nb)); cscore = to(zeros(T, nn * nb)); copt = to(fill(pts[1], nn * nb))
    cq = to(zeros(T, nn * nb)); cit = to(zeros(Int32, nn * nb))
    if parallel
        nearest_locations!(nc, bd, bj, V, dpts)
        evaluate_candidates!(cok, cscore, copt, cq, cit, bj, nc, V, dpts)
    else
        for k in 1:nb   # the same work items, one at a time
            r = ((k - 1) * nn + 1):(k * nn)
            nearest_locations!(view(nc, k:k), view(bd, :, k:k), view(bj, :, k:k), V, view(dpts, k:k))
            evaluate_candidates!(
                view(cok, r), view(cscore, r), view(copt, r), view(cq, r), view(cit, r), view(bj, :, k:k),
                view(nc, k:k), V, view(dpts, k:k)
            )
        end
    end
    hnc, hj = Array(nc), Array(bj)
    hok, hsc, ho, hq, hi = Array(cok), Array(cscore), Array(copt), Array(cq), Array(cit)
    # best candidate per location: minimum score, first in candidate order on ties
    resize!(pool.res, nb)
    res = pool.res
    for k in 1:nb
        best = Int32(0); bscore = -one(T); bopt = pts[k]; bq = zero(T); its = 0; nv = 0
        for c in 1:hnc[k]
            i = (k - 1) * nn + c
            its += hi[i]
            hok[i] || continue
            nv += 1
            if bscore < 0 || bscore > hsc[i]
                bscore = hsc[i]; best = hj[c, k]; bopt = ho[i]; bq = hq[i]
            end
        end
        res[k] = (best, bscore, bopt, bq, its, Int(hnc[k]), nv)
    end
    tc = time_ns()
    # 3. deterministic conflict-free selection in draw order
    thr = distance_threshold(m)
    h = max(T(4) * thr, T(1))
    empty!(pool.acc); empty!(pool.cells)
    for k in 1:nb
        s, _, _, _, its, nc, nv = res[k]
        st.candidates += nc; st.kamiya_iters += its; st.accepted += nv; st.addable_calls += nc
        s == 0 && continue
        _conflicts(E, pool, k, thr, h) && continue
        _accept!(E, pool, k, thr, h)
    end
    # 4. commit structure, grid and connectivity flags
    starts = Int32[]
    for k in pool.acc
        s, _, pOpt, qNew, _, _, _ = res[k]
        push!(E.isbase, E.isbase[s], false)
        push!(E.treeid, E.treeid[s], E.treeid[s])
        nl, nr = commit_structure!(E, s, pOpt, pts[k], qNew)
        grid_remove!(E.g, s)
        grid_insert!(E.g, m, s); grid_insert!(E.g, m, nl); grid_insert!(E.g, m, nr)
        push!(starts, nl, nr)
    end
    for k in pool.acc
        s = res[k][1]
        middle_update!(E, s, m.left[s], m.right[s])
    end
    # 5. hydraulics once per round
    tr0 = time_ns()
    D = _refresh_dirty!(E, pool, starts); ndirty = length(D)
    tr1 = time_ns()
    update_length_factor!(m)
    # radii: exact on this round's changed paths every round; a full refresh (O(N)) once the tree has grown by
    # RADII_REFRESH since the last one (the same rule as sequential growth, commit_fast!)
    E.F.since += length(pool.acc)
    if E.F.since >= radii_refresh_period(m)
        fast_refresh_radii!(E)
    else
        _path_radii!(E, D, n0)
    end
    st.t_refresh += (time_ns() - tr0) / 1.0e9
    if _round_log()
        @printf(
            stderr, "round %d nb %d acc %d dirty %d | loc %.0f eval %.0f select+commit %.0f dirty %.0f radii %.0f ms\n",
            pool.round, nb, length(pool.acc), ndirty, (tb - ta) / 1.0e6, (tc - tb) / 1.0e6, (tr0 - tc) / 1.0e6,
            (tr1 - tr0) / 1.0e6, (time_ns() - tr1) / 1.0e6
        )
    end
    st.rounds += 1
    na = length(pool.acc)
    st.steps += na; st.seeds_found += na
    td = time_ns()
    st.t_location += (tb - ta) / 1.0e9; st.t_candidates += (tc - tb) / 1.0e9; st.t_commit += (td - tc) / 1.0e9
    return na
end
