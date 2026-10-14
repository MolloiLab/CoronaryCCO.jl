# SPDX-License-Identifier: GPL-3.0-or-later
# Preorder layout of the tree for the whole-tree sweeps, and the per-step full hydraulic
# update (upstream updateLengthFactor + updateResistanceFromRoot(1) + updateRootRadius).
#
# Position k orders the real segments (upstream 1..N-1) in upstream's recursion order:
# parent, then left subtree, then right subtree. Parents precede children and every subtree
# is a contiguous range, so the top-down and bottom-up passes are plain array sweeps. The
# arithmetic per node is exactly upstream's; only the visiting order changes, and upstream's
# recursions have no cross-node dependencies other than parent/child.

mutable struct Layout{T}
    pre::Vector{Int32}    # position -> node
    ppos::Vector{Int32}   # node -> position (0 for the root special segment)
    sz::Vector{Int32}     # node -> subtree size
    par::Vector{Int32}    # position -> parent position (0 at the top)
    lc::Vector{Int32}     # position -> left child position (0 if none)
    rc::Vector{Int32}     # position -> right child position (0 if none)
    bin::Vector{Bool}     # position -> parent is a true bifurcation (radius uses beta)
    beta::Vector{T}       # position -> beta (base tree)
    len::Vector{T}        # position -> segment length * length factor (base tree)
    flow::Vector{T}       # position -> flow (scratch for the full update)
    res::Vector{T}        # position -> resistance (scratch for the full update)
    cum::Vector{T}        # position -> cumulative beta (scratch for the full update)
end

Layout{T}() where {T} = Layout{T}(
    Int32[], Int32[], Int32[], Int32[], Int32[], Int32[], Bool[],
    T[], T[], T[], T[], T[]
)

# Build `pre` and `sz` from scratch (once, after import).
function init_layout!(L::Layout, m::CCOModel)
    n = nsegments(m)
    empty!(L.pre)
    stack = Int32[2]
    while !isempty(stack)
        i = pop!(stack); push!(L.pre, i)
        m.right[i] != 0 && push!(stack, m.right[i])
        m.left[i] != 0 && push!(stack, m.left[i])
    end
    length(L.pre) == n - 1 || error("tree is not connected below upstream segment 1")
    L.sz = ones(Int32, n)
    for k in length(L.pre):-1:1
        i = L.pre[k]; p = m.parent[i]
        k > 1 && (L.sz[p] += L.sz[i])
    end
    refresh_layout!(L, m)
    return L
end

# Splice the two segments created by a commit at s (new left nl, new terminal nr) into the
# layout in place. Old positions after s's slot shift by one (inside s's old subtree, now
# under nl) or two (after it); nl goes right after s and nr right after nl's subtree.
function splice_layout!(L::Layout{T}, m::CCOModel, s, nl, nr) where {T}
    k = Int32(L.ppos[s]); old = L.sz[s]
    lim = k + old - Int32(1)
    @inline sh(q) = q <= k ? q : (q <= lim ? q + Int32(1) : q + Int32(2))
    olc = L.lc[k]; orc = L.rc[k]
    insert!(L.pre, k + 1, nl)
    insert!(L.pre, k + old + 1, nr)
    push!(L.sz, old, 1)       # sz[nl] = old size of s's subtree, sz[nr] = 1
    a = s
    while a != 1
        L.sz[a] += Int32(2)
        a = m.parent[a]
    end
    @inbounds for i in eachindex(L.ppos)
        q = L.ppos[i]; q != 0 && (L.ppos[i] = sh(q))
    end
    push!(L.ppos, k + Int32(1), k + old + Int32(1))
    @inbounds for v in (L.par, L.lc, L.rc)
        for j in eachindex(v)
            q = v[j]; q != 0 && (v[j] = sh(q))
        end
    end
    nlc = olc == 0 ? Int32(0) : sh(olc); nrc = orc == 0 ? Int32(0) : sh(orc)
    insert!(L.par, k + 1, k); insert!(L.lc, k + 1, nlc); insert!(L.rc, k + 1, nrc)
    insert!(L.par, k + old + 1, k); insert!(L.lc, k + old + 1, 0); insert!(L.rc, k + old + 1, 0)
    L.lc[k] = k + Int32(1); L.rc[k] = k + old + Int32(1)
    nlc != 0 && (L.par[nlc] = k + Int32(1))
    nrc != 0 && (L.par[nrc] = k + Int32(1))
    insert!(L.bin, k + 1, true); insert!(L.bin, k + old + 1, true)
    for v in (L.beta, L.len, L.flow, L.res, L.cum)
        insert!(v, k + 1, zero(T)); insert!(v, k + old + 1, zero(T))
    end
    return L
end

# Recompute every position-indexed array from `pre` and the node arrays.
function refresh_layout!(L::Layout{T}, m::CCOModel{T}) where {T}
    n = length(L.pre)
    resize!(L.ppos, nsegments(m)); fill!(L.ppos, 0)
    for k in 1:n
        L.ppos[L.pre[k]] = k
    end
    for v in (L.par, L.lc, L.rc)
        resize!(v, n)
    end
    for v in (L.beta, L.len, L.flow, L.res, L.cum)
        resize!(v, n)
    end
    resize!(L.bin, n)
    @inbounds for k in 1:n
        i = L.pre[k]
        p = m.parent[i]
        L.par[k] = k == 1 ? 0 : L.ppos[p]
        L.lc[k] = m.left[i] == 0 ? 0 : L.ppos[m.left[i]]
        L.rc[k] = m.right[i] == 0 ? 0 : L.ppos[m.right[i]]
        L.bin[k] = k == 1 ? false : (m.left[p] != 0 && m.right[p] != 0)
        L.beta[k] = m.beta[i]
        L.len[k] = m.glen[i] * m.lf
    end
    return L
end

# Bottom-up resistance over the position range lo:hi (a set of whole subtrees, or the
# upper part when the subtrees below have been done).
@inline function _resistance_range!(L::Layout{T}, gamma::T, lo, hi) where {T}
    len = L.len; lc = L.lc; rc = L.rc; beta = L.beta; res = L.res; flow = L.flow
    @inbounds for k in hi:-1:lo
        R = own_resistance(len[k])
        a = lc[k]; b = rc[k]
        if a != 0 && b == 0
            beta[a] = one(T); R += res[a]
        elseif b != 0 && a == 0
            beta[b] = one(T); R += res[b]
        elseif a != 0 && b != 0
            bL, bR, inv = resistance_combine(gamma, flow[a], flow[b], res[a], res[b])
            beta[a] = bL; beta[b] = bR
            R += inv
        end
        res[k] = R
    end
    return
end

"""
    full_update!(m, L)

Upstream's end-of-step update: length factor, every resistance and beta from the leaves up,
then the root radius and every radius from the root down. Subtrees are processed in
parallel (they are independent), then the few nodes above them serially.
"""
function full_update!(m::CCOModel{T}, L::Layout{T}) where {T}
    update_length_factor!(m)
    n = length(L.pre)
    pre = L.pre; lf = m.lf
    _gather!(L.len, m.glen, pre, lf)
    _gather!(L.flow, m.flow, pre)
    _gather!(L.beta, m.beta, pre)
    # Chunk roots: maximal subtrees no larger than `cap`.
    cap = max(Int32(256), Int32(cld(n, 8 * Threads.nthreads())))
    roots = Int32[]
    k = 1
    while k <= n
        if L.sz[L.pre[k]] <= cap
            push!(roots, k); k += L.sz[L.pre[k]]
        else
            k += 1
        end
    end
    AK.foreachindex(roots; max_tasks = Threads.nthreads(), min_elems = 1) do c
        r = roots[c]
        _resistance_range!(L, m.gamma, r, r + L.sz[L.pre[r]] - 1)
    end
    # Upper nodes (subtree larger than cap), deepest first: reverse position order.
    @inbounds for k in n:-1:1
        L.sz[L.pre[k]] > cap && _resistance_range!(L, m.gamma, k, k)
    end
    # Root radius and radii top-down.
    rroot = upow((L.res[1] * L.flow[1]) / m.pdrop, T(0.25))
    (m.root_radius_max > 0 && rroot > m.root_radius_max) && (rroot = m.root_radius_max)
    _radius_sweep!(L.cum, L.par, L.bin, L.beta, n)
    _scatter!(m.resistance, L.res, pre)
    _scatter!(m.beta, L.beta, pre)
    _scatter_radius!(m.radius, L.cum, pre, rroot)
    return m
end

function _gather!(dst, src, pre, s = nothing)
    return @inbounds for k in eachindex(pre)
        dst[k] = s === nothing ? src[pre[k]] : src[pre[k]] * s
    end
end
function _scatter!(dst, src, pre)
    return @inbounds for k in eachindex(pre)
        dst[pre[k]] = src[k]
    end
end
function _scatter_radius!(radius, cum, pre, rroot)
    @inbounds radius[pre[1]] = rroot
    return @inbounds for k in 2:length(pre)
        radius[pre[k]] = rroot * cum[k]
    end
end
function _radius_sweep!(cum::AbstractVector{T}, par, bin, beta, n) where {T}
    @inbounds cum[1] = one(T)
    return @inbounds for k in 2:n
        c = cum[par[k]]
        cum[k] = bin[k] ? c * beta[k] : c
    end
end
