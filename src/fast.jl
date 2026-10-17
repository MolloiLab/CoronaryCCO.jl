# SPDX-License-Identifier: GPL-3.0-or-later
# Fast mode: the same growth algorithm with an O(depth) hydraulic update.
#
# Upstream recomputes every resistance, radius and the total volume for every Kamiya
# iteration of every candidate, which makes each evaluation O(N). Two identities remove that:
#
#   * Every segment length is (geometric length) x (global length factor), so every subtree
#     resistance is proportional to the length factor and every radius ratio (beta) is
#     independent of it. We store the factor-free "baseline" resistance Rt (the HCT-CCO
#     paper's R-tilde); the length factor only enters the root radius.
#   * Radii are r_root x (product of betas from the root), so the volume of a subtree is its
#     top radius squared times W_j = l_j + sum over children c of beta_c^2 W_c (geometric
#     lengths). The total volume is r_root^2 x length factor x W_root.
#
# Adding a terminal changes Rt, W and the betas only on the path from the new bifurcation to
# the root, so candidate evaluation and commit are O(depth). Results equal exact mode to
# rounding (summation order differs), not bit for bit. Radii of segments off the updated path
# drift slightly between full refreshes (they only feed collision and spacing tests); the
# radii used in each candidate's own Kamiya solve are exact.
#
# `volume = :upstream` reproduces upstream's computeTotalVolume, which adds a node's children
# only at true bifurcations (so everything below a single-child node of an imported base
# tree is ignored); `volume = :full` counts the whole tree.


@inline own_baseline(lg::T) where {T} = (!(lg > zero(T)) || !isfinite(lg)) ? T(8) * T(MU) * T(KEPS) / T(π) :
    T(8) * T(MU) * lg / T(π)

# Recompute Rt[j], W[j] and the betas of j's children from the children's stored values.
@inline function fast_node!(E::Engine{T}, j) where {T}
    m = E.m; F = E.F; full = E.volume === :full
    lg = m.glen[j]
    R = own_baseline(lg); W = lg
    a = m.left[j]; b = m.right[j]
    if a != 0 && b == 0
        m.beta[a] = one(T); R += F.Rt[a]; full && (W += F.W[a])
    elseif b != 0 && a == 0
        m.beta[b] = one(T); R += F.Rt[b]; full && (W += F.W[b])
    elseif a != 0 && b != 0
        bA, bB, inv = resistance_combine(m.gamma, m.flow[a], m.flow[b], F.Rt[a], F.Rt[b])
        m.beta[a] = bA; m.beta[b] = bB
        R += inv
        W += bA * bA * F.W[a] + bB * bB * F.W[b]
    end
    F.Rt[j] = R; F.W[j] = W
    return
end

# Index of root segment c among E.roots (a forest has a few roots, a tree one).
@inline function _rootidx(E, c)
    roots = E.roots
    @inbounds for k in 1:length(roots)
        roots[k] == c && return k
    end
    return 1
end

# Pressure-law root radius of the tree whose root segment is c (no cap): every tree takes its
# root radius from the same pressure drop and the same terminal flow.
@inline function _root_unc(E, c)
    T = typeof(E.m.lf)
    m = E.m
    return upow((m.lf * E.F.Rt[c] * m.flow[c]) / m.pdrop, T(0.25))
end

# Root radius of tree k from its pressure-law radius r. Caps are upper bounds that act through
# one scale shared by the whole forest, s = min(1, min_j cap_j / r_j): scaling every radius
# by the same factor is the same as changing the perfusion pressure drop, so the trees stay
# at equal pressure, and the tree structure (volume ∝ r²) does not depend on the scale. `mk`
# is min over the other trees of cap_j / r_j (Inf for a single tree, which reduces this to
# min(r, cap)). Returns (root radius, s).
@inline function _scaled_root(cap::T, r::T, mk::T) where {T}
    own = cap > 0 ? cap / r : T(Inf)
    sc = min(one(T), own, mk)
    return (sc == own ? cap : r * sc), sc
end

@inline function fast_root_radius(E, c = E.roots[1])
    k = _rootidx(E, c)
    return _scaled_root(E.rcaps[k], _root_unc(E, c), E.rmin[k])[1]
end

# Refresh, for every tree k, the others' tightest cap ratio (rmin) and the others' unscaled
# volume (rest), from the current state.
function update_forest!(E)
    T = typeof(E.m.lf)
    m = E.m
    K = length(E.roots)
    if K == 1
        E.rmin[1] = T(Inf); E.rest[1] = zero(T)
        return E
    end
    for k in 1:K
        mk = T(Inf); A = zero(T)
        for j in 1:K
            j == k && continue
            c = E.roots[j]; r = _root_unc(E, c)
            E.rcaps[j] > 0 && (mk = min(mk, E.rcaps[j] / r))
            A += r * r * m.lf * E.F.W[c]
        end
        E.rmin[k] = mk; E.rest[k] = A
    end
    return E
end

# Exact radius of segment j: root radius times the betas on its path (unary links count 1).
@inline function fast_radius(E, j)
    T = typeof(E.m.lf)
    m = E.m
    B = one(T); c = j
    while m.parent[c] != 1
        par = m.parent[c]
        (m.left[par] != 0 && m.right[par] != 0) && (B *= m.beta[c])
        c = par
    end
    return fast_root_radius(E, c) * B
end

# Between full refreshes only the radii on changed paths are exact (they feed the candidate's own solve); the others
# drift slightly and feed only collision and spacing tests, whose grid bounds carry the drift. A full refresh is O(N), so
# it runs once the terminal count has grown by RADII_REFRESH since the last one, for sequential and batched growth alike.
const RADII_REFRESH = 0.05
radii_refresh_period(m) = max(1, round(Int, RADII_REFRESH * m.kterm_total))

# All radii from the root down, after which every stored radius is exact. A serial
# breadth-first pass down to a frontier of subtrees, then one depth-first pass per subtree
# in parallel. Each radius is the same product of betas in the same order either way.
function fast_refresh_radii!(E::Engine{T}) where {T}
    m = E.m
    update_forest!(E)
    rmax = zero(T)
    front = Tuple{Int32, T, T}[]          # (segment, product of betas, root radius of its tree)
    for c in E.roots
        r = fast_root_radius(E, c)
        m.radius[c] = r
        r > rmax && (rmax = r)
        push!(front, (c, one(T), r))
    end
    # small trees: one depth-first pass (task start-up would cost more than it saves)
    target = nsegments(m) < 262_144 ? 1 : 64 * Threads.nthreads()
    while 0 < length(front) < target
        nxt = Tuple{Int32, T, T}[]
        for (j, B, r) in front
            a = m.left[j]; b = m.right[j]
            bin = a != 0 && b != 0
            for c in (a, b)
                c == 0 && continue
                Bc = bin ? B * m.beta[c] : B
                m.radius[c] = r * Bc
                push!(nxt, (c, Bc, r))
            end
        end
        for (c, _, _) in nxt
            m.radius[c] > rmax && (rmax = m.radius[c])
        end
        length(nxt) == 0 && break
        front = nxt
    end
    nt = length(front)
    tmax = fill(rmax, nt)
    AK.foreachindex(1:nt; max_tasks = Threads.nthreads(), min_elems = 1) do f
        stack = Tuple{Int32, T}[(front[f][1], front[f][2])]
        r = front[f][3]
        rm = zero(T)
        @inbounds while !isempty(stack)
            j, B = pop!(stack)
            a = m.left[j]; b = m.right[j]
            bin = a != 0 && b != 0
            for c in (a, b)
                c == 0 && continue
                Bc = bin ? B * m.beta[c] : B
                rc = r * Bc
                m.radius[c] = rc
                rc > rm && (rm = rc)
                push!(stack, (c, Bc))
            end
        end
        tmax[f] = max(tmax[f], rm)
    end
    E.rmax = nt == 0 ? rmax : maximum(tmax)
    E.F.since = 0
    E.g.prune = true
    grid_refresh_radii!(E.g, m)
    return E
end

function init_fast!(E::Engine{T}) where {T}
    m = E.m; n = nsegments(m)
    E.F.Rt = zeros(T, n); E.F.W = zeros(T, n)
    for j in postorder(m)                 # every tree of a forest
        fast_node!(E, j)
    end
    fast_refresh_radii!(E)
    return E
end

"""
Fast-mode tentative tree for candidate segment `s` (bifurcation at `pOpt`, new terminal `p`):
O(depth) and allocation free (`ws` is unused). Returns (total volume, radius of s, radius of
new left, radius of new terminal).
"""
@inline function tentative_fast!(ws, E, s, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m; F = E.F; g = m.gamma; full = _fullvol(E)
    lg_s = norm3(sub3(pOpt, m.pos[m.parent[s]]))
    lg_nl = norm3(sub3(m.pos[s], pOpt))
    lg_nr = norm3(sub3(p, pOpt))
    a = m.left[s]; b = m.right[s]
    R_nl = own_baseline(lg_nl); W_nl = lg_nl
    if a != 0 && b == 0
        R_nl += F.Rt[a]; full && (W_nl += F.W[a])
    elseif b != 0 && a == 0
        R_nl += F.Rt[b]; full && (W_nl += F.W[b])
    elseif a != 0 && b != 0
        bA, bB, inv = resistance_combine(g, m.flow[a], m.flow[b], F.Rt[a], F.Rt[b])
        R_nl += inv
        W_nl += bA * bA * F.W[a] + bB * bB * F.W[b]
    end
    R_nr = own_baseline(lg_nr); W_nr = lg_nr
    b_nl, b_nr, inv = resistance_combine(g, m.flow[s], qNew, R_nl, R_nr)
    R_c = own_baseline(lg_s) + inv
    W_c = lg_s + b_nl * b_nl * W_nl + b_nr * b_nr * W_nr
    f_c = m.flow[s] + qNew
    B = one(T)              # product of the radius ratios on the path (the radius of s / root radius)
    c = s; par = m.parent[c]
    while par != 1
        la = m.left[par]; ra = m.right[par]
        lgp = m.glen[par]
        R_a = own_baseline(lgp); W_a = lgp
        if (la != 0) ⊻ (ra != 0)
            R_a += R_c; full && (W_a += W_c); bc = one(T)
        else
            fL = la == c ? f_c : m.flow[la]; RL = la == c ? R_c : F.Rt[la]; WL = la == c ? W_c : F.W[la]
            fR = ra == c ? f_c : m.flow[ra]; RR = ra == c ? R_c : F.Rt[ra]; WR = ra == c ? W_c : F.W[ra]
            bL, bR, inv2 = resistance_combine(g, fL, fR, RL, RR)
            R_a += inv2
            W_a += bL * bL * WL + bR * bR * WR
            bc = la == c ? bL : bR
        end
        B *= bc
        f_c = m.flow[par] + qNew; R_c = R_a; W_c = W_a
        c = par; par = m.parent[c]
    end
    lf = m.lf
    k = _rootidx(E, c)                    # c is now the root segment of s's tree
    rroot, sc = _scaled_root(E.rcaps[k], upow((lf * R_c * f_c) / m.pdrop, T(0.25)), E.rmin[k])
    r_s = rroot * B
    # total volume: this tree's tentative volume plus the other trees' at the forest's common
    # scale (zero for a single tree)
    return rroot * rroot * lf * W_c + sc * sc * E.rest[k], r_s, r_s * b_nl, r_s * b_nr
end

# Middle-stage flags, maintained incrementally: a segment qualifies if it is listed as a
# terminal or one of its descendants within three generations is (getConnectableSegments-
# ForMiddleStage expressed from the ancestor side).
function _middle_of(E::Engine, j)
    E.tcount[j] > 0 && return true
    m = E.m
    stack = Tuple{Int32, Int}[(Int32(j), 0)]
    while !isempty(stack)
        i, d = pop!(stack)
        for c in (m.left[i], m.right[i])
            c == 0 && continue
            E.tcount[c] > 0 && return true
            d + 1 < 3 && push!(stack, (c, d + 1))
        end
    end
    return false
end

function middle_init!(E::Engine)
    resize!(E.middle, nsegments(E.m))
    for j in eachindex(E.middle)
        E.middle[j] = j != 1 && _middle_of(E, j)
    end
    return E
end

# After a commit at s with new segments nl, nr: only s, nl, nr and s's two nearest
# ancestors can change membership.
function middle_update!(E::Engine, s, nl, nr)
    m = E.m
    resize!(E.middle, nsegments(m))
    for j in (nl, nr, s)
        E.middle[j] = _middle_of(E, j)
    end
    a = m.parent[s]
    for _ in 1:2
        a == 1 && break
        E.middle[a] = _middle_of(E, a)
        a = m.parent[a]
    end
    return E
end

# Topology part of a fast-mode commit (shared with batched growth): split s at pOpt, hang its
# old distal part (new left, nl) and the new terminal (nr) below it, and add qNew to the flow
# of every ancestor. Returns (nl, nr).
function commit_structure!(E::Engine{T}, s::Int32, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m; F = E.F
    nl = Int32(nsegments(m) + 1); nr = Int32(nl + 1)
    a = m.left[s]; b = m.right[s]
    push!(m.pos, m.pos[s], p)
    push!(m.radius, m.radius[s], m.radius[s])
    push!(m.kterm, m.kterm[s], Int32(1))
    push!(m.resistance, zero(T), zero(T))
    push!(m.flow, m.flow[s], qNew)
    push!(m.beta, one(T), one(T))
    push!(m.parent, s, s)
    push!(m.left, a, 0); push!(m.right, b, 0)
    push!(m.glen, zero(T), zero(T))
    push!(F.Rt, zero(T), zero(T)); push!(F.W, zero(T), zero(T))
    m.parent[a == 0 ? 1 : a] = nl
    m.parent[b == 0 ? 1 : b] = nl
    m.pos[s] = pOpt
    m.kterm[s] += Int32(1)
    m.flow[s] += qNew
    m.left[s] = nl; m.right[s] = nr
    push!(m.terminals, nr)
    push!(m.listed, false, true)
    push!(E.tcount, Int8(0), Int8(1))
    m.kterm_total += 1
    update_flow!(m, m.parent[s], qNew)
    refresh_glen!(m, s); refresh_glen!(m, nl); refresh_glen!(m, nr)
    return nl, nr
end

"Fast-mode commit: same topology update as `commit!`, then an O(depth) hydraulic update."
function commit_fast!(E::Engine{T}, s::Int32, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m; F = E.F
    nl, nr = commit_structure!(E, s, pOpt, p, qNew)
    fast_node!(E, nl); fast_node!(E, nr)
    j = s
    while j != 1
        fast_node!(E, j)
        j = m.parent[j]
    end
    update_length_factor!(m)
    # exact radii on the updated path and the new segments (as batched rounds do); full refresh periodically
    for j in (nl, nr)
        m.radius[j] = fast_radius(E, j)
        m.radius[j] > E.rmax && (E.rmax = m.radius[j])
    end
    j = s
    while j != 1
        r = fast_radius(E, j); m.radius[j] = r
        r > E.rmax && (E.rmax = r)
        j != s && length(E.g.rref) >= j && E.g.rref[j] > 0 && (E.g.growth = max(E.g.growth, r / E.g.rref[j]))
        j = m.parent[j]
    end
    F.since += 1
    every = F.every > 0 ? F.every : radii_refresh_period(m)
    F.since >= every && fast_refresh_radii!(E)
    grid_remove!(E.g, s)
    grid_insert!(E.g, m, s); grid_insert!(E.g, m, nl); grid_insert!(E.g, m, nr)
    return
end
# (grid_insert! raises the cell radius bounds for the three updated segments; every other
# stored radius changes only in fast_refresh_radii!, which recomputes all bounds.)
