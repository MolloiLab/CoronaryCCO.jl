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

@inline function fast_root_radius(E::Engine{T}) where {T}
    m = E.m
    r = upow((m.lf * E.F.Rt[2] * m.flow[2]) / m.pdrop, T(0.25))
    (m.root_radius_max > 0 && r > m.root_radius_max) && (r = m.root_radius_max)
    return r
end

# Exact radius of segment j: root radius times the betas on its path (unary links count 1).
function fast_radius(E::Engine{T}, j) where {T}
    m = E.m
    B = one(T); c = j
    while m.parent[c] != 1
        par = m.parent[c]
        (m.left[par] != 0 && m.right[par] != 0) && (B *= m.beta[c])
        c = par
    end
    return fast_root_radius(E) * B
end

# All radii from the root down (DFS), after which every stored radius is exact.
function fast_refresh_radii!(E::Engine{T}) where {T}
    m = E.m
    r = fast_root_radius(E)
    m.radius[2] = r
    stack = Tuple{Int32, T}[(Int32(2), one(T))]
    rmax = r
    @inbounds while !isempty(stack)
        j, B = pop!(stack)
        a = m.left[j]; b = m.right[j]
        bin = a != 0 && b != 0
        for c in (a, b)
            c == 0 && continue
            Bc = bin ? B * m.beta[c] : B
            rc = r * Bc
            m.radius[c] = rc
            rc > rmax && (rmax = rc)
            push!(stack, (c, Bc))
        end
    end
    E.rmax = rmax
    E.F.since = 0
    E.g.prune = true
    grid_refresh_radii!(E.g, m)
    return E
end

function init_fast!(E::Engine{T}) where {T}
    m = E.m; n = nsegments(m)
    E.F.Rt = zeros(T, n); E.F.W = zeros(T, n)
    for j in postorder(m)
        fast_node!(E, j)
    end
    fast_refresh_radii!(E)
    return E
end

"""
Fast-mode tentative tree for candidate segment `s` (bifurcation at `pOpt`, new terminal `p`):
O(depth). Returns (total volume, radius of s, radius of new left, radius of new terminal).
"""
function tentative_fast!(ws::Workspace{T}, E::Engine{T}, s, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m; F = E.F; g = m.gamma; full = E.volume === :full
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
    pb = ws.beta
    k = 0
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
        k += 1
        k > length(pb) && resize!(pb, 2k)
        pb[k] = bc
        f_c = m.flow[par] + qNew; R_c = R_a; W_c = W_a
        c = par; par = m.parent[c]
    end
    lf = m.lf
    rroot = upow((lf * R_c * f_c) / m.pdrop, T(0.25))
    (m.root_radius_max > 0 && rroot > m.root_radius_max) && (rroot = m.root_radius_max)
    B = one(T)
    for i in k:-1:1
        B *= pb[i]
    end
    r_s = rroot * B
    return rroot * rroot * lf * W_c, r_s, r_s * b_nl, r_s * b_nr
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

"Fast-mode commit: same topology update as `commit!`, then an O(depth) hydraulic update."
function commit_fast!(E::Engine{T}, s::Int32, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
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
    fast_node!(E, nl); fast_node!(E, nr)
    j = s
    while j != 1
        fast_node!(E, j)
        j = m.parent[j]
    end
    update_length_factor!(m)
    # exact radii on the updated path and the new segments; full refresh periodically
    for j in (s, nl, nr)
        m.radius[j] = fast_radius(E, j)
        m.radius[j] > E.rmax && (E.rmax = m.radius[j])
    end
    F.since += 1
    every = F.every > 0 ? F.every : max(1, nsegments(m) ÷ 16384)
    F.since >= every && fast_refresh_radii!(E)
    grid_remove!(E.g, s)
    grid_insert!(E.g, m, s); grid_insert!(E.g, m, nl); grid_insert!(E.g, m, nr)
    return
end
# (grid_insert! raises the cell radius bounds for the three updated segments; every other
# stored radius changes only in fast_refresh_radii!, which recomputes all bounds.)
