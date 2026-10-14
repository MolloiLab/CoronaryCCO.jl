# SPDX-License-Identifier: GPL-3.0-or-later
# Ports the growth loop of HCT-CCO commit 9bc141f: CoronaryGrowth::expandTree and the
# HctccoModel methods it drives (generateNewLocation, getN_NearestSegments,
# isIntersectingTree, isAddable, isDistanceValid, isValidConnection, isConnectable,
# getConnectableSegmentsForMiddleStage, terminalFlowWeight).
#
# Upstream evaluates each candidate bifurcation on a deep copy of the whole tree. Here every
# candidate reads the shared base tree and writes only its own workspace, so candidates run
# in parallel; the winner is picked in upstream's candidate order and committed once.

struct Workspace{T}
    cum::Vector{T}
    vol::Vector{T}
    beta::Vector{T}
end
Workspace{T}() where {T} = Workspace{T}(T[], T[], T[])

# Fast-mode state (see fast.jl).
mutable struct FastState{T}
    Rt::Vector{T}         # baseline resistance (resistance / length factor)
    W::Vector{T}          # subtree volume factor (see above)
    since::Int            # commits since the last full radius refresh
    every::Int            # refresh period (0: adaptive)
end
FastState{T}() where {T} = FastState{T}(T[], T[], 0, 0)

mutable struct GrowStats
    steps::Int
    locations::Int
    candidates::Int
    addable_calls::Int
    kamiya_iters::Int
    accepted::Int
    seeds_found::Int
    t_location::Float64
    t_nearest::Float64
    t_candidates::Float64
    t_commit::Float64
    t_view::Float64       # batched: snapshot and upload of the frozen tree
    t_refresh::Float64    # batched: hydraulic and radius refresh
    rounds::Int
end
GrowStats() = GrowStats(0, 0, 0, 0, 0, 0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)

mutable struct Engine{T}
    m::CCOModel{T}
    L::Layout{T}
    g::SegGrid{T}
    dom::MaskDomain{T}
    rng::CRand
    p::CCOParams
    seed::UInt64
    ws::Vector{Workspace{T}}
    rmax::T
    tcount::Vector{Int8}      # occurrences of each segment in the terminal list
    middle::Vector{Bool}      # scratch: middle-stage connectable set
    mode::Symbol              # :fast (production) or :exact (test-only bitwise reference)
    volume::Symbol            # fast mode only: :upstream or :full (see fast.jl)
    sampling::Symbol          # :voxel (upstream: integer voxel points) or :continuous
    fix_base::Bool            # keep the imaged base tree's centerlines where branches attach
    collision::Symbol         # :current (upstream) or :final radii in collision/spacing tests
    cscale::T                 # radius scale used by those tests (1 for :current)
    lscale::T                 # matching length scale for the degeneracy test (1 for :current)
    max_angle_cos::T          # cos of the largest allowed turn at a new bifurcation (-1: none)
    isbase::Vector{Bool}      # segment lies on an imaged base-tree vessel
    treeid::Vector{UInt8}     # tree of each segment (index into roots; 0 for the root sentinel)
    roots::Vector{Int32}      # root segment of each tree (one for a tree, several for a forest)
    rcaps::Vector{T}          # root radius cap of each tree (0: none)
    rest::Vector{T}           # per tree: unscaled volume of the other trees (update_forest!)
    rmin::Vector{T}           # per tree: tightest cap ratio cap/r of the other trees (update_forest!)
    nd::Matrix{T}             # nearest_k! buffers of the sequential step
    nj::Matrix{Int32}
    F::FastState{T}
    stats::GrowStats
end

# Tree of every segment: the index in `roots` of the root its path leads to.
function _tree_ids(m::CCOModel, roots)
    id = zeros(UInt8, nsegments(m))
    for j in 2:nsegments(m)
        k = j
        while m.parent[k] != 1
            k = m.parent[k]
        end
        id[j] = UInt8(something(findfirst(==(k), roots), 0))
    end
    return id
end

# Territories in a forest. A domain voxel value >= 0xC0 lists in its low six bits the trees a
# terminal there may join (bit k-1 for tree k); any other inside value (e.g. 255) allows all.
# This lets trees interdigitate across a boundary band instead of meeting at a hard seam.
@inline function _allowed_bits(d, p::P3{T}) where {T}
    q = (unsafe_trunc(Int32, floor(p[1] + T(0.5))), unsafe_trunc(Int32, floor(p[2] + T(0.5))), unsafe_trunc(Int32, floor(p[3] + T(0.5))))
    in_bounds(d, q) || return 0x3f
    v = img_at(d, q)
    return v >= 0xc0 ? (v & 0x3f) : 0x3f
end
@inline _tree_ok(treeid, j, bits::UInt8) = (t = @inbounds(treeid[j]); t == 0 || ((bits >> (t - 0x01)) & 0x01) == 0x01)

function Engine(
        m::CCOModel{T}, dom::MaskDomain{T}; seed::Integer = 0, cell::Real = 2.0,
        mode::Symbol = :fast, volume::Symbol = :upstream, sampling::Symbol = :voxel,
        fix_base::Bool = false, collision::Symbol = :current,
        max_branch_angle::Union{Nothing, Real} = nothing, root_caps = nothing
    ) where {T}
    collision in (:current, :final) || throw(ArgumentError("collision must be :current or :final"))
    sampling in (:voxel, :continuous) || throw(ArgumentError("sampling must be :voxel or :continuous"))
    mode in (:exact, :fast) || throw(ArgumentError("mode must be :exact or :fast"))
    volume in (:upstream, :full) || throw(ArgumentError("volume must be :upstream or :full"))
    (mode === :exact && volume !== :upstream) &&
        throw(ArgumentError("exact mode reproduces upstream, which uses volume = :upstream"))
    roots = tree_roots(m)
    (mode === :exact && length(roots) > 1) &&
        throw(ArgumentError("a forest of $(length(roots)) trees needs mode = :fast"))
    caps = root_caps === nothing ? fill(T(m.root_radius_max), length(roots)) :
        T.(root_caps isa Real ? fill(root_caps, length(roots)) : collect(root_caps))
    length(caps) == length(roots) ||
        throw(ArgumentError("root_caps has $(length(caps)) entries for $(length(roots)) trees"))
    L = mode === :exact ? init_layout!(Layout{T}(), m) : Layout{T}()   # preorder layout: exact mode only
    lo = (T(dom.lower[1]), T(dom.lower[2]), T(dom.lower[3]))
    hi = (T(dom.upper[1]), T(dom.upper[2]), T(dom.upper[3]))
    g = grid_build!(SegGrid{T}(lo, hi, cell), m)
    tc = zeros(Int8, nsegments(m))
    for t in m.terminals
        tc[t] += Int8(1)
    end
    E = Engine{T}(
        m, L, g, dom, CRand(seed), m.params, UInt64(seed),
        [Workspace{T}() for _ in 1:m.params.num_neighbors], zero(T), tc, Bool[],
        mode, volume, sampling, fix_base, collision, one(T), one(T),
        max_branch_angle === nothing ? -one(T) : T(cosd(max_branch_angle)),
        trues(nsegments(m)) |> Vector{Bool}, _tree_ids(m, roots), roots, caps, zeros(T, length(roots)),
        fill(T(Inf), length(roots)),
        zeros(T, m.params.num_neighbors, 1), zeros(Int32, m.params.num_neighbors, 1),
        FastState{T}(), GrowStats()
    )
    update_rmax!(E)
    middle_init!(E)
    mode === :fast && init_fast!(E)
    return E
end

update_rmax!(E::Engine) = (E.rmax = maximum(@view E.m.radius[2:end]); E)

# ---------------------------------------------------------------------------------------
# Scans over the tree through the spatial grid (same answers as upstream's full scans).

@inline segA(m, j) = @inbounds m.pos[j]
@inline segB(m, j) = @inbounds m.pos[m.parent[j]]

# HctccoModel::isIntersectingTree(ptA, ptB, r, idExcept) without the zero-index shortcut.
@inline function scan_hit(E, A::P3{T}, B::P3{T}, r::T, e1, e2, e3) where {T}
    m = E.m
    c = E.cscale; rq = r * c
    test(j) = j != e1 && j != e2 && j != e3 && is_intersecting(A, B, rq, segA(m, j), segB(m, j), m.radius[j] * c)
    return grid_any_r(test, E.g, A, B, r -> rq + r * c, E.rmax)
end

# Upstream's idExcept==0 shortcut: an excluded "index 0" (either "no child" or the root
# segment) makes the whole test report no intersection.
@inline _has_zero(m, idL, idR, idPar) = idL == 0 || idR == 0 || idPar == 1

# HctccoModel::isIntersectingTree(ptA, ptB, r, idSeg): the three segments of a bifurcation
# at ptB joining ptA to segment idSeg.
@inline function intersects_bifurcation(E, ptA::P3{T}, ptB::P3{T}, r::T, s) where {T}
    m = E.m
    par = m.parent[s]
    idL = m.left[par]; idR = m.right[par]
    if !_has_zero(m, idL, idR, par)
        scan_hit(E, ptB, m.pos[par], r, idL, idR, par) && return true
    end
    idL = m.left[s]; idR = m.right[s]
    if !_has_zero(m, idL, idR, s)
        scan_hit(E, ptB, m.pos[s], r, s, idL, idR) && return true
    end
    return scan_hit(E, ptB, ptA, r, s, 0, 0)
end

# HctccoModel::isDistanceValid(p, minPerfusionDistance)
@inline function distance_valid(E, p::P3{T}, mpd::T) where {T}
    m = E.m
    c = E.cscale
    seg_bad(j) = proj_distance(segA(m, j), segB(m, j), p) < cmax(mpd, m.radius[j] * c * T(2))
    term_bad(j) = E.tcount[j] > 0 && norm3(sub3(m.pos[j], p)) < mpd
    grid_any_r(seg_bad, E.g, p, p, r -> cmax(mpd, r * c * T(2)), E.rmax) && return false
    grid_any(term_bad, E.g, p, p, mpd) && return false
    return true
end

# ---------------------------------------------------------------------------------------
# Candidate locations (generateNewLocation / generateALocation).

function distance_threshold(m::CCOModel{T}) where {T}
    rp = m.rperf
    return upow((T(4) * T(π) * rp * rp * rp) / (T(3) * T(m.kterm_total)), one(T) / T(3))
end

function generate_a_location!(E::Engine{T}, thr::T, rng::CRand = E.rng) where {T}
    m = E.m
    q = random_point!(E.dom, rng)
    res = (T(q[1]), T(q[2]), T(q[3]))
    if E.sampling === :continuous
        # Upstream only proposes integer voxel points, which caps the terminal count near the
        # number of admissible voxels. Jitter within the voxel and re-check the wall clearance.
        u() = T(crand!(rng)) / T(2147483648) - T(0.5)
        res = (res[1] + u(), res[2] + u(), res[3] + u())
        interpolate_distance(E.dom, res) >= E.dom.min_dist || return res, false
    end
    nt = length(m.terminals)
    first = m.terminals[1]
    # terminals after the first list entry within thr of res
    tclose(j) = (E.tcount[j] - (j == first ? 1 : 0)) > 0 && !(norm3(sub3(m.pos[j], res)) > thr)
    grid_any(tclose, E.g, res, res, thr) && return res, false
    c = E.cscale
    sclose(j) = j > nt && !(proj_distance(segA(m, j), segB(m, j), res) > m.radius[j] * c)
    grid_any_r(sclose, E.g, res, res, r -> r * c, E.rmax) && return res, false
    return res, true
end

function generate_new_location!(E::Engine{T}, rng::CRand = E.rng) where {T}
    thr = distance_threshold(E.m)
    ntr = E.p.nb_try_candidate
    n = ntr
    while true
        n -= 1
        res, found = generate_a_location!(E, thr, rng)
        found && return res
        if n == 0
            n = ntr; thr *= T(0.9)
        end
    end
    return
end

# ---------------------------------------------------------------------------------------
# Connectivity rules.

function restriction_level(E::Engine)
    m = E.m; p = E.p
    prog = m.nterm == 0 ? 0.0 : clamp(Float64(m.kterm_total) / Float64(m.nterm), 0.0, 1.0)
    if p.extension_mode
        return prog < 0.25 ? 0.5 : 1.0
    end
    return prog < 0.1 ? 0.0 : prog < 0.35 ? 0.5 : 1.0
end

# Engine/TreeView accessors used by the shared candidate code.
@inline _isfast(E::Engine) = E.mode === :fast
@inline _fullvol(E::Engine) = E.volume === :full
@inline _level(E::Engine) = E.p.anywhere_mode ? 0.0 : restriction_level(E)

# HctccoModel::isConnectable: listed as a terminal, or a parent of a listed terminal.
@inline function is_connectable(E, j)
    m = E.m
    E.tcount[j] > 0 && return true
    L = m.left[j]; R = m.right[j]
    return (L != 0 && E.tcount[L] > 0) || (R != 0 && E.tcount[R] > 0)
end

# HctccoModel::getConnectableSegmentsForMiddleStage is kept incrementally in E.middle
# (middle_init!/middle_update! in fast.jl): terminals plus up to three generations of
# ancestors; upstream's 4-6 generation branch can never add anything.

# HctccoModel::getN_NearestSegments with stable ordering (ties by segment index), written
# into column `col` of the buffers `bd` (distances) and `bj` (segments), which need at least
# n rows. Returns the number found. Allocation free: a segment registered in several grid
# cells is counted only in the first cell (per axis) shared by its box and the query box.
@inline function nearest_k!(bd, bj, col, E, p::P3{T}, n::Int) where {T}
    m = E.m; g = E.g
    level = _level(E)
    bits = _allowed_bits(E.dom, p)
    R = g.h
    @inbounds while true
        k = 0
        bx = _cellbox(g, p, p, R + g.margin)
        covers_all = bx[1] == 1 && bx[3] == 1 && bx[5] == 1 &&
            bx[2] == g.dims[1] && bx[4] == g.dims[2] && bx[6] == g.dims[3]
        for kz in bx[5]:bx[6], ky in bx[3]:bx[4], kx in bx[1]:bx[2]
            for j in cell_items(g, _lin(g, kx, ky, kz))
                sb = g.box[j]
                (kx == max(sb[1], bx[1]) && ky == max(sb[3], bx[3]) && kz == max(sb[5], bx[5])) || continue
                keep = level >= 1 ? is_connectable(E, j) : level == 0.5 ? E.middle[j] : true
                keep || continue
                (bits == 0x3f || _tree_ok(E.treeid, j, bits)) || continue
                d = proj_distance(segA(m, j), segB(m, j), p)
                # insertion into the sorted top-n list
                if k < n
                    k += 1
                    i = k
                elseif d < bd[n, col] || (d == bd[n, col] && j < bj[n, col])
                    i = n
                else
                    continue
                end
                while i > 1 && (bd[i - 1, col] > d || (bd[i - 1, col] == d && bj[i - 1, col] > j))
                    bd[i, col] = bd[i - 1, col]; bj[i, col] = bj[i - 1, col]
                    i -= 1
                end
                bd[i, col] = d; bj[i, col] = j
            end
        end
        (covers_all || (k >= n && bd[n, col] < R - g.margin)) && return k
        R *= 2
    end
    return 0
end

# ---------------------------------------------------------------------------------------
# Candidate evaluation (HctccoModel::isAddable on a virtual copy of the tree).

# Direction change from u to v no larger than the angle whose cosine is cmin.
@inline function _turn_ok(u::P3{T}, v::P3{T}, cmin::T) where {T}
    nu = norm3(u); nv = norm3(v)
    (nu == 0 || nv == 0) && return true
    return dot3(u, v) >= cmin * nu * nv
end

# HctccoModel::terminalFlowWeight
# (Float64 runs hash and scale in Float64 like upstream; Float32 runs stay in Float32 so the
# same code compiles on GPUs without double precision.)
@inline function terminal_flow_weight(cv::T, p::P3{T}) where {T}
    cv <= 0 && return one(T)
    hashu(salt) = begin
        h = 0xcbf29ce484222325 ⊻ salt
        for i in 1:3
            q = round(Int64, p[i] * T(1000), RoundNearestTiesAway)
            h ⊻= reinterpret(UInt64, q)
            h *= 0x00000100000001b3
        end
        T((h >> 17) & 0x3FFFFFFF) / T(0x40000000)
    end
    u1 = hashu(0x9E3779B97F4A7C15); u2 = hashu(0xC2B2AE3D27D4EB4F)
    u1 < T(1.0e-12) && (u1 = T(1.0e-12))
    z = sqrt(T(-2) * log(u1)) * cos(T(2) * T(π) * u2)
    sigma = sqrt(log(one(T) + cv * cv))
    w = exp(sigma * z - T(0.5) * sigma * sigma)
    return clamp(w, T(0.05), T(20))
end

"""
Tentative tree for candidate segment `s` with bifurcation at `pOpt` and new terminal `p`:
path-only resistance update (every other node's inputs are unchanged), root radius, then a
full radius sweep and volume sweep with the two new segments as virtual nodes. Returns
(total volume, radius of s, radius of new left, radius of new terminal).
"""
function tentative!(ws::Workspace{T}, E::Engine{T}, s, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m; L = E.L; lf = m.lf; g = m.gamma
    ps = m.pos[s]; ppar = m.pos[m.parent[s]]
    l_s = norm3(sub3(pOpt, ppar)) * lf
    l_nl = norm3(sub3(ps, pOpt)) * lf
    l_nr = norm3(sub3(p, pOpt)) * lf
    a = m.left[s]; b = m.right[s]
    R_nl = own_resistance(l_nl)
    if a != 0 && b == 0
        R_nl += m.resistance[a]
    elseif b != 0 && a == 0
        R_nl += m.resistance[b]
    elseif a != 0 && b != 0
        _, _, inv = resistance_combine(g, m.flow[a], m.flow[b], m.resistance[a], m.resistance[b])
        R_nl += inv
    end
    R_nr = own_resistance(l_nr)
    b_nl, b_nr, inv = resistance_combine(g, m.flow[s], qNew, R_nl, R_nr)
    R_c = own_resistance(l_s) + inv
    f_c = m.flow[s] + qNew
    n = length(L.pre)
    beta = ws.beta
    copyto!(beta, 1, L.beta, 1, n)
    c = s
    par = m.parent[c]
    while par != 1
        R_a = own_resistance(seglen(m, par))
        la = m.left[par]; ra = m.right[par]
        if la != 0 && ra == 0
            beta[L.ppos[la]] = one(T); R_a += R_c
        elseif ra != 0 && la == 0
            beta[L.ppos[ra]] = one(T); R_a += R_c
        else
            fL = la == c ? f_c : m.flow[la]; RL = la == c ? R_c : m.resistance[la]
            fR = ra == c ? f_c : m.flow[ra]; RR = ra == c ? R_c : m.resistance[ra]
            bL, bR, inv2 = resistance_combine(g, fL, fR, RL, RR)
            beta[L.ppos[la]] = bL; beta[L.ppos[ra]] = bR
            R_a += inv2
        end
        f_c = m.flow[par] + qNew
        R_c = R_a
        c = par
        par = m.parent[c]
    end
    rroot = upow((R_c * f_c) / m.pdrop, T(0.25))
    (m.root_radius_max > 0 && rroot > m.root_radius_max) && (rroot = m.root_radius_max)
    # Top-down cumulative betas; the old children of s now hang below the new left node.
    ks = L.ppos[s]
    cum = ws.cum; vol = ws.vol
    cum_nl = zero(T); cum_nr = zero(T)
    cum[1] = one(T)
    if ks == 1
        cum_nl = b_nl; cum_nr = b_nr
    end
    @inbounds for k in 2:n
        pk = L.par[k]
        base = pk == ks ? cum_nl : cum[pk]
        cum[k] = L.bin[k] ? base * beta[k] : base
        if k == ks
            cum_nl = cum[k] * b_nl; cum_nr = cum[k] * b_nr
        end
    end
    # Bottom-up volume; at s, add the virtual new-left (with s's old children) and new terminal.
    s_bin = a != 0 && b != 0
    @inbounds for k in n:-1:1
        r = rroot * cum[k]
        lk = L.lc[k]; rk = L.rc[k]
        if k == ks
            r_nl = rroot * cum_nl; r_nr = rroot * cum_nr
            v_nl = (r_nl * r_nl) * l_nl
            s_bin && (v_nl = v_nl + (vol[lk] + vol[rk]))
            v_nr = (r_nr * r_nr) * l_nr
            vol[k] = (r * r) * l_s + (v_nl + v_nr)
        else
            vk = (r * r) * L.len[k]
            vol[k] = (lk != 0 && rk != 0) ? vk + (vol[lk] + vol[rk]) : vk
        end
    end
    return vol[1], rroot * cum[ks], rroot * cum_nl, rroot * cum_nr
end

"""
    try_candidate!(ws, E, p, s) -> (ok, vol, pOpt, qNew, iterations)

HctccoModel::isAddable(p, s, nbIter, tolerance) evaluated against the shared base tree.
"""
@inline function try_candidate!(ws, E, p::P3{T}, s::Int32) where {T}
    m = E.m; dom = E.dom; prm = E.p
    fail = (false, zero(T), p, zero(T), 0)
    (prm.anywhere_mode || is_connectable(E, s)) || return fail
    distance_valid(E, p, T(0.05)) || return fail
    # isValidConnection
    norm3(sub3(p, m.pos[s])) > T(25) && return fail
    rs = _isfast(E) ? fast_radius(E, s) : m.radius[s]
    rs * T(0.7) < T(0.005) && return fail
    ps = m.pos[s]
    par = m.parent[s]
    pParent = m.pos[par]
    newCenter = div3(add3(add3(p, ps), pParent), T(3))
    cDom = itrunc3(newCenter)
    is_inside(dom, cDom) || return fail
    check_no_intersect(dom, itrunc3(p), cDom) || return fail
    check_no_intersect(dom, itrunc3(ps), cDom) || return fail
    pParentDom = itrunc3(pParent)
    check_no_intersect(dom, pParentDom, cDom) || return fail
    qNew = m.qterm * terminal_flow_weight(T(prm.term_flow_cv), p)
    # sCurrent / sNewLeft / sNewRight state carried across iterations.
    rc = rs
    rL = rs; fL = m.flow[s]
    rR = rs; fR = qNew
    pCurrent = div3(add3(pParent, ps), T(2))
    lf = m.lf
    volCurr = zero(T)
    idL = m.left[par]; idR = m.right[par]
    sL = m.left[s]; sR = m.right[s]
    pLDom = itrunc3(ps); pRDom = itrunc3(p)
    for i in 0:(prm.kamiya_iters - 1)
        ok, pOpt, r0, r1, r2 = kamiya_optimization(m.gamma, lf, pCurrent, pParent, rc, ps, rL, fL, p, rR, fR)
        ok || return (false, zero(T), p, zero(T), i + 1)
        if E.fix_base && E.isbase[s]
            # Upstream moves the split segment's end to the optimal bifurcation point, which
            # bends imaged vessels. Keep the junction on the imaged centerline instead.
            ax = sub3(ps, pParent)
            t = clamp(dot3(sub3(pOpt, pParent), ax) / dot3(ax, ax), T(0.05), T(0.95))
            pOpt = add3(pParent, mul3(ax, t))
        end
        if E.max_angle_cos > -1
            # Minimum-work branching (Murray 1926) never sends a daughter backward relative to
            # the parent's flow; reject bifurcations whose new or moved segments turn more
            # than the configured angle.
            up = sub3(pOpt, pParent)
            _turn_ok(up, sub3(ps, pOpt), E.max_angle_cos) || return (false, zero(T), p, zero(T), i + 1)
            _turn_ok(up, sub3(p, pOpt), E.max_angle_cos) || return (false, zero(T), p, zero(T), i + 1)
            # s's existing children now hang below the new left segment (pOpt -> old end of s)
            dl = sub3(ps, pOpt)
            for ch in (m.left[s], m.right[s])
                ch == 0 && continue
                _turn_ok(dl, sub3(m.pos[ch], ps), E.max_angle_cos) || return (false, zero(T), p, zero(T), i + 1)
            end
            if par != 1 && m.parent[par] != 1
                _turn_ok(sub3(pParent, m.pos[m.parent[par]]), up, E.max_angle_cos) ||
                    return (false, zero(T), p, zero(T), i + 1)
            end
        end
        if !_has_zero(m, idL, idR, par)
            scan_hit(E, pOpt, pParent, r0, idL, idR, par) && return (false, zero(T), p, zero(T), i + 1)
        end
        if !_has_zero(m, sL, sR, s)
            scan_hit(E, pOpt, ps, r1, s, sL, sR) && return (false, zero(T), p, zero(T), i + 1)
        end
        scan_hit(E, pOpt, p, r2, s, 0, 0) && return (false, zero(T), p, zero(T), i + 1)
        oDom = itrunc3(pOpt)
        is_inside(dom, oDom) || return (false, zero(T), p, zero(T), i + 1)
        check_no_intersect(dom, pParentDom, oDom) || return (false, zero(T), p, zero(T), i + 1)
        check_no_intersect(dom, pLDom, oDom) || return (false, zero(T), p, zero(T), i + 1)
        check_no_intersect(dom, pRDom, oDom) || return (false, zero(T), p, zero(T), i + 1)
        vol, r_s, r_nl, r_nr = _isfast(E) ? tentative_fast!(ws, E, s, pOpt, p, qNew) :
            tentative!(ws, E, s, pOpt, p, qNew)
        if i == 0
            volCurr = vol
        else
            if abs(volCurr - vol) < T(prm.kamiya_tol)
                l0 = norm3(sub3(pOpt, pParent)) * lf
                l1 = norm3(sub3(pOpt, ps)) * lf
                l2 = norm3(sub3(pOpt, p)) * lf
                # degenerate check (2r <= l); with collision = :final both sides at their final scale
                c2 = T(2) * E.cscale; cl = E.lscale
                if (c2 * r0 <= cl * l0) && (c2 * r1 <= cl * l1) && (c2 * r2 <= cl * l2)
                    return (true, vol, pOpt, qNew, i + 1)
                end
            else
                volCurr = vol
            end
        end
        pCurrent = pOpt
        rc = r_s; rL = r_nl; rR = r_nr
    end
    return (false, zero(T), p, zero(T), prm.kamiya_iters)
end

# ---------------------------------------------------------------------------------------
# Commit the winning candidate and refresh the tree (end of an expandTree step).

function commit!(E::Engine{T}, s::Int32, pOpt::P3{T}, p::P3{T}, qNew::T) where {T}
    m = E.m
    nl = Int32(nsegments(m) + 1); nr = Int32(nl + 1)
    a = m.left[s]; b = m.right[s]
    # new left: the old distal part of s, inheriting its children
    push!(m.pos, m.pos[s], p)
    push!(m.radius, m.radius[s], m.radius[s])
    push!(m.kterm, m.kterm[s], Int32(1))
    push!(m.resistance, m.resistance[s], zero(T))
    push!(m.flow, m.flow[s], qNew)
    push!(m.beta, one(T), one(T))
    push!(m.parent, s, s)
    push!(m.left, a, 0); push!(m.right, b, 0)
    push!(m.glen, zero(T), zero(T))
    # upstream writes myVectParent[child] even when the child index is 0 (the root slot)
    m.parent[a == 0 ? 1 : a] = nl
    m.parent[b == 0 ? 1 : b] = nl
    m.pos[s] = pOpt
    m.kterm[s] = m.kterm[s] + Int32(1)
    m.flow[s] = m.flow[s] + qNew
    m.left[s] = nl; m.right[s] = nr
    push!(m.terminals, nr)
    push!(m.listed, false, true)
    push!(E.tcount, Int8(0), Int8(1))
    m.kterm_total += 1
    update_flow!(m, m.parent[s], qNew)
    refresh_glen!(m, s); refresh_glen!(m, nl); refresh_glen!(m, nr)
    splice_layout!(E.L, m, s, nl, nr)
    full_update!(m, E.L)
    update_rmax!(E)
    grid_remove!(E.g, s)
    grid_insert!(E.g, m, s); grid_insert!(E.g, m, nl); grid_insert!(E.g, m, nr)
    return
end

function _ensure_workspaces!(E::Engine{T}) where {T}
    n = length(E.L.pre)
    for ws in E.ws
        if length(ws.cum) < n + 2
            resize!(ws.cum, 2n + 64); resize!(ws.vol, 2n + 64); resize!(ws.beta, 2n + 64)
        end
    end
    return
end

# One line per step, same layout as cco_trace_step in oracle/cco_oracle.h.
function _trace_step(io, m, i, nbT, volopt)
    n = nsegments(m); s = m.parent[n]
    a = m.pos[n]; b = m.pos[s]
    return @printf(
        io, "%d %d %a %d %d %a %a %a %a %a %a\n", i, nbT, Float64(volopt), n, s - 1,
        Float64(a[1]), Float64(a[2]), Float64(a[3]), Float64(b[1]), Float64(b[2]), Float64(b[3])
    )
end

# Evaluate candidate segment s for terminal point pt (the bifurcation-intersection pre-check,
# then isAddable). Shared by sequential steps and batched rounds.
@inline function eval_candidate(E, ws, pt::P3{T}, s::Int32) where {T}
    m = E.m
    ptB = div3(add3(add3(pt, m.pos[s]), m.pos[m.parent[s]]), T(3))
    rs = _isfast(E) ? fast_radius(E, s) : m.radius[s]
    intersects_bifurcation(E, pt, ptB, rs, s) && return (false, zero(T), pt, zero(T), -1)
    return try_candidate!(ws, E, pt, s)
end

# Score of an accepted candidate (expandTree's side-branch bias).
@inline function candidate_score(E, s, vol::T) where {T}
    prm = E.p; m = E.m
    sb = T(prm.side_branch_bias)
    sb > 0 || return vol
    qt = m.qterm > 0 ? m.qterm : one(T)
    qr = m.flow[s] / qt
    qr < 1 && (qr = one(T))
    return vol * upow(qr, -sb)
end

# collision = :final: radii scale with the root radius for a fixed structure, so testing at
# cap / current root radius anticipates each vessel's size once all flow is in. The
# degeneracy test (2r <= l) compares a radius with a model length, l = geometric length x
# length factor, and the length factor grows as sqrt(terminals + 1); so it compares final
# radii with final lengths, scaled by sqrt((nterm + 1) / (terminals + 1)). Scaling only the
# radii would make every short segment look degenerate early in growth.
#
# For a forest, update_forest! refreshes what each tree needs to know about the others (their
# volume and cap ratios), so every candidate is scored by the total volume of the forest; the
# collision scale is the largest over the trees (conservative).
function update_cscale!(E::Engine{T}) where {T}
    m = E.m
    E.mode === :fast && update_forest!(E)
    if E.collision === :final
        cs = one(T)
        for (k, c) in enumerate(E.roots)
            cap = E.rcaps[k]
            cap > 0 || continue
            r_now = E.mode === :fast ? fast_root_radius(E, c) : m.radius[c]
            cs = max(cs, cap / r_now)
        end
        E.cscale = cs
        E.lscale = max(one(T), sqrt(T(m.nterm + 1) / T(m.kterm_total + 1)))
    end
    return E
end

# One expandTree step: sample locations until one has a valid candidate (or nb_max_search
# tries), commit the minimum-volume candidate.
function sequential_step!(E::Engine{T}, results, parallel::Bool) where {T}
    m = E.m; prm = E.p; st = E.stats
    st.steps += 1
    update_cscale!(E)
    nbsol = 0; nbT = 0
    volopt = -one(T)
    best_s = Int32(0); best_p = (zero(T), zero(T), zero(T)); best_opt = best_p; best_q = zero(T)
    while nbsol == 0 && nbT < prm.nb_max_search
        ta = time_ns()
        pt = generate_new_location!(E)
        tb = time_ns()
        nbT += 1; st.locations += 1
        nc = nearest_k!(E.nd, E.nj, 1, E, pt, prm.num_neighbors)
        vecN = @view E.nj[1:nc, 1]
        tc = time_ns()
        st.t_location += (tb - ta) / 1.0e9; st.t_nearest += (tc - tb) / 1.0e9
        st.candidates += nc
        E.mode === :exact && _ensure_workspaces!(E)
        if parallel && nc > 1
            AK.foreachindex(1:nc; max_tasks = nc, min_elems = 1) do c
                results[c] = eval_candidate(E, E.ws[c], pt, vecN[c])
            end
        else
            for c in 1:nc
                results[c] = eval_candidate(E, E.ws[c], pt, vecN[c])
            end
        end
        st.t_candidates += (time_ns() - tc) / 1.0e9
        for c in 1:nc
            ok, vol, pOpt, qNew, its = results[c]
            its >= 0 && (st.addable_calls += 1)
            st.kamiya_iters += max(its, 0)
            ok || continue
            st.accepted += 1
            score = candidate_score(E, vecN[c], vol)
            if volopt < 0 || volopt > score
                volopt = score; best_s = vecN[c]; best_p = pt; best_opt = pOpt; best_q = qNew
            end
            nbsol += 1
        end
    end
    nbT < prm.nb_max_search && (st.seeds_found += 1)
    td = time_ns()
    if best_s != 0
        push!(E.isbase, E.isbase[best_s], false)
        push!(E.treeid, E.treeid[best_s], E.treeid[best_s])
        E.mode === :fast ? commit_fast!(E, best_s, best_opt, best_p, best_q) :
            commit!(E, best_s, best_opt, best_p, best_q)
        middle_update!(E, best_s, Int32(nsegments(m) - 1), Int32(nsegments(m)))
    end
    st.t_commit += (time_ns() - td) / 1.0e9
    return nbT, volopt
end

"""
    grow!(E; parallel = true, progress = 0, growth = :batched, n_seq = 2_000,
          batch_frac = 0.01, batch_max = 16_384)

CoronaryGrowth::expandTree: add terminals until `nterm`. Candidates of each sampled location
are evaluated in parallel; the minimum-volume one (first in candidate order on ties) wins.

`growth = :batched` (fast mode only) grows the first `n_seq` terminals the same way, then
adds terminals in rounds of up to `batch_frac` × (current terminals) at once (at most
`batch_max`): every location of a round is evaluated in parallel against the same tree, a
conflict-free subset is committed, and the hydraulics are refreshed once per round. See
`docs/deviations.md`. `device` places the evaluation of each round: an array constructor
(`CuArray`, `ROCArray`, `MtlArray`, ...) or a KernelAbstractions backend; the default runs on
CPU threads.
"""
function grow!(
        E::Engine{T}; parallel::Bool = true, progress::Int = 0, trace::Union{Nothing, IO} = nothing,
        steps::Int = typemax(Int), growth::Symbol = E.mode === :fast ? :batched : :sequential, n_seq::Integer = 2_000,
        batch_frac::Real = 0.01, batch_max::Integer = 16_384, device = identity
    ) where {T}
    growth in (:sequential, :batched) || throw(ArgumentError("growth must be :sequential or :batched"))
    (growth === :batched && E.mode !== :fast) && throw(ArgumentError("growth = :batched needs mode = :fast"))
    m = E.m; st = E.stats
    reseed!(E.rng, E.seed)            # expandTree re-seeds rand()
    nseed = m.nterm
    results = Vector{Tuple{Bool, T, P3{T}, T, Int}}(undef, E.p.num_neighbors)
    t0 = time()
    total = min(nseed - 1, steps)
    i = 0
    nextp = progress
    pool = nothing
    while i < total
        nb = growth === :batched && i >= n_seq ?
            min(total - i, Int(batch_max), max(1, floor(Int, batch_frac * m.kterm_total))) : 1
        if nb == 1
            nbT, volopt = sequential_step!(E, results, parallel)
            i += 1
            trace === nothing || _trace_step(trace, m, i, nbT, volopt)
        else
            pool === nothing && (pool = BatchPool{T}(; to = _mover(device)))
            i += batched_round!(E, pool, nb; parallel = parallel)
        end
        if progress > 0 && i >= nextp
            @info "step $i/$(nseed - 1)" segments = nsegments(m) elapsed = round(time() - t0, digits = 1)
            nextp += progress
        end
    end
    growth === :batched && i > 0 && fast_refresh_radii!(E)
    st.seeds_found != nseed && @warn "All seeds not found due to too large domain constraints ($(st.seeds_found) over $nseed)"
    return E
end
