# SPDX-License-Identifier: GPL-3.0-or-later
# Ports the state and whole-tree hydraulics of HctccoModel (src/HctccoModel.h/.ih,
# HCT-CCO commit 9bc141f).
#
# Indexing: Julia index = upstream index + 1. The special zero-length root segment
# (upstream 0) is index 1 and the first real segment (upstream 1) is index 2. The child
# sentinel "no child" (upstream 0) is stored as 0. Upstream overloads index 0 as both
# "root" and "no child"; each port below says which one it means.

"""
    CCOModel{T}

Struct-of-arrays form of upstream `HctccoModel`: one entry per segment (the segment ends at
its distal point `pos[i]` and starts at `pos[parent[i]]`), plus global hydraulic state.
"""
mutable struct CCOModel{T}
    pos::Vector{P3{T}}        # myCoordinate
    radius::Vector{T}         # myRadius
    kterm::Vector{Int32}      # myKTerm
    resistance::Vector{T}     # myResistance
    flow::Vector{T}           # myFlow
    beta::Vector{T}           # myBeta
    parent::Vector{Int32}     # myVectParent
    left::Vector{Int32}       # myVectChildren[i].first
    right::Vector{Int32}      # myVectChildren[i].second
    terminals::Vector{Int32}  # myVectTerminals (never shrinks; see `listed`)
    listed::Vector{Bool}      # true if the segment appears in `terminals`
    glen::Vector{T}           # |pos[i] - pos[parent[i]]| (geometric length; 0 for the root)
    params::CCOParams
    nterm::Int                # my_NTerm
    kterm_total::Int          # myKTerm
    rsupp::T
    rperf::T
    qterm::T
    pperf::T
    pterm::T
    pdrop::T
    lf::T                     # myLengthFactor
    gamma::T
    root_radius_max::T
end

"""
    nsegments(m) -> Int

The number of segments of a tree.
"""
nsegments(m::CCOModel) = length(m.pos)

function CCOModel{T}(p::CCOParams) where {T}
    # 3D branch of the HctccoModel(aPerf, nTerm, ...) constructor.
    rsupp = upow(T(3) * T(p.aperf) / (T(4) * T(π) * T(p.nterm)), one(T) / T(3))
    rperf = upow(T(3) * T(p.aperf) / (T(4) * T(π)), one(T) / T(3))
    qterm = T(QPERF) / T(p.nterm)
    pterm = p.nterm > 250 ? T(7.98e3) : T(8.4e3)
    E = T[]
    return CCOModel{T}(
        P3{T}[], copy(E), Int32[], copy(E), copy(E), copy(E), Int32[], Int32[],
        Int32[], Int32[], Bool[], copy(E), p, p.nterm, 1, rsupp, rperf, qterm,
        T(PPERF), pterm, T(PPERF) - pterm, one(T), T(p.gamma),
        T(p.root_radius_max)
    )
end

const KEPS = 1.0e-12

@inline seglen(m::CCOModel{T}, j) where {T} = j == 1 ? zero(T) : m.glen[j] * m.lf
@inline refresh_glen!(m::CCOModel, j) = (m.glen[j] = j == 1 ? zero(eltype(m.glen)) : norm3(sub3(m.pos[j], m.pos[m.parent[j]])))

# HctccoModel::updateLengthFactor
update_length_factor!(m::CCOModel{T}) where {T} = (m.lf = sqrt(T(m.kterm_total + 1)) * m.rsupp / m.rperf; m)

# HctccoModel::updateFlow: walk from j to the root (upstream index 0) bumping kTerm and flow.
function update_flow!(m::CCOModel, j::Integer, q)
    while j != 1
        m.kterm[j] += Int32(1)
        m.flow[j] += q
        j = m.parent[j]
    end
    return m
end

# HctccoModel::applyResistanceUpdate, written as a pure function of the children's flow and
# resistance. Returns (betaL, betaR, 1/denom).
@inline function resistance_combine(gamma::T, segFlow::T, brotherFlow::T, segHydro::T, brotherHydro::T) where {T}
    keps = T(KEPS)
    (!isfinite(segFlow) || segFlow < keps) && (segFlow = keps)
    (!isfinite(brotherFlow) || brotherFlow < keps) && (brotherFlow = keps)
    (!isfinite(segHydro) || segHydro < keps) && (segHydro = keps)
    (!isfinite(brotherHydro) || brotherHydro < keps) && (brotherHydro = keps)
    ratioArgDen = brotherFlow * brotherHydro
    ratioArgNum = segFlow * segHydro
    (!isfinite(ratioArgNum) || ratioArgNum <= 0) && (ratioArgNum = keps)
    (!isfinite(ratioArgDen) || ratioArgDen <= 0) && (ratioArgDen = keps)
    ratioR = upow(ratioArgNum / ratioArgDen, T(0.25))
    isfinite(ratioR) || (ratioR = one(T))
    bL = upow(one(T) / (one(T) + upow(ratioR, -gamma)), one(T) / gamma)
    bR = upow(one(T) / (one(T) + upow(one(T) / ratioR, -gamma)), one(T) / gamma)
    isfinite(bL) || (bL = T(0.5))
    isfinite(bR) || (bR = T(0.5))
    denom = ((bL * bL * bL * bL) / segHydro) + ((bR * bR * bR * bR) / brotherHydro)
    (!isfinite(denom) || denom < keps) && (denom = keps)
    return bL, bR, one(T) / denom
end

# Float32 form of the same law with the powers shared: with ρ = ratioR = (num/den)^¼ and
# e = ρ^-γ, bL = (1 + e)^(-1/γ) and bR = (1 + 1/e)^(-1/γ). Two logs and four exps instead of
# five general pows; the Float64 method above stays the bitwise port of upstream.
@inline function resistance_combine(gamma::Float32, segFlow::Float32, brotherFlow::Float32, segHydro::Float32, brotherHydro::Float32)
    T = Float32
    keps = T(KEPS)
    (!isfinite(segFlow) || segFlow < keps) && (segFlow = keps)
    (!isfinite(brotherFlow) || brotherFlow < keps) && (brotherFlow = keps)
    (!isfinite(segHydro) || segHydro < keps) && (segHydro = keps)
    (!isfinite(brotherHydro) || brotherHydro < keps) && (brotherHydro = keps)
    ratioArgDen = brotherFlow * brotherHydro
    ratioArgNum = segFlow * segHydro
    (!isfinite(ratioArgNum) || ratioArgNum <= 0) && (ratioArgNum = keps)
    (!isfinite(ratioArgDen) || ratioArgDen <= 0) && (ratioArgDen = keps)
    lr = T(0.25) * @fastmath(log(ratioArgNum / ratioArgDen))
    isfinite(lr) || (lr = zero(T))
    ig = one(T) / gamma
    e = @fastmath exp(-gamma * lr)
    bL = @fastmath exp(-ig * log(one(T) + e))
    bR = @fastmath exp(-ig * log(one(T) + one(T) / e))
    isfinite(bL) || (bL = T(0.5))
    isfinite(bR) || (bR = T(0.5))
    denom = ((bL * bL * bL * bL) / segHydro) + ((bR * bR * bR * bR) / brotherHydro)
    (!isfinite(denom) || denom < keps) && (denom = keps)
    return bL, bR, one(T) / denom
end

# 8*mu*l/pi with upstream's guard against zero or non-finite length.
@inline function own_resistance(l::T) where {T}
    (!(l > zero(T)) || !isfinite(l)) && (l = T(KEPS))
    return T(8) * T(MU) * l / T(π)
end

# One node of HctccoModel::updateResistanceFromRoot, children already updated.
@inline function resistance_node!(m::CCOModel{T}, j) where {T}
    R = own_resistance(seglen(m, j))
    L = m.left[j]; Rr = m.right[j]
    if L != 0 && Rr == 0
        m.beta[L] = one(T); R += m.resistance[L]
    elseif Rr != 0 && L == 0
        m.beta[Rr] = one(T); R += m.resistance[Rr]
    elseif L != 0 && Rr != 0
        bL, bR, inv = resistance_combine(m.gamma, m.flow[L], m.flow[Rr], m.resistance[L], m.resistance[Rr])
        m.beta[L] = bL; m.beta[Rr] = bR
        R += inv
    end
    m.resistance[j] = R
    return
end

# Root segments: segments whose parent is the sentinel (index 1). One for a tree, several for
# a forest (for example left and right coronary arteries grown together).
tree_roots(m::CCOModel) = Int32[j for j in 2:nsegments(m) if m.parent[j] == 1]

# Postorder (children before parents) of the subtree under j, iteratively; with no j, of every
# tree in root order (for a single tree, the tree under segment 2).
function postorder(m::CCOModel, j = nothing)
    j === nothing && return reduce(vcat, (postorder(m, r) for r in tree_roots(m)); init = Int32[])
    out = Int32[]; stack = Int32[j]
    while !isempty(stack)
        i = pop!(stack); push!(out, i)
        m.left[i] != 0 && push!(stack, m.left[i])
        m.right[i] != 0 && push!(stack, m.right[i])
    end
    return reverse!(out)
end

"""
    tree_model(m, r) -> CCOModel

The tree under root segment `r` of a forest as a standalone model (segments renumbered, the
root at index 2), for per-tree morphometry and export.
"""
function tree_model(m::CCOModel{T}, r::Integer) where {T}
    segs = reverse(postorder(m, Int32(r)))          # parents before children
    idx = Dict{Int32, Int32}(Int32(1) => Int32(1))
    for (k, j) in enumerate(segs)
        idx[j] = Int32(k + 1)
    end
    t = deepcopy(m)
    sel = vcat(Int32(1), segs)
    for f in (:pos, :radius, :kterm, :resistance, :flow, :beta, :listed, :glen)
        setfield!(t, f, getfield(m, f)[sel])
    end
    t.parent = Int32[j == 1 ? Int32(2) : idx[m.parent[j]] for j in sel]
    t.left = Int32[m.left[j] == 0 || j == 1 ? Int32(0) : idx[m.left[j]] for j in sel]
    t.right = Int32[m.right[j] == 0 || j == 1 ? Int32(0) : idx[m.right[j]] for j in sel]
    t.terminals = Int32[idx[j] for j in m.terminals if haskey(idx, j) && j != 1]
    t.kterm_total = count(j -> m.left[j] == 0 && m.right[j] == 0, segs)
    return t
end

# HctccoModel::updateResistanceFromRoot(1)
function update_resistance_from_root!(m::CCOModel, order = postorder(m))
    for j in order
        resistance_node!(m, j)
    end
    return m
end

# HctccoModel::updateRootRadius + updateRadius(1, 1.0): radius = r_root * cumulative beta.
function update_root_radius!(m::CCOModel{T}, order = postorder(m)) where {T}
    Q = m.flow[2]
    r = upow((m.resistance[2] * Q) / m.pdrop, T(0.25))
    (m.root_radius_max > 0 && r > m.root_radius_max) && (r = m.root_radius_max)
    m.radius[2] = r
    cum = similar(m.radius); cum[2] = one(T)
    for k in length(order):-1:1          # reverse postorder visits parents before children
        j = order[k]
        L = m.left[j]; R = m.right[j]
        if L != 0 && R != 0
            cum[L] = cum[j] * m.beta[L]; cum[R] = cum[j] * m.beta[R]
            m.radius[L] = r * cum[L]; m.radius[R] = r * cum[R]
        elseif L != 0
            cum[L] = cum[j]; m.radius[L] = r * cum[L]
        elseif R != 0
            cum[R] = cum[j]; m.radius[R] = r * cum[R]
        end
    end
    return m
end
