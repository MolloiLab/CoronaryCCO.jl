# SPDX-License-Identifier: GPL-3.0-or-later
# Ports the tree representation and constructor of HctccoModel (src/HctccoModel.h,
# HCT-CCO commit 9bc141f).

"""
    CCOTree

Struct-of-arrays form of upstream `HctccoModel`'s `std::vector<Segment>` plus its parent,
children and terminal vectors. Arrays are plain `AbstractVector`s, so the same layout can
live on any KernelAbstractions backend.

Indexing: Julia index = upstream index + 1. The special zero-length root segment
(upstream 0) is index 1 and the first real segment (upstream 1) is index 2. The child
sentinel "no child" (upstream 0) is stored as 0.
"""
mutable struct CCOTree{VF<:AbstractVector{Float64},VI<:AbstractVector{Int32}}
    # Segment fields, one entry per segment (upstream struct Segment).
    x::VF                 # myCoordinate (distal point)
    y::VF
    z::VF
    radius::VF            # myRadius
    kterm::VI             # myKTerm
    resistance::VF        # myResistance
    flow::VF              # myFlow
    beta::VF              # myBeta
    # Topology.
    parent::VI            # myVectParent
    left::VI              # myVectChildren[i].first
    right::VI             # myVectChildren[i].second
    terminals::VI         # myVectTerminals
    # Global state.
    params::CCOParams
    kterm_total::Int      # myKTerm
    rsupp::Float64        # myRsupp
    rperf::Float64        # my_rPerf
    qterm::Float64        # my_qTerm
    pterm::Float64        # my_pTerm
    pdrop::Float64        # my_pDrop
    length_factor::Float64
end

nsegments(t::CCOTree) = length(t.x)

"""
    CCOTree(params, center; radius = 1.0)

Mirror of the 3D `HctccoModel(aPerf, nTerm, aDomCtr, aRadius, ...)` constructor: builds
the special root segment and the first segment ending at `center` (the domain center).
"""
function CCOTree(p::CCOParams, center::NTuple{3,Float64}; radius::Float64 = 1.0)
    p.nterm >= 1 || throw(ArgumentError("nterm must be >= 1"))
    rsupp = cpow(3.0 * p.aperf / (4.0 * π * p.nterm), 1.0 / 3.0)
    rperf = cpow(3.0 * p.aperf / (4.0 * π), 1.0 / 3.0)
    qterm = QPERF / p.nterm
    pterm = p.nterm > 250 ? 7.98e3 : 8.40e3
    F, I = Float64[], Int32[]
    t = CCOTree(copy(F), copy(F), copy(F), copy(F), copy(I), copy(F), copy(F), copy(F),
                copy(I), copy(I), copy(I), copy(I),
                p, 1, rsupp, rperf, qterm, pterm, PPERF - pterm, 1.0)

    # Special root segment (upstream index 0): proximal point of the first segment.
    cx, cy, cz = center
    push_segment!(t; x = cx, y = cy + rperf, z = cz, radius = radius, kterm = 0, parent = 1)
    update_length_factor!(t)

    # First segment (upstream index 1) ending at the tree center.
    len = norm3(cx - t.x[1], cy - t.y[1], cz - t.z[1]) * t.length_factor
    push_segment!(t; x = cx, y = cy, z = cz, radius = radius, kterm = 1, parent = 1,
                  resistance = 8.0 * MU * len / π, flow = qterm, beta = 1.0)
    push!(t.terminals, Int32(2))
    # Upstream calls updateRootRadius() here; it lands with the hydraulics port.
    return t
end

function push_segment!(t::CCOTree; x, y, z, radius, kterm, parent,
                       resistance = 0.0, flow = 0.0, beta = 1.0)
    push!(t.x, x); push!(t.y, y); push!(t.z, z)
    push!(t.radius, radius); push!(t.kterm, kterm)
    push!(t.resistance, resistance); push!(t.flow, flow); push!(t.beta, beta)
    push!(t.parent, parent); push!(t.left, 0); push!(t.right, 0)
    return nsegments(t)
end

# HctccoModel::updateLengthFactor
function update_length_factor!(t::CCOTree)
    t.length_factor = sqrt(t.kterm_total + 1) * t.rsupp / t.rperf
    return t
end
