# SPDX-License-Identifier: GPL-3.0-or-later
# Closed coronary circuit on a grown CCO tree (arteries → capillaries → veins → coronary sinus),
# steady flow, and exact plug-flow iodine transport with per-terminal capillary exchange.
#
# Design (sources cited inline):
# - Total intramyocardial blood volume is f_b × tissue volume, split 29 / 31 / 40% into arteries,
#   capillaries and veins (Kassab, Lin & Fung 1994, pig cast). These targets hold for any tree
#   size: whatever the grown tree does not resolve is lumped into its terminals.
# - Every arterial segment j has a companion vein that carries the same flow back, lies along the
#   same centerline, and holds 40/29 of the artery's volume (radius × √(40/29)). Copying only the
#   topology (not the geometry, as Kaimovitz et al. 2010 did) keeps venous volume on target.
# - Every terminal feeds a lumped bed in series: unresolved arterioles (plug flow), a capillary
#   bed, unresolved venules (plug flow), and an interstitial space v_e (its Voronoi territory ×
#   0.18). The unresolved volumes are the Kassab targets minus what the tree resolves, shared by
#   territory volume. Capillary exchange follows the adiabatic approximation to the tissue
#   homogeneity model (St Lawrence & Lee 1998) with PS ∝ F (Caldwell et al. 1994), one global
#   factor set so first-pass extraction is E = 0.33 (Canty et al. 1991).
# - Iodine moves with plasma. A segment's plasma transit is τ = V(1 − H_T)/(Q(1 − H_D)).
#   Plug flow through a fixed network is a pure delay per segment, so every iodine mass below is
#   a closed form of the cumulative input A(t) = ∫ C_p. Nothing is time-stepped.
#
# Cost: one O(N) pass for topology and flow, O(N) per frame for arteries and tissue, and
# O(Σ path length) per frame for veins, with early exit where iodine has not yet arrived.
# Everything is Float32.

module Flow

import AcceleratedKernels as AK
using Atomix: @atomic

export VTree, vtree, save_vtree, load_vtree, merge_vtrees, Physiology, Circuit, circuit, InputCurve,
    gamma_input, Transport, transport, frame_masses, sinus_outflow

const F32 = Float32

# ---------------------------------------------------------------------------------------------
# Tree
# ---------------------------------------------------------------------------------------------

"""
    VTree

Arterial tree as flat arrays. Segment `j` runs from `a[:, j]` to `b[:, j]` (mm) with radius
`r[j]` (mm). `parent[j] == 0` marks a tree root.
"""
struct VTree
    parent::Vector{Int32}
    a::Matrix{F32}
    b::Matrix{F32}
    r::Vector{F32}
    isbase::Vector{Bool}
end

nseg(t::VTree) = length(t.r)

"""
    vtree(E) -> VTree

Converts a Vasculature `Engine` (index 1 is upstream's zero-length root sentinel).
"""
function vtree(E)
    m = E.m
    n = length(m.pos) - 1
    parent = Int32[m.parent[j] <= 1 ? 0 : m.parent[j] - 1 for j in 2:(n + 1)]
    a = Matrix{F32}(undef, 3, n); b = similar(a)
    for s in 1:n
        j = s + 1
        pa = m.pos[m.parent[j]]; pb = m.pos[j]
        for d in 1:3
            a[d, s] = F32(pa[d]); b[d, s] = F32(pb[d])
        end
    end
    return VTree(parent, a, b, F32.(m.radius[2:end]), Bool.(E.isbase[2:end]))
end

"""
    merge_vtrees(trees...) -> VTree

One forest from several trees (e.g. LCA and RCA). Each keeps its own root (`parent == 0`).
"""
function merge_vtrees(trees::VTree...)
    off = 0; parent = Int32[]
    for t in trees
        append!(parent, [p == 0 ? Int32(0) : Int32(p + off) for p in t.parent]); off += nseg(t)
    end
    return VTree(
        parent, reduce(hcat, [t.a for t in trees]), reduce(hcat, [t.b for t in trees]),
        reduce(vcat, [t.r for t in trees]), reduce(vcat, [t.isbase for t in trees])
    )
end

function save_vtree(path, t::VTree)
    open(path, "w") do io
        write(io, Int64(nseg(t)), t.parent, t.a, t.b, t.r, UInt8.(t.isbase))
    end
    return path
end

function load_vtree(path)
    return open(path) do io
        n = read(io, Int64)
        parent = read!(io, Vector{Int32}(undef, n))
        a = read!(io, Matrix{F32}(undef, 3, n)); b = read!(io, Matrix{F32}(undef, 3, n))
        r = read!(io, Vector{F32}(undef, n)); ib = read!(io, Vector{UInt8}(undef, n))
        VTree(parent, a, b, r, ib .== 0x01)
    end
end

# ---------------------------------------------------------------------------------------------
# Physiology
# ---------------------------------------------------------------------------------------------

"""
    Physiology(; kwargs...)

Boundary conditions and tissue constants. Flow is not prescribed: it follows from the
geometry and the pressures (see [`circuit`](@ref)).

- `p_ao`: mean aortic pressure at the ostia, 93 mmHg (Chilian et al. 1986, cat; normal human
  mean arterial pressure is similar).
- `p_back`: effective outflow pressure, 13.6 mmHg: the zero-flow pressure of the maximally
  dilated beating heart (Kanatsuka et al. 1990). It lumps coronary-sinus pressure and the
  waterfall effect of intramyocardial pressure into one steady back pressure.
- `demand`: resting myocardial blood flow the arterioles autoregulate to, mL/min/g. 1.0
  (human PET at rest, Chareonthaitawee et al. 2001). Only used for the rest state.
- `hd`: discharge hematocrit 0.45. Tube hematocrit from the Fåhræus law (Pries et al. 1992).
- `fb`: intramyocardial blood volume fraction, mL/mL. Default 0.10 (human rest 0.09, Larsson
  et al. 1996; pig 0.078–0.10, Jerosch-Herold 1999; Kassab's dilated cast 0.122).
- `share_art`, `share_cap`, `share_ven`: 0.29 / 0.31 / 0.40 of that volume (Kassab et al. 1994).
- `cap_ht_ratio`: tube/discharge hematocrit 0.85 in unresolved microvessels (Bindschadler 2014).
- `cap_length_mm`, `cap_diameter_um`: functional capillary 0.5 mm × 6 µm (Kassab & Fung 1994),
  used for the capillary-bed resistance.
- `ve`: interstitial volume fraction 0.18 (ECV ≈ 0.25 minus plasma).
- `extraction`: first-pass extraction E = 0.33 (Canty et al. 1991). PS ∝ F_p keeps E fixed per
  terminal, so `PS = −ln(1 − E)·F_p`.
- `rho`: myocardial density 1.05 g/mL.
- `radius_cv`: biological diameter variability. CCO sets every radius exactly from its flow;
  real vessels of one order scatter by 10–25% (Kassab et al. 1993, Table 2). Each segment's
  radius is multiplied by a fixed lognormal factor with this coefficient of variation.
- `d_art_micro`, `d_ven_micro`, `r_micro_dilated`: the microcirculation (arteries below 170 µm,
  capillaries, veins below 150 µm, Chilian et al. 1989's classes) has a measured resistance at
  maximal vasodilation of 4 mmHg·min·g/mL (Chilian et al. 1989, dipyridamole, Table 2 / Fig. 5).
  Grown microvessels keep their Poiseuille resistances relative to each other, so flow
  heterogeneity still comes from the geometry, but one factor κ sets their absolute level to the
  measurement: CCO's smallest arterioles are about 2.4× longer than Kassab's (CCO runs each
  terminal to the centre of its territory, where real tissue uses the capillary network), so
  their raw Poiseuille resistance is far too high. Vessels above the cut-offs are pure Poiseuille.
"""
Base.@kwdef struct Physiology
    p_ao::F32 = 93.0f0
    p_back::F32 = 13.6f0
    demand::F32 = 1.0f0
    hd::F32 = 0.45f0
    fb::F32 = 0.1f0
    share_art::F32 = 0.29f0
    share_cap::F32 = 0.31f0
    share_ven::F32 = 0.4f0
    cap_length_mm::F32 = 0.5f0
    cap_diameter_um::F32 = 6.0f0
    cap_ht_ratio::F32 = 0.85f0
    ve::F32 = 0.18f0
    extraction::F32 = 0.33f0
    rho::F32 = 1.05f0
    plasma_viscosity_mPas::F32 = 1.2f0
    radius_cv::F32 = 0.12f0
    d_art_micro::F32 = 170.0f0
    d_ven_micro::F32 = 150.0f0
    r_micro_dilated::F32 = 4.0f0
end

# Tube/discharge hematocrit ratio (Pries, Neuhaus & Gaehtgens 1992), D in µm.
@inline fahraeus(d::F32, hd::F32) = hd + (1.0f0 - hd) * (1.0f0 + 1.7f0 * exp(-0.415f0 * d) - 0.6f0 * exp(-0.011f0 * d))

# In vivo relative apparent viscosity (Pries & Secomb 2005, as given in Secomb & Pries 2013), D in µm.
@inline function mu_rel(d::F32, hd::F32)
    d = max(d, 2.0f0)
    m45 = 6.0f0 * exp(-0.085f0 * d) + 3.2f0 - 2.44f0 * exp(-0.06f0 * d^0.645f0)
    s = 1.0f0 / (1.0f0 + 1.0f-11 * d^12)
    c = (0.8f0 + exp(-0.075f0 * d)) * (s - 1.0f0) + s
    w = (d / (d - 1.1f0))^2
    return (1.0f0 + (m45 - 1.0f0) * ((1.0f0 - hd)^c - 1.0f0) / (0.55f0^c - 1.0f0) * w) * w
end

# Poiseuille resistance in mmHg·s/mL for length L (mm), radius r (mm), viscosity μ (mPa·s).
@inline poiseuille(L::F32, r::F32, mu::F32) = 8.0f0 * (mu * 1.0f-3) * (L * 1.0f-3) / (F32(π) * (r * 1.0f-3)^4) * 1.0f-6 / 133.322f0

# ---------------------------------------------------------------------------------------------
# Circuit: topology, volumes, pressure-driven steady flow
# ---------------------------------------------------------------------------------------------

struct Circuit
    tree::VTree
    veins::VTree                 # venous tree (companion veins: arterial topology, radius × s)
    phys::Physiology
    state::Symbol                # :rest (autoregulated) or :stress (maximal vasodilation)
    order::Vector{Int32}         # arterial breadth-first order: every parent before its children
    terminals::Vector{Int32}     # arterial segment index of each terminal
    tidx::Vector{Int32}          # terminal number of an arterial segment (0 if internal)
    vleaf::Vector{Int32}         # venous segment that drains terminal i
    len::Vector{F32}             # arterial segment length, mm
    q::Vector{F32}               # arterial blood flow, mL/s
    qv::Vector{F32}              # venous blood flow, mL/s
    pa::Vector{F32}              # pressure at the distal end of each arterial segment, mmHg
    pv::Vector{F32}              # pressure at the upstream (distal) end of each vein, mmHg
    tau_a::Vector{F32}           # arterial plasma transit, s
    tau_v::Vector{F32}           # venous plasma transit, s
    t_end::Vector{F32}           # plasma arrival at the distal end of arterial segment j, s
    s_v::Vector{F32}             # venous delay from the inlet of vein k to the sinus, s
    vp_a::Vector{F32}            # arterial plasma volume, mL
    vp_v::Vector{F32}            # venous plasma volume, mL
    # per terminal
    fp::Vector{F32}              # plasma flow, mL/s
    ta::Vector{F32}              # unresolved-arteriole plasma transit, s
    tc::Vector{F32}              # capillary plasma transit, s
    tv::Vector{F32}              # unresolved-venule plasma transit, s
    vint::Vector{F32}            # interstitial volume, mL
    k::Vector{F32}               # ISF→plasma back-flux rate E·F_p/v_e, 1/s
    offset::Vector{F32}          # arrival delay to the sinus, without the capillary and exchange
    tone::Vector{F32}            # arteriolar tone resistance added at rest, mmHg·s/mL (0 = fully dilated)
    rcap::Vector{F32}            # capillary-bed resistance (dilated), mmHg·s/mL
    demand::Vector{F32}          # flow each terminal autoregulates to at rest, mL/s
    # summary
    venous_scale::F32
    mass_g::F32
    volumes::NamedTuple
    resist::NamedTuple
end

"""
    circuit(tree, territory_ml; phys = Physiology(), state = :rest, veins = nothing, vleaf = nothing) -> Circuit

Steady coronary circuit driven by pressures: aortic pressure `phys.p_ao` at every ostium,
`phys.p_back` at the coronary-sinus ostium. Every terminal `i` (perfusing `territory_ml[i]`)
feeds a capillary bed whose resistance follows its capillary count (Kassab capillary
dimensions), then drains into the venous tree.

- `state = :stress`: maximal vasodilation. No arteriolar tone; every flow and pressure follows
  from geometry and the two boundary pressures.
- `state = :rest`: autoregulation. Each terminal's arterioles add just the tone needed to
  deliver the metabolic demand `phys.demand × mass`; tone cannot be negative, so a terminal
  that cannot reach its demand even fully dilated (for example behind a stenosis) stays
  dilated and gets less (Gould et al. 1974).

The venous tree is either grown (`veins`, with `vleaf[i]` the venous segment receiving terminal
`i`), or, by default, companion veins: each artery's vein has the same topology and holds
`share_ven/share_art` of its volume. Arteries and veins are solved exactly as trees (Norton
equivalents bottom-up, pressures top-down) and coupled by iteration on the venous pressures.
"""
function circuit(
        t::VTree, territory_ml::AbstractVector; phys::Physiology = Physiology(), state::Symbol = :rest,
        veins::Union{Nothing, VTree} = nothing, vleaf = nothing, demand = nothing, kappa = nothing
    )
    state in (:rest, :stress) || throw(ArgumentError("state must be :rest or :stress"))
    cv(t, salt) = phys.radius_cv <= 0 ? t : (
            sig = sqrt(log(1 + Float64(phys.radius_cv)^2));
            VTree(t.parent, t.a, t.b, [t.r[j] * F32(exp(sig * _hashnormal(j + salt) - sig^2 / 2)) for j in 1:nseg(t)], t.isbase)
        )
    t = cv(t, 0)
    veins === nothing || (veins = cv(veins, 1 << 40))
    return _circuit(t, territory_ml; phys, state, veins, vleaf, demand, kappa)
end

# Fixed standard normal per segment index (SplitMix64 + Box–Muller), so results are reproducible.
function _hashnormal(j::Integer)
    sm(x) = (x += 0x9e3779b97f4a7c15; x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9; x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb; x ⊻ (x >> 31))
    a = sm(UInt64(j)); b = sm(a)
    u1 = (Float64(a >> 11) + 0.5) / 2.0^53; u2 = Float64(b >> 11) / 2.0^53
    return sqrt(-2log(u1)) * cospi(2u2)
end

# children (CSR) and breadth-first order of a tree
function _topology(parent::Vector{Int32})
    n = length(parent)
    nchild = zeros(Int32, n)
    for j in 1:n
        p = parent[j]; p == 0 || (nchild[p] += 1)
    end
    ptr = zeros(Int32, n + 1); ptr[1] = 1
    for j in 1:n
        ptr[j + 1] = ptr[j] + nchild[j]
    end
    kids = Vector{Int32}(undef, n); fill_ = copy(ptr)
    for j in 1:n
        p = parent[j]; p == 0 && continue
        kids[fill_[p]] = j; fill_[p] += 1
    end
    order = Int32[j for j in 1:n if parent[j] == 0]
    h = 1
    while h <= length(order)
        j = order[h]; h += 1
        for c in ptr[j]:(ptr[j + 1] - 1)
            push!(order, kids[c])
        end
    end
    return nchild, ptr, kids, order
end

function _lengths(t::VTree)
    n = nseg(t); len = Vector{F32}(undef, n); vol = similar(len)
    AK.foreachindex(1:n; min_elems = 4096) do j
        L = 0.0f0
        for d in 1:3
            L += (t.b[d, j] - t.a[d, j])^2
        end
        L = max(sqrt(L), 1.0f-3)
        len[j] = L; vol[j] = F32(π) * t.r[j]^2 * L * 1.0f-3
    end
    return len, vol
end

# Arterial tree + terminal beds + venous tree, solved for flows and pressures.
struct _Solver
    parent::Vector{Int32}; order::Vector{Int32}; ptr::Vector{Int32}; kids::Vector{Int32}
    terminals::Vector{Int32}; tidx::Vector{Int32}
    vparent::Vector{Int32}; vorder::Vector{Int32}; vleaf::Vector{Int32}
    ra::Vector{F32}; rv::Vector{F32}; rcap::Vector{F32}; ma::Vector{Bool}; mv::Vector{Bool}
    Pa::F32; Pb::F32
    G::Vector{F32}; E::Vector{F32}; I::Vector{F32}
    q::Vector{F32}; pa::Vector{F32}; qt::Vector{F32}; qv::Vector{F32}; pv::Vector{F32}
    pvt::Vector{F32}; tone::Vector{F32}; src::Vector{Bool}
    vE::Vector{F32}; vI::Vector{F32}; iters::Base.RefValue{Int}
end

function _Solver(parent, order, ptr, kids, terminals, tidx, vparent, vorder, vleaf, ra, rv, rcap, ma, mv, Pa, Pb)
    n = length(parent); nv = length(vparent); nt = length(terminals)
    z(k) = zeros(F32, k)
    return _Solver(
        parent, order, ptr, kids, terminals, tidx, vparent, vorder, vleaf, ra, rv, rcap, ma, mv, F32(Pa), F32(Pb),
        z(n), z(n), z(n), z(n), z(n), z(nt), z(nv), z(nv), fill(F32(Pb), nt), z(nt), fill(false, nt), z(nv), z(nv), Ref(0)
    )
end

# One steady solve. Microvessels (arteries `ma`, veins `mv`, capillary beds) are scaled by κ;
# `tone` is each terminal's added resistance. With `dem`, terminals autoregulate to those flows:
# a terminal is a flow source while its required tone is ≥ 0, else a fully dilated resistor.
#
# Unknowns are the terminal flows q. Arterial pressure drop to terminal i plus its own bed is
# (M_a q)_i (shared path resistances + bed), the venous rise from the sinus is (M_v q)_i, and
# continuity at every bed is (M_a + M_v) q = P_ao − P_back. Both matrices are symmetric positive
# definite and applied in O(N) by two tree passes, so the system is solved by conjugate gradients
# preconditioned with the exact arterial tree solve (M_a⁻¹ is the Norton solve of the arterial
# tree); its condition number is 1 + λ(M_a⁻¹M_v), a few, so it converges in tens of iterations.
function _solve!(S::_Solver, kappa::F32, tone::Vector{F32}, dem; ka::F32 = kappa)
    nt = length(S.terminals)
    S.tone .= tone
    if dem === nothing
        fill!(S.src, false)
    else
        fill!(S.src, true); fill!(S.tone, 0.0f0)
    end
    for outer in 1:50
        _pcg!(S, kappa, dem, ka)
        dem === nothing && break
        changed = false
        @inbounds for i in 1:nt
            S.src[i] || continue
            tn = (S.pa[S.terminals[i]] - S.pvt[i]) / dem[i] - ka * S.rcap[i]
            if tn < 0
                S.src[i] = false; S.tone[i] = 0.0f0; changed = true
            else
                S.tone[i] = tn
            end
        end
        changed || break
    end
    return S
end

# Terminal-path operator: out_i = Σ_j∈path(i) r_j Q_j (+ bed resistance × q_i for `bed`), where Q_j
# is the sum of q over the terminals below segment j. For the arteries the path runs root → i,
# for the veins from the sinus to the vein draining i.
function _apply_art!(out, S::_Solver, ka::F32, q, bed::Bool)
    n = length(S.parent); Q = S.I; D = S.E; ptr = S.ptr; kids = S.kids
    @inbounds for h in n:-1:1
        j = S.order[h]
        if ptr[j + 1] == ptr[j]
            Q[j] = q[S.tidx[j]]
        else
            s = 0.0f0
            for c in ptr[j]:(ptr[j + 1] - 1)
                s += Q[kids[c]]
            end
            Q[j] = s
        end
    end
    @inbounds for h in 1:n
        j = S.order[h]; p = S.parent[j]
        D[j] = (p == 0 ? 0.0f0 : D[p]) + Q[j] * (S.ma[j] ? ka * S.ra[j] : S.ra[j])
        if ptr[j + 1] == ptr[j]
            i = S.tidx[j]
            out[i] = D[j] + (bed ? (ka * S.rcap[i] + S.tone[i]) * q[i] : 0.0f0)
        end
    end
    return out
end

function _apply_ven!(out, S::_Solver, kappa::F32, q)
    nv = length(S.vparent); nt = length(S.terminals); Q = S.vI; D = S.vE
    fill!(Q, 0.0f0)
    @inbounds for i in 1:nt
        Q[S.vleaf[i]] += q[i]
    end
    @inbounds for h in nv:-1:1
        k = S.vorder[h]; p = S.vparent[k]
        p == 0 || (Q[p] += Q[k])
    end
    @inbounds for h in 1:nv
        k = S.vorder[h]; p = S.vparent[k]
        D[k] = (p == 0 ? 0.0f0 : D[p]) + Q[k] * (S.mv[k] ? kappa * S.rv[k] : S.rv[k])
    end
    @inbounds for i in 1:nt
        out[i] = D[S.vleaf[i]]
    end
    return out
end

# Preconditioner z = (M_a)_FF⁻¹ r over the free (non-source) terminals: the arterial tree with its
# root at 0 and free terminal i draining through its bed into pressure −r_i; sources are open.
function _precond!(z, S::_Solver, ka::F32, r)
    n = length(S.parent); G = S.G; E = S.E; ptr = S.ptr; kids = S.kids
    @inbounds for h in n:-1:1
        j = S.order[h]
        if ptr[j + 1] == ptr[j]
            i = S.tidx[j]
            if S.src[i]
                g = 0.0f0; e = 0.0f0
            else
                g = 1.0f0 / (ka * S.rcap[i] + S.tone[i]); e = -r[i]
            end
        else
            g = 0.0f0; ge = 0.0f0
            for c in ptr[j]:(ptr[j + 1] - 1)
                k = kids[c]; g += G[k]; ge += G[k] * E[k]
            end
            e = g > 0 ? ge / g : 0.0f0
        end
        rr = S.ma[j] ? ka * S.ra[j] : S.ra[j]
        G[j] = g / (1.0f0 + g * rr); E[j] = e
    end
    P = S.pa
    @inbounds for h in 1:n
        j = S.order[h]; p = S.parent[j]
        pin = p == 0 ? 0.0f0 : P[p]
        qj = G[j] * (pin - E[j])
        P[j] = pin - qj * (S.ma[j] ? ka * S.ra[j] : S.ra[j])
        if ptr[j + 1] == ptr[j]
            i = S.tidx[j]; z[i] = S.src[i] ? 0.0f0 : qj
        end
    end
    return z
end

function _pcg!(S::_Solver, kappa::F32, dem, ka::F32 = kappa)
    nt = length(S.terminals); dP = S.Pa - S.Pb
    q = S.qt; free = .!S.src
    # sources fixed at their demand; free terminals start from the solution with P_v = P_back
    @inbounds for i in 1:nt
        S.src[i] && (q[i] = dem[i])
    end
    t1 = zeros(F32, nt); t2 = zeros(F32, nt)
    Aop!(out, x) = (_apply_art!(t1, S, ka, x, true); _apply_ven!(t2, S, kappa, x); out .= t1 .+ t2; out)
    if any(free)
        qs = copy(q); qs[free] .= 0.0f0
        b = fill(dP, nt) .- Aop!(zeros(F32, nt), qs)       # right-hand side for the free set
        x = zeros(F32, nt); r = similar(x); z = similar(x); p = similar(x); Ap = similar(x)
        _precond!(x, S, ka, b)                              # initial guess: arteries alone
        Aop!(Ap, x); @. r = ifelse(free, b - Ap, 0.0f0)
        _precond!(z, S, ka, r); p .= z
        rz = sum(r .* z); bn = sqrt(sum(abs2, b[free]))
        for it in 1:500
            S.iters[] += 1
            Aop!(Ap, p); Ap[.!free] .= 0.0f0
            α = rz / sum(p .* Ap)
            @. x += α * p; @. r -= α * Ap
            sqrt(sum(abs2, r)) < 1.0f-6 * bn && break
            _precond!(z, S, ka, r)
            rz2 = sum(r .* z); β = rz2 / rz; rz = rz2
            @. p = z + β * p
        end
        @inbounds for i in 1:nt
            free[i] && (q[i] = x[i])
        end
    end
    # flows and pressures from the terminal flows
    _apply_art!(t1, S, ka, q, false)
    n = length(S.parent)
    @inbounds for j in 1:n
        S.q[j] = S.I[j]; S.pa[j] = S.Pa - S.E[j]
    end
    _apply_ven!(t2, S, kappa, q)
    @inbounds for k in eachindex(S.vparent)
        S.qv[k] = S.vI[k]; S.pv[k] = S.Pb + S.vE[k]
    end
    @inbounds for i in 1:nt
        S.pvt[i] = S.pv[S.vleaf[i]]
    end
    return S
end

# Microvascular resistance as Chilian et al. (1989) measured it: flow-weighted mean pressure where
# arterioles leave the larger arteries minus that where venules join the larger veins, divided by
# the flow per gram (mmHg·min·g/mL).
function _micro_resistance(S::_Solver, mass_g)
    sa = 0.0f0; wa = 0.0f0; sv = 0.0f0; wv = 0.0f0
    @inbounds for j in eachindex(S.parent)
        p = S.parent[j]
        if S.ma[j] && (p == 0 || !S.ma[p])
            pin = p == 0 ? S.Pa : S.pa[p]; sa += S.q[j] * pin; wa += S.q[j]
        end
    end
    @inbounds for k in eachindex(S.vparent)
        p = S.vparent[k]
        if S.mv[k] && (p == 0 || !S.mv[p])
            pout = p == 0 ? S.Pb : S.pv[p]; sv += S.qv[k] * pout; wv += S.qv[k]
        end
    end
    qtot = sum(S.qt)
    return (sa / wa - sv / wv) / (60.0f0 * qtot / mass_g), sa / wa, sv / wv
end

function _circuit(t::VTree, territory_ml::AbstractVector; phys::Physiology, state::Symbol, veins, vleaf, demand, kappa)
    n = nseg(t)
    nchild, ptr, kids, order = _topology(t.parent)
    terminals = Int32[j for j in 1:n if nchild[j] == 0]
    nt = length(terminals)
    length(territory_ml) == nt || throw(ArgumentError("territory_ml needs one entry per terminal ($nt)"))
    tidx = zeros(Int32, n)
    for (i, j) in enumerate(terminals)
        tidx[j] = i
    end
    vol_tissue = F32(sum(territory_ml)); mass_g = vol_tissue * phys.rho; hd = phys.hd; mup = phys.plasma_viscosity_mPas
    len, vol = _lengths(t); vart = sum(vol)
    vblood = phys.fb * vol_tissue
    # venous tree (companion veins share the arterial topology)
    s = sqrt(phys.share_ven / phys.share_art)
    if veins === nothing
        vt = VTree(t.parent, t.a, t.b, s .* t.r, t.isbase); vl = terminals
        vorder = order; vlen = len; vvol = (s * s) .* vol
    else
        length(vleaf) == nt || throw(ArgumentError("vleaf needs one venous segment per terminal"))
        vl = Int32.(vleaf)
        # A grown venous tree is sized by CCO's arterial law; scale its radii uniformly so it holds
        # share_ven/share_art of the resolved arterial volume (Kassab et al. 1994), as companion
        # veins do. Veins are wider and lower-resistance than their arteries.
        vlen, vvol = _lengths(veins)
        fr = sqrt(s * s * vart / sum(vvol))
        vt = VTree(veins.parent, veins.a, veins.b, veins.r .* F32(fr), veins.isbase)
        vorder = _topology(vt.parent)[4]; vlen, vvol = _lengths(vt)
    end
    nv = nseg(vt); vven = sum(vvol)
    va_unres = max(0.0f0, phys.share_art * vblood - vart)
    vv_unres = max(0.0f0, phys.share_ven * vblood - vven)
    vcap_total = phys.share_cap * vblood
    ra = similar(len); rv = Vector{F32}(undef, nv)
    AK.foreachindex(1:n; min_elems = 4096) do j
        ra[j] = poiseuille(len[j], t.r[j], mup * mu_rel(2000.0f0 * t.r[j], hd))
    end
    AK.foreachindex(1:nv; min_elems = 4096) do k
        rv[k] = poiseuille(vlen[k], vt.r[k], mup * mu_rel(2000.0f0 * vt.r[k], hd))
    end
    # capillary bed of terminal i: its share of the capillary volume as parallel functional capillaries
    rc = phys.cap_diameter_um * 0.5f-3
    v1 = F32(π) * rc^2 * phys.cap_length_mm * 1.0f-3
    r1 = poiseuille(phys.cap_length_mm, rc, mup * mu_rel(phys.cap_diameter_um, hd))
    w = F32.(territory_ml) ./ vol_tissue
    rcap = [r1 / (vcap_total * w[i] / v1) for i in 1:nt]
    # microcirculation: arteries < d_art_micro, veins < d_ven_micro and the capillary beds
    ma = [2000.0f0 * t.r[j] < phys.d_art_micro for j in 1:n]
    mv = [2000.0f0 * vt.r[k] < phys.d_ven_micro for k in 1:nv]
    S = _Solver(t.parent, order, ptr, kids, terminals, tidx, vt.parent, vorder, vl, ra, rv, rcap, ma, mv, phys.p_ao, phys.p_back)

    # ---- microvascular scale κ (maximal dilation): Chilian-style resistance between the inlets
    # of the arterioles and the outlets of the venules equals the measured value
    if kappa === nothing
        kappa = 1.0f0
        for _ in 1:4
            _solve!(S, kappa, zeros(F32, nt), nothing)
            rm = _micro_resistance(S, mass_g)[1]
            kappa *= phys.r_micro_dilated / rm
        end
    end
    kappa = F32(kappa)
    # ---- the state
    tone = zeros(F32, nt); tau = 0.0f0
    if state === :stress
        _solve!(S, kappa, tone, nothing)
        dem = copy(S.qt)
    elseif demand === nothing
        # healthy rest: resting tone constricts every arteriole (arteries below d_art_micro) and
        # bed by the same factor (1 + τ) until total flow meets the metabolic demand, so the
        # network keeps its own flow heterogeneity; each terminal's flow is then its demand
        target = phys.demand * mass_g / 60.0f0
        f(τ) = (_solve!(S, kappa, tone, nothing; ka = kappa * (1 + τ)); sum(S.qt) - target)
        lo, hi = 0.0f0, 1.0f0
        while f(hi) > 0
            lo = hi; hi *= 2
        end
        for _ in 1:40
            mid = (lo + hi) / 2
            f(mid) > 0 ? (lo = mid) : (hi = mid)
            (hi - lo) < 1.0f-4 * hi && break
        end
        tau = (lo + hi) / 2; f(tau)
        dem = copy(S.qt)
    else
        # given demand (e.g. the healthy heart's, with a stenosis added): each terminal adjusts its
        # own tone to keep its demand; tone cannot go below zero (Gould et al. 1974)
        dem = F32.(demand)
        _solve!(S, kappa, tone, dem)
    end
    q = S.q; pa = S.pa; qt = S.qt; qv = S.qv; pv = S.pv
    nneg = count(<=(0), qt)
    nneg > 0 && @warn "$nneg terminals with reversed or zero flow (venous pressure above arterial); their beds are treated as stagnant"
    @inbounds for i in 1:nt
        qt[i] = max(qt[i], 1.0f-9)
    end
    tone = S.tone

    # ---- plasma volumes and transit times
    vp_a = similar(len); tau_a = similar(len)
    AK.foreachindex(1:n; min_elems = 4096) do j
        vp_a[j] = vol[j] * (1.0f0 - fahraeus(2000.0f0 * t.r[j], hd) * hd)
        tau_a[j] = vp_a[j] / (q[j] * (1.0f0 - hd))
    end
    vp_v = Vector{F32}(undef, nv); tau_v = similar(vp_v)
    AK.foreachindex(1:nv; min_elems = 4096) do k
        vp_v[k] = vvol[k] * (1.0f0 - fahraeus(2000.0f0 * vt.r[k], hd) * hd)
        tau_v[k] = vp_v[k] / max(qv[k] * (1.0f0 - hd), 1.0f-12)
    end
    t_end = similar(len)
    for h in 1:n
        j = order[h]; p = t.parent[j]
        t_end[j] = (p == 0 ? 0.0f0 : t_end[p]) + tau_a[j]
    end
    s_v = similar(vp_v)
    for h in 1:nv
        k = vorder[h]; p = vt.parent[k]
        s_v[k] = (p == 0 ? 0.0f0 : s_v[p]) + tau_v[k]
    end

    # ---- terminals: capillary bed and interstitium
    fp = Vector{F32}(undef, nt); ta = similar(fp); tc = similar(fp); tv = similar(fp)
    vint = similar(fp); kk = similar(fp); off = similar(fp)
    pf = 1.0f0 - phys.cap_ht_ratio * hd      # plasma fraction of microvascular blood
    Ex = phys.extraction
    AK.foreachindex(1:nt; min_elems = 4096) do i
        j = terminals[i]
        fp[i] = qt[i] * (1.0f0 - hd)
        ta[i] = w[i] * va_unres * pf / fp[i]
        tc[i] = w[i] * vcap_total * pf / fp[i]
        tv[i] = w[i] * vv_unres * pf / fp[i]
        vint[i] = phys.ve * F32(territory_ml[i])
        kk[i] = Ex * fp[i] / vint[i]
        off[i] = t_end[j] + ta[i] + tv[i] + s_v[vl[i]]
    end

    # ---- pressure-drop shares (flow-weighted) and summary flows
    ka = kappa * (1 + tau)          # arterioles and beds carry the resting tone
    sa = [ma[j] ? ka : 1.0f0 for j in 1:n]; sv = [mv[k] ? kappa : 1.0f0 for k in 1:nv]
    dpa = sum(q .* ra .* sa); dpv = sum(qv .* rv .* sv); dpc = sum(qt .* rcap) * ka; dpt = sum(qt .* tone)
    roots = Int32[j for j in 1:n if t.parent[j] == 0]
    qtot = sum(q[roots])
    resist = (;
        dp_total = phys.p_ao - phys.p_back, art_share = dpa / (dpa + dpv + dpc + dpt), tone_share = dpt / (dpa + dpv + dpc + dpt),
        cap_share = dpc / (dpa + dpv + dpc + dpt), ven_share = dpv / (dpa + dpv + dpc + dpt),
        mbf = 60 * qtot / mass_g, root_flow = q[roots], n_dilated = count(i -> tone[i] == 0, 1:nt), kappa, tau,
        r_micro = (rmc = _micro_resistance(S, mass_g); rmc[1]), sweeps = S.iters[],
        chilian_shares = ((phys.p_ao - rmc[2]), (rmc[2] - rmc[3]), (rmc[3] - phys.p_back)) ./ (phys.p_ao - phys.p_back),
    )
    return Circuit(
        t, vt, phys, state, order, terminals, tidx, vl, len, q, qv, pa, pv, tau_a, tau_v, t_end, s_v, vp_a, vp_v,
        fp, ta, tc, tv, vint, kk, off, tone, rcap .* kappa, dem, s, mass_g,
        (;
            blood = vblood, art_resolved = vart, ven_resolved = vven, art_unresolved = va_unres,
            ven_unresolved = vv_unres, capillary = vcap_total,
        ), resist
    )
end

# ---------------------------------------------------------------------------------------------
# Input and exchange tables
# ---------------------------------------------------------------------------------------------

"""
    InputCurve

Plasma iodine concentration at the coronary ostia, sampled every `dt` seconds from t = 0.
Stores the cumulative integral A(t) and, on a log-spaced grid of back-flux rates k, the
exchange integrals X_k(u) = ∫₀ᵘ C(s) e^{−k(u−s)} ds and Y_k(u) = ∫₀ᵘ X_k.
"""
struct InputCurve
    dt::F32
    c::Vector{F32}
    A::Vector{F32}
    logk::Vector{F32}
    X::Matrix{F32}     # (nt, nk)
    Y::Matrix{F32}
end

function InputCurve(c::Vector{F32}, dt::Real; kmin = 1.0f-4, kmax = 10.0f0, nk = 160)
    dt = F32(dt); nt = length(c)
    A = zeros(F32, nt)
    for i in 2:nt
        A[i] = A[i - 1] + 0.5f0 * dt * (c[i - 1] + c[i])
    end
    logk = collect(F32, range(log(F32(kmin)), log(F32(kmax)), length = nk))
    X = zeros(F32, nt, nk); Y = zeros(F32, nt, nk)
    AK.foreachindex(1:nk; min_elems = 1) do m
        k = exp(logk[m]); e = exp(-k * dt)
        # exact for piecewise-linear C
        a = (1.0f0 - e) / k; b = (dt - a) / (k * dt)   # weights of C[i] and C[i]-C[i-1] slope part
        for i in 2:nt
            ci0 = c[i - 1]; ci1 = c[i]
            X[i, m] = e * X[i - 1, m] + ci1 * a - (ci1 - ci0) * b
            Y[i, m] = Y[i - 1, m] + 0.5f0 * dt * (X[i - 1, m] + X[i, m])
        end
    end
    return InputCurve(dt, c, A, logk, X, Y)
end

"""
    gamma_input(; t0, tpeak, peak_blood, alpha, hd, dt, tmax) -> InputCurve

Gamma-variate whole-blood iodine curve (mg/mL) arriving at `t0`, peaking at `tpeak`, converted
to plasma concentration C_p = C_b/(1 − H_D). Times in seconds from the start of the simulation.
"""
function gamma_input(;
        t0 = 2.0f0, tpeak = 10.0f0, peak_blood = 12.0f0, alpha = 3.0f0, hd = 0.45f0,
        dt = 0.01f0, tmax = 40.0f0
    )
    ts = collect(F32, 0:dt:tmax)
    b = (tpeak - t0) / alpha
    c = map(ts) do t
        t <= t0 ? 0.0f0 : peak_blood / (1.0f0 - hd) * ((t - t0) / (tpeak - t0))^alpha * exp(alpha - (t - t0) / b)
    end
    return InputCurve(c, dt)
end

@inline function _lin(v::AbstractVector{F32}, u::F32, dt::F32)
    u <= 0.0f0 && return 0.0f0
    x = u / dt
    i = unsafe_trunc(Int, x)
    n = length(v)
    i >= n - 1 && return v[n]
    f = x - F32(i)
    return @inbounds v[i + 1] + f * (v[i + 2] - v[i + 1])
end

@inline function _bilin(M::AbstractMatrix{F32}, logk::AbstractVector{F32}, k::F32, u::F32, dt::F32)
    u <= 0.0f0 && return 0.0f0
    nk = length(logk)
    y = (log(k) - logk[1]) / (logk[2] - logk[1])
    y = clamp(y, 0.0f0, F32(nk - 1) - 1.0f-3)
    m = unsafe_trunc(Int, y); g = y - F32(m)
    x = u / dt; nt = size(M, 1)
    i = unsafe_trunc(Int, x)
    if i >= nt - 1
        return @inbounds (1.0f0 - g) * M[nt, m + 1] + g * M[nt, m + 2]
    end
    f = x - F32(i)
    @inbounds begin
        lo = M[i + 1, m + 1] + f * (M[i + 2, m + 1] - M[i + 1, m + 1])
        hi = M[i + 1, m + 2] + f * (M[i + 2, m + 2] - M[i + 1, m + 2])
    end
    return (1.0f0 - g) * lo + g * hi
end

# Cumulative plasma iodine that has left terminal i's capillary bed by local time u (mg):
# G(u) = F_p [(1 − E) A(u − T_c) + E k Y_k(u − T_c)]   (AATH, St Lawrence & Lee 1998)
@inline function _gout(A, Y, logk, dt, fp, tc, k, E, u)
    w = u - tc
    w <= 0.0f0 && return 0.0f0
    return fp * ((1.0f0 - E) * _lin(A, w, dt) + E * k * _bilin(Y, logk, k, w, dt))
end

# Iodine in terminal i's capillary bed + interstitium at local time u (mg):
# q(u) = F_p [A(u) − A(u − T_c)] + F_p E X_k(u − T_c)
@inline function _gcap(A, X, logk, dt, fp, tc, k, E, u)
    u <= 0.0f0 && return 0.0f0
    return fp * (_lin(A, u, dt) - _lin(A, u - tc, dt) + E * _bilin(X, logk, k, u - tc, dt))
end

# Whole lumped bed of terminal i at local time u (from the terminal segment's distal end):
# unresolved arterioles (plug, T_a) → capillary + interstitium → unresolved venules (plug, T_v).
@inline function _gtissue(A, X, Y, logk, dt, fp, ta, tc, tv, k, E, u)
    u <= 0.0f0 && return 0.0f0
    art = fp * (_lin(A, u, dt) - _lin(A, u - ta, dt))
    ven = _gout(A, Y, logk, dt, fp, tc, k, E, u - ta) - _gout(A, Y, logk, dt, fp, tc, k, E, u - ta - tv)
    return art + _gcap(A, X, logk, dt, fp, tc, k, E, u - ta) + ven
end

# ---------------------------------------------------------------------------------------------
# Iodine masses at chosen times
# ---------------------------------------------------------------------------------------------

"""
    frame_masses(c::Circuit, inp::InputCurve, times; to = identity) -> (artery, vein, tissue)

Iodine mass (mg) in every arterial segment, every vein of `c.veins` and every terminal's lumped bed
(unresolved arterioles, capillaries, interstitium, unresolved venules) at each time in `times`.
Matrices are (ntimes, n) so one segment's frames are contiguous.

`to` moves arrays to a device (for example `to = MtlArray` or `CuArray`); the kernels then run
there and the results stay on that device. Each mass is a closed form of the cumulative input,
so the cost per frame is O(segments) for arteries and terminals and O(Σ path length) for veins,
with early exit where iodine has not arrived yet.
"""
function frame_masses(c::Circuit, inp::InputCurve, times::AbstractVector; to = identity)
    tm = to(F32.(collect(times))); nf = length(times)
    n = nseg(c.tree); nt = length(c.terminals)
    E = c.phys.extraction; hd = c.phys.hd; dt = inp.dt
    parent = to(c.tree.parent); terminals = to(c.terminals); vparent = to(c.veins.parent); vleaf = to(c.vleaf)
    t_end = to(c.t_end); s_v = to(c.s_v); q = to(c.q); nv = nseg(c.veins)
    fp = to(c.fp); ta = to(c.ta); tc = to(c.tc); tv = to(c.tv); kk = to(c.k); off = to(c.offset)
    A = to(inp.A); X = to(inp.X); Y = to(inp.Y); logk = to(inp.logk)
    art = to(zeros(F32, nf, n)); ven = to(zeros(F32, nf, nv)); tis = to(zeros(F32, nf, nt))
    # arteries: plasma inside segment j entered the root between t − t_end[j] and t − t_beg[j]
    AK.foreachindex(q; min_elems = 4096) do j
        p = parent[j]
        tb = p == 0 ? 0.0f0 : t_end[p]; te = t_end[j]
        qp = q[j] * (1.0f0 - hd)
        for f in 1:nf
            art[f, j] = qp * (_lin(A, tm[f] - tb, dt) - _lin(A, tm[f] - te, dt))
        end
    end
    # lumped bed of every terminal
    AK.foreachindex(fp; min_elems = 1024) do i
        j = terminals[i]
        for f in 1:nf
            tis[f, i] = _gtissue(A, X, Y, logk, dt, fp[i], ta[i], tc[i], tv[i], kk[i], E, tm[f] - t_end[j])
        end
    end
    # veins: walk each terminal's path to the sinus; the mass of its outflow inside vein j is
    # G(t − O_i + S[j]) − G(t − O_i + S[parent j]). Stop where its iodine has not yet arrived.
    AK.foreachindex(fp; min_elems = 256) do i
        j0 = vleaf[i]
        f_ = fp[i]; tcap = tc[i]; k = kk[i]; base = off[i]
        for f in 1:nf
            j = j0
            ghi = _gout(A, Y, logk, dt, f_, tcap, k, E, tm[f] - base + s_v[j])
            while ghi > 0.0f0
                p = vparent[j]
                glo = _gout(A, Y, logk, dt, f_, tcap, k, E, tm[f] - base + (p == 0 ? 0.0f0 : s_v[p]))
                @atomic ven[f, j] += ghi - glo
                p == 0 && break
                ghi = glo; j = p
            end
        end
    end
    return art, ven, tis
end

"""
    sinus_outflow(c, inp, times; to = identity) -> Vector

Cumulative iodine (mg) that has left through the coronary sinus by each time.
"""
function sinus_outflow(c::Circuit, inp::InputCurve, times::AbstractVector; to = identity)
    tm = to(F32.(collect(times))); nf = length(times); nt = length(c.terminals)
    E = c.phys.extraction; dt = inp.dt
    fp = to(c.fp); tc = to(c.tc); kk = to(c.k); off = to(c.offset)
    A = to(inp.A); Y = to(inp.Y); logk = to(inp.logk)
    # terminals in chunks, so memory stays (frames × chunk) however many terminals there are
    chunk = max(1, min(nt, 1 << 24 ÷ max(nf, 1)))
    out = to(zeros(F32, nf, chunk)); total = zeros(F32, nf)
    for c0 in 1:chunk:nt
        m = min(chunk, nt - c0 + 1)
        fill!(out, 0.0f0)
        AK.foreachindex(1:m; min_elems = 1024) do ii
            i = c0 + ii - 1
            for f in 1:nf
                out[f, ii] = _gout(A, Y, logk, dt, fp[i], tc[i], kk[i], E, tm[f] - off[i])
            end
        end
        total .+= vec(Array(sum(out; dims = 2)))
    end
    return total
end

end # module
