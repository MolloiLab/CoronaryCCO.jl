# SPDX-License-Identifier: GPL-3.0-or-later
# Closed coronary circuit on a grown CCO tree (arteries → capillaries → veins → coronary sinus),
# steady flow, and exact plug-flow iodine transport with per-terminal capillary exchange.
#
# Design (sources cited inline):
# - Total intramyocardial blood volume is f_b × tissue volume, split 29 / 31 / 40% into arteries,
#   capillaries and veins (Kassab, Lin & Fung 1994, pig cast). These targets hold for any tree
#   size: whatever the grown tree does not resolve is lumped into its terminals.
# - Flow is driven by two pressures (aortic, coronary back pressure) through arteries, capillary
#   beds and veins. Arteries below 150 µm are walls in force balance whose smooth-muscle tone
#   follows wall tension (myogenic) and the O₂ state of the beds they feed (ATP from desaturating
#   red cells, conducted upstream), so resting flow, its distribution, the pressure–flow curve and
#   the response to a stenosis all follow from the physics (see `Physiology`, `circuit`).
# - Veins are a grown venous tree (each bed drains into a venule), or, by default, companion veins:
#   every arterial segment j has a vein that carries the same flow back, lies along the same
#   centerline, and holds 40/29 of the artery's volume (radius × √(40/29)). Copying only the
#   topology (not the geometry, as Kaimovitz et al. 2010 did) keeps venous volume on target.
# - Every terminal feeds a lumped bed in series: unresolved arterioles (plug flow), a capillary
#   bed, unresolved venules (plug flow), and an interstitial space v_e (its Voronoi territory ×
#   0.18). The unresolved volumes are the Kassab targets minus what the tree resolves, shared by
#   territory volume. Capillary exchange follows the adiabatic approximation to the tissue
#   homogeneity model (St Lawrence & Lee 1998), with first-pass extraction from a constant
#   permeability (Renkin–Crone, E = 1 − e^{−PS/F}; PS from Daghini et al. 2007).
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

export VTree, vtree, save_vtree, load_vtree, merge_vtrees, Physiology, Circuit, circuit, with_extraction, InputCurve,
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

Boundary conditions, tissue constants and vessel-wall mechanics. Flow is not prescribed: it
follows from the geometry, the two pressures and the arteriolar diameters (see [`circuit`](@ref)).

Boundary conditions and blood
- `p_ao`: mean aortic pressure at the ostia, 93 mmHg (Chilian et al. 1986, cat; normal human
  mean arterial pressure is similar).
- `p_back`: effective outflow pressure, 13.6 mmHg: the zero-flow pressure of the maximally
  dilated beating heart (Kanatsuka et al. 1990). It lumps coronary-sinus pressure and the
  waterfall effect of intramyocardial pressure into one steady back pressure.
- `hd`: discharge hematocrit 0.45. Tube hematocrit from the Fåhræus law (Pries et al. 1992).

Volumes and exchange
- `fb`: intramyocardial blood volume fraction, mL/mL. Default 0.10 (human rest 0.09, Larsson
  et al. 1996; pig 0.078–0.10, Jerosch-Herold 1999; Kassab's dilated cast 0.122).
- `share_art`, `share_cap`, `share_ven`: 0.29 / 0.31 / 0.40 of that volume (Kassab et al. 1994).
- `cap_ht_ratio`: tube/discharge hematocrit 0.85 in unresolved microvessels (Bindschadler 2014).
- `cap_length_mm`, `cap_diameter_um`: functional capillary 0.5 mm × 6 µm (Kassab & Fung 1994),
  used for the capillary-bed resistance.
- `ve`: interstitial volume fraction 0.18 (ECV ≈ 0.25 minus plasma).
- `ps`: capillary permeability–surface product for iodinated contrast, 0.37 mL/min per g: the
  Patlak influx constant of iopamidol in normal pig myocardium, the same by 64-slice CT and by
  electron-beam CT, and "similar" during adenosine (Daghini et al. 2007). Each bed's first-pass
  extraction follows Renkin–Crone, E = 1 − exp(−PS·m/F_p), so E ≈ 0.5 at rest and ≈ 0.12 at
  maximal dilation. A Patlak constant is E·F_p, a lower bound on PS. Other measurements spread
  around it: 0.13–0.19 for iohexol under coronary adenosine (Ritman 1994), and E = 0.33 at full
  vasodilation in dogs (Canty et al. 1991, whose per-gram flow is unknown). These are checks.
- `rho`: myocardial density 1.05 g/mL.
- `radius_cv`: biological diameter variability. CCO sets every radius exactly from its flow;
  real vessels of one order scatter by 10–25% (Kassab et al. 1993, Table 2). Each segment's
  radius is multiplied by a fixed lognormal factor with this coefficient of variation.

Microvascular resistance at maximal dilation
- `d_art_micro`, `d_ven_micro`, `r_micro_dilated`: the microcirculation (arteries below 170 µm,
  capillaries, veins below 150 µm, Chilian et al. 1989's classes) has a measured resistance at
  maximal vasodilation of 4 mmHg·min·g/mL (Chilian et al. 1989, dipyridamole, Table 2 / Fig. 5).
  Grown microvessels keep their Poiseuille resistances relative to each other, so flow
  heterogeneity still comes from the geometry, but one factor κ sets their absolute level to the
  measurement: CCO's smallest arterioles are about 2.4× longer than Kassab's (CCO runs each
  terminal to the centre of its territory, where real tissue uses the capillary network), so
  their raw Poiseuille resistance is far too high. Vessels above the cut-offs are pure Poiseuille.
  This anchors the scaled quantity itself to its own measurement; stress flow and flow reserve
  stay predictions.

Arteriolar wall mechanics (autoregulation; Carlson, Arciero & Secomb 2008, Arciero et al. 2008)
- `d_active`: arteries whose passive diameter is below 150 µm carry smooth-muscle tone. In the
  beating dog heart adenosine and dipyridamole dilate only arterial microvessels below 150 µm
  (Kanatsuka et al. 1989), so maximal dilation does not reach the larger ones. Larger arteries
  are rigid conduits.
- Wall tension at equilibrium: P·D/2 = T_pass(D) + A·T_act,max(D), with the diameter-dependent
  fits of Carlson & Secomb (2005) as tabulated by Fry et al. (2013): C_pass = 6.666·D₀ dyn/cm,
  C′_pass = 12.52 − 0.027·D₀, C_act = 1.30·D₀^1.48 dyn/cm, C′_act = 0.773 − 0.00059·D₀,
  C″_act = 0.415 − 0.00080·D₀ (D₀ in µm, the passive diameter at 100 mmHg). At its own pressure
  each artery can therefore be no wider than its relaxed wall (A = 0) and no narrower than its
  fully contracted wall (A = 1); its diameter lies between them as 1 − A. This is its capacity.
- Activation A = 1/(1 + e^{−S}), S = C_myo·T − `g_meta`·M + C_tone, with C_myo = 1.369/D₀ cm/dyn
  (Carlson & Secomb 2005): the myogenic and metabolic responses. The shear (flow-induced)
  response of Carlson et al. (2008) is left out: a CCO twig carries a real arteriole's flow in a
  vessel that is too long and too thin for it (see κ above), so its wall shear is not
  physiological, and shear plays a secondary role in coronary autoregulation (Namani et al. 2018).
- Metabolic signal M: red cells release ATP at a rate ∝ 1 − `atp_r1`·SO₂ as they desaturate
  (Arciero et al. 2008, from Ellsworth's measurements). M is a bed's release relative to its
  release at the reference saturation, minus 1, conducted upstream along the vessel wall with
  length constant `atp_l0` = 1 cm; each artery takes the mean over the beds it feeds (Namani et
  al. 2018). `g_meta` = 15 is the metabolic sensitivity: the size of the metabolic stimulus
  C_meta·S_CR in the control state of Carlson et al. (2008) (C_meta = 30 µM⁻¹cm⁻¹, S_CR ≈ 0.5
  µM·cm). It is not fitted to any flow.
- Oxygen: `mvo2` 0.116 mL O₂/min/g (5.2 µmol/g/min in healthy volunteers, McCommis et al. 2010), arterial
  saturation `sao2`, hemoglobin `hb` g/dL, 1.34 mL O₂ per g; venous saturation of each bed follows
  from its own flow (Fick).
- `oef`: the reference. Coronary flow at rest is matched to O₂ demand, and the regulated quantity
  is coronary venous oxygen (Feigl 1983; Tune et al. 2004). In the healthy heart at rest the
  single tone constant C_tone is set so that the coronary-sinus blood has the measured resting
  O₂ extraction, 0.73 (McCommis et al. 2010). No flow is prescribed: each arteriole finds its
  own diameter from its wall mechanics, its pressure and the O₂ state of its own beds, so rest
  flow, its distribution and every response to pressure or a stenosis follow. Pass the healthy
  rest circuit as `reference` to keep C_tone fixed (then the extraction is free to change).
"""
Base.@kwdef struct Physiology
    p_ao::F32 = 93.0f0
    p_back::F32 = 13.6f0
    hd::F32 = 0.45f0
    fb::F32 = 0.1f0
    share_art::F32 = 0.29f0
    share_cap::F32 = 0.31f0
    share_ven::F32 = 0.4f0
    cap_length_mm::F32 = 0.5f0
    cap_diameter_um::F32 = 6.0f0
    cap_ht_ratio::F32 = 0.85f0
    ve::F32 = 0.18f0
    ps::F32 = 0.37f0
    rho::F32 = 1.05f0
    plasma_viscosity_mPas::F32 = 1.2f0
    radius_cv::F32 = 0.12f0
    d_art_micro::F32 = 170.0f0
    d_ven_micro::F32 = 150.0f0
    r_micro_dilated::F32 = 4.0f0
    # autoregulation
    d_active::F32 = 150.0f0
    oef::F32 = 0.73f0
    g_meta::F32 = 15.0f0
    mvo2::F32 = 0.116f0
    sao2::F32 = 0.97f0
    hb::F32 = 14.5f0
    atp_r1::F32 = 1.026f0
    atp_l0::F32 = 1.0f0
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

const DYN_PER_MMHG = 1333.22f0

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
    ex::Vector{F32}              # first-pass extraction E of each bed (Renkin–Crone)
    k::Vector{F32}               # ISF→plasma back-flux rate E·F_p/v_e, 1/s
    offset::Vector{F32}          # arrival delay to the sinus, without the capillary and exchange
    rcap::Vector{F32}            # capillary-bed resistance, mmHg·s/mL
    so2::Vector{F32}             # venous O₂ saturation leaving each bed
    demand::Vector{F32}          # relative O₂ demand of each terminal's tissue (mean 1)
    # arteriolar state, per arterial segment
    xd::Vector{F32}              # diameter / passive diameter at 100 mmHg (1 for rigid arteries)
    act::Vector{F32}             # smooth-muscle activation A (0 = fully dilated; 0 for rigid arteries)
    ctone::Vector{F32}           # [C_tone], the tone constant of the healthy rest state
    # summary
    venous_scale::F32
    mass_g::F32
    volumes::NamedTuple
    resist::NamedTuple
end

"""
    circuit(tree, territory_ml; phys = Physiology(), state = :rest, veins = nothing, vleaf = nothing,
            reference = nothing, kappa = nothing, stenosis = nothing) -> Circuit

Steady coronary circuit driven by pressures: aortic pressure `phys.p_ao` at every ostium,
`phys.p_back` at the coronary-sinus ostium. Every terminal `i` (perfusing `territory_ml[i]`)
feeds a capillary bed whose resistance follows its capillary count (Kassab capillary
dimensions), then drains into the venous tree.

Every artery with a passive diameter ≤ `phys.d_active` is a vessel wall in force balance: its
diameter is where the wall tension P·D/2 equals passive plus active smooth-muscle tension (see
[`Physiology`](@ref)). Flows, pressures, oxygen saturations and the conducted metabolic signal
are solved together with all those diameters until they agree.

- `state = :stress`: maximal vasodilation (adenosine): activation 0 everywhere, so vessels sit
  on their passive pressure–diameter curves.
- `state = :rest`: autoregulation. Activation follows the myogenic and metabolic stimuli.
  Without `reference`, this is the healthy reference state: the one tone constant is set so the
  coronary-sinus O₂ extraction equals `phys.oef`. With `reference` (a rest `Circuit` of the
  healthy heart), that constant is kept, so a changed geometry (a stenosis) or changed pressures
  give the autoregulated response; flow is limited only when arterioles run out of dilation.

The venous tree is either grown (`veins`, with `vleaf[i]` the venous segment receiving terminal
`i`), or, by default, companion veins: each artery's vein has the same topology and holds
`share_ven/share_art` of its volume. `kappa` reuses the microvascular scale of another circuit.

`demand`: relative O₂ demand of each terminal's tissue (mean 1; default uniform). Local O₂ use
varies across the wall and local flow follows it (Bassingthwaighte, Beard & Li 2001; Deussen
2001); see [`demand_field`](@ref). A `reference` passes its own on.

`stenosis`: a list of `(seg, area, length_mm)`: arterial segment `seg` narrows to the fraction
`area` of its lumen area over `length_mm`. Its pressure loss is ΔP = f·Q + s·Q² (Young & Tsai
1973; Gould 1978): viscous f = 8πμL_s/A_s², separation s = ρ·k_e/(2A₀²)·(A₀/A_s − 1)² with
k_e = 1.52 (counted once, at the narrowing's exit). It is solved with the rest of the network, so the flow through it (and its loss) is
whatever the pressures and the distal arterioles make it.
"""
function circuit(
        t::VTree, territory_ml::AbstractVector; phys::Physiology = Physiology(), state::Symbol = :rest,
        veins::Union{Nothing, VTree} = nothing, vleaf = nothing, reference::Union{Nothing, Circuit} = nothing,
        kappa = nothing, stenosis = nothing, demand = nothing, verbose::Bool = false
    )
    state in (:rest, :stress) || throw(ArgumentError("state must be :rest or :stress"))
    cv(t, salt) = phys.radius_cv <= 0 ? t : (
            sig = sqrt(log(1 + Float64(phys.radius_cv)^2));
            VTree(t.parent, t.a, t.b, [t.r[j] * F32(exp(sig * _hashnormal(j + salt) - sig^2 / 2)) for j in 1:nseg(t)], t.isbase)
        )
    t = cv(t, 0)
    veins === nothing || (veins = cv(veins, 1 << 40))
    return _circuit(t, territory_ml; phys, state, veins, vleaf, reference, kappa, stenosis, demand, verbose)
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
    fa::Vector{F32}              # arterial resistance factor from the current diameters, (D₀/D)⁴·μ(D)/μ(D₀)
    Pa::F32; Pb::F32
    G::Vector{F32}; E::Vector{F32}; I::Vector{F32}
    q::Vector{F32}; pa::Vector{F32}; qt::Vector{F32}; qv::Vector{F32}; pv::Vector{F32}
    pvt::Vector{F32}; vE::Vector{F32}; vI::Vector{F32}; iters::Base.RefValue{Int}
    lev::Vector{Int32}           # start of each depth level in `order` (breadth-first), plus n + 1
    rr::Vector{F32}              # effective arterial resistances of the current solve
end

function _Solver(parent, order, ptr, kids, terminals, tidx, vparent, vorder, vleaf, ra, rv, rcap, ma, mv, Pa, Pb)
    n = length(parent); nv = length(vparent); nt = length(terminals)
    z(k) = zeros(F32, k)
    # depth levels are contiguous in the breadth-first order
    dep = zeros(Int32, n); lev = Int32[1]
    for h in 1:n
        j = order[h]; p = parent[j]; dep[j] = p == 0 ? 1 : dep[p] + 1
        h > 1 && dep[j] != dep[order[h - 1]] && push!(lev, h)
    end
    push!(lev, n + 1)
    return _Solver(
        parent, order, ptr, kids, terminals, tidx, vparent, vorder, vleaf, ra, rv, rcap, ma, mv, ones(F32, n), F32(Pa), F32(Pb),
        z(n), z(n), z(n), z(n), z(n), z(nt), z(nv), z(nv), fill(F32(Pb), nt), z(nv), z(nv), Ref(0), lev, z(n)
    )
end

# Level-synchronous tree pass: every segment of one depth level at once (a 1 M-terminal CCO
# tree is only ~150 levels deep, and most segments sit on a few dozen wide levels). `up` runs
# from the deepest level to the roots (children before parents), otherwise roots first.
@inline function _levels(f::F, S::_Solver, up::Bool) where {F}
    nl = length(S.lev) - 1
    for l in (up ? (nl:-1:1) : (1:nl))
        lo = S.lev[l]; m = S.lev[l + 1] - lo
        AK.foreachindex(1:m; min_elems = 2048) do k
            f(lo + k - 1)
        end
    end
    return nothing
end

@inline _ra(S::_Solver, kappa::F32, j) = (S.ma[j] ? kappa * S.ra[j] : S.ra[j]) * S.fa[j]

# Terminal-path operator: out_i = Σ_j∈path(i) r_j Q_j (+ bed resistance × q_i for `bed`), where Q_j
# is the sum of q over the terminals below segment j. For the arteries the path runs root → i,
# for the veins from the sinus to the vein draining i.
function _apply_art!(out, S::_Solver, kappa::F32, q, bed::Bool)
    Q = S.I; D = S.E; ptr = S.ptr; kids = S.kids
    _levels(S, true) do h
        @inbounds begin
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
    end
    _levels(S, false) do h
        @inbounds begin
            j = S.order[h]; p = S.parent[j]
            D[j] = (p == 0 ? 0.0f0 : D[p]) + Q[j] * S.rr[j]
            if ptr[j + 1] == ptr[j]
                i = S.tidx[j]
                out[i] = D[j] + (bed ? kappa * S.rcap[i] * q[i] : 0.0f0)
            end
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

# Preconditioner z = M_a⁻¹ r: the arterial tree with its root at 0 and terminal i draining
# through its bed into pressure −r_i (exact Norton solve).
function _precond!(z, S::_Solver, kappa::F32, r)
    G = S.G; E = S.E; ptr = S.ptr; kids = S.kids
    _levels(S, true) do h
        @inbounds begin
            j = S.order[h]
            if ptr[j + 1] == ptr[j]
                i = S.tidx[j]
                g = 1.0f0 / (kappa * S.rcap[i]); e = -r[i]
            else
                g = 0.0f0; ge = 0.0f0
                for c in ptr[j]:(ptr[j + 1] - 1)
                    k = kids[c]; g += G[k]; ge += G[k] * E[k]
                end
                e = g > 0 ? ge / g : 0.0f0
            end
            G[j] = g / (1.0f0 + g * S.rr[j]); E[j] = e
        end
    end
    P = S.pa
    _levels(S, false) do h
        @inbounds begin
            j = S.order[h]; p = S.parent[j]
            pin = p == 0 ? 0.0f0 : P[p]
            qj = G[j] * (pin - E[j])
            P[j] = pin - qj * S.rr[j]
            ptr[j + 1] == ptr[j] && (z[S.tidx[j]] = qj)
        end
    end
    return z
end

# One steady solve for the current resistances.
#
# Unknowns are the terminal flows q. Arterial pressure drop to terminal i plus its own bed is
# (M_a q)_i (shared path resistances + bed), the venous rise from the sinus is (M_v q)_i, and
# continuity at every bed is (M_a + M_v) q = P_ao − P_back. Both matrices are symmetric positive
# definite and applied in O(N) by two tree passes, so the system is solved by conjugate gradients
# preconditioned with the exact arterial tree solve (M_a⁻¹ is the Norton solve of the arterial
# tree); its condition number is 1 + λ(M_a⁻¹M_v), a few, so it converges in tens of iterations.
# `warm` starts from the previous terminal flows.
function _solve!(S::_Solver, kappa::F32; warm::Bool = false, tol::F32 = 1.0f-6)
    nt = length(S.terminals); dP = S.Pa - S.Pb
    rr = S.rr
    AK.foreachindex(rr; min_elems = 4096) do j
        rr[j] = _ra(S, kappa, j)
    end
    t1 = zeros(F32, nt); t2 = zeros(F32, nt)
    Aop!(out, x) = (_apply_art!(t1, S, kappa, x, true); _apply_ven!(t2, S, kappa, x); out .= t1 .+ t2; out)
    b = fill(dP, nt)
    x = warm ? copy(S.qt) : _precond!(zeros(F32, nt), S, kappa, b)     # cold start: arteries alone
    r = similar(x); z = similar(x); p = similar(x); Ap = similar(x)
    Aop!(Ap, x); @. r = b - Ap
    _precond!(z, S, kappa, r); p .= z
    rz = sum(r .* z); bn = sqrt(sum(abs2, b))
    for it in 1:500
        sqrt(sum(abs2, r)) < tol * bn && break
        S.iters[] += 1
        Aop!(Ap, p)
        α = rz / sum(p .* Ap)
        @. x += α * p; @. r -= α * Ap
        _precond!(z, S, kappa, r)
        rz2 = sum(r .* z); β = rz2 / rz; rz = rz2
        @. p = z + β * p
    end
    S.qt .= x
    # flows and pressures from the terminal flows
    _apply_art!(t1, S, kappa, x, false)
    @inbounds for j in eachindex(S.parent)
        S.q[j] = S.I[j]; S.pa[j] = S.Pa - S.E[j]
    end
    _apply_ven!(t2, S, kappa, x)
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

# ---------------------------------------------------------------------------------------------
# Arteriolar wall mechanics and the metabolic signal
# ---------------------------------------------------------------------------------------------

# Vasoactive arteries. Wall constants follow D₀ (µm; tensions dyn/cm): the Carlson & Secomb
# (2005) fits as tabulated by Fry et al. (2013, Table 1).
struct _Wall
    seg::Vector{Int32}           # arterial segment
    d0::Vector{F32}              # passive diameter at 100 mmHg, µm
    xd::Vector{F32}              # D/D₀
    act::Vector{F32}             # activation A
    xhi::Vector{F32}             # D/D₀ of the relaxed wall (A = 0) at the current pressure
    xlo::Vector{F32}             # D/D₀ of the fully contracted wall (A = 1) at the current pressure
end
_Wall(seg, d0) = (n = length(seg); _Wall(seg, d0, ones(F32, n), zeros(F32, n), ones(F32, n), fill(0.5f0, n)))

# The relaxed and fully contracted equilibria depend only on D₀ and pressure, so they are
# tabulated once (log D₀ × P, followed along each branch as pressure rises from 1 mmHg) and
# interpolated bilinearly.
struct _Limits
    ld0::F32; dld::F32; nd::Int      # log D₀ grid
    dp::F32; np::Int                 # pressure grid from 1 mmHg
    hi::Matrix{F32}; lo::Matrix{F32}
end
function _Limits(dmin::F32, dmax::F32, pmax::F32)
    nd = 96; ld0 = log(dmin) - 0.01f0; dld = (log(dmax) + 0.01f0 - ld0) / (nd - 1)
    dp = 0.25f0; np = ceil(Int, (pmax + 5.0f0 - 1.0f0) / dp) + 1
    hi = Matrix{F32}(undef, np, nd); lo = similar(hi)
    Threads.@threads for k in 1:nd
        d0 = exp(ld0 + (k - 1) * dld); xh = 1.0f0; xl = 0.5f0
        for i in 1:np
            P = 1.0f0 + (i - 1) * dp
            xh = _equilibrium(xh, 0.0f0, P, d0); xl = _equilibrium(xl, 1.0f0, P, d0)
            hi[i, k] = xh; lo[i, k] = xl
        end
    end
    return _Limits(ld0, dld, nd, dp, np, hi, lo)
end
@inline function _limits(L::_Limits, d0::F32, P::F32)
    u = clamp((log(d0) - L.ld0) / L.dld, 0.0f0, F32(L.nd - 1) - 1.0f-4); k = floor(Int, u); fu = u - k
    v = clamp((P - 1.0f0) / L.dp, 0.0f0, F32(L.np - 1) - 1.0f-4); i = floor(Int, v); fv = v - i
    @inbounds begin
        h = (1 - fu) * ((1 - fv) * L.hi[i + 1, k + 1] + fv * L.hi[i + 2, k + 1]) + fu * ((1 - fv) * L.hi[i + 1, k + 2] + fv * L.hi[i + 2, k + 2])
        l = (1 - fu) * ((1 - fv) * L.lo[i + 1, k + 1] + fv * L.lo[i + 2, k + 1]) + fu * ((1 - fv) * L.lo[i + 1, k + 2] + fv * L.lo[i + 2, k + 2])
    end
    return h, l
end

@inline function _tension(x::F32, A::F32, d0::F32)
    cp = 6.666f0 * d0; cp1 = 12.52f0 - 0.027f0 * d0
    ca = 1.3f0 * d0^1.48f0; ca1 = 0.773f0 - 0.00059f0 * d0; ca2 = 0.415f0 - 0.0008f0 * d0
    return cp * exp(cp1 * (x - 1.0f0)) + A * ca * exp(-((x - ca1) / ca2)^2)
end

# Stable equilibrium diameter ratio near x0 at activation A and transmural pressure P (mmHg):
# wall tension T_pass + A·T_act,max minus the load P·D/2 crosses zero upward there (a smaller
# vessel yields, a larger one contracts back).
@inline function _equilibrium(x0::F32, A::F32, P::F32, d0::F32)
    load = P * DYN_PER_MMHG * d0 * 1.0f-4 / 2.0f0
    g(x) = _tension(x, A, d0) - load * x
    lo = 0.05f0; hi = 2.0f0; h = 0.02f0
    x0 = clamp(x0, lo, hi)
    if g(x0) > 0
        b = x0; a = b - h
        while a > lo && g(a) > 0
            b = a; a -= h
        end
        a <= lo && g(lo) > 0 && return lo
        a = max(a, lo)
    else
        a = x0; b = a + h
        while b < hi && g(b) <= 0
            a = b; b += h
        end
        b = min(b, hi)
    end
    for _ in 1:20                     # bracket 0.02 wide → 2e-8
        m = (a + b) / 2
        g(m) > 0 ? (b = m) : (a = m)
    end
    return (a + b) / 2
end

# Metabolic signal at the midpoint of every vasoactive artery. Each bed's venous saturation
# follows from its flow and O₂ consumption (Fick). Its red cells release ATP at a rate
# ∝ 1 − R₁·SO₂ (Arciero et al. 2008, from Ellsworth's measurements); the signal is that release
# relative to its value at the reference saturation, minus 1 (0 at the reference, > 0 when the
# bed is short of O₂). It is conducted upstream along the vessel wall with length constant L₀,
# and each artery takes the mean over the beds it feeds, each decayed by its distance to them
# (Namani et al. 2018).
function _metabolic!(sig::Vector{F32}, so2::Vector{F32}, S::_Solver, phys::Physiology, terr_g::Vector{F32}, len::Vector{F32})
    nt = length(S.terminals); n = length(S.parent)
    cO2 = 1.34f0 * phys.hb / 100.0f0; r1 = phys.atp_r1; l0 = 10.0f0 * phys.atp_l0      # L₀ in mm
    rref = 1.0f0 - r1 * phys.sao2 * (1 - phys.oef)
    AK.foreachindex(so2; min_elems = 4096) do i
        qi = max(S.qt[i], 1.0f-12)
        so2[i] = clamp(phys.sao2 - phys.mvo2 * terr_g[i] / 60.0f0 / (qi * cO2), 0.0f0, phys.sao2)
    end
    sumv = S.G; cnt = S.E                         # scratch (overwritten by the next solve)
    _levels(S, true) do h
        @inbounds begin
            j = S.order[h]
            if S.ptr[j + 1] == S.ptr[j]
                sumv[j] = (1.0f0 - r1 * so2[S.tidx[j]]) / rref - 1.0f0; cnt[j] = 1.0f0
            else
                s = 0.0f0; c = 0.0f0
                for cc in S.ptr[j]:(S.ptr[j + 1] - 1)
                    k = S.kids[cc]
                    s += exp(-len[k] / l0) * sumv[k]; c += cnt[k]
                end
                sumv[j] = s; cnt[j] = c
            end
            sig[j] = sumv[j] * exp(-0.5f0 * len[j] / l0) / cnt[j]
        end
    end
    return sig
end

# Autoregulated (or maximally dilated) steady state: flows and diameters solved together.
# Every vasoactive artery's diameter lies between its relaxed (A = 0) and fully contracted
# (A = 1) wall equilibria at its own pressure, D = D_A=1 + (D_A=0 − D_A=1)·(1 − A). Activation
# A = 1/(1 + e^{−S}) with S = C_myo·P·D/2 − `g_meta`·(metabolic signal) + C_tone relaxes toward
# its stimulus (τ_a ≫ τ_d, Carlson et al. 2008); the network is re-solved after every step.
# Healthy reference (`fixtone = false`): C_tone is integrated until the coronary-sinus saturation
# equals the measured resting extraction (`phys.oef`).
function _autoregulate!(
        S::_Solver, W::_Wall, kappa::F32, phys::Physiology, state::Symbol, ctone::Vector{F32}, fixtone::Bool,
        terr_g, len, mu0, sten, lim, verbose
    )
    n = length(S.parent); nw = length(W.seg); hd = phys.hd
    sig = zeros(F32, n); so2 = zeros(F32, length(S.terminals))
    sref = phys.sao2 * (1 - phys.oef); gain = 2.0f0; lasterr = 0.0f0
    # Activation relaxes toward its stimulus with a step per vessel: the loop through a terminal
    # arteriole's own bed has a gain of ~20, so a fixed explicit step either crawls or makes some
    # vessels flip back and forth. A vessel's step halves whenever its activation reverses
    # direction and grows by 10% otherwise (bounded 0.005–0.3).
    ωw = fill(0.1f0, nw); dA = zeros(F32, nw)
    function setfa!()
        @inbounds for w in 1:nw
            j = W.seg[w]; d = W.d0[w] * W.xd[w]
            S.fa[j] = mu_rel(d, hd) / mu0[j] / W.xd[w]^4
        end
        # stenoses: unchanged lumen outside the narrowing plus the f·Q + s·Q² loss of the throat
        for (j, frac, f, s2) in sten
            S.fa[j] = (1 - frac) + (f + s2 * abs(S.q[j])) / S.ra[j]
        end
        return nothing
    end
    setfa!(); _solve!(S, kappa; warm = sum(S.qt) > 0)   # warm when called again (κ loop)
    qprev = sum(S.qt); iters = 0; dx = 0.0f0; ssin = 0.0f0; qhist = fill(qprev, 10); dq1 = 1.0f0
    maxit = state === :stress ? 60 : 600
    stim = zeros(F32, nw)
    for it in 1:maxit
        iters = it
        if state === :rest
            _metabolic!(sig, so2, S, phys, terr_g, len)
            ssin = sum(max.(S.qt, 0.0f0) .* so2) / sum(max.(S.qt, 0.0f0))
        end
        Threads.@threads for w in 1:nw
            @inbounds begin
                j = W.seg[w]; p = S.parent[j]
                P = max(((p == 0 ? S.Pa : S.pa[p]) + S.pa[j]) / 2, 1.0f0)
                stim[w] = 0.09128f0 * P * W.xd[w] - phys.g_meta * sig[j]
            end
        end
        if state === :rest && !fixtone
            if it == 1
                ctone[1] = -sort(stim)[(nw + 1) ÷ 2]          # start from half activation on average
            else
                err = ssin - sref                                # too much flow → more tone
                err * lasterr < 0 && (gain *= 0.7f0)
                ctone[1] += clamp(gain * err, -0.5f0, 0.5f0); lasterr = err
            end
        end
        xold = copy(W.xd)
        Threads.@threads for w in 1:nw
            @inbounds begin
                j = W.seg[w]; p = S.parent[j]
                P = max(((p == 0 ? S.Pa : S.pa[p]) + S.pa[j]) / 2, 1.0f0)
                if state === :stress
                    W.act[w] = 0.0f0
                else
                    d = 1.0f0 / (1.0f0 + exp(-(stim[w] + ctone[1]))) - W.act[w]
                    ωw[w] = d * dA[w] < 0 ? max(0.005f0, 0.5f0 * ωw[w]) : min(0.3f0, 1.1f0 * ωw[w])
                    dA[w] = d; W.act[w] += ωw[w] * d
                end
                W.xhi[w], W.xlo[w] = _limits(lim, W.d0[w], P)
                W.xd[w] = W.xlo[w] + (W.xhi[w] - W.xlo[w]) * (1.0f0 - W.act[w])
            end
        end
        # a few arterioles may keep cycling (vasomotion), so convergence is judged on total flow
        # (and on the diameters' 99th-percentile step at maximal dilation, which has no cycling)
        dx = nw > 0 ? (dd = abs.(W.xd .- xold); partialsort!(dd, max(1, ceil(Int, 0.99 * nw)))) : 0.0f0
        # inexact network solves while the diameters are still moving; exact once they settle
        setfa!(); _solve!(S, kappa; warm = true, tol = clamp(0.1f0 * dq1, 1.0f-6, 1.0f-3))
        qt = sum(S.qt); dq = abs(qt - qhist[mod1(it, 10)]) / qt; qhist[mod1(it, 10)] = qt
        dq1 = abs(qt - qprev) / qt; qprev = qt
        verbose && (it % 10 == 0 || it < 5) && (println("  autoregulation $it: flow $(round(60qt, digits = 3)) mL/min, sinus SO₂ $(round(ssin, digits = 3)), C_tone $(round(ctone[1], digits = 3)), p99 |ΔD/D₀| $(round(dx, sigdigits = 2)), mean A $(round(sum(W.act) / max(nw, 1), digits = 3))"); flush(stdout))
        # maximal dilation has no vasomotion: consecutive steps decide; at rest, ten steps
        if state === :stress
            it > 2 && dq1 < 1.0f-5 && dx < 1.0f-3 && break
        else
            it > 10 && dq < 1.0f-4 && (fixtone || abs(ssin - sref) < 1.0f-3) && break
        end
    end
    _solve!(S, kappa; warm = true)
    state === :rest && _metabolic!(sig, so2, S, phys, terr_g, len)
    return so2, iters, dx
end

function _circuit(t::VTree, territory_ml::AbstractVector; phys::Physiology, state::Symbol, veins, vleaf, reference, kappa, stenosis, demand, verbose)
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
    vv_unres = max(0.0f0, phys.share_ven * vblood - vven)
    vcap_total = phys.share_cap * vblood
    ra = similar(len); rv = Vector{F32}(undef, nv); mu0 = similar(len)
    AK.foreachindex(1:n; min_elems = 4096) do j
        mu0[j] = mu_rel(2000.0f0 * t.r[j], hd)
        ra[j] = poiseuille(len[j], t.r[j], mup * mu0[j])
    end
    AK.foreachindex(1:nv; min_elems = 4096) do k
        rv[k] = poiseuille(vlen[k], vt.r[k], mup * mu_rel(2000.0f0 * vt.r[k], hd))
    end
    # capillary bed of terminal i: its share of the capillary volume as parallel functional capillaries
    rc = phys.cap_diameter_um * 0.5f-3
    v1 = F32(π) * rc^2 * phys.cap_length_mm * 1.0f-3
    r1 = poiseuille(phys.cap_length_mm, rc, mup * mu_rel(phys.cap_diameter_um, hd))
    w = F32.(territory_ml) ./ vol_tissue
    ncap = [vcap_total * w[i] / v1 for i in 1:nt]
    rcap = [r1 / ncap[i] for i in 1:nt]
    # microcirculation: arteries < d_art_micro, veins < d_ven_micro and the capillary beds
    ma = [2000.0f0 * t.r[j] < phys.d_art_micro for j in 1:n]
    mv = [2000.0f0 * vt.r[k] < phys.d_ven_micro for k in 1:nv]
    S = _Solver(t.parent, order, ptr, kids, terminals, tidx, vt.parent, vorder, vl, ra, rv, rcap, ma, mv, phys.p_ao, phys.p_back)
    # vasoactive arteries
    act = Int32[j for j in 1:n if 2000.0f0 * t.r[j] < phys.d_active]
    d0 = F32[2000.0f0 * t.r[j] for j in act]
    W = _Wall(act, d0)
    lim = _Limits(isempty(d0) ? 1.0f0 : minimum(d0), isempty(d0) ? 2.0f0 : maximum(d0), max(phys.p_ao, 1.0f0))
    if reference !== nothing
        rx = reference.xd; ra_ = reference.act
        for (w, j) in enumerate(act)
            W.xd[w] = rx[j]; W.act[w] = ra_[j]
        end
    end
    # O₂ use of each terminal's tissue: its mass × the local relative demand (mean 1)
    demand === nothing && reference !== nothing && (demand = reference.demand)
    dem = demand === nothing ? ones(F32, nt) : F32.(demand)
    terr_g = F32.(territory_ml) .* phys.rho .* dem
    # stenoses: (segment, stenosed fraction of its length, f, s) with f in mmHg·s/mL and s in
    # mmHg·s²/mL² (blood density 1060 kg/m³, the segment's own apparent viscosity)
    # A narrowing that runs over several segments loses its separation (exit) pressure once, at its
    # most downstream segment.
    sten = Tuple{Int32, F32, F32, F32}[]
    stens = stenosis === nothing ? () : stenosis
    segs = Set(Int32(e[1]) for e in stens)
    for (j, area, Ls) in stens
        A0 = F32(π) * (t.r[j] * 1.0f-3)^2; As = F32(area) * A0
        mu = mup * 1.0f-3 * mu_rel(2000.0f0 * t.r[j], hd)
        Ls = min(F32(Ls), len[j])
        f = 8.0f0 * F32(π) * mu * Ls * 1.0f-3 / As^2 * 1.0f-6 / 133.322f0
        exit = !any(k -> t.parent[k] == j, segs)
        s2 = exit ? 1060.0f0 * 1.52f0 / (2A0^2) * (1 / F32(area) - 1)^2 * 1.0f-12 / 133.322f0 : 0.0f0
        push!(sten, (Int32(j), Ls / len[j], f, s2))
    end

    # ---- microvascular scale κ at maximal dilation (vessels on their passive curves): the
    # Chilian-style resistance between arteriolar inlets and venular outlets equals the measurement
    ctone = reference === nothing ? F32[-5.0f0] : copy(reference.ctone)
    kappa === nothing && reference !== nothing && (kappa = reference.resist.kappa)
    if kappa === nothing
        kappa = 1.0f0
        Wd = _Wall(act, d0)
        for _ in 1:12
            _autoregulate!(S, Wd, F32(kappa), phys, :stress, ctone, true, terr_g, len, mu0, empty(sten), lim, false)   # healthy geometry
            ratio = phys.r_micro_dilated / _micro_resistance(S, mass_g)[1]
            kappa *= ratio
            abs(ratio - 1) < 1.0f-3 && break
        end
    end
    kappa = F32(kappa)
    # ---- the state
    so2, iters, dxlast = _autoregulate!(S, W, kappa, phys, state, ctone, reference !== nothing, terr_g, len, mu0, sten, lim, verbose)
    q = S.q; pa = S.pa; qt = S.qt; qv = S.qv; pv = S.pv
    nneg = count(<=(0), qt)
    nneg > 0 && @warn "$nneg terminals with reversed or zero flow (venous pressure above arterial); their beds are treated as stagnant"
    @inbounds for i in 1:nt
        qt[i] = max(qt[i], 1.0f-9)
    end
    xd = ones(F32, n); av = zeros(F32, n)
    for (w, j) in enumerate(act)
        xd[j] = W.xd[w]; av[j] = W.act[w]
    end
    # constricted arterioles hold less blood; the unresolved arterial volume is what the dilated
    # tree leaves of the Kassab target, constricted with the terminal's own arteriole
    vol .*= xd .^ 2; vart_now = sum(vol)
    va_unres = max(0.0f0, phys.share_art * vblood - vart)

    # ---- plasma volumes and transit times
    vp_a = similar(len); tau_a = similar(len)
    AK.foreachindex(1:n; min_elems = 4096) do j
        vp_a[j] = vol[j] * (1.0f0 - fahraeus(2000.0f0 * t.r[j] * xd[j], hd) * hd)
        tau_a[j] = vp_a[j] / (max(q[j], 1.0f-12) * (1.0f0 - hd))      # reversed segments: stagnant
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
    vint = similar(fp); kk = similar(fp); off = similar(fp); ex = similar(fp)
    pf = 1.0f0 - phys.cap_ht_ratio * hd      # plasma fraction of microvascular blood
    AK.foreachindex(1:nt; min_elems = 4096) do i
        j = terminals[i]
        fp[i] = qt[i] * (1.0f0 - hd)
        ta[i] = w[i] * va_unres * xd[j]^2 * pf / fp[i]
        tc[i] = w[i] * vcap_total * pf / fp[i]
        tv[i] = w[i] * vv_unres * pf / fp[i]
        vint[i] = phys.ve * F32(territory_ml[i])
        ex[i] = 1.0f0 - exp(-phys.ps / 60.0f0 * F32(territory_ml[i]) * phys.rho / fp[i])
        kk[i] = ex[i] * fp[i] / vint[i]
        off[i] = t_end[j] + ta[i] + tv[i] + s_v[vl[i]]
    end

    # ---- pressure-drop shares (flow-weighted) and summary flows
    dpa = sum(q[j] * _ra(S, kappa, j) * q[j] for j in 1:n)
    dpv = sum(qv[k]^2 * (mv[k] ? kappa : 1.0f0) * rv[k] for k in 1:nv); dpc = sum(qt .^ 2 .* rcap) * kappa
    roots = Int32[j for j in 1:n if t.parent[j] == 0]
    qtot = sum(q[roots])
    nw = length(act)
    resist = (;
        dp_total = phys.p_ao - phys.p_back, art_share = dpa / (dpa + dpv + dpc),
        cap_share = dpc / (dpa + dpv + dpc), ven_share = dpv / (dpa + dpv + dpc),
        mbf = 60 * qtot / mass_g, root_flow = q[roots], kappa, iters, sweeps = S.iters[], residual_dx = dxlast,
        r_micro = (rmc = _micro_resistance(S, mass_g); rmc[1]),
        chilian_shares = ((phys.p_ao - rmc[2]), (rmc[2] - rmc[3]), (rmc[3] - phys.p_back)) ./ (phys.p_ao - phys.p_back),
        so2_sinus = sum(qt .* so2) / sum(qt), starved = count(<(0.05f0), so2) / nt, mean_act = nw > 0 ? sum(W.act) / nw : 0.0f0,
        n_dilated = count(<(0.01f0), W.act), n_max_tone = count(>(0.99f0), W.act),
        art_volume_ml = vart_now,
    )
    return Circuit(
        t, vt, phys, state, order, terminals, tidx, vl, len, q, qv, pa, pv, tau_a, tau_v, t_end, s_v, vp_a, vp_v,
        fp, ta, tc, tv, vint, ex, kk, off, rcap .* kappa, so2, dem, xd, av, ctone, s, mass_g,
        (;
            blood = vblood, art_resolved = vart_now, ven_resolved = vven, art_unresolved = va_unres,
            ven_unresolved = vv_unres, capillary = vcap_total,
        ), resist
    )
end

"""
    with_extraction(c::Circuit, E) -> Circuit

The same steady circuit (flows, diameters, volumes) with first-pass extraction `E` in every bed;
`E = 0` is an intravascular tracer. Flow does not depend on the tracer, so nothing is re-solved.
"""
function with_extraction(c::Circuit, E::Real)
    ex = fill(F32(E), length(c.ex)); k = ex .* c.fp ./ c.vint
    return Circuit((f === :ex ? ex : f === :k ? k : getfield(c, f) for f in fieldnames(Circuit))...)
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
    hd = c.phys.hd; dt = inp.dt; ex = to(c.ex)
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
            tis[f, i] = _gtissue(A, X, Y, logk, dt, fp[i], ta[i], tc[i], tv[i], kk[i], ex[i], tm[f] - t_end[j])
        end
    end
    # veins: walk each terminal's path to the sinus; the mass of its outflow inside vein j is
    # G(t − O_i + S[j]) − G(t − O_i + S[parent j]). Stop where its iodine has not yet arrived.
    AK.foreachindex(fp; min_elems = 256) do i
        j0 = vleaf[i]
        f_ = fp[i]; tcap = tc[i]; k = kk[i]; base = off[i]; E = ex[i]
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
    dt = inp.dt; ex = to(c.ex)
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
                out[f, ii] = _gout(A, Y, logk, dt, fp[i], tc[i], kk[i], ex[i], tm[f] - off[i])
            end
        end
        total .+= vec(Array(sum(out; dims = 2)))
    end
    return total
end

end # module
