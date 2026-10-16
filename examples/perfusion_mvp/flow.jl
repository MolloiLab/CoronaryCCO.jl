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
using TOML

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
    Physiology(; organ, kwargs...)

Boundary conditions, tissue constants and vessel-wall mechanics of one organ's vascular bed.
Flow is not prescribed: it follows from the geometry, the two pressures and the arteriolar
diameters (see [`circuit`](@ref)). Values are data, not code: they are read from
`data/organs/<organ>.toml` (or from a path), where each one names its sources in
`data/reference/`; keyword arguments override single values. The values below are the heart's.

Boundary conditions and blood
- `p_in`: mean inflow pressure. Heart: aortic pressure at the ostia, 90 mmHg (Chilian et al.
  1989, cat, the preparation of the resistance anchor below; normal human mean arterial pressure
  is similar).
- `p_out`: effective outflow pressure. Heart: 13.6 mmHg, the zero-flow pressure of the maximally
  dilated beating heart (Kanatsuka et al. 1990). It lumps coronary-sinus pressure and the
  waterfall effect of intramyocardial pressure into one steady back pressure.
- `regulation`: `:o2` (arteriolar tone follows the O₂ state of each bed, below) or `:passive`
  (no active tone: every state is the relaxed one).
- `venous_scale`: `:volume_share` (a grown venous tree is scaled to `share_ven/share_art` of the
  arterial volume) or `:grown` (its radii stay as grown, e.g. anchored at an imaged root).
- `capillary`: `:tube` (parallel functional capillaries, below) or `:sheet` (Fung & Sobin sheet
  flow, for alveolar walls): a bed of sheet area A between arterial pressure P_a and venous
  pressure P_v passes Q = S·A/(4µ·kf·L̄²·α)·[h(P_a)⁴ − h(P_v)⁴], with sheet thickness
  h(P) = max(h₀ + α·(P − `p_ext`), 0) (`p_ext`: alveolar pressure). `sheet_h0_um`, `sheet_alpha_um_cmH2O`,
  `sheet_s` (vascular-space-tissue ratio), `sheet_area_m2_per_ml` (sheet area per mL of tissue
  domain), `sheet_mu_cP` (apparent viscosity), `sheet_kf` (k·f) and `sheet_path_um` (mean path
  length L̄) set it; the sheet stops thickening above a transmural pressure of `sheet_pmax_cmH2O`.
  Each bed's resistance is then (P_a − P_v)/Q, updated with the pressures.
- `hd`: discharge hematocrit 0.45. Tube hematocrit from the Fåhræus law (Pries et al. 1992).

Volumes and exchange
- `fb`: intramyocardial blood volume fraction, mL/mL. Heart 0.09 (healthy human at rest by CMR,
  Nickander et al. 2019).
- `share_art`, `share_cap`, `share_ven`: 0.288 / 0.317 / 0.395 of that volume (Kassab et al.
  1994 cast volumes).
- `cap_ht_ratio`: tube/discharge hematocrit 0.85 in unresolved microvessels (Bindschadler 2014).
- `cap_length_mm`, `cap_diameter_um`: functional capillary 0.5 mm × 6 µm (Kassab & Fung 1994),
  used for the capillary-bed resistance.
- `ve`: interstitial volume fraction. Heart 0.224 (ECV 0.273 minus plasma, Nickander et al. 2019).
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
struct Physiology
    p_in::F32
    p_out::F32
    regulation::Symbol
    capillary::Symbol
    venous_scale::Symbol
    hd::F32
    fb::F32
    share_art::F32
    share_cap::F32
    share_ven::F32
    cap_length_mm::F32
    cap_diameter_um::F32
    cap_ht_ratio::F32
    ve::F32
    ps::F32
    rho::F32
    plasma_viscosity_mPas::F32
    radius_cv::F32
    d_art_micro::F32
    d_ven_micro::F32
    r_micro_dilated::F32
    # autoregulation
    d_active::F32
    oef::F32
    g_meta::F32
    mvo2::F32
    sao2::F32
    hb::F32
    atp_r1::F32
    atp_l0::F32
    # sheet-flow capillary bed (capillary = :sheet; Fung & Sobin)
    p_ext::F32
    sheet_h0_um::F32
    sheet_alpha_um_cmH2O::F32
    sheet_s::F32
    sheet_area_m2_per_ml::F32
    sheet_mu_cP::F32
    sheet_kf::F32
    sheet_path_um::F32
    sheet_pmax_cmH2O::F32
end

const ORGAN_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "organs"))
const REFERENCE_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "reference"))
const _REGULATION = (:o2, :passive)
const _WARNED = Set{String}()
const _CAPILLARY = (:tube, :sheet)
const _VENOUS_SCALE = (:volume_share, :grown)

# One organ's parameter set: name → value, every cited reference file checked to exist.
function organ_parameters(organ::AbstractString)
    path = isfile(organ) ? organ : joinpath(ORGAN_DIR, organ * ".toml")
    isfile(path) || throw(ArgumentError("no organ parameter file $path"))
    d = TOML.parsefile(path)["parameters"]
    for (k, v) in d, id in get(v, "sources", String[])
        isfile(joinpath(REFERENCE_DIR, id * ".toml")) || throw(ArgumentError("$path: $k cites missing reference $id"))
    end
    # values the model uses that cite no reference file in the repo: say so, once per file
    loose = sort([k for (k, v) in d if isempty(get(v, "sources", String[])) && !any(startswith(get(v, "note", ""), w) for w in ("Not used", "Definition"))])
    isempty(loose) || path in _WARNED || (push!(_WARNED, path); @warn "$(basename(path)): no reference file in data/reference for $(join(loose, ", "))")
    return Dict(Symbol(k) => v["value"] for (k, v) in d)
end

function Physiology(; organ::AbstractString, kwargs...)
    d = organ_parameters(organ)
    bad = setdiff(keys(kwargs), fieldnames(Physiology))
    isempty(bad) || throw(ArgumentError("unknown Physiology fields $(collect(bad))"))
    vals = map(fieldnames(Physiology)) do f
        v = haskey(kwargs, f) ? kwargs[f] : haskey(d, f) ? d[f] : throw(ArgumentError("organ $organ: no value for $f"))
        fieldtype(Physiology, f) === Symbol ? Symbol(v) : F32(v)
    end
    p = Physiology(vals...)
    p.regulation in _REGULATION || throw(ArgumentError("regulation must be one of $_REGULATION"))
    p.capillary in _CAPILLARY || throw(ArgumentError("capillary must be one of $_CAPILLARY"))
    p.venous_scale in _VENOUS_SCALE || throw(ArgumentError("venous_scale must be one of $_VENOUS_SCALE"))
    return p
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
    circuit(tree, territory_ml; phys, state = :rest, veins = nothing, vleaf = nothing,
            reference = nothing, kappa = nothing, kappa_init = 1, stenosis = nothing, demand = nothing, to = identity) -> Circuit

Steady coronary circuit driven by pressures: inflow pressure `phys.p_in` at every ostium,
outflow pressure `phys.p_out` at the coronary-sinus ostium. Every terminal `i` (perfusing `territory_ml[i]`)
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
`share_ven/share_art` of its volume. `kappa` reuses the microvascular scale of another circuit. Without it, κ is found by fixed-point
iteration from `kappa_init` (a nearby circuit's κ makes that one or two solves).

`demand`: relative O₂ demand of each terminal's tissue (mean 1; default uniform). Local O₂ use
varies across the wall and local flow follows it (Bassingthwaighte, Beard & Li 2001; Deussen
2001); see [`demand_field`](@ref). A `reference` passes its own on.

`stenosis`: a list of `(seg, area, length_mm)`: arterial segment `seg` narrows to the fraction
`area` of its lumen area over `length_mm`. Its pressure loss is ΔP = f·Q + s·Q² (Young & Tsai
1973; Gould 1978): viscous f = 8πμL_s/A_s², separation s = ρ·k_e/(2A₀²)·(A₀/A_s − 1)² with
k_e = 1.52 (counted once, at the narrowing's exit). It is solved with the rest of the network, so the flow through it (and its loss) is
whatever the pressures and the distal arterioles make it.

`to` moves the solver's arrays to a device (for example `to = MtlArray` or `CuArray`), where the
whole iteration then runs; the returned `Circuit` holds host arrays either way.
"""
function circuit(
        t::VTree, territory_ml::AbstractVector; phys::Physiology, state::Symbol = :rest,
        veins::Union{Nothing, VTree} = nothing, vleaf = nothing, reference::Union{Nothing, Circuit} = nothing,
        kappa = nothing, kappa_init = 1, stenosis = nothing, demand = nothing, verbose::Bool = false, to = identity
    )
    state in (:rest, :stress) || throw(ArgumentError("state must be :rest or :stress"))
    cv(t, salt) = phys.radius_cv <= 0 ? t : (
            sig = sqrt(log(1.0f0 + phys.radius_cv^2));
            VTree(t.parent, t.a, t.b, [t.r[j] * F32(exp(sig * _hashnormal(j + salt) - sig^2 / 2)) for j in 1:nseg(t)], t.isbase)
        )
    t = cv(t, 0)
    veins === nothing || (veins = cv(veins, 1 << 40))
    return _circuit(t, territory_ml; phys, state, veins, vleaf, reference, kappa, kappa_init, stenosis, demand, verbose, to)
end

# Fixed standard normal per segment index (SplitMix64 + Box–Muller), so results are reproducible.
function _hashnormal(j::Integer)
    sm(x) = (x += 0x9e3779b97f4a7c15; x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9; x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb; x ⊻ (x >> 31))
    a = sm(UInt64(j)); b = sm(a)
    u1 = (F32(a >> 40) + 0.5f0) / 2.0f0^24; u2 = F32(b >> 40) / 2.0f0^24            # 24-bit uniforms, Float32
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

# ---------------------------------------------------------------------------------------------
# Network solver: tree passes as prefix scans
# ---------------------------------------------------------------------------------------------
#
# The solver renumbers both trees in depth-first preorder, so that every subtree is a contiguous
# range. Arterial segments with children come first (1..ni, in preorder), then the terminal
# segments (ni + i is terminal i, also in preorder), so the terminals below any segment are a
# contiguous range of terminal numbers as well. Every tree pass is then a fixed number of kernels,
# however deep the tree:
# - a subtree sum (the flow through a segment is the sum of its terminals' flows) is the
#   difference of one prefix sum over the terminals at the two ends of the subtree's range;
# - a path sum from the root (the pressure drop down to a segment) is one prefix sum over the
#   Euler tour of the tree: entering segment k adds its value, leaving it subtracts it again. In
#   preorder, k is entered at position 2k − depth(k) and left 2·size(k) − 1 positions later.
# Both scans run in 64-bit integers, on a fixed-point scale (a power of two chosen from a bound on
# the result, so every result is below 2⁶¹), which makes them exact and independent of the order
# in which the scan adds: a difference of two prefix sums loses nothing, where in Float32 a
# terminal's flow (~10⁻⁶ of the total) would drown in the rounding of the total. Every other
# value is Float32.
#
# All per-iteration arrays live on the device of `to` (plain `Array`s and Julia threads by default).

# Power of two σ with σ·a·b < 2⁶¹ (a, b ≥ 0 and finite)
@inline _expo(a::F32) = Int32((reinterpret(UInt32, a) >> 23) & 0xff) - Int32(127)
@inline function _scale(a::F32, b::F32 = 1.0f0)
    k = clamp(Int32(59) - _expo(a) - _expo(b), Int32(-126), Int32(126))
    return reinterpret(F32, UInt32(k + Int32(127)) << 23)
end
@inline _fix(x::F32, σ::F32) = unsafe_trunc(Int64, x * σ)
@inline _unfix(v::Int64, σ::F32) = F32(v) / σ

const _CHUNK = 256      # largest subtree (segments with children) that one thread reduces alone

_get(r) = Array(r)[1]                               # a device scalar on the host
# Reduction of f(1), …, f(m) with `op` into the device scalar R[1]: partial results in `part`,
# each over runs of 16 consecutive elements spaced 16·length(part) apart (neighbouring threads
# read neighbouring runs), then one small reduction. Nothing returns to the host.
function _reduce!(R, part, op, init::F32, m::Int, f::F) where {F}
    B = length(part)
    AK.foreachindex(part; min_elems = 64) do b
        @inbounds begin
            acc = init
            for i0 in (16 * (b - 1)):(16 * B):(m - 1)
                for i in (i0 + 1):min(i0 + 16, m)
                    acc = op(acc, f(i))
                end
            end
            part[b] = acc
        end
    end
    AK.mapreducedim!(identity, op, R, part; init)
    return R
end
_dot!(R, S, a, b) = _reduce!(R, S.part, +, 0.0f0, length(a), i -> a[i] * b[i])

# Depth-first preorder of a forest (children in index order)
function _preorder(parent::Vector{Int32}, ptr, kids)
    n = length(parent); pre = Vector{Int32}(undef, n); h = 0
    stack = Int32[j for j in n:-1:1 if parent[j] == 0]
    while !isempty(stack)
        j = pop!(stack); h += 1; pre[h] = j
        for c in (ptr[j + 1] - 1):-1:ptr[j]
            push!(stack, kids[c])
        end
    end
    h == n || throw(ArgumentError("parent pointers do not form a forest"))
    return pre
end

# Arterial tree + terminal beds + venous tree, solved for flows and pressures, in the solver's own
# numbering (see above). Venous segments are numbered in preorder; `vsorted` lists the terminals
# by the preorder number of the vein they drain into, those of vein k at vcum[k]+1 .. vcum[k+1].
struct _Solver{VF <: AbstractVector{F32}, VI <: AbstractVector{Int32}, VB <: AbstractVector{Bool}, VL <: AbstractVector{Int64}, MF, BE}
    n::Int; ni::Int; nt::Int; nv::Int
    be::BE                        # backend of the arrays
    lev::Vector{Int32}            # start of each depth level in `lord`, plus its length + 1 (host)
    # arterial layout
    par::VI                       # parent (0 at a root), all segments
    sz::VI; dep::VI               # subtree size (segments with children) and depth, segments 1..ni
    tb::VI; te::VI                # terminals below segment k are tb[k]+1 .. te[k]
    cptr::VI; ckid::VI            # children of segments 1..ni (CSR)
    croot::VI                     # roots of the maximal subtrees of ≤ _CHUNK segments with children
    lord::VI                      # the larger segments sorted by depth
    # venous layout
    vpar::VI; vsz::VI; vdep::VI; vcum::VI
    vsorted::VI                   # terminals sorted by their draining vein
    vdrain::VI                    # draining vein of terminal i
    # resistances
    ra::VF; rv::VF; rcap::VF; ma::VB; mv::VB
    fa::VF                        # arterial resistance factor from the current diameters, (D₀/D)⁴·μ(D)/μ(D₀)
    rr::VF; rvv::VF               # effective arterial and venous resistances of the current solve
    Pa::F32; Pb::F32
    # scan buffers
    X::VL                         # 0, terminal flows, the flows draining into each vein
    Xm::VL                        # 0, metabolic signal of the terminals, by band
    Eu::VL                        # Euler tours: arterial (segments 1..ni), then venous
    Xp::VL; Ep::VL; Xv::VL; Ev::VL
    Qa::VF; Qv::VF; Dv::VF
    # preconditioner (Norton factors of the arterial tree)
    G::VF                         # conductance through each segment, G = g/(1 + g·r)
    lc::VF                        # log c, c = 1/(1 + g·r)
    gi::VF                        # conductance g below each segment with children
    lpi::VF                       # log π, π the path product of c (segments with children)
    lF::VF; lFt::VF               # log F, F = π·g, at the segments with children and the terminals
    Ek::VF; Pk::VF                # Norton sources and pressures of one application
    nbu::Base.RefValue{Int}; nbd::Base.RefValue{Int}   # bands of the two passes (see `_precond!`)
    # solution
    q::VF; pa::VF; qt::VF; qv::VF; pv::VF
    x::VF; r::VF; z::VF; p::VF; Ap::VF; Mp::VF
    # metabolic signal (static weights, see `_metabolic!`)
    mlf::VF; mbw::VF; mnf::VF; mband::VI; nband::Int
    val::VF; sig::VF; so2::VF
    # wall-mechanics tables on the device, and the apparent viscosity μ_rel(D) on a log D grid
    lhi::MF; llo::MF
    mut::VF; mua::F32; mudl::F32
    # device scalars
    part::VF                      # partial results of reductions
    s_x::VF; s_rz::VF; s_rz2::VF; s_pap::VF; s_rn::VF; s_Ra::VF; s_Rv::VF; s_lF::VF; s_lFm::VF; s_lpm::VF; s_w0::VF; s_wn::VF; s_b::VF; s_v::VF
    iters::Base.RefValue{Int}
    solved::Base.RefValue{Bool}
    kv::Base.RefValue{F32}        # κ of the current venous resistances
    # host copies for the summaries
    hpar::Vector{Int32}; hvpar::Vector{Int32}; hma::Vector{Bool}; hmv::Vector{Bool}
end

# Builds the solver in its own numbering. Returns it and the permutations back to the caller's
# numbering: `perm[k]` is the caller's arterial segment k of the solver, `tperm[i]` the caller's
# terminal i, `vperm[k]` the caller's vein k.
function _solver(
        parent, nchild, ptr, kids, vparent, vptr, vkids, vleaf, tidx, ra, rv, rcap, ma, mv, len, lim, Pa, Pb, l0, hd, to
    )
    n = length(parent); nv = length(vparent)
    # arteries: segments with children first, then the terminals, both in preorder
    ni = count(>(0), nchild); nt = n - ni
    perm = Vector{Int32}(undef, n); tb = zeros(Int32, ni); a = 0; b = ni
    for j in _preorder(parent, ptr, kids)
        if nchild[j] > 0
            a += 1; perm[a] = j; tb[a] = b - ni
        else
            b += 1; perm[b] = j
        end
    end
    iperm = similar(perm); iperm[perm] = 1:n
    par = Int32[parent[perm[k]] == 0 ? 0 : iperm[parent[perm[k]]] for k in 1:n]
    sz = ones(Int32, ni); nl = zeros(Int32, ni)
    for k in n:-1:(ni + 1)
        p = par[k]; p > 0 && (nl[p] += 1)
    end
    for k in ni:-1:1
        p = par[k]; p > 0 && (sz[p] += sz[k]; nl[p] += nl[k])
    end
    te = tb .+ nl
    dep = zeros(Int32, ni)
    for k in 1:ni
        p = par[k]; dep[k] = p == 0 ? 1 : dep[p] + 1
    end
    # bottom-up passes: small subtrees one per thread, the segments above them by depth level
    croot = Int32[k for k in 1:ni if sz[k] <= _CHUNK && (par[k] == 0 || sz[par[k]] > _CHUNK)]
    big = Int32[k for k in 1:ni if sz[k] > _CHUNK]
    nd = maximum(view(dep, big); init = Int32(0))
    lev = zeros(Int32, nd + 1); lev[1] = 1
    for k in big
        lev[dep[k] + 1] += 1
    end
    lev = accumulate(+, lev); fillp = lev[1:nd] .+ 0; lord = Vector{Int32}(undef, length(big))
    for k in big
        lord[fillp[dep[k]]] = k; fillp[dep[k]] += 1
    end
    cptr = zeros(Int32, ni + 1); cptr[1] = 1
    for k in 1:n
        p = par[k]; p > 0 && (cptr[p + 1] += 1)
    end
    cptr = accumulate(+, cptr); fillp = cptr[1:ni] .+ 0; ckid = Vector{Int32}(undef, n)
    for k in 1:n
        p = par[k]; p > 0 && (ckid[fillp[p]] = k; fillp[p] += 1)
    end
    resize!(ckid, cptr[end] - 1)
    tperm = Int32[tidx[perm[ni + i]] for i in 1:nt]
    # veins: preorder
    vperm = _preorder(vparent, vptr, vkids)
    viperm = similar(vperm); viperm[vperm] = 1:nv
    vpar = Int32[vparent[vperm[k]] == 0 ? 0 : viperm[vparent[vperm[k]]] for k in 1:nv]
    vsz = ones(Int32, nv)
    for k in nv:-1:1
        p = vpar[k]; p > 0 && (vsz[p] += vsz[k])
    end
    vdep = zeros(Int32, nv)
    for k in 1:nv
        p = vpar[k]; vdep[k] = p == 0 ? 1 : vdep[p] + 1
    end
    vdrain = Int32[viperm[vleaf[tperm[i]]] for i in 1:nt]
    vsorted = Int32.(sortperm(vdrain))
    vcum = zeros(Int32, nv + 1)
    for i in 1:nt
        vcum[vdrain[i] + 1] += 1
    end
    vcum = accumulate(+, vcum)
    # path lengths ℓ for the conducted metabolic signal, whose subtree sums run in bands of ℓ
    lp = len[perm]; Lp = zeros(F32, n)
    for k in 1:n
        p = par[k]; Lp[k] = (p == 0 ? 0.0f0 : Lp[p]) + lp[k]
    end
    band = 12.0f0 * l0
    nband = floor(Int, maximum(view(Lp, 1:ni); init = 0.0f0) / band) + 1
    mlf = Vector{F32}(undef, nt); mbw = Vector{F32}(undef, nt * nband)
    mnf = Vector{F32}(undef, ni); mband = Vector{Int32}(undef, ni)
    AK.foreachindex(mlf; min_elems = 4096) do i
        @inbounds begin
            L = Lp[ni + i]; mlf[i] = exp(-0.5f0 * lp[ni + i] / l0)
            for b in 0:(nband - 1)
                mbw[b * nt + i] = L >= b * band ? exp(-(L - b * band) / l0) : 0.0f0
            end
        end
    end
    AK.foreachindex(mnf; min_elems = 4096) do k
        @inbounds begin
            bk = floor(Int32, Lp[k] / band); mband[k] = bk
            mnf[k] = exp((Lp[k] - bk * band) / l0) * exp(-0.5f0 * lp[k] / l0) / F32(te[k] - tb[k])
        end
    end
    # μ_rel on a log-spaced grid of diameters from 2 µm (where `mu_rel` clamps) to 400 µm
    nmu = 16384; mua = log(2.0f0); mudl = (log(400.0f0) - mua) / (nmu - 1)
    mut = F32[mu_rel(exp(mua + (m - 1) * mudl), hd) for m in 1:nmu]
    z(k) = to(zeros(F32, k)); zl(k) = to(zeros(Int64, k)); s() = to(zeros(F32, 1))
    S = _Solver(
        n, ni, nt, nv, AK.get_backend(to(zeros(F32, 1))), lev,
        to(par), to(sz), to(dep), to(tb), to(te), to(cptr), to(ckid), to(croot), to(lord),
        to(vpar), to(vsz), to(vdep), to(vcum), to(vsorted), to(vdrain),
        to(ra[perm]), to(rv[vperm]), to(F32.(rcap[tperm])), to(ma[perm]), to(mv[vperm]),
        to(ones(F32, n)), z(n), z(nv), F32(Pa), F32(Pb),
        zl(1 + nt + nv), zl(1 + nband * nt), zl(2ni + 2nv), zl(1 + nt), zl(2ni), zl(1 + nv), zl(2nv),
        z(ni), z(nv), z(nv),
        z(n), z(n), z(ni), z(ni), z(ni), z(nt), z(ni), z(ni), Ref(1), Ref(1),
        z(n), z(n), z(nt), z(nv), z(nv),
        z(nt), z(nt), z(nt), z(nt), z(nt), z(nt),
        to(mlf), to(mbw), to(mnf), to(mband), nband, z(nt), z(n), z(nt),
        to(lim.hi), to(lim.lo), to(mut), mua, mudl,
        z(4096), s(), s(), s(), s(), s(), s(), s(), s(), s(), s(), s(), to(F32[3nt]), s(), s(),
        Ref(0), Ref(false), Ref(-1.0f0),
        par, vpar, ma[perm], mv[vperm]
    )
    return S, perm, iperm, tperm, vperm
end

# Path sums D[k] = Σ w over the path from the root to segment k (inclusive), for the first
# length(D) segments of a preorder tree with subtree sizes `sz` and depths `dep`; E is the
# Euler-tour buffer (2·length(D)) and `bnd` a device scalar for the bound Σ|w|.
function _pathsum!(D, w, bnd, part, E, sz, dep)
    _reduce!(bnd, part, +, 0.0f0, length(w), i -> abs(w[i]))
    AK.foreachindex(D) do k
        @inbounds begin
            v = _fix(w[k], _scale(bnd[1])); e = 2k - dep[k]
            E[e] = v; E[e + 2sz[k] - 1] = -v
        end
    end
    AK.accumulate!(+, E)
    AK.foreachindex(D) do k
        @inbounds begin
            D[k] = _unfix(E[2k - dep[k]], _scale(bnd[1]))
        end
    end
    return D
end

# Resistances of the current diameters, the bounds for the fixed-point scale of the operator,
# and the factors of the preconditioner. (Keeping the factors of an earlier solve for a few solves
# was tried: at rest the arterioles move enough between steps that CG then needs ~25× more
# iterations.)
function _resistances!(S::_Solver, kappa::F32)
    ni = S.ni; rr = S.rr; ma = S.ma; ra = S.ra; fa = S.fa
    G = S.G; lc = S.lc; gi = S.gi; rcap = S.rcap; cptr = S.cptr; ckid = S.ckid; sz = S.sz
    # resistances; and at the terminals the start of the Norton reduction below
    AK.foreachindex(rr) do k
        @inbounds begin
            r = (ma[k] ? kappa * ra[k] : ra[k]) * fa[k]; rr[k] = r
            if k > ni
                g = 1.0f0 / (kappa * rcap[k - ni])
                G[k] = g / (1.0f0 + g * r); lc[k] = -log1p(g * r)
            end
        end
    end
    if S.kv[] != kappa
        rvv = S.rvv; mv = S.mv; rv = S.rv
        AK.foreachindex(rvv) do k
            @inbounds begin
                rvv[k] = mv[k] ? kappa * rv[k] : rv[k]
            end
        end
        _pathsum!(S.Dv, rvv, S.s_b, S.part, S.Ev, S.vsz, S.vdep)
        Dv = S.Dv; _reduce!(S.s_Rv, S.part, max, 0.0f0, S.nv, k -> Dv[k])
        S.kv[] = kappa
    end
    # |path sum of r·Q| ≤ Σ|x| · Σ r over the segments with children (a loose bound costs bits,
    # not accuracy: the scale keeps ~2⁶⁰ quanta for it)
    _reduce!(S.s_Ra, S.part, +, 0.0f0, ni, k -> rr[k])
    # Norton reduction of the arterial tree with every terminal draining through its bed to 0:
    # conductance below segment k, g, and through it, G = g/(1 + g·r). It is nonlinear, so it runs
    # bottom-up: every small subtree on one thread (children follow their parent in preorder),
    # then the few larger segments above them level by level.
    lpi = S.lpi
    @inline function norton!(k)
        @inbounds begin
            g = 0.0f0
            for c in cptr[k]:(cptr[k + 1] - 1)
                g += G[ckid[c]]
            end
            G[k] = g / (1.0f0 + g * rr[k]); gi[k] = g
            l = -log1p(g * rr[k]); lc[k] = l; lpi[k] = l
        end
        return nothing
    end
    croot = S.croot
    AK.foreachindex(croot; min_elems = 64) do c
        @inbounds begin
            k0 = croot[c]
            for k in (k0 + sz[k0] - Int32(1)):Int32(-1):k0
                norton!(k)
            end
        end
    end
    lev = S.lev; lord = S.lord
    for l in (length(lev) - 1):-1:1
        lo = lev[l]
        AK.foreachindex(1:(lev[l + 1] - lo); backend = S.be, min_elems = 2048) do h
            @inbounds begin
                norton!(lord[lo + h - 1])
            end
        end
    end
    # π[k]: product of c = 1/(1 + g·r) from the root down to k, the pressure at the distal end of k
    # as a fraction of the inlet pressure when the beds drain to 0, and F = π·g, the flow through
    # k in that problem. Both span many decades when beds are narrow and arteries wide, so they
    # are kept as logarithms and the preconditioner works in bands of them.
    lF = S.lF; lFt = S.lFt; par = S.par
    _pathsum!(lpi, lpi, S.s_b, S.part, S.Ep, S.sz, S.dep)        # log π
    AK.foreachindex(rr) do k
        @inbounds begin
            if k <= ni
                lF[k] = lpi[k] + log(gi[k])
            else
                p = par[k]; i = k - ni
                lFt[i] = (p == 0 ? 0.0f0 : lpi[p]) + lc[k] - log(kappa * rcap[i])
            end
        end
    end
    _reduce!(S.s_lF, S.part, max, -Inf32, ni, k -> lF[k])
    _reduce!(S.s_lFm, S.part, min, Inf32, ni, k -> lF[k])
    _reduce!(S.s_lpm, S.part, max, 0.0f0, ni, k -> -lpi[k])
    s_lF = S.s_lF
    _reduce!(S.s_w0, S.part, +, 0.0f0, S.nt, i -> exp(min(lFt[i] - s_lF[1], 1.0f0)))
    span = _get(s_lF) - _get(S.s_lFm)
    S.nbu[] = span <= _BAND_UP ? 1 : 2 + floor(Int, (span - _BAND_UP) / _BAND_UP2)
    S.nbd[] = max(1, ceil(Int, _get(S.s_lpm) / _BAND_DOWN))
    return S
end

# Terminal-path operator: out_i = Σ_j∈path(i) r_j Q_j (+ bed resistance × x_i for `bed`), where Q_j
# is the sum of x over the terminals below segment j; arterial path from the root to terminal i
# plus venous path from the sinus to the vein draining i. With `final`, the flows and pressures of
# every segment and vein are stored.
function _apply!(out, S::_Solver, kappa::F32, x, bed::Bool, final::Bool = false)
    ni = S.ni; nt = S.nt; n = S.n; nv = S.nv
    s_x = S.s_x; s_Ra = S.s_Ra; s_Rv = S.s_Rv
    _reduce!(s_x, S.part, +, 0.0f0, nt, i -> abs(x[i]))
    X = S.X; vsorted = S.vsorted; vcum = S.vcum
    AK.foreachindex(1:(nt + nv); backend = S.be) do i
        @inbounds begin
            σ = _scale(s_x[1])
            if i <= nt
                X[i + 1] = _fix(x[i], σ)
            else
                k = i - nt; s = 0.0f0
                for m in (vcum[k] + 1):vcum[k + 1]
                    s += x[vsorted[m]]
                end
                X[i + 1] = _fix(s, σ)
            end
        end
    end
    AK.accumulate!(+, X)
    # flow through every segment with children and every vein, times its resistance, on the Euler tours
    Eu = S.Eu; tb = S.tb; te = S.te; sz = S.sz; dep = S.dep; rr = S.rr; Qa = S.Qa
    vsz = S.vsz; vdep = S.vdep; rvv = S.rvv; Qv = S.Qv
    AK.foreachindex(1:(ni + nv); backend = S.be) do k
        @inbounds begin
            σ = _scale(s_x[1]); σw = _scale(s_x[1], s_Ra[1] + s_Rv[1])
            if k <= ni
                Q = _unfix(X[te[k] + 1] - X[tb[k] + 1], σ); Qa[k] = Q
                v = _fix(rr[k] * Q, σw); e = 2k - dep[k]
                Eu[e] = v; Eu[e + 2sz[k] - 1] = -v
            else
                m = k - ni
                Q = _unfix(X[nt + m + vsz[m]] - X[nt + m], σ); Qv[m] = Q
                v = _fix(rvv[m] * Q, σw); e = 2ni + 2m - vdep[m]
                Eu[e] = v; Eu[e + 2vsz[m] - 1] = -v
            end
        end
    end
    AK.accumulate!(+, Eu)
    par = S.par; vdrain = S.vdrain; rcap = S.rcap
    AK.foreachindex(out) do i
        @inbounds begin
            σw = _scale(s_x[1], s_Ra[1] + s_Rv[1])
            L = ni + i; p = par[L]; m = vdrain[i]
            D = (p == 0 ? Int64(0) : Eu[2p - dep[p]]) + Eu[2ni + 2m - vdep[m]]
            out[i] = _unfix(D, σw) + rr[L] * x[i] + (bed ? kappa * rcap[i] * x[i] : 0.0f0)
        end
    end
    if final
        q = S.q; pa = S.pa; qv = S.qv; pv = S.pv; Pa = S.Pa; Pb = S.Pb
        AK.foreachindex(1:(n + nv); backend = S.be) do k
            @inbounds begin
                σw = _scale(s_x[1], s_Ra[1] + s_Rv[1])
                if k <= ni
                    q[k] = Qa[k]; pa[k] = Pa - _unfix(Eu[2k - dep[k]], σw)
                elseif k <= n
                    p = par[k]; i = k - ni
                    Dp = p == 0 ? 0.0f0 : _unfix(Eu[2p - dep[p]], σw)
                    q[k] = x[i]; pa[k] = Pa - (Dp + rr[k] * x[i])
                else
                    m = k - n
                    qv[m] = Qv[m]; pv[m] = Pb + _unfix(Eu[2ni + 2m - vdep[m]], σw)
                end
            end
        end
    end
    return out
end

# Preconditioner z = M_a⁻¹ r: the arterial tree with its root at 0 and terminal i draining
# through its bed into pressure −r_i (exact Norton solve). With the conductances G fixed by
# `_resistances!`, the two tree passes are linear in r and become scans:
# - upward, the Norton source of segment k is the flow-weighted mean E_k = Σ_{i below k} F_i (−r_i) / F_k;
# - downward, the pressure at its distal end is P_k = Σ_{a on the path to k} (π_k/π_a)(1 − c_a) E_a.
# F and π can span many decades, more than one fixed-point scale resolves, so each pass runs in
# bands of its log weight. Upward, a band of segments (log F within `_BAND_UP` below its top) scans
# the terminals weighted by F/F_top ≲ 1, leaving out the terminals above the band (no subtree of
# the band holds them) and those `_FAR` below it (they carry < e^-_FAR of any band segment's
# flow). Downward, a band (log π within `_BAND_DOWN` below its top) weighs ancestor a by
# (1 − c_a) π_top/π_a = G_a r_a π_top/π_a: these weights stay below e^_BAND_DOWN and so does their
# sum along any path; ancestors more than `_FAR` above the band add < e^-_FAR·max|r| and are left
# out. With the
# usual dynamic range each pass is a single band.
const _BAND_UP = 20.0f0      # first upward band; the further ones are `_BAND_UP2` wide
const _BAND_UP2 = 8.0f0
const _BAND_DOWN = 20.0f0
const _FAR = 40.0f0
# upward band of a segment from log(F/F_max) ≤ 0, and the top of band b
@inline _bandup(λ0::F32, nb::Int) = λ0 > -_BAND_UP ? 0 : clamp(1 + floor(Int, (-_BAND_UP - λ0) / _BAND_UP2), 0, nb - 1)
@inline _topup(b::Int) = b == 0 ? 0.0f0 : -_BAND_UP - (b - 1) * _BAND_UP2
function _precond!(z, S::_Solver, r)
    ni = S.ni
    s_r = S.s_rn; s_lF = S.s_lF          # max|r| ≤ ‖r‖ = √(r·r), which the caller formed
    Xp = S.Xp; lFt = S.lFt; lF = S.lF; tb = S.tb; te = S.te; Ek = S.Ek; s_w0 = S.s_w0
    nbu = S.nbu[]; nbd = S.nbd[]; one = nbu == 1      # one upward band: E is formed where it is used
    for b in 0:(nbu - 1)
        # weights F_i/F_top: their sum is below 3 Σ F_i/F_max (first band, about 3 × the number of
        # roots) or 3·(number of terminals)
        top = _topup(b); wid = b == 0 ? _BAND_UP : _BAND_UP2
        bnd = b == 0 ? s_w0 : S.s_wn; bf = b == 0 ? 3.0f0 : 1.0f0
        AK.foreachindex(r) do i
            @inbounds begin
                λ = lFt[i] - s_lF[1] - top
                w = -wid - 1 - _FAR < λ <= 1 ? exp(λ) * -r[i] : 0.0f0
                Xp[i + 1] = _fix(w, _scale(sqrt(s_r[1]), bf * bnd[1]))
            end
        end
        AK.accumulate!(+, Xp)
        one && break
        AK.foreachindex(Ek) do k
            @inbounds begin
                λ0 = lF[k] - s_lF[1]
                if _bandup(λ0, nbu) == b
                    Ek[k] = _unfix(Xp[te[k] + 1] - Xp[tb[k] + 1], _scale(sqrt(s_r[1]), bf * bnd[1])) * exp(top - λ0)
                end
            end
        end
    end
    lpi = S.lpi; G = S.G; rr = S.rr; Pk = S.Pk; Ep = S.Ep; dep = S.dep; sz = S.sz
    eB = exp(_BAND_DOWN + 1)
    for b in 0:(nbd - 1)
        top = -b * _BAND_DOWN
        AK.foreachindex(Ek) do a
            @inbounds begin
                l = lpi[a]
                E = one ? _unfix(Xp[te[a] + 1] - Xp[tb[a] + 1], _scale(sqrt(s_r[1]), 3.0f0 * s_w0[1])) * exp(s_lF[1] - lF[a]) : Ek[a]
                w = top - _BAND_DOWN - 1 < l <= top + _FAR ? G[a] * rr[a] * exp(top - l) * E : 0.0f0
                v = _fix(w, _scale(sqrt(s_r[1]), eB)); e = 2a - dep[a]
                Ep[e] = v; Ep[e + 2sz[a] - 1] = -v
            end
        end
        AK.accumulate!(+, Ep)
        nbd == 1 && break
        AK.foreachindex(Pk) do k
            @inbounds begin
                l = lpi[k]
                if clamp(floor(Int, -l / _BAND_DOWN), 0, nbd - 1) == b
                    Pk[k] = exp(l - top) * _unfix(Ep[2k - dep[k]], _scale(sqrt(s_r[1]), eB))
                end
            end
        end
    end
    par = S.par; oned = nbd == 1                      # one downward band: P is read off the scan
    AK.foreachindex(z) do i
        @inbounds begin
            p = par[ni + i]
            P = p == 0 ? 0.0f0 : oned ? exp(lpi[p]) * _unfix(Ep[2p - dep[p]], _scale(sqrt(s_r[1]), eB)) : Pk[p]
            z[i] = G[ni + i] * (P + r[i])
        end
    end
    return z
end

# Venous path operator alone: out_i = Σ_k∈path(i) r_k Q_k from the sinus to the vein draining
# terminal i, Q_k the sum of x over the terminals draining into the venous subtree of k.
function _venous!(out, S::_Solver, x)
    nt = S.nt
    s_x = S.s_x; s_Rv = S.s_Rv; Xv = S.Xv; Ev = S.Ev
    vsorted = S.vsorted; vcum = S.vcum; vsz = S.vsz; vdep = S.vdep; rvv = S.rvv; vdrain = S.vdrain
    _reduce!(s_x, S.part, +, 0.0f0, nt, i -> abs(x[i]))
    AK.foreachindex(rvv) do k
        @inbounds begin
            s = 0.0f0
            for m in (vcum[k] + 1):vcum[k + 1]
                s += x[vsorted[m]]
            end
            Xv[k + 1] = _fix(s, _scale(s_x[1]))
        end
    end
    AK.accumulate!(+, Xv)
    AK.foreachindex(rvv) do k
        @inbounds begin
            Q = _unfix(Xv[k + vsz[k]] - Xv[k], _scale(s_x[1]))
            v = _fix(rvv[k] * Q, _scale(s_x[1], s_Rv[1])); e = 2k - vdep[k]
            Ev[e] = v; Ev[e + 2vsz[k] - 1] = -v
        end
    end
    AK.accumulate!(+, Ev)
    AK.foreachindex(out) do i
        @inbounds begin
            m = vdrain[i]
            out[i] = _unfix(Ev[2m - vdep[m]], _scale(s_x[1], s_Rv[1]))
        end
    end
    return out
end

# One steady solve for the current resistances.
#
# Unknowns are the terminal flows q. Arterial pressure drop to terminal i plus its own bed is
# (M_a q)_i (shared path resistances + bed), the venous rise from the sinus is (M_v q)_i, and
# continuity at every bed is (M_a + M_v) q = P_ao − P_back. Both matrices are symmetric positive
# definite and applied in O(N) by tree passes, so the system is solved by conjugate gradients
# preconditioned with the exact arterial tree solve (M_a⁻¹ is the Norton solve of the arterial
# tree); its condition number is 1 + λ(M_a⁻¹M_v), a few, so it converges in tens of iterations.
# As z = M_a⁻¹ r, the search directions p = z + β p satisfy M_a p = r + β M_a p, so an iteration
# needs the preconditioner and the venous operator only; the full operator runs once for the
# starting residual and once at the end for the flows and pressures of every segment.
# `warm` starts from the previous terminal flows.
function _solve!(S::_Solver, kappa::F32; warm::Bool = false, tol::F32 = 1.0f-6)
    nt = S.nt; dP = S.Pa - S.Pb
    _resistances!(S, kappa)
    x = S.x; r = S.r; z = S.z; p = S.p; Ap = S.Ap; Mp = S.Mp
    s_rz = S.s_rz; s_rz2 = S.s_rz2; s_pap = S.s_pap; s_rn = S.s_rn
    if warm
        copyto!(x, S.qt)
    else                                                # cold start: arteries alone
        fill!(r, dP); _dot!(s_rn, S, r, r); _precond!(x, S, r)
    end
    _apply!(Ap, S, kappa, x, true)
    AK.foreachindex(r) do i
        @inbounds begin
            r[i] = dP - Ap[i]
        end
    end
    _dot!(s_rn, S, r, r)
    bn = sqrt(F32(nt)) * dP
    if sqrt(_get(s_rn)) >= tol * bn
        _precond!(z, S, r); copyto!(p, z); copyto!(Mp, r)
        _dot!(s_rz, S, r, z)
        for it in 1:500
            S.iters[] += 1
            _venous!(Ap, S, p)                          # A p = M_a p + M_v p
            _reduce!(s_pap, S.part, +, 0.0f0, nt, i -> p[i] * (Mp[i] + Ap[i]))
            AK.foreachindex(x) do i
                @inbounds begin
                    α = s_rz[1] / s_pap[1]
                    x[i] += α * p[i]; r[i] -= α * (Mp[i] + Ap[i])
                end
            end
            _dot!(s_rn, S, r, r)
            (sqrt(_get(s_rn)) < tol * bn || it == 500) && break
            _precond!(z, S, r)
            _dot!(s_rz2, S, r, z)
            AK.foreachindex(p) do i
                @inbounds begin
                    β = s_rz2[1] / s_rz[1]
                    p[i] = z[i] + β * p[i]; Mp[i] = r[i] + β * Mp[i]
                end
            end
            copyto!(s_rz, s_rz2)
        end
    end
    copyto!(S.qt, x)
    _apply!(Ap, S, kappa, x, false, true)               # flows and pressures of every segment
    S.solved[] = true
    return S
end

# Microvascular resistance as Chilian et al. (1989) measured it: flow-weighted mean pressure where
# arterioles leave the larger arteries minus that where venules join the larger veins, divided by
# the flow per gram (mmHg·min·g/mL).
function _micro_resistance(S::_Solver, mass_g)
    q = Array(S.q); pa = Array(S.pa); qv = Array(S.qv); pv = Array(S.pv)
    par = S.hpar; ma = S.hma; vpar = S.hvpar; mv = S.hmv
    sa = 0.0f0; wa = 0.0f0; sv = 0.0f0; wv = 0.0f0
    @inbounds for j in eachindex(par)
        p = par[j]
        if ma[j] && (p == 0 || !ma[p])
            pin = p == 0 ? S.Pa : pa[p]; sa += q[j] * pin; wa += q[j]
        end
    end
    @inbounds for k in eachindex(vpar)
        p = vpar[k]
        if mv[k] && (p == 0 || !mv[p])
            pout = p == 0 ? S.Pb : pv[p]; sv += qv[k] * pout; wv += qv[k]
        end
    end
    qtot = AK.reduce(+, S.qt; init = 0.0f0)
    return (sa / wa - sv / wv) / (60.0f0 * qtot / mass_g), sa / wa, sv / wv
end


# ---------------------------------------------------------------------------------------------
# Arteriolar wall mechanics and the metabolic signal
# ---------------------------------------------------------------------------------------------

# Vasoactive arteries. Wall constants follow D₀ (µm; tensions dyn/cm): the Carlson & Secomb
# (2005) fits as tabulated by Fry et al. (2013, Table 1).
struct _Wall{VI <: AbstractVector{Int32}, VF <: AbstractVector{F32}}
    seg::VI                      # arterial segment (solver numbering)
    d0::VF                       # passive diameter at 100 mmHg, µm
    u::VF                        # its coordinate on the D₀ grid of the `_Limits` tables
    xd::VF                       # D/D₀
    act::VF                      # activation A
end

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
    AK.foreachindex(1:nd; min_elems = 1) do k
        d0 = exp(ld0 + (k - 1) * dld); xh = 1.0f0; xl = 0.5f0
        for i in 1:np
            P = 1.0f0 + (i - 1) * dp
            xh = _equilibrium(xh, 0.0f0, P, d0); xl = _equilibrium(xl, 1.0f0, P, d0)
            hi[i, k] = xh; lo[i, k] = xl
        end
    end
    return _Limits(ld0, dld, nd, dp, np, hi, lo)
end
function _Wall(seg, d0::Vector{F32}, L::_Limits, to = identity)
    n = length(d0); u = F32[clamp((log(d) - L.ld0) / L.dld, 0.0f0, F32(L.nd - 1) - 1.0f-4) for d in d0]
    return _Wall(to(seg), to(d0), to(u), to(ones(F32, n)), to(zeros(F32, n)))
end
# (tables `hi`, `lo` of a `_Limits`, passed apart so that they can live on a device; `u` is the
# vessel's coordinate on the D₀ grid)
@inline function _limits(hi, lo, u::F32, dp::F32, np::Int, P::F32)
    k = floor(Int, u); fu = u - k
    v = clamp((P - 1.0f0) / dp, 0.0f0, F32(np - 1) - 1.0f-4); i = floor(Int, v); fv = v - i
    @inbounds begin
        h = (1 - fu) * ((1 - fv) * hi[i + 1, k + 1] + fv * hi[i + 2, k + 1]) + fu * ((1 - fv) * hi[i + 1, k + 2] + fv * hi[i + 2, k + 2])
        l = (1 - fu) * ((1 - fv) * lo[i + 1, k + 1] + fv * lo[i + 2, k + 1]) + fu * ((1 - fv) * lo[i + 1, k + 2] + fv * lo[i + 2, k + 2])
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
#
# With path lengths ℓ from the root, the decayed sum over the beds below segment k is
# e^{ℓ_k/L₀} Σ_i s_i e^{−ℓ_i/L₀}, a subtree sum. Path lengths span ~20 L₀, more than one
# fixed-point scale holds, so the sums run in bands of 12 L₀ of ℓ_k: a band's scan weighs bed i by
# e^{−(ℓ_i − ℓ_band)/L₀} ≤ 1 and leaves out the beds above the band, which no subtree of the band
# contains. The weights are fixed (`_solver`).
function _metabolic!(S::_Solver, phys::Physiology, terr_g)
    ni = S.ni; nt = S.nt
    cO2 = 1.34f0 * phys.hb / 100.0f0; r1 = phys.atp_r1
    rref = 1.0f0 - r1 * phys.sao2 * (1 - phys.oef)
    sao2 = phys.sao2; mvo2 = phys.mvo2
    qt = S.qt; so2 = S.so2; val = S.val; sig = S.sig; mlf = S.mlf
    AK.foreachindex(so2) do i
        @inbounds begin
            qi = max(qt[i], 1.0f-12)
            s = clamp(sao2 - mvo2 * terr_g[i] / 60.0f0 / (qi * cO2), 0.0f0, sao2)
            so2[i] = s
            v = (1.0f0 - r1 * s) / rref - 1.0f0; val[i] = v
            sig[ni + i] = v * mlf[i]
        end
    end
    s_v = S.s_v; Xm = S.Xm; tb = S.tb; te = S.te; ntf = F32(nt); nb = S.nband
    mbw = S.mbw; mnf = S.mnf; mband = S.mband
    _reduce!(s_v, S.part, max, 0.0f0, nt, i -> abs(val[i]))
    # all bands in one scan, one after the other (each band's sums are differences within it)
    AK.foreachindex(val) do i
        @inbounds begin
            σ = _scale(s_v[1], ntf)
            for b in 0:(nb - 1)
                Xm[b * nt + i + 1] = _fix(val[i] * mbw[b * nt + i], σ)
            end
        end
    end
    AK.accumulate!(+, Xm)
    AK.foreachindex(mnf) do k
        @inbounds begin
            o = mband[k] * nt + 1
            sig[k] = _unfix(Xm[o + te[k]] - Xm[o + tb[k]], _scale(s_v[1], ntf)) * mnf[k]
        end
    end
    return S
end

# Autoregulated (or maximally dilated) steady state: flows and diameters solved together.
# Every vasoactive artery's diameter lies between its relaxed (A = 0) and fully contracted
# (A = 1) wall equilibria at its own pressure, D = D_A=1 + (D_A=0 − D_A=1)·(1 − A). Activation
# A = 1/(1 + e^{−S}) with S = C_myo·P·D/2 − `g_meta`·(metabolic signal) + C_tone relaxes toward
# its stimulus (τ_a ≫ τ_d, Carlson et al. 2008); the network is re-solved after every step.
# Healthy reference (`fixtone = false`): C_tone is integrated until the coronary-sinus saturation
# equals the measured resting extraction (`phys.oef`).
# Sheet-flow capillary beds (Fung & Sobin): K per bed in mL/s per µm⁴, h₀ in µm, α in µm/mmHg,
# p_ext in mmHg. `_sheet!` sets each bed's resistance from its current end pressures.
struct _Sheet{VF <: AbstractVector{F32}}
    K::VF; h0::F32; al::F32; pe::F32; hmax::F32
end
function _sheet!(S::_Solver, sh::_Sheet, kappa::F32)
    ni = S.ni; pa = S.pa; pv = S.pv; vd = S.vdrain; rcap = S.rcap; K = sh.K; h0 = sh.h0; al = sh.al; pe = sh.pe; hm = sh.hmax
    AK.foreachindex(rcap) do i
        @inbounds begin
            Pa = pa[ni + i]; Pv = pv[vd[i]]
            ha = clamp(h0 + al * (Pa - pe), 0.0f0, hm); hv = clamp(h0 + al * (Pv - pe), 0.0f0, hm)
            q = K[i] * (ha^4 - hv^4)
            r = Pa - Pv > 1.0f-3 && q > 0 ? (Pa - Pv) / q : 1.0f0 / (4.0f0 * K[i] * al * max(ha, 1.0f-3)^3)
            rcap[i] = sqrt(rcap[i] * r / kappa)          # geometric relaxation toward the new resistance
        end
    end
    return nothing
end

function _autoregulate!(
        S::_Solver, W::_Wall, kappa::F32, phys::Physiology, state::Symbol, ctone::Vector{F32}, fixtone::Bool,
        terr_g, mu0, sten, lim::_Limits, verbose, sheet = nothing
    )
    nw = length(W.seg); hd = phys.hd; stress = state === :stress
    sref = phys.sao2 * (1 - phys.oef); gain = 2.0f0; lasterr = 0.0f0
    # Activation relaxes toward its stimulus with a step per vessel: the loop through a terminal
    # arteriole's own bed has a gain of ~20, so a fixed explicit step either crawls or makes some
    # vessels flip back and forth. A vessel's step halves whenever its activation reverses
    # direction and grows by 10% otherwise (bounded 0.005–0.3).
    seg = W.seg; d0 = W.d0; uw = W.u; xd = W.xd; act = W.act
    ωw = similar(xd); fill!(ωw, 0.1f0); dA = similar(xd); fill!(dA, 0.0f0)
    stim = similar(xd); xold = similar(xd)
    fa = S.fa; ra = S.ra; q = S.q; par = S.par; pa = S.pa; Pa = S.Pa; sig = S.sig
    lhi = S.lhi; llo = S.llo; dp = lim.dp; np = lim.np
    mut = S.mut; mua = S.mua; mudl = S.mudl; numax = F32(length(mut) - 1) - 1.0f-3
    gmeta = phys.g_meta
    function setfa!()
        # μ_rel(D) interpolated in log D (relative error < 10⁻⁷)
        AK.foreachindex(xd) do w
            @inbounds begin
                j = seg[w]; u = clamp((log(d0[w] * xd[w]) - mua) / mudl, 0.0f0, numax)
                m = unsafe_trunc(Int, u); f = u - m
                fa[j] = (mut[m + 1] + f * (mut[m + 2] - mut[m + 1])) / mu0[j] / xd[w]^4
            end
        end
        # stenoses: unchanged lumen outside the narrowing plus the f·Q + s·Q² loss of the throat
        if length(sten) > 0
            AK.foreachindex(sten) do m
                @inbounds begin
                    j, frac, f, s2 = sten[m]
                    fa[j] = (1 - frac) + (f + s2 * abs(q[j])) / ra[j]
                end
            end
        end
        return nothing
    end
    qsum() = AK.reduce(+, S.qt; init = 0.0f0)
    k99 = max(1, ceil(Int, 0.99 * nw))
    # 99th percentile of the last step |ΔD/D₀| (exact, on the host)
    p99() = nw > 0 ? partialsort!(abs.(Array(xd) .- Array(xold)), k99) : 0.0f0
    setfa!(); _solve!(S, kappa; warm = S.solved[])   # warm when called again (κ loop)
    qprev = qsum(); iters = 0; dxok = true; ssin = 0.0f0; qhist = fill(qprev, 10); dq1 = 1.0f0
    maxit = stress ? 60 : 600
    for it in 1:maxit
        iters = it
        if !stress
            _metabolic!(S, phys, terr_g)
            ssin = AK.mapreduce((a, b) -> max(a, 0.0f0) * b, +, S.qt, S.so2; init = 0.0f0) /
                AK.mapreduce(a -> max(a, 0.0f0), +, S.qt; init = 0.0f0)
        end
        if !stress && !fixtone
            if it == 1                                           # start from half activation on average
                AK.foreachindex(stim) do w
                    @inbounds begin
                        j = seg[w]; p = par[j]
                        P = max(((p == 0 ? Pa : pa[p]) + pa[j]) / 2, 1.0f0)
                        stim[w] = 0.09128f0 * P * xd[w] - gmeta * sig[j]
                    end
                end
                ctone[1] = -sort!(Array(stim))[(nw + 1) ÷ 2]
            else
                err = ssin - sref                                # too much flow → more tone
                err * lasterr < 0 && (gain *= 0.7f0)
                ctone[1] += clamp(gain * err, -0.5f0, 0.5f0); lasterr = err
            end
        end
        copyto!(xold, xd)
        ct = ctone[1]
        AK.foreachindex(xd) do w
            @inbounds begin
                j = seg[w]; p = par[j]
                P = max(((p == 0 ? Pa : pa[p]) + pa[j]) / 2, 1.0f0)
                if stress
                    act[w] = 0.0f0
                else
                    sw = 0.09128f0 * P * xd[w] - gmeta * sig[j]
                    d = 1.0f0 / (1.0f0 + exp(-(sw + ct))) - act[w]
                    ωw[w] = d * dA[w] < 0 ? max(0.005f0, 0.5f0 * ωw[w]) : min(0.3f0, 1.1f0 * ωw[w])
                    dA[w] = d; act[w] += ωw[w] * d
                end
                h, l = _limits(lhi, llo, uw[w], dp, np, P)
                xd[w] = l + (h - l) * (1.0f0 - act[w])
            end
        end
        # a few arterioles may keep cycling (vasomotion), so convergence is judged on total flow
        # (and on the diameters' 99th-percentile step at maximal dilation, which has no cycling):
        # the 99th percentile is below 10⁻³ when at least 99% of the steps are
        stress && (dxok = nw == 0 || AK.mapreduce((a, b) -> abs(a - b) < 1.0f-3 ? 1 : 0, +, xd, xold; init = 0) >= k99)
        # inexact network solves while the diameters are still moving; exact once they settle
        sheet === nothing || _sheet!(S, sheet, kappa)
        setfa!(); _solve!(S, kappa; warm = true, tol = clamp(0.1f0 * dq1, 1.0f-6, 1.0f-3))
        qt = qsum(); dq = abs(qt - qhist[mod1(it, 10)]) / qt; qhist[mod1(it, 10)] = qt
        dq1 = abs(qt - qprev) / qt; qprev = qt
        verbose && (it % 10 == 0 || it < 5) && (println("  autoregulation $it: flow $(round(60qt, digits = 3)) mL/min, sinus SO₂ $(round(ssin, digits = 3)), C_tone $(round(ctone[1], digits = 3)), p99 |ΔD/D₀| $(round(p99(), sigdigits = 2)), mean A $(round(AK.reduce(+, act; init = 0.0f0) / max(nw, 1), digits = 3))"); flush(stdout))
        # maximal dilation has no vasomotion: consecutive steps decide; at rest, ten steps
        if stress
            it > 2 && dq1 < 1.0f-5 && dxok && break
        else
            it > 10 && dq < 1.0f-4 && (fixtone || abs(ssin - sref) < 1.0f-3) && break
        end
    end
    _solve!(S, kappa; warm = true)
    stress || _metabolic!(S, phys, terr_g)
    return iters, p99()
end

function _circuit(t::VTree, territory_ml::AbstractVector; phys::Physiology, state::Symbol, veins, vleaf, reference, kappa, kappa_init, stenosis, demand, verbose, to)
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
        vptr = ptr; vkids = kids; vorder = order; vlen = len; vvol = (s * s) .* vol
    else
        length(vleaf) == nt || throw(ArgumentError("vleaf needs one venous segment per terminal"))
        vl = Int32.(vleaf)
        # A grown venous tree is sized by CCO's arterial law. With venous_scale = :volume_share its
        # radii are scaled uniformly so it holds share_ven/share_art of the resolved arterial volume
        # (Kassab et al. 1994), as companion veins do; with :grown its own (anatomically anchored)
        # radii stay.
        vlen, vvol = _lengths(veins)
        fr = phys.venous_scale === :grown ? 1.0f0 : sqrt(s * s * vart / sum(vvol))
        vt = VTree(veins.parent, veins.a, veins.b, veins.r .* F32(fr), veins.isbase)
        _, vptr, vkids, vorder = _topology(vt.parent); vlen, vvol = _lengths(vt)
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
    # vasoactive arteries
    d0c = F32[2000.0f0 * t.r[j] for j in 1:n if 2000.0f0 * t.r[j] < phys.d_active]
    lim = _Limits(isempty(d0c) ? 1.0f0 : minimum(d0c), isempty(d0c) ? 2.0f0 : maximum(d0c), max(phys.p_in, 1.0f0))
    S, perm, iperm, tperm, vperm = _solver(
        t.parent, nchild, ptr, kids, vt.parent, vptr, vkids, vl, tidx, ra, rv, rcap, ma, mv, len, lim,
        phys.p_in, phys.p_out, 10.0f0 * phys.atp_l0, hd, to
    )
    sheet = nothing
    if phys.capillary === :sheet
        al = phys.sheet_alpha_um_cmH2O * 1.35951f0                   # µm per mmHg
        alsi = phys.sheet_alpha_um_cmH2O * 1.0f-6 / 98.0665f0          # m per Pa
        k1 = phys.sheet_s / (4.0f0 * phys.sheet_mu_cP * 1.0f-3 * phys.sheet_kf * (phys.sheet_path_um * 1.0f-6)^2 * alsi)
        Kt = F32[k1 * phys.sheet_area_m2_per_ml * F32(territory_ml[i]) * 1.0f-18 for i in 1:nt]   # mL/s per µm⁴
        h = clamp(phys.sheet_h0_um + al * ((phys.p_in + phys.p_out) / 2 - phys.p_ext), 0.1f0, phys.sheet_h0_um + al * phys.sheet_pmax_cmH2O / 1.35951f0)
        copyto!(S.rcap, to(F32[1.0f0 / (4.0f0 * Kt[i] * al * h^3) for i in tperm]))
        sheet = _Sheet(to(Kt[tperm]), phys.sheet_h0_um, al, phys.p_ext, phys.sheet_h0_um + al * phys.sheet_pmax_cmH2O / 1.35951f0)
    end
    seg = Int32[k for k in 1:n if 2000.0f0 * t.r[perm[k]] < phys.d_active]   # solver numbering
    act = perm[seg]                                                          # caller numbering
    d0 = F32[2000.0f0 * t.r[j] for j in act]
    W = _Wall(seg, d0, lim, to)
    if reference !== nothing
        copyto!(W.xd, to(reference.xd[act])); copyto!(W.act, to(reference.act[act]))
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
        push!(sten, (iperm[j], Ls / len[j], f, s2))
    end
    tg = to(terr_g[tperm]); mu0d = to(mu0[perm]); stend = to(sten)

    # ---- microvascular scale κ at maximal dilation (vessels on their passive curves): the
    # Chilian-style resistance between arteriolar inlets and venular outlets equals the measurement
    ctone = reference === nothing ? F32[-5.0f0] : copy(reference.ctone)
    kappa === nothing && reference !== nothing && (kappa = reference.resist.kappa)
    kappa === nothing && phys.r_micro_dilated <= 0 && (kappa = 1.0f0)       # no anchor: pure Poiseuille
    if kappa === nothing
        kappa = F32(kappa_init)
        Wd = _Wall(seg, d0, lim, to)
        for _ in 1:12
            _autoregulate!(S, Wd, F32(kappa), phys, :stress, ctone, true, tg, mu0d, to(empty(sten)), lim, false, sheet)   # healthy geometry
            ratio = phys.r_micro_dilated / _micro_resistance(S, mass_g)[1]
            kappa *= ratio
            abs(ratio - 1) < 1.0f-3 && break
        end
    end
    kappa = F32(kappa)
    # ---- the state (a bed without active tone is always in its relaxed state)
    phys.regulation === :passive && (state = :stress)
    iters, dxlast = _autoregulate!(S, W, kappa, phys, state, ctone, reference !== nothing, tg, mu0d, stend, lim, verbose, sheet)
    sheet === nothing || (rcap[tperm] = Array(S.rcap))     # the sheet beds' final resistances
    rmc = _micro_resistance(S, mass_g)
    # back to the caller's numbering
    q = Vector{F32}(undef, n); pa = similar(q); rr = similar(q); qt = Vector{F32}(undef, nt); so2 = similar(qt)
    qv = Vector{F32}(undef, nv); pv = similar(qv)
    q[perm] = Array(S.q); pa[perm] = Array(S.pa); rr[perm] = Array(S.rr)
    qt[tperm] = Array(S.qt); so2[tperm] = Array(S.so2); qv[vperm] = Array(S.qv); pv[vperm] = Array(S.pv)
    nneg = count(<=(0), qt)
    nneg > 0 && @warn "$nneg terminals with reversed or zero flow (venous pressure above arterial); their beds are treated as stagnant"
    @inbounds for i in 1:nt
        qt[i] = max(qt[i], 1.0f-9)
    end
    wx = Array(W.xd); wa = Array(W.act)
    xd = ones(F32, n); av = zeros(F32, n)
    xd[act] = wx; av[act] = wa
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
    dpa = sum(q[j] * rr[j] * q[j] for j in 1:n)
    dpv = sum(qv[k]^2 * (mv[k] ? kappa : 1.0f0) * rv[k] for k in 1:nv); dpc = sum(qt .^ 2 .* rcap) * kappa
    roots = Int32[j for j in 1:n if t.parent[j] == 0]
    qtot = sum(q[roots])
    nw = length(act)
    resist = (;
        dp_total = phys.p_in - phys.p_out, art_share = dpa / (dpa + dpv + dpc),
        cap_share = dpc / (dpa + dpv + dpc), ven_share = dpv / (dpa + dpv + dpc),
        mbf = 60 * qtot / mass_g, root_flow = q[roots], kappa, iters, sweeps = S.iters[], residual_dx = dxlast,
        r_micro = rmc[1],
        chilian_shares = ((phys.p_in - rmc[2]), (rmc[2] - rmc[3]), (rmc[3] - phys.p_out)) ./ (phys.p_in - phys.p_out),
        so2_sinus = sum(qt .* so2) / sum(qt), starved = count(<(0.05f0), so2) / nt, mean_act = nw > 0 ? sum(wa) / nw : 0.0f0,
        n_dilated = count(<(0.01f0), wa), n_max_tone = count(>(0.99f0), wa),
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
        AK.foreachindex(1:m; backend = AK.get_backend(out), min_elems = 1024) do ii
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
