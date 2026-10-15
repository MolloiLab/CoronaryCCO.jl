# SPDX-License-Identifier: GPL-3.0-or-later
# Whole-body iodine circulation for a clinical IV injection: arm vein → SVC → RA → RV → pulmonary
# trunk → lungs → LA → LV → aorta → organs → veins → SVC / IVC → RA, with the coronary sinus
# draining into the RA. Iodine leaves the body only through the kidneys (GFR), so what leaves the
# coronary sinus recirculates and is diluted and cleared instead of accumulating.
#
# Design (sources in the parameter docstrings):
# - A physiologically based compartment model of the Bae, Heiken & Brink (1998) type, with the
#   compartment set of Sahbaee et al. (2016), who built theirs for XCAT phantoms. Heart chamber
#   volumes and the stroke volume come from the XCAT phantom itself.
# - Every compartment is well mixed and holds plasma (iodine does not enter red cells). Organs
#   exchange with their extracellular space through a permeability-surface product PS.
# - Pipes with plug flow (arm veins, lungs) are chains of small well-mixed stages, so a varying
#   flow (the injection stops, the saline chaser pushes the last contrast out) is handled exactly.
# - Explicit integration with a step well below every stage's mean transit time. Float32.

module Body

export Injection, BodyParams, simulate

const F32 = Float32

"""
    Injection(; conc, volume, rate, saline, saline_rate)

Contrast `conc` mgI/mL, `volume` mL at `rate` mL/s from t = 0, then `saline` mL at
`saline_rate`. Default: iohexol 350, 50 mL at 5 mL/s, 40 mL saline chaser (Lee et al. 2022
multicentre adenosine CT perfusion protocol; chaser per Mushtaq et al. 2025), right antecubital
vein (avoids reflux through the left brachiocephalic vein, Tseng et al. 2007).
"""
Base.@kwdef struct Injection
    conc::F32 = 350.0f0
    volume::F32 = 50.0f0
    rate::F32 = 5.0f0
    saline::F32 = 40.0f0
    saline_rate::F32 = 5.0f0
end

inj_rate(j::Injection, t) = t < j.volume / j.rate ? j.rate : t < j.volume / j.rate + j.saline / j.saline_rate ? j.saline_rate : 0.0f0
inj_conc(j::Injection, t) = t < j.volume / j.rate ? j.conc : 0.0f0

"""
    BodyParams(; kwargs...)

Blood volumes in mL, flows as fractions of cardiac output.

- `co`: cardiac output, mL/s. Default from XCAT vmale_50: stroke volume 132.2 − 50.6 = 81.6 mL
  at 60 beats/min = 81.6 mL/s (4.9 L/min; CT perfusion cohorts sit near 5 L/min).
- `hct`: large-vessel hematocrit 0.43 (ICRP 89 reference male).
- `v_ra`, `v_rv`, `v_la`, `v_lv`: XCAT chamber volumes at this frame (52.0, 85.4, 37.6, 95.2 mL).
- `v_pulm`: pulmonary blood 550 mL (ICRP 89: 10.5% of 5.3 L), split into the pulmonary trunk
  (`v_pa`) and `n_lung` mixing stages. Four stages give the relative dispersion of lung transit
  Knopp et al. (1976) quote for the dog lung (0.46 ≈ 1/√4), and the mean PA → LA transit
  (550 − v_pa)/CO ≈ 6 s matches MRI pulmonary transit 5.7 ± 1.4 s (Nelsson et al. 2021).
- organ flow fractions (Williams & Leggett 1989 / ICRP 89 adult male): brain 0.12, kidneys 0.19,
  hepatic artery 0.065, gut + spleen (portal) 0.185, muscle 0.17, skin 0.05, fat 0.05, bone 0.05,
  coronary from the coronary model, the rest "other".
- `upper`: fraction of the peripheral tissues (muscle, skin, fat, bone, other) draining to the SVC
  (arms, head and neck, upper thorax). With the brain, SVC flow is ≈ 0.28 of the venous return.
- organ blood volumes (ICRP 89, % of 5.3 L): brain 1.2, kidneys 2, liver 10, gut + spleen 4.7,
  muscle 14, skin 3, fat 5, bone 7; large arteries 6% and large veins 18% are transit volumes.
- extracellular volumes and PS (mL, mL/s): muscle 4000 / 7, skin 1300 / 1, fat 2000 / 2, bone
  1500 / 1, other 1500 / 3; gut, liver and kidney exchange is flow-limited (PS = 2 F_p); the
  brain has none (blood–brain barrier). Total ≈ 11 L of the 0.23–0.27 L/kg volume of
  distribution measured for iohexol (Olsson et al. 1983; Edelson et al. 1984). PS values are
  scaled from EDTA (Levitt 2003) and are the least certain numbers here; they shape only the
  recirculation and washout, not the first pass.
- `gfr`: 2.0 mL/s of plasma (120 mL/min; iohexol renal clearance 120 ± 19 mL/min, Edelson 1984).
- `arm_flow`: native flow of the injected arm vein, 0.02 CO. `v_arm`: 24 mL of basilic, axillary,
  subclavian and brachiocephalic vein (about 45 cm) as 30 stages (near plug flow), giving arm → SVC arrival ≈ 3 s at 4–5 mL/s
  (Borges et al. 2022: 3.0 ± 2.4 s).
- `beat`: cardiac cycle, 1 s (XCAT). Blood entering a ventricle in diastole leaves at the next
  systole, so each ventricle adds half a beat of pure delay before mixing.
- `lung_plug`: a quarter of the pulmonary transit is plug flow (no tracer leaves the lung faster
  than the shortest capillary path); the rest is `n_lung` mixing stages.
"""
Base.@kwdef struct BodyParams
    co::F32 = 81.6f0
    hct::F32 = 0.43f0
    v_ra::F32 = 52.0f0
    v_rv::F32 = 85.4f0
    v_la::F32 = 37.6f0
    v_lv::F32 = 95.2f0
    v_pa::F32 = 40.0f0
    v_pulm::F32 = 550.0f0
    n_lung::Int = 4
    v_aorta::F32 = 60.0f0
    v_svc::F32 = 15.0f0
    v_ivc::F32 = 25.0f0
    v_arm::F32 = 24.0f0
    n_arm::Int = 30
    beat::F32 = 1.0f0
    lung_plug::F32 = 0.25f0
    arm_flow::F32 = 0.02f0
    upper::F32 = 0.4f0
    f_brain::F32 = 0.12f0
    f_kidney::F32 = 0.19f0
    f_hepart::F32 = 0.065f0
    f_portal::F32 = 0.185f0
    f_muscle::F32 = 0.17f0
    f_skin::F32 = 0.05f0
    f_fat::F32 = 0.05f0
    f_bone::F32 = 0.05f0
    f_cor::F32 = 0.033f0
    v_brain::F32 = 64.0f0
    v_kidney::F32 = 106.0f0
    v_liver::F32 = 530.0f0
    v_portal::F32 = 250.0f0
    v_muscle::F32 = 742.0f0
    v_skin::F32 = 160.0f0
    v_fat::F32 = 265.0f0
    v_bone::F32 = 370.0f0
    v_other::F32 = 300.0f0
    v_cor::F32 = 16.0f0
    v_art::F32 = 320.0f0
    v_vein_upper::F32 = 300.0f0
    v_vein_lower::F32 = 650.0f0
    e_muscle::F32 = 4000.0f0; ps_muscle::F32 = 7.0f0
    e_skin::F32 = 1300.0f0; ps_skin::F32 = 1.0f0
    e_fat::F32 = 2000.0f0; ps_fat::F32 = 2.0f0
    e_bone::F32 = 1500.0f0; ps_bone::F32 = 1.0f0
    e_other::F32 = 1500.0f0; ps_other::F32 = 3.0f0
    e_portal::F32 = 600.0f0
    e_liver::F32 = 500.0f0
    e_kidney::F32 = 100.0f0
    gfr::F32 = 2.0f0
end

# A well-mixed plasma compartment, optionally exchanging with an extracellular space.
mutable struct Tank
    v::F32      # plasma volume, mL
    c::F32      # plasma concentration, mgI/mL
    ve::F32     # extracellular volume, mL (0 = none)
    ce::F32
    ps::F32     # mL/s
    cl::F32     # clearance from plasma, mL/s
end
Tank(v; ve = 0.0f0, ps = 0.0f0, cl = 0.0f0) = Tank(F32(v), 0.0f0, F32(ve), 0.0f0, F32(ps), F32(cl))

# Pure delay of `d` seconds on a concentration signal (ring buffer at the integration step).
mutable struct Delay
    buf::Vector{F32}
    k::Int
end
Delay(d, dt) = Delay(zeros(F32, max(1, round(Int, d / dt))), 1)
@inline function push!(D::Delay, c::F32)
    out = D.buf[D.k]; D.buf[D.k] = c; D.k = D.k == length(D.buf) ? 1 : D.k + 1
    return out
end

# advance one step with plasma inflow f (mL/s) at concentration cin
@inline function step!(T::Tank, f::F32, cin::F32, dt::F32)
    x = T.ps * (T.c - T.ce)
    T.c += dt * (f * (cin - T.c) - x - T.cl * T.c) / T.v
    T.ve > 0 && (T.ce += dt * x / T.ve)
    return T.c
end

"""
    simulate(; inj = Injection(), bp = BodyParams(), tmax = 90, dt = 0.005, every = 0.1,
             coronary = nothing) -> NamedTuple

Integrates the circulation. Returns sampled times and, for every compartment, the plasma iodine
concentration (mgI/mL): `arm`, `svc`, `ivc`, `ra`, `rv`, `pa`, `lung` (each stage), `la`, `lv`,
`aorta` (outflow of the ascending aorta = coronary input), `art` (systemic arteries feeding the
organs), and organ blood and extracellular curves. `coronary`, if given, is a function
`t -> plasma concentration leaving the coronary sinus` from the detailed coronary model; by
default the coronary bed is a single exchanging compartment.
"""
function simulate(;
        inj::Injection = Injection(), bp::BodyParams = BodyParams(), tmax::Real = 90, dt::Real = 0.005,
        every::Real = 0.1, coronary = nothing
    )
    dt = F32(dt); pf = 1.0f0 - bp.hct; Q = bp.co * pf       # plasma cardiac output, mL/s
    P(v) = v * pf
    arm = [Tank(P(bp.v_arm / bp.n_arm)) for _ in 1:bp.n_arm]
    svc = Tank(P(bp.v_svc)); ivc = Tank(P(bp.v_ivc)); ra = Tank(P(bp.v_ra)); rv = Tank(P(bp.v_rv)); pa = Tank(P(bp.v_pa))
    nl = bp.n_lung; vlung = bp.v_pulm - bp.v_pa
    lung = [Tank(P((1 - bp.lung_plug) * vlung / nl)) for _ in 1:nl]
    dlung = Delay(bp.lung_plug * vlung / bp.co, dt); drv = Delay(bp.beat / 2, dt); dlv = Delay(bp.beat / 2, dt)
    la = Tank(P(bp.v_la)); lv = Tank(P(bp.v_lv)); ao = Tank(P(bp.v_aorta)); art = Tank(P(bp.v_art - bp.v_aorta))
    # organ shares fill what the coronary flow leaves (at stress the heart takes a larger share)
    f0 = bp.f_brain + bp.f_kidney + bp.f_hepart + bp.f_portal + bp.f_muscle + bp.f_skin + bp.f_fat + bp.f_bone
    f_other0 = 1.0f0 - 0.033f0 - f0
    ks = (1.0f0 - bp.f_cor) / (1.0f0 - 0.033f0)
    f_other = ks * f_other0
    fq(f) = (f === bp.f_cor ? f : ks * f) * Q
    brain = Tank(P(bp.v_brain)); kid = Tank(P(bp.v_kidney); ve = bp.e_kidney, ps = 2fq(bp.f_kidney), cl = bp.gfr)
    portal = Tank(P(bp.v_portal); ve = bp.e_portal, ps = 2fq(bp.f_portal)); liver = Tank(P(bp.v_liver); ve = bp.e_liver, ps = 2fq(bp.f_hepart + bp.f_portal))
    mus = Tank(P(bp.v_muscle); ve = bp.e_muscle, ps = bp.ps_muscle); skin = Tank(P(bp.v_skin); ve = bp.e_skin, ps = bp.ps_skin)
    fat = Tank(P(bp.v_fat); ve = bp.e_fat, ps = bp.ps_fat); bone = Tank(P(bp.v_bone); ve = bp.e_bone, ps = bp.ps_bone)
    oth = Tank(P(bp.v_other); ve = bp.e_other, ps = bp.ps_other); cor = Tank(P(bp.v_cor); ve = 9.0f0, ps = 0.4f0 * fq(bp.f_cor))
    vup = Tank(P(bp.v_vein_upper)); vlo = Tank(P(bp.v_vein_lower))
    per = bp.f_muscle + bp.f_skin + bp.f_fat + bp.f_bone + f_other0
    q_arm = bp.arm_flow * Q
    nsteps = round(Int, tmax / dt); stride = max(1, round(Int, every / dt))
    names = (
        :arm, :svc, :ivc, :ra, :rv, :pa, :la, :lv, :aorta, :art, :brain, :kidney, :portal, :liver, :muscle, :skin, :fat, :bone, :other,
        :cor, :vein_upper, :vein_lower, :muscle_e, :liver_e, :portal_e, :kidney_e, :fat_e, :skin_e, :bone_e,
    )
    ns = nsteps ÷ stride + 1
    out = Dict{Symbol, Vector{F32}}(k => zeros(F32, ns) for k in names)
    lungout = zeros(F32, ns, nl); tt = zeros(F32, ns)
    excreted = 0.0f0; injected = 0.0f0
    for s in 0:nsteps
        t = F32(s * dt)
        if s % stride == 0
            k = s ÷ stride + 1; tt[k] = t
            for (nm, T) in (
                    (:arm, arm[end]), (:svc, svc), (:ivc, ivc), (:ra, ra), (:rv, rv), (:pa, pa), (:la, la), (:lv, lv), (:aorta, ao), (:art, art),
                    (:brain, brain), (:kidney, kid), (:portal, portal), (:liver, liver), (:muscle, mus), (:skin, skin), (:fat, fat), (:bone, bone),
                    (:other, oth), (:cor, cor), (:vein_upper, vup), (:vein_lower, vlo),
                )
                out[nm][k] = T.c
            end
            out[:muscle_e][k] = mus.ce; out[:liver_e][k] = liver.ce; out[:portal_e][k] = portal.ce; out[:kidney_e][k] = kid.ce
            out[:fat_e][k] = fat.ce; out[:skin_e][k] = skin.ce; out[:bone_e][k] = bone.ce
            for l in 1:nl
                lungout[k, l] = lung[l].c
            end
        end
        # injected arm: native arm plasma flow plus the injector
        ri = F32(inj_rate(inj, t)); ci = F32(inj_conc(inj, t)); injected += dt * ri * ci
        fa = q_arm + ri
        cin = (q_arm * vup.c + ri * ci) / fa
        for T in arm
            cin = (old = T.c; step!(T, fa, cin, dt); old)
        end
        carm = arm[end].c
        # venous returns (organs at their current outflow concentrations)
        q_up = fq(bp.f_brain) + bp.upper * fq(per) - q_arm
        c_up_in = (fq(bp.f_brain) * brain.c + bp.upper * (fq(bp.f_muscle) * mus.c + fq(bp.f_skin) * skin.c + fq(bp.f_fat) * fat.c + fq(bp.f_bone) * bone.c + fq(f_other0) * oth.c)) /
            (fq(bp.f_brain) + bp.upper * fq(per))
        q_lo = (1 - bp.upper) * fq(per) + fq(bp.f_kidney) + fq(bp.f_hepart + bp.f_portal)
        c_lo_in = (
            (1 - bp.upper) * (fq(bp.f_muscle) * mus.c + fq(bp.f_skin) * skin.c + fq(bp.f_fat) * fat.c + fq(bp.f_bone) * bone.c + fq(f_other0) * oth.c) +
                fq(bp.f_kidney) * kid.c + fq(bp.f_hepart + bp.f_portal) * liver.c
        ) / q_lo
        step!(vup, q_up + q_arm, c_up_in, dt); step!(vlo, q_lo, c_lo_in, dt)
        fsvc = q_up + fa
        step!(svc, fsvc, (q_up * vup.c + fa * carm) / fsvc, dt)
        step!(ivc, q_lo, vlo.c, dt)
        ccs = coronary === nothing ? cor.c : F32(coronary(t))
        fra = fsvc + q_lo + fq(bp.f_cor)
        step!(ra, fra, (fsvc * svc.c + q_lo * ivc.c + fq(bp.f_cor) * ccs) / fra, dt)
        step!(rv, fra, push!(drv, ra.c), dt); step!(pa, fra, rv.c, dt)
        c = push!(dlung, pa.c)
        for T in lung
            c = (old = T.c; step!(T, fra, c, dt); old)
        end
        step!(la, fra, lung[end].c, dt); step!(lv, fra, push!(dlv, la.c), dt); step!(ao, fra, lv.c, dt)
        fart = fra - fq(bp.f_cor)
        step!(art, fart, ao.c, dt)
        step!(cor, fq(bp.f_cor), ao.c, dt)
        ca = art.c
        step!(brain, fq(bp.f_brain), ca, dt); step!(kid, fq(bp.f_kidney), ca, dt)
        excreted += dt * bp.gfr * kid.c
        step!(portal, fq(bp.f_portal), ca, dt)
        step!(liver, fq(bp.f_hepart + bp.f_portal), (fq(bp.f_hepart) * ca + fq(bp.f_portal) * portal.c) / fq(bp.f_hepart + bp.f_portal), dt)
        step!(mus, fq(bp.f_muscle), ca, dt); step!(skin, fq(bp.f_skin), ca, dt); step!(fat, fq(bp.f_fat), ca, dt)
        step!(bone, fq(bp.f_bone), ca, dt); step!(oth, fq(f_other0), ca, dt)
    end
    return (; t = tt, out..., lung = lungout, injected, excreted, hct = bp.hct, params = bp, inj)
end

end # module
