# SPDX-License-Identifier: GPL-3.0-or-later
# Contrast bolus through the circulation as a graph of named vessel segments.
#
# - A `Circulation` is a list of `Segment`s joined at named nodes. The injection site and every
#   measurement site (an AIF in the aortic root for coronary CTA, in the pulmonary trunk for
#   pulmonary CTA, in an internal carotid for the brain) are just names in that graph.
# - Each segment is a lumped transfer element: a volume carried by plug flow (exact, parcels
#   in a queue, so a changing flow such as the injector stopping is handled exactly), an optional
#   pure delay in seconds (a ventricle ejects at the next beat), then well-mixed stages whose
#   number follows from the segment's relative dispersion. Beds can exchange iodine with an
#   extracellular space (permeability-surface product PS) and clear it (kidney).
# - Any segment can be replaced by a detailed model (the grown coronary tree now; grown lung or
#   brain trees later) through `detailed`, a map from segment name to t -> outflow concentration.
# - The graph need not be a whole body. A study of one organ can feed that organ from a single
#   path segment with a delay and a dispersion (`organ_only`), with no heart or lungs at all.
# - Concentrations are mgI per mL of whole blood (what CT sees). Iodine stays in plasma, so
#   exchange uses the plasma concentration c / (1 − Hct). Float32 throughout.
#
# Volumes and flow: flows are fractions of cardiac output given for the beds (segments from an
# arterial node to a venous node); artery and vein flows follow from conservation. The injected
# volume raises venous return and, through the Frank-Starling mechanism, cardiac output while
# the injector runs; the systemic beds store it (venous compliance), so total iodine and volume
# are conserved exactly.

module Body

export Segment, Circulation, Injection, human_circulation, adenosine, organ_only, simulate, curve, transit_table, aif

const F32 = Float32

"""
    Segment(name, from, to; kind, q, v, plug = 0, delay = 0, rd = 1, ve = 0, ps = 0, cl = 0)

- `kind`: `:bed` (arterial node → venous node, has its own flow fraction `q` of cardiac
  output), `:artery`, `:vein` or `:central` (heart and lungs; carries the whole output).
  Artery and vein flows are derived from the beds unless `q` is given (a vessel that feeds
  only part of a bed, such as the hepatic artery).
- `v`: blood volume, mL. `plug`: fraction of `v` moved by plug flow (first). `delay`: pure
  delay, s. The rest of `v` is `n` well-mixed stages with `n` set so that the segment's transit
  time has relative dispersion `rd` (SD / mean) at its nominal flow; `rd = 1` with no plug is
  one well-mixed tank.
- `ve`, `ps`: extracellular volume (mL) and permeability-surface product (mL/s) of a bed;
  `cl`: plasma clearance, mL/s (kidney: GFR).
- `c`, `r`: venous compliance (mL/mmHg) and outflow resistance (mmHg·s/mL) of a vein that
  receives an injection. Extra inflow distends it: the stored volume ΔV raises its pressure by
  ΔV / c and its outflow by ΔV / (r c), so the injected volume leaves with time constant r·c
  and keeps leaving after the injector stops. 0 = rigid.
"""
Base.@kwdef struct Segment
    name::Symbol
    from::Symbol
    to::Symbol
    kind::Symbol = :bed
    q::F32 = 0.0f0
    v::F32
    plug::F32 = 0.0f0
    delay::F32 = 0.0f0
    rd::F32 = 1.0f0
    ve::F32 = 0.0f0
    ps::F32 = 0.0f0
    cl::F32 = 0.0f0
    c::F32 = 0.0f0
    r::F32 = 0.0f0
end
Segment(name, from, to; kw...) = Segment(; name, from, to, kw...)

"""
    Circulation(segments; co, hct, inlet = :ra_in)

`co` cardiac output (mL/s of blood), `hct` large-vessel hematocrit. `inlet` is the node where
the venous return meets the heart; the injected volume travels with the venous flow up to it.
"""
struct Circulation
    segs::Vector{Segment}
    co::F32
    hct::F32
    inlet::Symbol
    q::Vector{F32}          # nominal flow of each segment, fraction of cardiac output
end

function Circulation(segs::Vector{Segment}; co, hct, inlet = :ra_in)
    q = _flows(segs)
    return Circulation(segs, F32(co), F32(hct), inlet, q)
end

# Flow fractions: beds as given, central = 1, arteries = sum of what leaves their distal node,
# veins = sum of what enters their proximal node.
function _flows(segs)
    n = length(segs); q = fill(F32(NaN), n)
    for (i, s) in enumerate(segs)
        (s.kind === :bed || s.q > 0) && (q[i] = s.q); s.kind === :central && (q[i] = 1.0f0)
    end
    out_of(node) = [i for (i, s) in enumerate(segs) if s.from === node]
    into(node) = [i for (i, s) in enumerate(segs) if s.to === node]
    for _ in 1:n, (i, s) in enumerate(segs)
        isnan(q[i]) || continue
        ks = s.kind === :artery ? out_of(s.to) : s.kind === :vein ? into(s.from) : Int[]
        (!isempty(ks) && all(!isnan, q[ks])) && (q[i] = sum(q[ks]))
    end
    bad = [segs[i].name for i in 1:n if isnan(q[i])]
    isempty(bad) || error("cannot derive flows for $(bad): every artery must end in beds, every vein start from them")
    return q
end

"""
    Injection(; site, phases)
    Injection(; site, conc, volume, rate, saline, saline_rate)

`site` is a node of the circulation (`:antecubital_r`, `:antecubital_l`, `:jugular_r`,
`:jugular_l`, `:femoral_r`, `:femoral_l`, `:svc_top` for a central line, ...). `phases` is a
list of (volume mL, rate mL/s, concentration mgI/mL) given in order from t = 0, so a test
bolus, a biphasic protocol or a saline chaser are all phases. `dead_space` is the volume of the
connecting tube and cannula between syringe and vein (Medrad low-pressure line 2.78 mL plus
about 1.2 mL syringe dead space, [cowan2017]; long extension sets about 12 mL), plug flow,
primed with contrast (`prime = :contrast`, the usual air purge) or saline. Default: iohexol 350, 50 mL at
5 mL/s, then 40 mL saline (Lee et al. 2022 adenosine CT perfusion protocol), right antecubital
vein (avoids the left brachiocephalic vein, Tseng et al. 2007).
"""
struct Injection
    site::Symbol
    phases::Vector{NTuple{3, F32}}
    dead_space::F32
    prime::F32
end
function Injection(;
        site = :antecubital_r, phases = nothing, conc = 350, volume = 50, rate = 5, saline = 40, saline_rate = rate,
        dead_space = 0, prime = :contrast
    )
    ph = phases !== nothing ? [F32.(p) for p in phases] : NTuple{3, F32}[(volume, rate, conc)]
    phases === nothing && saline > 0 && push!(ph, (saline, saline_rate, 0))
    return Injection(site, ph, F32(dead_space), prime === :contrast ? ph[1][3] : 0.0f0)
end
function _inj(j::Injection, t)
    t0 = 0.0f0
    for (v, r, c) in j.phases
        d = v / r
        t < t0 + d && return (r, c)
        t0 += d
    end
    return (0.0f0, 0.0f0)
end

# ---- runtime elements
mutable struct Pipe                      # exact plug flow: a queue of fluid parcels
    vol::Vector{F32}
    c::Vector{F32}
    head::Int
    V::F32
end
function Pipe(V)
    return Pipe(F32[V], F32[0], 1, F32(V))
end
function pass!(p::Pipe, dv::F32, cin::F32)
    p.V <= 0 && return cin
    push!(p.vol, dv); push!(p.c, cin)
    need = dv; m = 0.0f0
    while need > 0 && p.head <= length(p.vol)
        take = min(need, p.vol[p.head]); m += take * p.c[p.head]
        p.vol[p.head] -= take; need -= take
        p.vol[p.head] <= 1.0f-9 && (p.head += 1)
    end
    if p.head > 4096                     # compact the queue now and then
        deleteat!(p.vol, 1:(p.head - 1)); deleteat!(p.c, 1:(p.head - 1)); p.head = 1
    end
    return m / dv
end
pipe_mass(p::Pipe) = sum(p.vol[k] * p.c[k] for k in p.head:length(p.vol); init = 0.0f0)

mutable struct Delay                     # pure delay in seconds (ring buffer at the step)
    buf::Vector{F32}
    k::Int
end
Delay(d, dt) = Delay(zeros(F32, max(1, round(Int, d / dt))), 1)
function pass!(D::Delay, c::F32)
    out = D.buf[D.k]; D.buf[D.k] = c; D.k = D.k == length(D.buf) ? 1 : D.k + 1
    return out
end

mutable struct Run
    pipe::Pipe
    delay::Union{Nothing, Delay}
    c::Vector{F32}                       # stage concentrations
    vs::Vector{F32}                      # stage volumes (the last one stores infused volume)
    ce::F32
    out::F32                             # outflow concentration during the last step
end

# mixing volume: what plug flow and the pure delay (blood waiting in a chamber) leave of `v`
_vmix(s::Segment, qn, co) = max(0.0f0, (1 - s.plug) * s.v - qn * co * s.delay)
function _stages(s::Segment, qn::F32, co::F32)
    vm = _vmix(s, qn, co)
    vm <= 0 && return 0
    T = s.v / (qn * co); τm = vm / (qn * co)
    return max(1, round(Int, (τm / (s.rd * T))^2))
end

"""
    transit_table(C) -> Vector of NamedTuples

Nominal flow (mL/s), mean transit time (s), mixing stages and relative dispersion of every
segment, for checking a circulation against measured transit times.
"""
function transit_table(C::Circulation)
    return [
        (;
            name = s.name, flow = C.q[i] * C.co, mtt = s.v / (C.q[i] * C.co), n = _stages(s, C.q[i], C.co),
            rd = (n = _stages(s, C.q[i], C.co); n == 0 ? 0.0f0 : _vmix(s, C.q[i], C.co) / sqrt(F32(n)) / s.v),
        )
            for (i, s) in enumerate(C.segs)
    ]
end

"""
    simulate(C, inj; tmax = 90, dt = 0.002, every = 0.1, detailed = Dict()) -> NamedTuple

Integrates the bolus. Returns `t` and, per segment name, `out` (outflow blood concentration,
mgI/mL), `mean` (mean blood concentration inside the segment) and `ecv` (extracellular
concentration of a bed), plus `injected`, `excreted` and `held` (in the body at `tmax`) iodine, mg. `detailed[name]` replaces
segment `name` by a model: a function of t giving its outflow blood concentration.
"""
function simulate(
        C::Circulation, inj::Injection = Injection(); tmax::Real = 90, dt::Real = 0.002, every::Real = 0.1,
        detailed::AbstractDict = Dict{Symbol, Any}()
    )
    dt = F32(dt); segs = C.segs; ns = length(segs); co = C.co; pf = 1.0f0 - C.hct
    nodes = unique(vcat([s.from for s in segs], [s.to for s in segs]))
    inj.site in nodes || error("injection site $(inj.site) is not a node; nodes: $(nodes)")
    nid = Dict(n => k for (k, n) in enumerate(nodes))
    from = [nid[s.from] for s in segs]; to = [nid[s.to] for s in segs]
    # venous path carrying the injected volume: from the site, follow veins to the heart inlet
    path = Int[]; node = inj.site
    while node !== C.inlet
        k = findfirst(i -> segs[i].from === node && segs[i].kind === :vein, 1:ns)
        k === nothing && break
        push!(path, k); node = segs[k].to
    end
    onpath = falses(ns); onpath[path] .= true
    dvol = zeros(F32, ns)                                 # volume stored in compliant veins, mL
    line = Pipe(inj.dead_space); line.c[1] = inj.prime
    R = Run[]
    for (i, s) in enumerate(segs)
        n = _stages(s, C.q[i], co); vm = _vmix(s, C.q[i], co)
        push!(R, Run(Pipe(s.plug * s.v), s.delay > 0 ? Delay(s.delay, dt) : nothing, zeros(F32, n), fill(F32(n == 0 ? 0 : vm / n), n), 0.0f0, 0.0f0))
    end
    names = [s.name for s in segs]
    nt = round(Int, tmax / dt); stride = max(1, round(Int, every / dt)); nk = nt ÷ stride + 1
    rec_out = zeros(F32, nk, ns); rec_mean = zeros(F32, nk, ns); rec_ecv = zeros(F32, nk, ns); tt = zeros(F32, nk)
    det = [get(detailed, s.name, nothing) for s in segs]
    nodemass = zeros(F32, length(nodes)); nodeflow = zeros(F32, length(nodes)); cnode = zeros(F32, length(nodes))
    qin = zeros(F32, ns); qout = zeros(F32, ns)
    # a bed takes its share of whatever reaches its inlet node (the liver: hepatic artery + portal vein)
    beds = [i for i in 1:ns if segs[i].kind === :bed]; others = [i for i in 1:ns if segs[i].kind !== :bed]
    inc = [[k for k in 1:ns if to[k] == from[i]] for i in 1:ns]
    share = [C.q[i] / sum(C.q[j] for j in 1:ns if from[j] == from[i]) for i in 1:ns]
    injected = 0.0f0; excreted = 0.0f0
    for step in 0:nt
        t = F32(step * dt)
        ri, ci = _inj(inj, t); injected += dt * ri * ci
        ri > 0 && (ci = pass!(line, ri * dt, ci))           # through the connecting tube
        # the injected flow moves down the venous path; a compliant vein stores part of it
        extra = ri
        for i in path
            s = segs[i]; qin[i] = C.q[i] * co + extra
            if s.c > 0
                out = dvol[i] / (s.r * s.c); dvol[i] += dt * (extra - out); extra = out
            end
            qout[i] = C.q[i] * co + extra
        end
        coe = co + extra                                  # Frank-Starling: output follows return
        for i in others
            onpath[i] && continue
            qin[i] = segs[i].kind === :vein ? C.q[i] * co : C.q[i] * coe
            qout[i] = qin[i]
        end
        for i in beds                                     # beds store the infused volume
            qin[i] = share[i] * sum(qout[k] for k in inc[i]; init = 0.0f0); qout[i] = C.q[i] * co
        end
        if step % stride == 0
            kk = step ÷ stride + 1; tt[kk] = t
            for i in 1:ns
                r = R[i]; m = pipe_mass(r.pipe) + sum(r.c .* r.vs; init = 0.0f0); v = r.pipe.V + sum(r.vs; init = 0.0f0)
                rec_out[kk, i] = r.out; rec_mean[kk, i] = v > 0 ? m / v : r.out; rec_ecv[kk, i] = r.ce
            end
        end
        # node concentrations from last step's outflows
        fill!(nodemass, 0.0f0); fill!(nodeflow, 0.0f0)
        for i in 1:ns
            nodemass[to[i]] += qout[i] * R[i].out; nodeflow[to[i]] += qout[i]
        end
        s0 = nid[inj.site]; nodemass[s0] += ri * ci; nodeflow[s0] += ri
        for k in eachindex(cnode)
            cnode[k] = nodeflow[k] > 0 ? nodemass[k] / nodeflow[k] : 0.0f0
        end
        for i in 1:ns
            s = segs[i]; r = R[i]; cin = cnode[from[i]]; f = qin[i]
            if det[i] !== nothing
                r.out = F32(det[i](t)); continue
            end
            c = pass!(r.pipe, f * dt, cin)
            r.delay === nothing || (c = pass!(r.delay, c))
            n = length(r.c)
            if n == 0
                r.out = c; continue
            end
            psn = s.ps / n; cln = s.cl / n; xs = 0.0f0
            for k in 1:n
                old = r.c[k]; x = psn * (old / pf - r.ce)
                fo = k == n ? qout[i] : f
                # stage with inflow f, outflow fo: d(cV)/dt = f cin − fo c − exchange − clearance
                r.c[k] += dt * (f * (c - old) - x - cln * old / pf) / r.vs[k]
                k == n && (r.vs[k] += dt * (f - fo))
                excreted += dt * cln * old / pf; xs += x
                c = old
            end
            s.ve > 0 && (r.ce += dt * xs / s.ve)
            r.out = c
        end
    end
    held = sum(
        pipe_mass(R[i].pipe) + sum(R[i].c .* R[i].vs; init = 0.0f0) + R[i].ce * segs[i].ve +
            (R[i].delay === nothing ? 0.0f0 : sum(R[i].delay.buf) * qin[i] * dt) for i in 1:ns
    )
    res = Dict{Symbol, Any}(:t => tt, :injected => injected, :excreted => excreted, :held => held, :names => names, :circ => C, :inj => inj)
    res[:out] = Dict(names[i] => rec_out[:, i] for i in 1:ns)
    res[:mean] = Dict(names[i] => rec_mean[:, i] for i in 1:ns)
    res[:ecv] = Dict(names[i] => rec_ecv[:, i] for i in 1:ns)
    return (; (Symbol(k) => v for (k, v) in res)...)
end

"""
    curve(r, name; what = :out) -> Vector{Float32}

Blood concentration (mgI/mL) leaving segment `name` (`what = :out`), inside it (`:mean`), or
its extracellular concentration (`:ecv`).
"""
curve(r, name::Symbol; what::Symbol = :out) = getfield(r, what)[name]

"""
    aif(r, site) -> (t, c)

Arterial input function at a named site: the outflow of that segment.
"""
aif(r, site::Symbol) = (r.t, curve(r, site))

# ---------------------------------------------------------------------------------------------
"""
    adenosine(; co = 81.6, beat = 1.0, hr = 96/64, sv = 51.0/44.8) -> (; co, beat)

Cardiac output (mL/s) and beat length (s) during i.v. adenosine 140 µg/kg/min, from the resting
values. In 25 healthy volunteers heart rate rose from 64 ± 8 to 96 ± 13 /min and stroke volume
index from 44.8 ± 10.1 to 51.0 ± 11.8 mL/m², so cardiac output rises × 1.71 (Thomas et al. 2017,
Table 2, [thomas2017]). Blood pressure barely moves (systolic −5, diastolic 0 mmHg), so the coronary
model keeps its aortic pressure. Other series bracket the rise: +28% (abstract of the same group)
and × 2.2–2.4 by gated Rb-82 PET, whose resting volumes read low [lassen2023].
Use as `human_circulation(; adenosine()..., f_cor = q_cor / adenosine().co)`.
"""
adenosine(; co = 81.6f0, beat = 1.0f0, hr = 96.0f0 / 64.0f0, sv = 51.0f0 / 44.8f0) = (co = F32(co * hr * sv), beat = F32(beat / hr))

"""
    human_circulation(; kwargs...) -> Circulation

Adult male (ICRP 89 reference, 5.3 L blood) with a named venous tree for every common injection
site, the four heart chambers, the lungs, the coronary bed, and the systemic organs. All numbers
are keyword arguments; sources are in `docs` of each block below.

Coronary bed: `cor = (; ml, fb, ve, ps, rho)` (required): myocardial tissue volume and the heart's
blood fraction, interstitial fraction, contrast PS and density (`Flow.Physiology(organ = "heart")`).

Heart (from the XCAT phantom): `co` 81.6 mL/s (stroke volume 132.2 − 50.6 mL at 60/min),
chamber volumes `v_ra` 52.0, `v_rv` 85.4, `v_la` 37.6, `v_lv` 95.2 mL, `beat` 1 s. Blood that
enters a ventricle in diastole leaves at the next systole, so each ventricle adds half a beat of
pure delay before its mixing volume.

Lungs: `v_pa` 40 mL pulmonary trunk and main branches (plug), `v_lung` 510 mL (ICRP 89: pulmonary
blood 10.5% of 5.3 L ≈ 550 mL in all), `lung_plug` and `lung_rd` give its transit-time shape.

Veins (lumen volume = area × length, supine; data/reference ids in brackets):
- arm (antecubital → axilla): basilic 4.4 mm [salinaro2024] over 25 cm (no in vivo length
  found) = 3.8 mL; axillary 8 mm (between basilic and subclavian; no in vivo series) × 12 cm =
  6.0 mL; subclavian 9.9 mm [berk2026] × 4 cm = 3.1 mL.
- brachiocephalic: right 10.9 mm × 38 mm = 3.5 mL, left 12.1 mm × 63 mm = 7.3 mL (cadaver
  circumference and length [wang2024]; the left crosses the mediastinum).
- internal jugular: right 96 mm², left 65 mm² [kawalecrutkowska2025] × 12 cm = 11.5 / 7.8 mL.
- SVC: 70.7 mm long [mahlon2007] × 2.18 cm² (geometric mean of the 5th–95th percentile range
  1.08–4.42 cm² [lin2009]; cadaver diameter 1.67 cm gives 2.19 cm² [wang2024]) = 15.4 mL.
- femoral + external iliac, common iliac, IVC (infrarenal, renal → hepatic, hepatic → RA):
  standard anatomical sizes, no series read yet.
Large veins are near plug flow (95% plug, relative dispersion 0.04): with the arm compliance
below, this is what keeps a hand-pushed radionuclide bolus as narrow in the SVC as measured
(0.60 ± 0.21 s at half height, Dymond et al. 1982, [dymond1982]; the model gives 0.76 s).
The injected arm's upper-arm veins are compliant (the valves keep the forearm out of it):
- `arm_c` 0.88 mL/mmHg: forearm venous compliance at ~10 mmHg, β1 + 2β2·P = 0.108 − 2·0.00102·10
  ≈ 0.088, and 0.071 in the second group, so 0.08 mL per 100 mL per mmHg [halliwill1999; Oue et al.
  2019 give 0.08–0.11 at 20 mmHg, oue2019], times the upper-arm tissue volume 1.1 L (arm from wrist
  1.92 L [buffa2015] minus forearm 0.82 L [harvey2014]).
- `arm_r` 1.76 mmHg·s/mL: arm-vein minus central-vein pressure 2.2 mmHg (spontaneous breathing,
  [amar2001]; 2.0 mmHg in [coskun2025], 3 mmHg anaesthetised in [munis2001]) over resting brachial
  flow 74.9 mL/min (304 adults supine, [zhang2025]).
So a 5 mL/s injection raises arm venous pressure by about 9 mmHg and the arm holds part of the
dose with τ = r·c ≈ 1.5 s. Set `arm_c = 0` for a rigid path.

Beds (Williams & Leggett 1989 / ICRP 89 flow fractions and blood volumes): brain 0.12 (half per
side, draining to each internal jugular), kidneys 0.19 (GFR 2 mL/s plasma, Edelson et al. 1984),
hepatic artery 0.065 and gut + spleen 0.185 (portal → liver), coronary `f_cor` 0.033, and
peripheral tissue (muscle 0.17, skin 0.05, fat 0.05, bone 0.05, other rest) split by region:
each arm 5%, each leg 15%, head-neck-thorax 30% (to the brachiocephalic veins and SVC), pelvis and
lower trunk 30% (to the IVC). Extracellular volumes and PS: muscle 4000 mL / 7 mL/s, skin 1300 /
1, fat 2000 / 2, bone 1500 / 1, other 1500 / 3 (scaled from EDTA, Levitt 2003); gut, liver and
kidney exchange is flow-limited; the brain has none. Systemic arteries beyond the aorta (160 mL)
are each bed's plug part.
"""
function human_circulation(;
        co = 81.6f0, hct = 0.43f0, beat = 1.0f0,
        v_ra = 52.0f0, v_rv = 85.4f0, v_la = 37.6f0, v_lv = 95.2f0,
        v_pa = 40.0f0, pa_plug = 0.5f0, pa_rd = 0.2f0, v_lung = 510.0f0, lung_plug = 0.25f0, lung_rd = 0.46f0,
        v_aortic_root = 25.0f0, v_arch = 35.0f0, v_desc = 100.0f0, v_sys_art = 160.0f0,
        vein_plug = 0.95f0, vein_rd = 0.04f0, arm_c = 0.88f0, arm_r = 1.76f0,
        v_arm = 3.8f0, v_axillary = 6.0f0, v_subclavian = 3.1f0, v_bc_r = 3.5f0, v_bc_l = 7.3f0, v_ij_r = 11.5f0, v_ij_l = 7.8f0, v_svc = 15.4f0,
        v_leg = 17.0f0, v_iliac = 7.7f0, v_ivc_infra = 35.0f0, v_ivc_renal = 27.0f0, v_ivc_hep = 21.0f0,
        f_cor = 0.033f0, cor,
        f_brain = 0.12f0, v_brain = 64.0f0, f_kidney = 0.19f0, v_kidney = 106.0f0, e_kidney = 100.0f0, gfr = 2.0f0,
        f_hepart = 0.065f0, f_portal = 0.185f0, v_portal = 250.0f0, e_portal = 600.0f0, v_liver = 530.0f0, e_liver = 500.0f0,
        tissues = (
            (:muscle, 0.17f0, 742.0f0, 4000.0f0, 7.0f0), (:skin, 0.05f0, 160.0f0, 1300.0f0, 1.0f0), (:fat, 0.05f0, 265.0f0, 2000.0f0, 2.0f0),
            (:bone, 0.05f0, 370.0f0, 1500.0f0, 1.0f0), (:other, NaN32, 300.0f0, 1500.0f0, 3.0f0),
        ),
        regions = ((:arm_r, :antecubital_r, 0.05f0), (:arm_l, :antecubital_l, 0.05f0), (:leg_r, :femoral_r, 0.15f0), (:leg_l, :femoral_l, 0.15f0), (:upper, :svc_top, 0.3f0), (:lower, :ivc_low, 0.3f0)),
        v_veins_extra = 750.0f0,
    )
    F = F32
    S = Segment[]
    vein(name, a, b, v) = push!(S, Segment(name, a, b; kind = :vein, v = F(v), plug = F(vein_plug), rd = F(vein_rd)))
    # heart and lungs
    push!(S, Segment(:ra, :ra_in, :tv; kind = :central, v = F(v_ra)))
    push!(S, Segment(:rv, :tv, :pv; kind = :central, v = F(v_rv), delay = F(beat / 2)))
    push!(S, Segment(:pulm_trunk, :pv, :pa_br; kind = :central, v = F(v_pa), plug = F(pa_plug), rd = F(pa_rd)))
    push!(S, Segment(:lung, :pa_br, :pulm_vein; kind = :central, v = F(v_lung), plug = F(lung_plug), rd = F(lung_rd)))
    push!(S, Segment(:la, :pulm_vein, :mv; kind = :central, v = F(v_la)))
    push!(S, Segment(:lv, :mv, :av; kind = :central, v = F(v_lv), delay = F(beat / 2)))
    # aorta and beds; the flow fractions of the beds other than the coronary sum to 1 − f_cor
    fo = f_brain + f_kidney + f_hepart + f_portal + sum(x[2] for x in tissues if !isnan(x[2]))
    f_other = 1.0f0 - 0.033f0 - fo
    ks = (1.0f0 - f_cor) / (1.0f0 - 0.033f0)          # at stress the heart takes a larger share
    push!(S, Segment(:aortic_root, :av, :root; kind = :artery, v = F(v_aortic_root), plug = 0.8f0, rd = 0.1f0))
    # lumped coronary bed from the heart's own physiology: `cor = (; ml, fb, ve, ps, rho)` with the
    # myocardial tissue volume (mL) and the blood fraction, interstitial fraction, contrast PS
    # (mL/min/g) and density of `Flow.Physiology(organ = "heart")`
    push!(S, Segment(:coronary, :root, :cs; kind = :bed, q = F(f_cor), v = F(cor.fb * cor.ml), ve = F(cor.ve * cor.ml), ps = F(cor.ps / 60 * cor.ml * cor.rho)))
    push!(S, Segment(:coronary_sinus, :cs, :ra_in; kind = :vein, v = 3.0f0, plug = 0.8f0, rd = 0.1f0))
    push!(S, Segment(:arch, :root, :arch_out; kind = :artery, v = F(v_arch), plug = 0.8f0, rd = 0.1f0))
    push!(S, Segment(:desc_aorta, :arch_out, :abd; kind = :artery, v = F(v_desc), plug = 0.8f0, rd = 0.1f0))
    ftot = ks * (fo + f_other)                         # all non-coronary beds
    va(f) = F(v_sys_art * f / ftot)                    # each bed's share of the systemic arteries
    fres = ks * (f_portal + sum(x[2] for x in tissues if !isnan(x[2])) + f_other)
    vx(f) = F(v_veins_extra * f / fres)                # unnamed systemic veins: limbs, trunk, splanchnic
    function bed(name, a, b, f, vb; ve = 0, ps = 0, cl = 0, reservoir = false)
        v = va(f) + F(vb) + (reservoir ? vx(f) : 0.0f0)
        return push!(S, Segment(name, a, b; kind = :bed, q = F(f), v, plug = va(f) / v, ve = F(ve), ps = F(ps), cl = F(cl)))
    end
    bed(:brain_r, :arch_out, :jugular_r, ks * f_brain / 2, v_brain / 2)
    bed(:brain_l, :arch_out, :jugular_l, ks * f_brain / 2, v_brain / 2)
    bed(:kidney, :abd, :ivc_renal, ks * f_kidney, v_kidney; ve = e_kidney, ps = 2 * ks * f_kidney * co * (1 - hct), cl = gfr)
    bed(:gut, :abd, :portal_in, ks * f_portal, v_portal; ve = e_portal, ps = 2 * ks * f_portal * co * (1 - hct), reservoir = true)
    push!(S, Segment(:portal_vein, :portal_in, :liver_in; kind = :vein, v = 15.0f0, plug = 0.5f0, rd = 0.3f0))
    push!(S, Segment(:hepatic_artery, :abd, :liver_in; kind = :artery, q = F(ks * f_hepart), v = 3.0f0, plug = 0.8f0, rd = 0.1f0))
    v = F(v_liver)
    push!(S, Segment(:liver, :liver_in, :ivc_hepatic; kind = :bed, q = F(ks * (f_hepart + f_portal)), v, ve = F(e_liver), ps = F(2 * ks * (f_hepart + f_portal) * co * (1 - hct))))
    for (region, node, share) in regions, (tname, ft, vt, et, pst) in tissues
        f = ks * (isnan(ft) ? f_other : ft) * share
        bed(Symbol(tname, :_, region), region in (:arm_r, :arm_l, :upper) ? :arch_out : :abd, node, f, vt * share; ve = et * share, ps = pst * share, reservoir = true)
    end
    # named veins
    armv(name, a, b) = push!(S, Segment(name, a, b; kind = :vein, v = F(v_arm), plug = F(vein_plug), rd = F(vein_rd), c = F(arm_c), r = F(arm_r)))
    armv(:arm_vein_r, :antecubital_r, :axillary_r_in); vein(:axillary_r, :axillary_r_in, :subclavian_r_in, v_axillary)
    vein(:subclavian_r, :subclavian_r_in, :bc_r, v_subclavian); vein(:ij_r, :jugular_r, :bc_r, v_ij_r); vein(:brachiocephalic_r, :bc_r, :svc_top, v_bc_r)
    armv(:arm_vein_l, :antecubital_l, :axillary_l_in); vein(:axillary_l, :axillary_l_in, :subclavian_l_in, v_axillary)
    vein(:subclavian_l, :subclavian_l_in, :bc_l, v_subclavian); vein(:ij_l, :jugular_l, :bc_l, v_ij_l); vein(:brachiocephalic_l, :bc_l, :svc_top, v_bc_l)
    vein(:svc, :svc_top, :ra_in, v_svc)
    vein(:femoral_vein_r, :femoral_r, :iliac_r_in, v_leg); vein(:common_iliac_r, :iliac_r_in, :ivc_low, v_iliac)
    vein(:femoral_vein_l, :femoral_l, :iliac_l_in, v_leg); vein(:common_iliac_l, :iliac_l_in, :ivc_low, v_iliac)
    vein(:ivc_infrarenal, :ivc_low, :ivc_renal, v_ivc_infra); vein(:ivc_renal, :ivc_renal, :ivc_hepatic, v_ivc_renal); vein(:ivc_hepatic, :ivc_hepatic, :ra_in, v_ivc_hep)
    return Circulation(S; co, hct)
end

"""
    organ_only(organ::Segment; site = :inlet, delay, mtt, rd, co, hct) -> Circulation

A circulation with no heart or lungs: the injection site feeds one lumped path (pure `delay`
plus a transit of mean `mtt` and relative dispersion `rd`) that ends at the organ's inlet node,
and the organ drains to the outlet. Use it for a single-organ study when a whole-body model is not
wanted; `co` here is the organ's own flow.
"""
function organ_only(organ::Segment; site = :inlet, delay = 0.0f0, mtt = 5.0f0, rd = 0.4f0, co, hct = 0.43f0)
    path = Segment(:path, site, organ.from; kind = :central, v = F32((delay + mtt) * co), delay = F32(delay), rd = F32(rd))
    o = Segment(organ.name, organ.from, organ.to; kind = :central, v = organ.v, plug = organ.plug, delay = organ.delay, rd = organ.rd, ve = organ.ve, ps = organ.ps, cl = organ.cl)
    return Circulation(Segment[path, o], F32(co), F32(hct), site, F32[1, 1])
end

end # module
