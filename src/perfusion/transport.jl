# SPDX-License-Identifier: GPL-3.0-or-later
# The one transport entry point: a contrast injection through the whole-body circulation, with the organ's own bed.

using .Body: Injection, human_circulation, simulate, ORGAN_SLOTS

"""
    TransportSpec(; phys, tissue_ml, bed, input, co, beat, f_coronary, tmax, dt, every, coronary)
    TransportSpec(spec::TransportSpec; kwargs...)

Everything that sets where an injected bolus goes. Every field is a keyword and every one is required: an organ is an
instance (`HEART_TRANSPORT`, `LUNG_TRANSPORT`, organs.jl), and a run varies one with the second form. `phys` is the
organ's [`FlowSpec`](@ref): its blood hematocrit (`hd`, used for the whole body), blood and interstitial fractions,
contrast permeability and tissue density. `tissue_ml` is the organ's tissue volume (mL). `bed` is the body segment the
organ is (`:coronary` for the heart, `:lung` for the lungs) and `input` the segment whose outflow feeds it (`:lv` for an
organ fed by the aorta, `:rv` for the lungs). `co` is the cardiac output (mL/s; rest 81.6, adenosine stress
`Body.adenosine()`), `beat` the heart period (s), and `f_coronary` the coronary bed's share of the cardiac output (for
the heart: its root flow / `co`). When the organ is not the heart, `coronary` gives the body's own lumped coronary bed as
`(; ml, fb, ve, ps, rho)` (the heart's FlowSpec values and myocardial volume); otherwise `nothing`. `tmax`, `dt` and
`every` are the simulated time, the step and the output spacing (s).
"""
struct TransportSpec
    phys::FlowSpec
    tissue_ml::Float32
    bed::Symbol
    input::Symbol
    co::Float32
    beat::Float32
    f_coronary::Float32
    tmax::Float32
    dt::Float32
    every::Float32
    coronary::Union{Nothing, NamedTuple}
end
function TransportSpec(; kwargs...)
    bad = setdiff(keys(kwargs), fieldnames(TransportSpec))
    isempty(bad) || throw(ArgumentError("unknown TransportSpec fields $(collect(bad))"))
    miss = setdiff(fieldnames(TransportSpec), keys(kwargs))
    isempty(miss) || throw(ArgumentError("TransportSpec: no value for $(collect(miss))"))
    return TransportSpec((kwargs[f] for f in fieldnames(TransportSpec))...)
end
TransportSpec(s::TransportSpec; kwargs...) = TransportSpec(; merge((; (f => getfield(s, f) for f in fieldnames(TransportSpec))...), values(kwargs))...)

"""
    circulation(spec::TransportSpec) -> Body.Circulation

The whole-body circulation (`Body.human_circulation`) with the organ's own values in its body segment `spec.bed` and its
blood everywhere.
"""
function circulation(s::TransportSpec)
    s.bed in ORGAN_SLOTS || throw(ArgumentError("TransportSpec: bed must be a body organ slot $(ORGAN_SLOTS), got :$(s.bed)"))
    (s.bed === :coronary && s.coronary !== nothing) && throw(ArgumentError("TransportSpec: the organ is the coronary bed; do not also pass `coronary`"))
    s.coronary === nothing && s.bed !== :coronary && throw(ArgumentError("TransportSpec: the body needs its coronary bed: pass `coronary` when the organ is not the heart"))
    organ = (; ml = s.tissue_ml, s.phys.fb, s.phys.ve, s.phys.ps, s.phys.rho)
    beds = merge((; coronary = s.coronary), NamedTuple{(s.bed,)}((organ,))) # the organ takes its own segment
    BC = human_circulation(; co = s.co, beat = s.beat, hct = s.phys.hd, f_cor = s.f_coronary, beds...)
    any(sg -> sg.name === s.input, BC.segs) || throw(ArgumentError("TransportSpec: input :$(s.input) is not a body segment"))
    return BC
end

"""
    transport(spec::TransportSpec, injection = Injection()) -> result of `Body.simulate`

Iodine everywhere in the body over time for one injection (dose, rate, saline flush, site): `t`, the outflow
concentration of every segment (`out`), the mean blood and interstitial concentrations of every bed (`mean`, `ecv`), and
the iodine injected, held and excreted by the kidneys.
"""
transport(s::TransportSpec, inj::Injection = Injection()) = simulate(circulation(s), inj; tmax = s.tmax, dt = s.dt, every = s.every)

using .Flow: Circuit, InputCurve, frame_masses, sinus_outflow

# Linear interpolation of samples y at times ts, held flat outside.
_interp(t, ts, y) = t <= ts[1] ? y[1] : t >= ts[end] ? y[end] :
    (k = searchsortedlast(ts, t); w = (t - ts[k]) / (ts[k + 1] - ts[k]); y[k] + w * (y[k + 1] - y[k]))

"""
    transport(spec::TransportSpec, C::Circuit, injection = Injection(); passes = 10, tol = 1e-3, fine = 0.02, dt_out = 0.1,
              times = 0:60, to = identity)
        -> (; body, input, t_out, outflow, artery, vein, tissue, passes, change)

The injection through the whole body and the organ's own solved circuit `C` (from [`flow`](@ref)). The organ's arterial
input is the plasma concentration leaving body segment `spec.input`. The organ's venous outflow is fed back into the
body in place of its lumped segment `spec.bed`, so recirculation carries the organ's own transit; that changes the
organ's input, so the feedback is repeated until the outflow changes by at most `tol` of its peak (at most `passes`
times; `passes = 0` is the plain chain with the lumped segment). Every curve returned is then from the same final body
run: the organ's input, its outflow and the body's segment that receives it agree. Returns the body result, the organ's
plasma input curve, its venous plasma outflow concentration every `dt_out` s, the iodine (mg) in every artery, vein and
capillary bed at `times` (host arrays), the passes run and the last pass's relative change.
"""
function transport(
        s::TransportSpec, C::Circuit, inj::Injection = Injection(); passes::Integer = 10, tol::Real = 1.0f-3, fine::Real = 0.02f0, dt_out::Real = 0.1f0,
        times = 0.0f0:1.0f0:60.0f0, to = identity
    )
    BC = circulation(s); pf = 1 - s.phys.hd
    fine = Float32(fine); dt_out = Float32(dt_out)
    r = simulate(BC, inj; tmax = s.tmax, dt = s.dt, every = s.every)
    curve_in(r) = InputCurve(Float32[_interp(t, r.t, r.out[s.input]) / pf for t in collect(0.0f0:fine:s.tmax)], fine) # plasma
    tc = collect(0.0f0:dt_out:s.tmax); fpt = sum(C.fp)
    outflow(inp) = (cum = sinus_outflow(C, inp, tc; to); [k == 1 ? 0.0f0 : (cum[k] - cum[k - 1]) / dt_out / fpt for k in eachindex(tc)])
    cs = outflow(curve_in(r)); npass = 0; change = 0.0f0
    while npass < passes
        fb = cs # the outflow fed back in this pass
        r = simulate(BC, inj; tmax = s.tmax, dt = s.dt, every = s.every, detailed = Dict(s.bed => t -> pf * _interp(Float32(t), tc, fb)))
        cs = outflow(curve_in(r)); npass += 1
        change = maximum(abs.(cs .- fb)) / max(maximum(abs.(cs)), floatmin(Float32))
        change <= tol && break
    end
    inp = curve_in(r)
    am, vm, tm = Array.(frame_masses(C, inp, collect(Float32, times); to))
    return (; body = r, input = inp, t_out = tc, outflow = cs, artery = am, vein = vm, tissue = tm, passes = npass, change)
end
