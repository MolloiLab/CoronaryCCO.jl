# SPDX-License-Identifier: GPL-3.0-or-later
# The one transport entry point: a contrast injection through the whole-body circulation, with the organ's own bed.

using .Body: Injection, human_circulation, simulate

"""
    TransportSpec(; phys, tissue_ml, co = 81.6, beat = 1, f_organ = 0.033, tmax = 70, dt = 0.002, every = 0.1)

Everything that sets where an injected bolus goes. `phys` is the organ's [`FlowSpec`](@ref): its blood hematocrit (`hd`,
used for the whole body), blood and interstitial fractions, contrast permeability and tissue density. `tissue_ml` is the
organ's tissue volume (mL). `co` is the cardiac output (mL/s; rest 81.6, adenosine stress `Body.adenosine()`), `beat` the
heart period (s), and `f_organ` the organ bed's share of the cardiac output (for the heart: its root flow / `co`).
`tmax`, `dt` and `every` are the simulated time, the step and the output spacing (s).
"""
Base.@kwdef struct TransportSpec
    phys::FlowSpec
    tissue_ml::Float32
    co::Float32 = 81.6f0
    beat::Float32 = 1.0f0
    f_organ::Float32 = 0.033f0
    tmax::Float32 = 70.0f0
    dt::Float32 = 0.002f0
    every::Float32 = 0.1f0
end

"""
    circulation(spec::TransportSpec) -> Body.Circulation

The whole-body circulation (`Body.human_circulation`) with the organ's lumped bed and its blood.
"""
circulation(s::TransportSpec) = human_circulation(;
    co = s.co, beat = s.beat, hct = s.phys.hd, f_cor = s.f_organ,
    cor = (; ml = s.tissue_ml, s.phys.fb, s.phys.ve, s.phys.ps, s.phys.rho),
)

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
    transport(spec::TransportSpec, C::Circuit, injection = Injection(); input = :lv, bed = :coronary, passes = 1,
              fine = 0.02, dt_out = 0.1, times = 0:60, to = identity)
        -> (; body, input, t_out, outflow, artery, vein, tissue)

The injection through the whole body and the organ's own solved circuit `C` (from [`flow`](@ref)). The organ's arterial
input is the plasma concentration leaving body segment `input` (`:lv` for an organ fed by the aorta, `:rv` for the lungs).
`passes` times, the organ's venous outflow is fed back into the body in place of its lumped bed `bed`, so recirculation
carries the organ's own transit. Returns the body result, the organ's plasma input curve, its venous plasma outflow
concentration every `dt_out` s, and the iodine (mg) in every artery, vein and capillary bed at `times` (host arrays).
"""
function transport(
        s::TransportSpec, C::Circuit, inj::Injection = Injection(); input::Symbol = :lv, bed::Symbol = :coronary,
        passes::Integer = 1, fine::Real = 0.02f0, dt_out::Real = 0.1f0, times = 0.0f0:1.0f0:60.0f0, to = identity
    )
    BC = circulation(s); pf = 1 - s.phys.hd
    fine = Float32(fine); dt_out = Float32(dt_out)
    r = simulate(BC, inj; tmax = s.tmax, dt = s.dt, every = s.every)
    curve_in(r) = InputCurve(Float32[_interp(t, r.t, r.out[input]) / pf for t in collect(0.0f0:fine:s.tmax)], fine) # plasma
    tc = collect(0.0f0:dt_out:s.tmax); fpt = sum(C.fp)
    outflow(inp) = (cum = sinus_outflow(C, inp, tc; to); [k == 1 ? 0.0f0 : (cum[k] - cum[k - 1]) / dt_out / fpt for k in eachindex(tc)])
    for _ in 1:passes
        cs = outflow(curve_in(r))
        r = simulate(BC, inj; tmax = s.tmax, dt = s.dt, every = s.every, detailed = Dict(bed => t -> pf * _interp(Float32(t), tc, cs)))
    end
    inp = curve_in(r); cs = outflow(inp)
    am, vm, tm = Array.(frame_masses(C, inp, collect(Float32, times); to))
    return (; body = r, input = inp, t_out = tc, outflow = cs, artery = am, vein = vm, tissue = tm)
end
