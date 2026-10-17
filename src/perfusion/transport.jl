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
