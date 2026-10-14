# SPDX-License-Identifier: GPL-3.0-or-later
# Times the growth of an n-terminal tree on HCT-CCO's case0 sample.
#   julia -t auto --project=bench bench/run.jl 40000 [Float32|Float64] [upstream|full]
using CoronaryCCO
using Downloads

const UPSTREAM = "https://raw.githubusercontent.com/Sujal-vajire/HCT-CCO/9bc141fe3438312f84d5e878cc9830b166ad63fd/samples"
data(name) = (p = joinpath(@__DIR__, "..", "test", "data", name); isfile(p) || (mkpath(dirname(p)); Downloads.download("$UPSTREAM/$name", p)); p)

n = parse(Int, get(ARGS, 1, "40000"))
T = get(ARGS, 2, "Float32") == "Float64" ? Float64 : Float32
volume = Symbol(get(ARGS, 3, "upstream"))
run(n) = generate(
    data("case0.vol"), data("case0.xml"); nterm = n, aperf = 10000, min_dist = 1,
    root = (-26, 41, 13), seed = 1, T = T, volume = volume
)
run(100)                                   # compile
t = @elapsed E = run(n)
st = E.stats
println(
    "n=$n T=$T volume=$volume threads=$(Threads.nthreads()): $(round(t, digits = 2)) s, ",
    "$(nsegments(E.m)) segments, root radius $(E.m.radius[2])"
)
println(
    "  location $(round(st.t_location, digits = 2)) s, nearest $(round(st.t_nearest, digits = 2)) s, ",
    "candidates $(round(st.t_candidates, digits = 2)) s, commit $(round(st.t_commit, digits = 2)) s"
)
