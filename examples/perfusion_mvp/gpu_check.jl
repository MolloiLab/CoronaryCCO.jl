# SPDX-License-Identifier: GPL-3.0-or-later
# Runs the transport kernels on the CPU and on a Metal GPU and compares them.
using Vasculature, Metal, Printf
@eval(import Vasculature: Flow); @eval(import Vasculature: Territory)
using .Flow, .Territory

n = parse(Int, get(ARGS, 1, "20000"))
tree = load_vtree(joinpath(@__DIR__, "out", "tree_$n.bin"))
nchild = zeros(Int, length(tree.r)); foreach(p -> p > 0 && (nchild[p] += 1), tree.parent)
terms = findall(==(0), nchild)
terr = fill(255.0f0 / length(terms), length(terms))
C = circuit(tree, terr; state = :stress, phys = HEART_FLOW)
inp = gamma_input()
times = 0.0f0:0.5f0:30.0f0
frame_masses(C, inp, times[1:2]); frame_masses(C, inp, times[1:2]; to = MtlArray)   # compile
tc = @elapsed (a, v, t) = frame_masses(C, inp, times)
tg = @elapsed begin
    ag, vg, tg_ = frame_masses(C, inp, times; to = MtlArray)
    Metal.synchronize()
end
ag, vg, tg_ = Array(ag), Array(vg), Array(tg_)
rel(x, y) = maximum(abs.(x .- y)) / maximum(abs.(x))
@printf("%d segments, %d terminals, %d frames\n", length(tree.r), length(terms), length(times))
@printf("CPU (%d threads) %.3f s | Metal %.3f s\n", Threads.nthreads(), tc, tg)
@printf(
    "max rel diff: artery %.1e, vein %.1e, tissue %.1e; totals CPU %.2f vs GPU %.2f mg at t=15 s\n",
    rel(a, ag), rel(v, vg), rel(t, tg_), sum(a[31, :]) + sum(v[31, :]) + sum(t[31, :]), sum(ag[31, :]) + sum(vg[31, :]) + sum(tg_[31, :])
)
