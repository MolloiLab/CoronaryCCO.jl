# SPDX-License-Identifier: GPL-3.0-or-later
# Smoke test: grow a small tree on case0, build the closed circuit, transport a bolus, and check
# that iodine is conserved (inflow = in network + out through the sinus).
using Vasculature, Printf
include("flow.jl"); include("territory.jl")
using .Flow, .Territory

const S = joinpath(@__DIR__, "..", "..", "test", "data")
n = parse(Int, get(ARGS, 1, "20000"))
treefile = joinpath(@__DIR__, "out", "tree_$n.bin")
mkpath(dirname(treefile))
if isfile(treefile)
    tree = load_vtree(treefile)
else
    t = @elapsed E = generate(
        joinpath(S, "case0.vol"), joinpath(S, "case0.xml"); nterm = n, aperf = 10000, min_dist = 1,
        root = (-26, 41, 13), seed = 1, T = Float32, gamma = 2.33, volume = :full, metric = :euclidean,
        sampling = :continuous, fix_base = true, collision = :final, max_branch_angle = 90
    )
    @printf("grew %d terminals in %.1f s\n", n, t)
    tree = vtree(E); save_vtree(treefile, tree)
end

# territories: nearest terminal for every 1 mm mask voxel
img, lower = read_vol(joinpath(S, "case0.vol"))
idx = findall(>(128), img)
q = Matrix{Float32}(undef, 3, length(idx))
for (v, I) in enumerate(idx), d in 1:3
    q[d, v] = lower[d] + I[d] - 1
end
nchild = zeros(Int, length(tree.r)); foreach(p -> p > 0 && (nchild[p] += 1), tree.parent)
terms = findall(==(0), nchild)
tpts = tree.b[:, terms]
tv = @elapsed owner = nearest_terminal(tpts, q)
terr = zeros(Float32, length(terms))
for o in owner
    terr[o] += 1.0f-3      # 1 mm³ = 1e-3 mL
end
@printf("territories: %d voxels → %d terminals in %.2f s (empty %d)\n", length(idx), length(terms), tv, count(==(0), terr))
terr .= max.(terr, 1.0f-4)

for label in ("rest", "stress")
    tc = @elapsed C = circuit(tree, terr; state = Symbol(label), phys = Physiology(organ = "heart"))
    inp = gamma_input()
    times = 0.0f0:1.0f0:30.0f0
    tf = @elapsed art, ven, tis = frame_masses(C, inp, times)
    tf2 = @elapsed frame_masses(C, inp, times)
    out = sinus_outflow(C, inp, times)
    qp = sum(C.q[j] for j in 1:length(C.q) if C.tree.parent[j] == 0) * (1 - C.phys.hd)
    inflow = [qp * Flow._lin(inp.A, Float32(t), inp.dt) for t in times]
    inside = vec(sum(art; dims = 2) .+ sum(ven; dims = 2) .+ sum(tis; dims = 2))
    err = maximum(abs.(inflow .- inside .- out)) / maximum(inflow)
    mtt = (sum(C.vp_a) + sum(C.vp_v) + sum(C.fp .* (C.ta .+ C.tc .+ C.tv))) / qp
    vo = C.volumes; @printf(
        "  blood %.1f mL: art %.1f resolved + %.1f lumped, cap %.1f, ven %.1f resolved + %.1f lumped\n",
        vo.blood, vo.art_resolved, vo.art_unresolved, vo.capillary, vo.ven_resolved, vo.ven_unresolved
    )
    @printf(
        "%s: mass %.1f g, Q %.2f mL/s (MBF %.2f mL/min/g) (art %.0f%% cap %.0f%% ven %.0f%%), venous scale %.2f\n",
        label, C.mass_g, qp / (1 - C.phys.hd), C.resist.mbf, 100C.resist.art_share, 100C.resist.cap_share, 100C.resist.ven_share, C.venous_scale
    )
    ta = [C.t_end[j] for j in C.terminals]; sa = [C.offset[i] - C.t_end[C.terminals[i]] - C.ta[i] for i in eachindex(C.terminals)]
    @printf(
        "  arterial transit to terminals: median %.2f s (5–95%% %.2f–%.2f); capillary %.2f s; venous %.2f s; plasma MTT (vascular) %.2f s\n",
        sort(ta)[end ÷ 2], sort(ta)[end ÷ 20], sort(ta)[19end ÷ 20], sort(C.tc)[end ÷ 2], sort(sa)[end ÷ 2], mtt
    )
    first_out = times[findfirst(>(1.0f-3 * maximum(inflow)), out)]
    @printf(
        "  circuit %.2f s, frames %.2f s (warm %.3f s) for %d times; conservation error %.1e; first sinus outflow at %.0f s\n",
        tc, tf, tf2, length(times), err, first_out
    )
    for t in (10, 15, 20)
        f = findfirst(==(t), times)
        @printf(
            "  t=%2d s: in %.1f mg | artery %.1f vein %.1f tissue %.1f | left %.1f (%.1f%%)\n", t, inflow[f], sum(art[f, :]),
            sum(ven[f, :]), sum(tis[f, :]), out[f], 100out[f] / inflow[f]
        )
    end
end
