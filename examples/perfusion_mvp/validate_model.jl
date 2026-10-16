# SPDX-License-Identifier: GPL-3.0-or-later
# Model-side validation metrics for the closed circuit, to plot against experimental data
# (data/reference/). Writes one JSON per state.
#
#   julia -t 4 --project=. validate_model.jl <tree.bin> <out.json> [mask] [veins.bin]
# mask: omitted → case0 sample; or a directory holding heart_labels_1mm.u8 + grid.txt (XCAT prep),
# whose labels 1 (LV myocardium) and 2 (RV myocardium) form the perfused domain.
# veins.bin: a grown venous tree (each terminal drains into its nearest venous terminal); omitted →
# companion veins. Flow is pressure-driven: stress = maximal vasodilation, rest = autoregulated
# arteriolar diameters (myogenic, shear and metabolic tone).
using Vasculature, Printf, Serialization
(isdefined(Main, :Flow) || include("flow.jl")); (isdefined(Main, :Territory) || include("territory.jl")); include("coupling.jl")
using .Flow, .Territory

const DATA = joinpath(@__DIR__, "..", "..", "test", "data")
treefile = get(ARGS, 1, joinpath(@__DIR__, "out", "tree_20000.bin")); outjson = get(ARGS, 2, joinpath(@__DIR__, "out", "validate_20000.json"))
tree = load_vtree(treefile); ns = length(tree.r); n = ns
nchild = zeros(Int, ns); foreach(p -> p > 0 && (nchild[p] += 1), tree.parent)
terms = findall(==(0), nchild)
if length(ARGS) >= 3
    gl = split.(readlines(joinpath(ARGS[3], "grid.txt"))); lower = parse.(Int, gl[1]); gd = Tuple(parse.(Int, gl[2]))
    lab = Array{UInt8}(undef, gd); read!(joinpath(ARGS[3], "heart_labels_1mm.u8"), lab); mask = (lab .== 0x01) .| (lab .== 0x02)
else
    img, lower = read_vol(joinpath(DATA, "case0.vol")); mask = img .> 128
end
vidx = findall(mask); q = Matrix{Float32}(undef, 3, length(vidx))
for (v, I) in enumerate(vidx), d in 1:3
    q[d, v] = lower[d] + I[d] - 1
end
terr, tshare = territory_shares(tree.b[:, terms], q)
terr .= max.(terr, 1.0f-7)
veins = length(ARGS) >= 4 ? load_vtree(ARGS[4]) : nothing
vleaf = veins === nothing ? nothing : haskey(ENV, "CIRC") ? deserialize(ENV["CIRC"]).Cr.vleaf : (
        vn = zeros(Int, length(veins.r)); foreach(p -> p > 0 && (vn[p] += 1), veins.parent);
        couple_veins(tree, terr, veins, terms, findall(==(0), vn))[1]
    )
# coverage: distance from every myocardial voxel to its nearest terminal, and terminal density
own1 = nearest_terminal(tree.b[:, terms], q)
gapd = [sqrt(sum(abs2, q[:, v] .- tree.b[:, terms[own1[v]]])) for v in 1:size(q, 2)]
sg = sort(gapd)
dens = Dict{NTuple{3, Int}, Int}(); for j in terms
    k = Tuple(floor.(Int, tree.b[:, j] ./ 4)); dens[k] = get(dens, k, 0) + 1
end
vcnt = Dict{NTuple{3, Int}, Int}(); for v in 1:size(q, 2)
    k = Tuple(floor.(Int, q[:, v] ./ 4)); vcnt[k] = get(vcnt, k, 0) + 1
end
dv = [get(dens, k, 0) / c for (k, c) in vcnt if c >= 48]          # terminals per mm³ in well-filled 4 mm cubes
coverage = Dict(
    "gap_p99_mm" => sg[ceil(Int, 0.99 * end)], "gap_max_mm" => sg[end], "frac_voxels_gt_1p5mm" => count(>(1.5), gapd) / length(gapd),
    "density_cv_4mm" => (m_ = sum(dv) / length(dv); sqrt(sum((dv .- m_) .^ 2) / length(dv)) / m_), "terminals_per_mm3" => length(terms) / size(q, 2)
)
@printf(
    "coverage: nearest-terminal distance p99 %.2f mm, max %.2f mm, voxels > 1.5 mm %.3f%%; terminal density CV (4 mm cubes) %.2f\n",
    coverage["gap_p99_mm"], coverage["gap_max_mm"], 100coverage["frac_voxels_gt_1p5mm"], coverage["density_cv_4mm"]
)

const P_AO = Physiology().p_in     # mean aortic pressure, mmHg (pressures are reported relative to it)

js(x) = x isa AbstractString ? "\"" * x * "\"" : x isa AbstractDict ? "{" * join(["\"$k\":" * js(v) for (k, v) in x], ",") * "}" :
    x isa Union{AbstractVector, Tuple} ? "[" * join(js.(collect(x)), ",") * "]" : x isa AbstractFloat ? (isfinite(x) ? @sprintf("%.5g", x) : "null") : string(x)

logbins(lo, hi, k) = exp10.(range(log10(lo), log10(hi), length = k + 1))
function binstats(x, y, edges; w = ones(length(x)))
    out = []
    for b in 1:(length(edges) - 1)
        sel = findall(i -> edges[b] <= x[i] < edges[b + 1], eachindex(x))
        isempty(sel) && continue
        ys = sort(y[sel])
        push!(
            out, Dict(
                "x" => sqrt(edges[b] * edges[b + 1]), "n" => length(sel), "median" => ys[(end + 1) ÷ 2],
                "p10" => ys[max(1, ceil(Int, 0.1 * end))], "p90" => ys[max(1, ceil(Int, 0.9 * end))]
            )
        )
    end
    return out
end

res = Dict{String, Any}("n" => n, "segments" => ns, "terminals" => length(terms), "coverage" => coverage)
# maximal dilation first: its microvascular scale κ is shared by the rest state
Cstress = haskey(ENV, "CIRC") ? nothing : circuit(tree, terr; state = :stress, veins, vleaf, phys = Physiology(radius_cv = parse(Float32, get(ENV, "RADIUS_CV", "0.12"))))
for (st, state) in (("rest", :rest), ("stress", :stress))
    rcv = parse(Float32, get(ENV, "RADIUS_CV", "0.12"))
    C = if haskey(ENV, "CIRC")                 # solved circuits (; Cs, Cr) from an earlier run on the same inputs
        (cc = deserialize(ENV["CIRC"]); state === :rest ? cc.Cr : cc.Cs)
    else
        state === :rest ? circuit(tree, terr; state, veins, vleaf, phys = Physiology(radius_cv = rcv), kappa = Cstress.resist.kappa) : Cstress
    end
    hd = C.phys.hd
    d_um = 2000 .* tree.r
    pd = C.pa                                                     # mmHg at each segment's distal end
    vel = C.q ./ (Float32(π) .* (tree.r .* C.xd ./ 10) .^ 2)      # cm/s, at the solved diameters
    e = logbins(8, 4000, 18)
    # regional flow heterogeneity: flow per gram in cubes of side L (Bassingthwaighte et al. 1989)
    fvox = zeros(Float64, size(q, 2))                                # flow per voxel (mL/s per mm³)
    dvox = zeros(Float64, size(q, 2))                                # relative O₂ demand per voxel
    for (v, i, fr) in tshare
        fvox[v] += fr * C.q[C.terminals[i]] / (terr[i] * 1.0f3)
        dvox[v] += fr * C.demand[i] / (terr[i] * 1.0f3)
    end
    # pieces: 1 mm voxels, then Voronoi cells of K random voxel seeds (equal mean volume; they
    # follow the thin wall, unlike cubes). RD is mass-weighted.
    rd = []; rdd = []
    nv = length(vidx)
    for target_mm3 in (1, 8, 27, 64, 216, 1000, 4000, 16000)
        if target_mm3 == 1
            lab = collect(1:nv); K = nv
        else
            K = max(4, nv ÷ target_mm3)
            seeds = q[:, round.(Int, range(1, nv, length = K))] .+ 0.1f0 .* randn(Float32, 3, K)
            lab = nearest_terminal(seeds, q; h = 4.0f0)
        end
        for (vox, out) in ((fvox, rd), (dvox, rdd))
            fsum = zeros(Float64, K); msum = zeros(Float64, K)
            for v in 1:nv
                fsum[lab[v]] += vox[v]; msum[lab[v]] += 1
            end
            ok = msum .> 0; f = fsum[ok] ./ msum[ok]; m = msum[ok]
            fm = sum(f .* m) / sum(m)
            push!(out, Dict("mass_g" => 1.05e-3 * sum(m) / length(m), "rd" => sqrt(sum(m .* (f .- fm) .^ 2) / sum(m)) / fm, "pieces" => length(m)))
        end
    end
    # impulse response at the coronary sinus: a 0.02 s pulse at the ostia
    pulse = Float32[t < 0.02f0 ? 50.0f0 : 0.0f0 for t in 0.0f0:0.005f0:150.0f0]
    inp = Flow.InputCurve(pulse, 0.005f0)
    tt = 0.0f0:0.05f0:150.0f0
    cum = sinus_outflow(C, inp, tt)
    tot = cum[end]; h = diff(cum) ./ step(tt) ./ tot; tm = collect(tt[1:(end - 1)]) .+ step(tt) / 2
    mtt_out = sum(tm .* h) * step(tt)
    appear = tm[findfirst(>(0.05 * maximum(h)), h)]
    tpeak = tm[argmax(h)]
    # intravascular tracer (no extraction), as in labelled-red-cell / albumin indicator dilution
    C0 = with_extraction(C, 0)
    cum0 = sinus_outflow(C0, inp, tt); h0 = diff(cum0) ./ step(tt) ./ cum0[end]
    mtt0 = sum(tm .* h0) * step(tt); appear0 = tm[findfirst(>(0.05 * maximum(h0)), h0)]
    # vascular volumes and the central-volume MTT
    qp = sum(C.q[j] for j in 1:ns if tree.parent[j] == 0) * (1 - hd)
    vblood = C.volumes.blood
    # arrival at the most distal imaged (base-tree) point: TIMI-like ostium → distal landmark time
    base = findall(tree.isbase); tdist = maximum(C.t_end[base])
    # per root: arrival at its imaged distal tips (base segments with no base child)
    rootof = zeros(Int32, ns)
    for h in 1:ns
        j = C.order[h]; p = tree.parent[j]; rootof[j] = p == 0 ? j : rootof[p]
    end
    hasbasechild = falses(ns); for j in base
        p = tree.parent[j]; p > 0 && (hasbasechild[p] = true)
    end
    tips = [j for j in base if !hasbasechild[j]]
    tip_arrival = Dict(string(r) => sort([C.t_end[j] for j in tips if rootof[j] == r]) for r in unique(rootof[tips]))
    root = findfirst(j -> tree.parent[j] == 0, 1:ns)
    res[st] = Dict(
        "mbf" => C.resist.mbf, "kappa" => C.resist.kappa, "r_micro" => C.resist.r_micro, "chilian_shares" => collect(C.resist.chilian_shares), "flow_ml_min" => qp / (1 - hd) * 60, "mass_g" => C.mass_g,
        "pressure_vs_diameter" => binstats(d_um, pd ./ P_AO, e),
        "flow_vs_diameter" => binstats(d_um, C.q .* 60, e),
        "velocity_vs_diameter" => binstats(d_um, vel, e),
        "flow_scatter" => [[d_um[j], C.q[j] * 60] for j in 1:max(1, ns ÷ 4000):ns],
        "rd_vs_size" => rd, "rd_demand_vs_size" => rdd,
        "impulse" => Dict(
            "t" => tm[1:4:min(end, 600)], "h" => h[1:4:min(end, 600)], "mtt_s" => mtt_out, "appear_s" => appear,
            "peak_s" => tpeak, "appear_over_mtt" => appear / mtt_out,
            "h_intravascular" => h0[1:4:min(end, 600)], "mtt_intravascular_s" => mtt0, "appear_intravascular_s" => appear0,
            "appear_over_mtt_intravascular" => appear0 / mtt0, "central_volume_mtt_s" => (sum(C.vp_a) + sum(C.vp_v) + sum(C.fp .* (C.ta .+ C.tc .+ C.tv))) / qp
        ),
        "volumes_ml_per_100g" => Dict(k => v / C.mass_g * 100 for (k, v) in pairs(C.volumes)),
        "blood_fraction" => vblood / (C.mass_g / 1.05f0),
        "ostium_to_distal_base_s" => tdist, "tip_arrival_by_root" => tip_arrival, "root_velocity_cm_s" => vel[root], "root_diameter_mm" => 2tree.r[root],
        "resistance_share" => Dict("art" => C.resist.art_share, "cap" => C.resist.cap_share, "ven" => C.resist.ven_share),
        "so2_sinus" => C.resist.so2_sinus, "mean_activation" => C.resist.mean_act, "ctone" => C.ctone,
        "diameter_ratio_vs_d0" => binstats(d_um[C.xd .< 1], C.xd[C.xd .< 1], logbins(5, 260, 12)),
        "capillary_transit_s" => Dict("median" => sort(C.tc)[end ÷ 2], "p10" => sort(C.tc)[end ÷ 10], "p90" => sort(C.tc)[9end ÷ 10]),
        "terminal_flow_cv" => (f = C.q[C.terminals] ./ terr; sqrt(sum((f .- sum(f) / length(f)) .^ 2) / length(f)) / (sum(f) / length(f))),
    )
    @printf(
        "%s: intravascular MTT %.2f s (central volume %.2f), appearance %.2f s (%.2f × MTT) | with exchange MTT %.2f s; blood %.1f mL/100g, distal base arrival %.2f s, root v %.0f cm/s, RD %.2f (1 mg) → %.2f (%.1f g)\n",
        st, mtt0, res_cv = (sum(C.vp_a) + sum(C.vp_v) + sum(C.fp .* (C.ta .+ C.tc .+ C.tv))) / qp, appear0, appear0 / mtt0, mtt_out, vblood / C.mass_g * 100, tdist, vel[root], rd[1]["rd"], rd[end]["rd"], rd[end]["mass_g"]
    )
end
open(io -> print(io, js(res)), outjson, "w")
println("wrote ", outjson)
