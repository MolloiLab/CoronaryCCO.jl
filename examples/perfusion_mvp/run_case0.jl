# SPDX-License-Identifier: GPL-3.0-or-later
# End-to-end ground truth on the HCT-CCO case0 sample: grown tree → closed circuit → iodine
# transport → voxel phantoms (labels + iodine mg/mL) at chosen times, plus the exact truth
# (flow, myocardial mass, iodine inside the mask, sinus outflow, input integral).
#
#   julia -t 8 --project=. run_case0.jl <nterm> [h_mm]
using Vasculature, Printf, TOML, Unitful
include("flow.jl"); include("territory.jl"); include("phantom.jl")
using .Flow, .Territory, .Phantom

const DATA = joinpath(@__DIR__, "..", "..", "test", "data")
n = parse(Int, get(ARGS, 1, "20000")); h = parse(Float32, get(ARGS, 2, "0.5"))
tree = load_vtree(joinpath(@__DIR__, "out", "tree_$n.bin"))
nchild = zeros(Int, length(tree.r)); foreach(p -> p > 0 && (nchild[p] += 1), tree.parent)
terms = findall(==(0), nchild)
img, lower = read_vol(joinpath(DATA, "case0.vol"))
mask = img .> 128
tg = @elapsed g = phantom_grid(mask, lower, tree.b[:, terms]; h, nearest = nearest_terminal)
vv = h^3 / 1000
terr = g.count .* Float32(vv)
@printf(
    "grid %s at %.2f mm (%.1f M voxels) in %.1f s; %d myocardial voxels, %d cavity voxels, empty territories %d\n",
    g.dims, h, prod(g.dims) / 1.0e6, tg, count(g.myo), count(==(0x03), g.labels), count(==(0), g.count)
)
terr .= max.(terr, 1.0f-6)

# segments whose midpoint lies in the myocardium (inside the myocardial mask)
inside = map(1:length(tree.r)) do j
    i, k, l = Phantom._vox(g, (tree.a[1, j] + tree.b[1, j]) / 2, (tree.a[2, j] + tree.b[2, j]) / 2, (tree.a[3, j] + tree.b[3, j]) / 2)
    all((1, 1, 1) .<= (i, k, l) .<= g.dims) && g.myo[i, k, l]
end
mark_lumen!(g, tree, sqrt(FlowSpec(organ = "heart").share_ven / FlowSpec(organ = "heart").share_art))
@printf("lumen: %d myocardial voxels are resolved vessel lumen\n", count(g.myo .& (g.labels .== 0x03)))

inp = gamma_input(; t0 = 2.0f0, tpeak = 10.0f0, peak_blood = 12.0f0)
hd = 0.45f0
cblood(t) = (1 - hd) * Flow._lin(inp.c, Float32(t) + 0.3f0, inp.dt) / 1     # LV cavity leads the ostia by ~0.3 s
Cb(t) = (1 - hd) * (t <= 0 ? 0.0f0 : inp.c[min(length(inp.c), Int(floor(t / inp.dt)) + 1)])
scan_times = Float32[0, 6, 8, 10, 12, 14]
truth_times = 0.0f0:0.25f0:30.0f0

for state in ("rest", "stress")
    dir = joinpath(@__DIR__, "out", "case0_$(n)_$(state)"); mkpath(dir)
    C = circuit(tree, terr; state = Symbol(state), phys = FlowSpec(organ = "heart"))
    t1 = @elapsed art, ven, tis = frame_masses(C, inp, truth_times)
    out = sinus_outflow(C, inp, truth_times)
    m_in = [sum(@view tis[f, :]) + sum(@view(art[f, inside])) + sum(@view(ven[f, inside])) for f in eachindex(truth_times)]
    m_all = [sum(@view tis[f, :]) + sum(@view art[f, :]) + sum(@view ven[f, :]) for f in eachindex(truth_times)]
    qp = sum(C.q[j] for j in eachindex(C.q) if tree.parent[j] == 0) * (1 - hd)
    infl = [qp * Flow._lin(inp.A, t, inp.dt) for t in truth_times]
    # voxelize the scan frames
    files = String[]; conc = zeros(Float32, g.dims); m_vox = Float64[]; m_mat = Float64[]; nmat = Int[]; wmyo = Float32[]
    t2 = @elapsed for t in scan_times
        f = findfirst(==(t), truth_times)
        deposit!(conc, g, tree, art, ven, tis, f, C.venous_scale; blood_conc = cblood(t))
        name = @sprintf("iodine_t%03d.f32", round(Int, 10t)); push!(files, name)
        write_frame(joinpath(dir, name), conc)
        push!(m_vox, sum(Float64(conc[v]) for v in eachindex(conc) if g.myo[v]) * vv)
        idx, mats, widths = materialize(g.labels, conc; seed = round(Int, 10t) + 1)
        tag = @sprintf("t%03d", round(Int, 10t))
        write(joinpath(dir, "index_$tag.u16"), idx); write_materials(joinpath(dir, "materials_$tag.toml"), mats)
        # iodine carried by the materialized myocardium (density excess × voxel volume)
        mq = zeros(Float64, length(mats))
        for (k, m) in enumerate(mats)
            mq[k] = get(m.composition, 53, 0.0) * ustrip(u"g/cm^3", m.density)
        end
        push!(m_mat, sum(mq[idx[v]] for v in eachindex(idx) if g.myo[v]) * vv * 1000)
        push!(nmat, length(mats)); push!(wmyo, widths[0x02])
    end
    write_meta(dir, g, scan_times, files)
    truth = Dict(
        "state" => state, "mbf_ml_min_g" => Float64(C.resist.mbf), "myocardial_mass_g" => Float64(C.mass_g),
        "flow_ml_s" => Float64(qp / (1 - hd)), "hematocrit" => Float64(hd),
        "vascular_mtt_s" => Float64((sum(C.vp_a) + sum(C.vp_v) + sum(C.fp .* (C.ta .+ C.tc .+ C.tv))) / qp),
        "input" => Dict("kind" => "gamma variate at the coronary ostia", "t0_s" => 2.0, "tpeak_s" => 10.0, "peak_blood_mg_ml" => 12.0),
        "curves" => Dict(
            "t_s" => collect(Float64.(truth_times)),
            "inflow_mg" => Float64.(infl), "network_mg" => Float64.(m_all), "mask_mg" => Float64.(m_in),
            "sinus_outflow_mg" => Float64.(out), "blood_conc_lv_mg_ml" => Float64.(cblood.(truth_times)),
        ),
        "scan_frames" => Dict(
            "t_s" => collect(Float64.(scan_times)), "mask_mg_voxelized" => m_vox,
            "mask_mg_materialized" => m_mat, "n_materials" => nmat, "myocardium_bin_mg_ml" => Float64.(wmyo)
        ),
    )
    open(io -> TOML.print(io, truth), joinpath(dir, "truth.toml"), "w")
    @printf("%s: transport %.2f s (%d times), voxelize+write %.1f s (%d frames) → %s\n", state, t1, length(truth_times), t2, length(scan_times), dir)
    for t in (6, 10, 14)
        f = findfirst(==(Float32(t)), truth_times); k = findfirst(==(Float32(t)), scan_times)
        @printf(
            "  t=%2d s: inflow %.1f mg, in mask %.1f (voxelized %.1f, materialized %.1f; %d materials, myo bin %.3f mg/mL), left via sinus %.1f (%.1f%%)\n",
            t, infl[f], m_in[f], m_vox[k], m_mat[k], nmat[k], wmyo[k], out[f], 100out[f] / max(infl[f], 1.0f-9)
        )
    end
end
