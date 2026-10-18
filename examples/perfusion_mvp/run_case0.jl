# SPDX-License-Identifier: GPL-3.0-or-later
# The four stages end to end on HCT-CCO's public case0 sample (no licensed data): grow -> flow -> transport -> phantom,
# at rest and at stress. Writes a DynamicPhantom bundle (compartment volume + label_healthy.csv) with its ground truth.
#
#   julia -t 8 --project=. run_case0.jl [nterm] [h_mm] [out_dir]
using Vasculature, Printf, Downloads
import Vasculature as V

n = parse(Int, get(ARGS, 1, "2000")); h = parse(Float64, get(ARGS, 2, "0.5"))
out = get(ARGS, 3, joinpath(@__DIR__, "out", "case0_phantom_$n"))
const UPSTREAM = "https://raw.githubusercontent.com/Sujal-vajire/HCT-CCO/9bc141fe3438312f84d5e878cc9830b166ad63fd/samples"
data(name) = (p = joinpath(@__DIR__, "..", "..", "test", "data", name); isfile(p) || (mkpath(dirname(p)); Downloads.download("$UPSTREAM/$name", p)); p)

# 1. grow: the heart's instance, varied to this sample's size
spec = GrowthSpec(HEART_GROWTH; nterm = n, tissue_ml = 100.0, wall_clearance = 1.0, root_caps = nothing, symmetry_floor = Tuple{Float64, Float64}[])
t = @elapsed tree = grow(spec, data("case0.vol"), [data("case0.xml")])
vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
@printf("grow: %d segments, %d terminals in %.1f s\n", length(vt.r), count(==(0), nch), t)

# the tissue: the sample's domain voxels (1 mm), each shared among the terminals that serve it
img, lower = read_vol(data("case0.vol")); vox = findall(>(0x80), img)
q = Float32[lower[d] + I[d] - 1 for d in 1:3, I in vox]
terr, pairs = V.Territory.territory_shares(Float32.(vt.b[:, nch .== 0]), q)

# 2. flow, from pressures, at rest and at maximal dilation
Cr = V.Flow.circuit(vt, terr; phys = HEART_FLOW, state = :rest); Cs = V.Flow.circuit(vt, terr; phys = HEART_FLOW, state = :stress)
@printf("flow: MBF %.2f mL/min/g at rest, %.2f at stress\n", Cr.resist.mbf, Cs.resist.mbf)

# 3. transport: a clinical injection through the whole body into this tree
ph = PhantomSpec(
    HEART_PHANTOM; tissue_classes = 16, vessel_classes = 8,
    pools = [Pool("LV blood", "ncat_blood", :lv, :mean), Pool("aorta", "ncat_blood", :lv, :out)]
)
scans = [(st, C, transport(TransportSpec(HEART_TRANSPORT; tissue_ml = sum(terr)), C; times = ph.times)) for (st, C) in (("rest", Cr), ("stress", Cs))]

# 4. phantom: the classifier is the only subject-specific step (here: domain = tissue, a slab beside it = aorta, else air)
inside = Set(Tuple(I) for I in vox); row = Dict(Tuple(I) => Int32(k) for (k, I) in enumerate(vox))
cell(p) = Tuple(round.(Int, p .- lower) .+ 1)
lo = lower .- 3.0; hi = lower .+ size(img) .+ 3.0
classify(p) = (
    c = cell(p); c in inside ? (VOX_TISSUE, row[c], 0.0f0) :
        p[1] < lo[1] + 2 ? (VOX_POOL, Int32(2), 0.0f0) : (VOX_STATIC, Int32(0), 0.0f0)
)
grid = V.VoxelGrid(Tuple(lo), h, Tuple(ceil.(Int, (hi .- lo) ./ h)))
t = @elapsed B = phantom(ph, scans, classify, grid; pairs, terr, nvox = length(vox), out, name = "case0")
@printf("phantom: %d compartments in %s (%.1f s)\n", length(B.ids), out, t)
for (st, f) in B.flow
    @printf("  ground truth %s: %.3f mL/min per mL of tissue\n", st, f)
end
