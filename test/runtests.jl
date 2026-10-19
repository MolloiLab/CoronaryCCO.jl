# SPDX-License-Identifier: GPL-3.0-or-later
using Vasculature
using Downloads
using Test

const C = Vasculature

# HCT-CCO's sample case (the HCT-CCO authors' data, under the HCT-CCO license), fetched from
# the pinned upstream commit rather than redistributed here.
const UPSTREAM = "https://raw.githubusercontent.com/Sujal-vajire/HCT-CCO/9bc141fe3438312f84d5e878cc9830b166ad63fd/samples"
const DATA = joinpath(@__DIR__, "data")

function sample(name)
    mkpath(DATA)
    path = joinpath(DATA, name)
    isfile(path) || Downloads.download("$UPSTREAM/$name", path)
    return path
end

grow_case0(n; kw...) = Vasculature.generate(
    sample("case0.vol"), sample("case0.xml"); nterm = n, aperf = 10000,
    min_dist = 1, root = (-26, 41, 13), seed = 1, kw...
)

same_tree(a, b) = a.pos == b.pos && a.radius == b.radius && a.parent == b.parent &&
    a.left == b.left && a.right == b.right && a.flow == b.flow &&
    a.resistance == b.resistance && a.terminals == b.terminals

@testset "Vasculature" begin
    @testset "defaults match upstream HctCco CLI" begin
        p = CCOParams()
        @test p.nterm == 500
        @test p.aperf == 20000.0
        @test p.gamma == 3.0
        @test p.root_radius_max == 1.588
        @test p.min_distance_to_border == 5.0
        @test p.num_neighbors == 20
    end

    @testset "portable rand() matches cco_rand()" begin
        # first draws of SplitMix64 >> 33 from seed 1 (computed independently)
        r = C.CRand(1)
        @test [C.crand!(r) for _ in 1:3] == Int32[1216681718, 1601554128, 2085212535]
        @test all(0 <= C.crand!(r) <= typemax(Int32) for _ in 1:1000)
    end

    @testset "exact Lp distance transform (p = $p)" for p in (2, 3)
        # p = 3 is upstream's metric (dimension used as the exponent); p = 2 is Euclidean
        img = fill(0xff, 7, 6, 5); img[1, 1, 1] = 0x00; img[7, 6, 5] = 0x00; img[4, 2, 3] = 0x80
        sites = [(1, 1, 1), (7, 6, 5), (4, 2, 3)]
        d = C.distance_map(Float64, img; p = p)
        ok = true
        for k in 1:5, j in 1:6, i in 1:7
            raw = minimum(abs(i - a)^p + abs(j - b)^p + abs(k - c)^p for (a, b, c) in sites)
            ok &= d[i, j, k] == C.cpow(Float64(raw), 1.0 / p)
        end
        @test ok
    end

    @testset "case0 domain and base tree" begin
        img, lower = read_vol(sample("case0.vol"))
        @test size(img) == (131, 144, 105)
        @test lower == (Int32(-66), Int32(-70), Int32(-55))      # DGtal centred domain
        m = import_gxl(Float64, sample("case0.xml"), CCOParams(nterm = 500, aperf = 10000))
        @test m.terminals[1] == 2                                 # upstream's stale entry
        @test m.parent[1] == 2                                    # myVectParent[0] = 1
        @test m.kterm_total == length(m.terminals)
        @test m.radius[2] <= 1.588
    end

    @testset "fast mode: parallel = serial, deterministic, tracks exact mode ($T)" for T in (Float64, Float32)
        ser = grow_case0(300; T = T, parallel = false)
        par = grow_case0(300; T = T, parallel = true)
        @test same_tree(ser.m, par.m)           # parallel candidates pick what serial picks
        @test same_tree(par.m, grow_case0(300; T = T).m)   # fixed seed, fixed tree
        ex = grow_case0(300; T = T, mode = :exact)
        @test nsegments(ex.m) == nsegments(par.m)
        # same algorithm up to rounding: work counts agree closely
        @test abs(par.stats.kamiya_iters - ex.stats.kamiya_iters) <= 0.02 * ex.stats.kamiya_iters
        @test abs(par.m.radius[2] - ex.m.radius[2]) <= 1.0e-3 * ex.m.radius[2]
    end

    @testset "corrected growth options" begin
        E = grow_case0(
            2000; T = Float32, gamma = 2.33, volume = :full, metric = :euclidean,
            sampling = :continuous, fix_base = true, collision = :final, max_branch_angle = 90
        )
        q = tree_qc(E; rmin = 0.1)
        @test q.grown_turns_90 == 0          # no daughter turns back past 90°
        @test q.base_turns_135 == 0          # imaged vessels are not folded
        R = dds_orders(E.m)
        @test R.stats[1].count > R.stats[end].count
        @test all(s.mean > 0 for s in R.stats)
    end

    @testset "batched growth" begin
        kw = (;
            T = Float32, gamma = 2.33, volume = :full, metric = :euclidean, sampling = :continuous,
            fix_base = true, collision = :final, max_branch_angle = 90,
        )
        seq = grow_case0(1200; kw...)
        # all terminals before n_seq, or batches of one, are exactly sequential growth
        @test same_tree(seq.m, grow_case0(1200; kw..., growth = :batched, n_seq = 1200).m)
        @test same_tree(seq.m, grow_case0(1200; kw..., growth = :batched, n_seq = 200, batch_frac = 1.0e-9).m)
        bat = grow_case0(3000; kw..., growth = :batched, n_seq = 500)
        # deterministic, and independent of threading
        @test same_tree(bat.m, grow_case0(3000; kw..., growth = :batched, n_seq = 500).m)
        @test same_tree(bat.m, grow_case0(3000; kw..., growth = :batched, n_seq = 500, parallel = false).m)
        # the device argument: an array constructor or a KernelAbstractions backend
        @test same_tree(bat.m, grow_case0(3000; kw..., growth = :batched, n_seq = 500, device = Array).m)
        @test same_tree(bat.m, grow_case0(3000; kw..., growth = :batched, n_seq = 500, device = C.KA.CPU()).m)
        # a unified-memory device takes the tree in place (here the host itself): the same tree
        @test same_tree(bat.m, grow_case0(3000; kw..., growth = :batched, n_seq = 500, device = C.Unified(x -> x, copy)).m)
        ref = grow_case0(3000; kw...)
        @test length(bat.m.terminals) == length(ref.m.terminals)
        vol(E) = sum(Float64(E.m.radius[j])^2 * Float64(E.m.glen[j]) for j in 2:nsegments(E.m))
        @test abs(vol(bat) - vol(ref)) <= 0.02 * vol(ref)
        q = tree_qc(bat; rmin = 0.1)
        @test q.grown_turns_90 == 0
        @test q.base_turns_135 == 0
        @test_throws ArgumentError grow_case0(100; T = Float64, mode = :exact, growth = :batched)
    end

    @testset "finer grid cell, same tree" begin
        kw = (;
            T = Float32, gamma = 2.33, volume = :full, metric = :euclidean, sampling = :continuous,
            fix_base = true, collision = :final, max_branch_angle = 90,
        )
        for growth in (:sequential, :batched)
            E1 = grow_case0(1500; kw..., growth)
            E2 = deepcopy(E1)
            v = E2.vox; E2.vox = 1; C.regrid!(E2); E2.vox = v # force one halving of the cell
            @test E2.g.h == E1.g.h / 2
            grow!(E1; steps = 400, growth, n_seq = 0); grow!(E2; steps = 400, growth, n_seq = 0)
            @test same_tree(E1.m, E2.m)
        end
    end

    @testset "root_fixed, upow, ordering reference" begin
        # the root sits exactly at its cap after growth, and every radius scales by the same factor
        kw = (; T = Float32, volume = :full, sampling = :continuous, collision = :final, root_radius_max = 3.0)
        a = grow_case0(600; kw...); b = grow_case0(600; kw..., root_fixed = true)
        ra = a.m.radius[a.roots[1]]; rb = b.m.radius[b.roots[1]]
        @test rb ≈ 3.0f0
        @test a.m.pos == b.m.pos
        @test all(isapprox.(b.m.radius[2:end] ./ a.m.radius[2:end], rb / ra; rtol = 1.0e-5))
        # C pow semantics for Float32: negative base, non-integer exponent → NaN (no throw)
        @test isnan(Vasculature.upow(-1.0f0, 0.5f0))
        @test Vasculature.upow(-2.0f0, 2.0f0) == 4.0f0
        # the default ordering reference is upstream's Kassab arrays
        @test Vasculature._literature_boundaries(Vasculature.LITERATURE_MU, Vasculature.LITERATURE_SD) == Vasculature._LITERATURE_DB
    end

    @testset "README calls run" begin
        # the README's growth step on the package's sample case: GrowthSpec, grow, vtree, mesh_domain
        tree = grow(GrowthSpec(HEART_GROWTH; nterm = 300, tissue_ml = 100.0, wall_clearance = 1.0, root_caps = nothing, symmetry_floor = Tuple{Float64, Float64}[]), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        @test count(==(0), vt.parent) >= 1 && count(==(0), nch) >= 300                  # a root, and the 300 grown terminals
        X = Float64[0 20 0 0; 0 0 20 0; 0 0 0 20]; T = [1 1 1 2; 3 2 4 3; 2 4 3 4]   # a tetrahedron, 20 mm edges
        img, _ = read_vol(mesh_domain(tempname() * ".vol", X, T))
        @test 1000 < count(>(0x80), img) < 1600                                      # 1333 mm³
    end

    @testset "transport is the whole-body bolus with the organ's bed" begin
        heart = HEART_FLOW
        spec = TransportSpec(HEART_TRANSPORT; tissue_ml = 148.1, tmax = 20, dt = 0.005, every = 1)
        r = transport(spec)
        circ = C.Body.human_circulation(; hct = heart.hd, coronary = (; ml = 148.1f0, heart.fb, heart.ve, heart.ps, heart.rho))
        r2 = C.Body.simulate(circ, C.Body.Injection(); tmax = 20, dt = 0.005, every = 1)
        @test r.t == r2.t && r.out[:lv] == r2.out[:lv] && r.mean[:coronary] == r2.mean[:coronary]
        @test r.injected > 0 && maximum(r.out[:lv]) > 0
        @test circulation(spec).hct == heart.hd
        # with the organ's own circuit: passes = 0 is the plain chain (body -> plasma input -> masses); one pass feeds the
        # organ's venous outflow back into the body in place of its lumped bed
        tree = grow(GrowthSpec(HEART_GROWTH; nterm = 300, tissue_ml = 100.0, wall_clearance = 1.0, root_caps = nothing, symmetry_floor = Tuple{Float64, Float64}[]), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        Cr = C.Flow.circuit(vt, fill(0.1f0, count(==(0), nch)); phys = heart, state = :rest)
        T0 = transport(spec, Cr; passes = 0, times = 0.0f0:2.0f0:20.0f0)
        # the heart's body takes its coronary share from the grown tree's root flow
        qroot = sum(Cr.q[j] for j in eachindex(vt.parent) if vt.parent[j] == 0)
        @test coronary_share(Cr, spec.co) == Float32(qroot / spec.co)
        r = transport(TransportSpec(spec; f_coronary = coronary_share(Cr, spec.co)))
        pf = 1 - heart.hd; tf = collect(0.0f0:0.02f0:20.0f0)
        inp = C.Flow.InputCurve(Float32[C._interp(t, r.t, r.out[:lv]) / pf for t in tf], 0.02f0)
        @test T0.body.out[:lv] == r.out[:lv] && T0.input.c == inp.c
        @test transport(TransportSpec(spec; f_coronary = 0.5f0), Cr; passes = 0, times = 0.0f0:2.0f0:20.0f0).body.out[:lv] == r.out[:lv]
        am, vm, tm = Array.(C.Flow.frame_masses(Cr, inp, collect(0.0f0:2.0f0:20.0f0)))
        @test T0.artery == am && T0.vein == vm && T0.tissue == tm
        T1 = transport(spec, Cr; passes = 1, times = 0.0f0:2.0f0:20.0f0)
        @test T1.body.out[:coronary] != T0.body.out[:coronary] && all(>=(0), T1.tissue)
        # the lungs as the organ: the body's lung segment takes the lung's own blood, interstitium and PS; the coronary bed
        # stays the reference heart's
        lungp = LUNG_FLOW
        hb = (; ml = 148.1f0, heart.fb, heart.ve, heart.ps, heart.rho)
        BL = circulation(TransportSpec(LUNG_TRANSPORT; tissue_ml = 2813, coronary = hb))
        sl = BL.segs[findfirst(x -> x.name === :lung, BL.segs)]; sc = BL.segs[findfirst(x -> x.name === :coronary, BL.segs)]
        @test sl.v ≈ lungp.fb * 2813 && sl.ve ≈ lungp.ve * 2813 && sl.ps ≈ lungp.ps / 60 * 2813 * lungp.rho && sc.v ≈ heart.fb * 148.1f0
        @test_throws ArgumentError circulation(TransportSpec(LUNG_TRANSPORT; tissue_ml = 2813, coronary = nothing)) # no coronary bed
        @test_throws ArgumentError circulation(TransportSpec(LUNG_TRANSPORT; tissue_ml = 1, bed = :brain, input = :lv, coronary = hb))
        @test_throws ArgumentError circulation(TransportSpec(HEART_TRANSPORT; tissue_ml = 1, coronary = hb))
        @test_throws ArgumentError circulation(TransportSpec(HEART_TRANSPORT; tissue_ml = 1, input = :nowhere))
        # an organ in the lung slot fed back once: its outflow replaces the body's lung segment
        sL = TransportSpec(LUNG_TRANSPORT; tissue_ml = 2813, coronary = hb, tmax = 20, dt = 0.005, every = 1)
        L0 = transport(sL, Cr; passes = 0, times = 0.0f0:2.0f0:20.0f0); L1 = transport(sL, Cr; passes = 1, times = 0.0f0:2.0f0:20.0f0)
        @test L1.body.out[:lung] != L0.body.out[:lung] && all(>=(0), L1.tissue)
        # by default the feedback repeats until the organ's outflow stops changing, so the curves come from one body run
        LC = transport(sL, Cr; times = 0.0f0:2.0f0:20.0f0)
        @test L1.passes == 1 && LC.passes >= 1 && LC.change <= 1.0e-3
    end

    @testset "phantom stage: spec rules and a small bundle that loads" begin
        @test_throws ArgumentError PhantomSpec(; times = Float32[0]) # every field required
        @test_throws ArgumentError PhantomSpec(HEART_PHANTOM; colour = 1)
        @test_throws ArgumentError PhantomSpec(HEART_PHANTOM; materials = Dict(100 => "air"))
        @test_throws ArgumentError PhantomSpec(HEART_PHANTOM; pools = [Pool("x", "ncat_blood", :lv, :peak)])
        @test_throws ArgumentError PhantomSpec(HEART_PHANTOM; pools = [Pool("x", "ncat_liver", :liver, :organ)])        # no ml
        @test_throws ArgumentError PhantomSpec(HEART_PHANTOM; tissue_classes = 1000)                                    # ids 1001-1999
        @test LUNG_PHANTOM.pools == HEART_PHANTOM.pools && LUNG_PHANTOM.tissue_material == "ncat_lung"
        # a 300-terminal tree on HCT-CCO's sample, its terminals' voxels as the tissue, everything else LV outflow or air
        tree = grow(GrowthSpec(HEART_GROWTH; nterm = 300, tissue_ml = 100.0, wall_clearance = 1.0, root_caps = nothing, symmetry_floor = Tuple{Float64, Float64}[]), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        h = 2.0; tp = Float32.(vt.b[:, nch .== 0]); key(p) = Tuple(round.(Int, p ./ h))
        cells = unique(key(tp[:, i]) for i in axes(tp, 2)); vox = Dict(c => Int32(k) for (k, c) in enumerate(cells))
        q = Float32.(reduce(hcat, [collect(c) .* h for c in cells])) ./ Float32(h) # unit grid for territory_shares
        terr, pairs = C.Territory.territory_shares(tp ./ Float32(h), q; voxel_ml = h^3 / 1000)
        Cr = C.Flow.circuit(vt, terr; phys = HEART_FLOW, state = :rest)
        spec = PhantomSpec(
            HEART_PHANTOM; times = collect(0.0f0:2.0f0:20.0f0), tissue_classes = 4, vessel_classes = 3,
            pools = [Pool("LV outflow", "ncat_blood", :lv, :out), Pool("tissue mean", "ncat_heart", :coronary, :tissue)]
        )
        T = transport(TransportSpec(HEART_TRANSPORT; tissue_ml = sum(terr), tmax = 20, dt = 0.005, every = 1), Cr; times = spec.times)
        lo = minimum(vt.b; dims = 2) .- 4; hi = maximum(vt.b; dims = 2) .+ 4
        grid = C.VoxelGrid(Tuple(Float64.(vec(lo))), h, Tuple(ceil.(Int, vec(hi .- lo) ./ h) .+ 1))
        classify(p) = (
            c = key(collect(p)); haskey(vox, c) ? (VOX_TISSUE, vox[c], 0.0f0) :
                p[1] < (lo[1] + hi[1]) / 2 ? (VOX_POOL, Int32(1), Float32(0.1 * (p[3] - lo[3]) / (hi[3] - lo[3]))) : (VOX_STATIC, Int32(0), 0.0f0)
        )
        out = mktempdir()
        half = [c[1] <= cells[1][1] for c in cells] # a region of the tissue
        R = phantom(
            spec, [("rest", Cr, T)], classify, grid; pairs, terr, nvox = length(cells), out, name = "case0",
            regions = ["half" => half, "all" => fill(true, length(cells))], rois = ["outflow" => (lo[1] + 0.5, lo[2] + 2.0, lo[3] + 2.0)]
        )
        @test sum(values(R.counts)) == prod(grid.dims) && 0 in R.ids && any(i -> 100 <= i < 1000, R.ids)
        @test any(i -> 1001 <= i < 2000, R.ids) && any(i -> i > 2000, R.ids) # tissue and drawn vessels
        @test R.flow[1][2] ≈ sum(Cr.q[j] for j in eachindex(vt.parent) if vt.parent[j] == 0) * 60 / sum(terr)
        ph = C.DynamicPhantom.load(R.csv; volume = true)
        @test ph.counts == R.counts && ph.facts["states"] == "rest" && size(ph.labels) == grid.dims
        # each tissue class carries its ground-truth flow; the voxel-weighted mean over classes is the organ's flow per gram
        A = C.DynamicPhantom.readdlm(R.csv, ','; header = false, quotes = true); hdr = string.(A[1, :]); rw(n) = A[findfirst(==(n), string.(A[:, 1])), :]
        tcols = [c for c in 3:length(hdr) if all(isdigit, hdr[c]) && 1000 < parse(Int, hdr[c]) < 2000]
        mbf = Float64.(rw("mbf_rest_ml_min_g")[tcols]); nv = Float64.(rw("voxels")[tcols])
        @test all(>(0), mbf) && isapprox(sum(mbf .* nv) / sum(nv), R.flow[1][2] / HEART_FLOW.rho; rtol = 0.1)
        # regions carry the same ground truth; the whole tissue as a region is the organ's flow; the ROI is in ΔHU
        ac = findfirst(==("region:all"), hdr)
        @test isapprox(Float64(rw("mbf_rest_ml_min_g")[ac]), R.flow[1][2] / HEART_FLOW.rho; rtol = 0.02)
        @test haskey(ph.reference, ("rest", "roi:outflow")) && maximum(ph.reference[("rest", "roi:outflow")]) > 80
        st = C.DynamicPhantom.scan_times(ph; state = "rest"); @test st.arrival <= st.peak
        ps = C.DynamicPhantom.load(R.csv; zrange = 2:3) # a z slab of a large phantom
        @test ps.labels == ph.labels[:, :, 2:3] && sum(values(ps.counts)) == length(ps.labels) && ps.facts["zrange"] == "2:3"
        @test_throws ArgumentError C.DynamicPhantom.load(R.csv; zrange = 0:2)
        tc = sort([i for i in keys(ph.counts) if 1001 <= i < 2000]) # classes rank by enhancement
        @test issorted([sum(ph.curves[("rest", i)]) for i in tc])
        m = C.DynamicPhantom.materials_at(ph, "rest", 10.0)
        @test length(m.materials) <= 64 && all(l -> haskey(m.lut, l), keys(ph.counts))
    end

    @testset "balanced pairing of beds and venules" begin
        # beds and venules on two interleaved grids: nearest pairing would load some venules twice; balanced is one-to-one
        P = Float32.(hcat([[i, 0, 0] for i in 0:9]...)); Q = Float32.(hcat([[i + 0.6f0 * (i % 2), 0.1f0, 0] for i in 0:9]...))
        b = C.Territory.balanced_terminal(P, Q); n = C.Territory.nearest_terminal(P, Q)
        @test length(unique(n)) < 10                                   # nearest pairing loads some venules twice
        @test sort(b) == collect(1:10) && b == C.Territory.balanced_terminal(P, Q)
        @test all(i -> abs(P[1, b[i]] - Q[1, i]) <= 1, 1:10)            # every bed drains into a neighbouring venule
        # one bed more than venules, two beds at each even venule with the odd venules free beside them: capacity
        # max(1, 11 ÷ 10) = 1 sends each second bed to a free neighbour and leaves at most one venule with two beds (a bed
        # with no free venule within reach falls back to its nearest), where ⌈11 / 10⌉ = 2 would load the five even ones
        Q2 = Float32.(hcat([[j + 0.1f0 * k, 0.1f0, 0] for j in 0:2:8 for k in 1:2]..., [9.5f0, 0.1f0, 0])); b2 = C.Territory.balanced_terminal(P, Q2)
        @test count(j -> count(==(j), b2) == 2, 1:10) <= 1 && all(in(1:10), b2)
    end

    @testset "flow summaries after the solver's work arrays are freed (CPU, > 4096 segments)" begin
        # on the CPU a solver array can be its own host copy (par === hpar): freeing the work arrays must keep it (audit 29)
        tree = grow(GrowthSpec(HEART_GROWTH; nterm = 2500, tissue_ml = 100.0, wall_clearance = 1.0, root_caps = nothing, symmetry_floor = Tuple{Float64, Float64}[]), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        @test length(vt.r) > 4096
        Cr = C.Flow.circuit(vt, fill(0.04f0, count(==(0), nch)); phys = HEART_FLOW, state = :rest)
        @test isfinite(Cr.resist.r_micro) && all(isfinite, Cr.resist.chilian_shares)
    end

    @testset "sub-terminal vessels in series with a sheet bed" begin
        F = C.Flow; z = zeros(Float32, 3, 1)
        K = 0.02f0; h0 = 3.5f0; al = 0.17f0; pe = 0.0f0; hm = 20.0f0
        # no sub-terminal orders: the Fung sheet alone
        q, _ = F._bed(17.0f0, 11.0f0, K, h0, al, pe, hm, 1.0f0, Int32(0), Int32(0), z, z, 0.0f0)
        @test q ≈ K * ((h0 + al * 17)^4 - (h0 + al * 11)^4)
        # a sheet far wider than the vessels: the bed is the vessels' series resistance (rigid, beta = 0), weight w
        sa = Float32[2.0 3.0; 0 0; 0 0]; sv = Float32[5.0; 0; 0;;]; w = 0.5f0
        q, r = F._bed(17.0f0, 11.0f0, 1.0f6, h0, al, pe, hm, w, Int32(2), Int32(1), sa, sv, 0.0f0)
        @test r ≈ (2 + 3 + 5) / w rtol = 1.0e-3
        @test q ≈ 6 / r rtol = 1.0e-3
        # reversed pressures reverse the flow
        q2, _ = F._bed(11.0f0, 17.0f0, 1.0f6, h0, al, pe, hm, w, Int32(1), Int32(2), sv[:, [1, 1]], sa, 0.0f0)
        @test q2 < 0
        @test F.HEART_FLOW.sub_art === nothing
        # sub-terminal vessels are solved with sheet beds only (a tube bed's κ would scale them; audit A42-2)
        T1 = F.OrderTable([1.0], [0.01], [0.1], [0.0])
        @test_throws ArgumentError F.FlowSpec(F.HEART_FLOW; sub_art = T1, sub_ven = T1)
        @test F.FlowSpec(F.LUNG_FLOW; sub_art = T1, sub_ven = T1).sub_art === T1
        @test F.FlowSpec(F.LUNG_FLOW; p_ref_cmH2O = 10).p_ref_cmH2O == 10
        # the lung carries its intra-acinar vessels per order (Huang 2011: 16 arterial, 15 venous orders)
        @test length(F.LUNG_FLOW.sub_art.N) == 16 && length(F.LUNG_FLOW.sub_ven.N) == 15
    end

    @testset "organ presets cite existing reference data" begin
        # the presets are code (organs.jl); every [id] in their comments must be a measurement file in data/reference
        src = read(joinpath(pkgdir(C), "src", "perfusion", "organs.jl"), String)
        ids = unique(vcat([split(m.captures[1], r",\s*") for m in eachmatch(r"\[([a-z][a-z0-9_]*(?:,\s*[a-z][a-z0-9_]*)*)\]", src)]...))
        @test length(ids) > 20
        @test all(id -> isfile(joinpath(pkgdir(C), "data", "reference", id * ".toml")), ids)
    end

    @testset "terminal count from organ flow" begin
        # the reference heart's flow gives its validated count; any organ's count scales with its flow
        @test terminal_count(GrowthSpec(HEART_GROWTH; nterm = nothing, flow_ml_min = 135.6)) == 957_045
        @test terminal_count(GrowthSpec(HEART_GROWTH; nterm = nothing, flow_ml_min = 2 * 135.6, tissue_ml = 300.0)) == 2 * 957_045
        @test terminal_count(GrowthSpec(HEART_GROWTH; nterm = 500, tissue_ml = 1.0)) == 500
        @test_throws ArgumentError terminal_count(GrowthSpec(HEART_GROWTH; nterm = nothing, tissue_ml = 1.0))
        @test_throws ArgumentError terminal_count(GrowthSpec(HEART_GROWTH; nterm = 1, flow_ml_min = 1.0, tissue_ml = 1.0))
    end

    @testset "export round trip" begin
        E = grow_case0(100; T = Float64)
        path = tempname() * ".xml"
        export_gxl(E.m, path)
        m2 = import_gxl(Float64, path, CCOParams(nterm = 100))
        @test nsegments(m2) == nsegments(E.m)
        @test all(isapprox.(collect(Iterators.flatten(m2.pos)), collect(Iterators.flatten(E.m.pos)); atol = 1.0e-3))
    end

    include("test_kamiya.jl")
    include("test_parity.jl")
    include("test_postprocess.jl")
    include("test_supply_gaps.jl")
end
