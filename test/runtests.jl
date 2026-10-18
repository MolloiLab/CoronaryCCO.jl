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
        tree = grow(GrowthSpec(nterm = 300, tissue_ml = 100.0, wall_clearance = 1.0), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        @test count(==(0), vt.parent) >= 1 && count(==(0), nch) >= 300                  # a root, and the 300 grown terminals
        X = Float64[0 20 0 0; 0 0 20 0; 0 0 0 20]; T = [1 1 1 2; 3 2 4 3; 2 4 3 4]   # a tetrahedron, 20 mm edges
        img, _ = read_vol(mesh_domain(tempname() * ".vol", X, T))
        @test 1000 < count(>(0x80), img) < 1600                                      # 1333 mm³
    end

    @testset "transport is the whole-body bolus with the organ's bed" begin
        heart = FlowSpec(organ = "heart")
        spec = TransportSpec(; phys = heart, tissue_ml = 148.1, bed = :coronary, input = :lv, tmax = 20, dt = 0.005, every = 1)
        r = transport(spec)
        circ = C.Body.human_circulation(; hct = heart.hd, coronary = (; ml = 148.1f0, heart.fb, heart.ve, heart.ps, heart.rho))
        r2 = C.Body.simulate(circ, C.Body.Injection(); tmax = 20, dt = 0.005, every = 1)
        @test r.t == r2.t && r.out[:lv] == r2.out[:lv] && r.mean[:coronary] == r2.mean[:coronary]
        @test r.injected > 0 && maximum(r.out[:lv]) > 0
        @test circulation(spec).hct == heart.hd
        # with the organ's own circuit: passes = 0 is the plain chain (body -> plasma input -> masses); one pass feeds the
        # organ's venous outflow back into the body in place of its lumped bed
        tree = grow(GrowthSpec(nterm = 300, tissue_ml = 100.0, wall_clearance = 1.0), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        Cr = C.Flow.circuit(vt, fill(0.1f0, count(==(0), nch)); phys = heart, state = :rest)
        T0 = transport(spec, Cr; passes = 0, times = 0.0f0:2.0f0:20.0f0)
        pf = 1 - heart.hd; tf = collect(0.0f0:0.02f0:20.0f0)
        inp = C.Flow.InputCurve(Float32[C._interp(t, r.t, r.out[:lv]) / pf for t in tf], 0.02f0)
        @test T0.body.out[:lv] == r.out[:lv] && T0.input.c == inp.c
        am, vm, tm = Array.(C.Flow.frame_masses(Cr, inp, collect(0.0f0:2.0f0:20.0f0)))
        @test T0.artery == am && T0.vein == vm && T0.tissue == tm
        T1 = transport(spec, Cr; passes = 1, times = 0.0f0:2.0f0:20.0f0)
        @test T1.body.out[:coronary] != T0.body.out[:coronary] && all(>=(0), T1.tissue)
        # the lungs as the organ: the body's lung segment takes the lung's own blood, interstitium and PS; the coronary bed
        # stays the reference heart's
        lungp = FlowSpec(organ = "lung")
        hb = (; ml = 148.1f0, heart.fb, heart.ve, heart.ps, heart.rho)
        BL = circulation(TransportSpec(; phys = lungp, tissue_ml = 2813, bed = :lung, input = :rv, coronary = hb))
        sl = BL.segs[findfirst(x -> x.name === :lung, BL.segs)]; sc = BL.segs[findfirst(x -> x.name === :coronary, BL.segs)]
        @test sl.v ≈ lungp.fb * 2813 && sl.ve ≈ lungp.ve * 2813 && sl.ps ≈ lungp.ps / 60 * 2813 * lungp.rho && sc.v ≈ heart.fb * 148.1f0
        @test_throws ArgumentError circulation(TransportSpec(; phys = lungp, tissue_ml = 2813, bed = :lung, input = :rv)) # no coronary bed
        @test_throws ArgumentError circulation(TransportSpec(; phys = lungp, tissue_ml = 1, bed = :brain, input = :lv, coronary = hb))
        @test_throws ArgumentError circulation(TransportSpec(; phys = heart, tissue_ml = 1, bed = :coronary, input = :lv, coronary = hb))
        @test_throws ArgumentError circulation(TransportSpec(; phys = heart, tissue_ml = 1, bed = :coronary, input = :nowhere))
        # an organ in the lung slot fed back once: its outflow replaces the body's lung segment
        sL = TransportSpec(; phys = lungp, tissue_ml = 2813, bed = :lung, input = :rv, coronary = hb, tmax = 20, dt = 0.005, every = 1)
        L0 = transport(sL, Cr; passes = 0, times = 0.0f0:2.0f0:20.0f0); L1 = transport(sL, Cr; passes = 1, times = 0.0f0:2.0f0:20.0f0)
        @test L1.body.out[:lung] != L0.body.out[:lung] && all(>=(0), L1.tissue)
        # by default the feedback repeats until the organ's outflow stops changing, so the curves come from one body run
        LC = transport(sL, Cr; times = 0.0f0:2.0f0:20.0f0)
        @test L1.passes == 1 && LC.passes >= 1 && LC.change <= 1.0e-3
    end

    @testset "flow summaries after the solver's work arrays are freed (CPU, > 4096 segments)" begin
        # on the CPU a solver array can be its own host copy (par === hpar): freeing the work arrays must keep it (audit 29)
        tree = grow(GrowthSpec(nterm = 2500, tissue_ml = 100.0, wall_clearance = 1.0), sample("case0.vol"), [sample("case0.xml")])
        vt = vtree(tree); nch = zeros(Int, length(vt.r)); foreach(p -> p > 0 && (nch[p] += 1), vt.parent)
        @test length(vt.r) > 4096
        Cr = C.Flow.circuit(vt, fill(0.04f0, count(==(0), nch)); phys = FlowSpec(organ = "heart"), state = :rest)
        @test isfinite(Cr.resist.r_micro) && all(isfinite, Cr.resist.chilian_shares)
    end

    @testset "terminal count from organ flow" begin
        # the reference heart's flow gives its validated count; any organ's count scales with its flow
        @test terminal_count(GrowthSpec(flow_ml_min = 135.6, tissue_ml = 148.1)) == 957_045
        @test terminal_count(GrowthSpec(flow_ml_min = 2 * 135.6, tissue_ml = 300.0)) == 2 * 957_045
        @test terminal_count(GrowthSpec(nterm = 500, tissue_ml = 1.0)) == 500
        @test_throws ArgumentError terminal_count(GrowthSpec(tissue_ml = 1.0))
        @test_throws ArgumentError terminal_count(GrowthSpec(nterm = 1, flow_ml_min = 1.0, tissue_ml = 1.0))
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
