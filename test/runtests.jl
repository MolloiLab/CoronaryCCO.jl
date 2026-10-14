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

grow_case0(n; kw...) = generate(
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
end
