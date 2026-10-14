# SPDX-License-Identifier: GPL-3.0-or-later
using CoronaryCCO
using Test

@testset "CoronaryCCO" begin
    @testset "defaults match upstream HctCco CLI" begin
        p = CCOParams()
        @test p.nterm == 500
        @test p.aperf == 20000.0
        @test p.gamma == 3.0
        @test p.root_radius_max == 1.588
        @test p.min_distance_to_border == 5.0
        @test p.num_neighbors == 20
    end

    @testset "constructor mirrors HctccoModel(aPerf, nTerm, ...)" begin
        p = CCOParams(nterm = 40_000, aperf = 10_000.0)
        t = CCOTree(p, (-26.0, 41.0, 13.0))
        @test nsegments(t) == 2
        @test t.parent == Int32[1, 1]
        @test t.terminals == Int32[2]
        @test t.pterm == 7.98e3                      # nTerm > 250 branch
        @test t.qterm == 8.33e3 / 40_000
        @test t.rperf ≈ (3 * 10_000 / (4π))^(1 / 3)
        @test t.y[1] == 41.0 + t.rperf               # root lies rPerf above the center
    end
end
