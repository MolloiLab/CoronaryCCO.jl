# Supply gaps of an imaged tree (src/xcat/supply_gaps.jl): a tree that reaches all tissue is left alone; a truncated
# tree is completed until no tissue is farther than max_gap_mm from a vessel.
module SupplyGapsTest
using Test
import AcceleratedKernels as AK
include(joinpath(@__DIR__, "..", "src", "xcat", "supply_gaps.jl"))
@testset "supply gaps" begin
    dom = trues(40, 40, 40)                                  # a 40 mm cube of tissue, 1 mm voxels
    lower = (0.5, 0.5, 0.5)
    # a "full" tree: a trunk along z with 8 branches reaching every octant
    par = [0, 1]; X = [20.0 20.0; 20.0 20.0; 0.0 20.0]; r = [3.0, 3.0]
    for (dx, dy, dz) in Iterators.product((-10, 10), (-10, 10), (-10, 10))
        push!(par, 2); X = hcat(X, [20 + dx, 20 + dy, 20 + dz]); push!(r, 1.0)
    end
    full = fill_supply_gaps(par, X, r, dom, lower; r_resolved = 1.0, max_gap_mm = 18, min_gap_ml = 0.5)    # cube corners are 17 mm from the branch tips
    @test isempty(full.added) && full.parent == par && full.X == X
    # truncated: only the first 10 mm of the trunk
    t = fill_supply_gaps([0, 1], [20.0 20.0; 20.0 20.0; 0.0 8.0], [3.0, 3.0], dom, lower; r_resolved = 1.0, max_gap_mm = 10, min_gap_ml = 0.5)
    @test !isempty(t.added)
    @test maximum(t.gaps[end]; init = 0.0) < 0.5               # the last round finds no gap of 0.5 mL or more
    @test all(t.r[a.from] >= a.radius_mm for a in t.added)    # every branch starts from a vessel at least as thick as itself
    @test all(count(==(k), t.parent) <= 2 for k in eachindex(t.parent))    # binary
    # blockers: a wall across the cube (forbidden) is never crossed, and branches stay in their target's region
    wall = falses(40, 40, 40); wall[:, 20, 12:40] .= true                  # a slab at y = 20 above z = 12
    reg = ones(Int, 40, 40, 40); reg[:, 21:40, :] .= 2
    dom2 = dom .& .!wall
    b = fill_supply_gaps([0, 1], [20.0 20.0; 10.0 10.0; 0.0 8.0], [3.0, 3.0], dom2, lower; r_resolved = 1.0, max_gap_mm = 10, min_gap_ml = 0.5, forbidden = wall, region = reg)
    crosses(a, c) = any(wall[clamp.(round.(Int, (a .+ s .* (c .- a)) .- 0.5) .+ 1, 1, 40)...] for s in 0:0.01:1)
    @test !isempty(b.added) && all(!crosses(b.X[:, b.parent[k]], b.X[:, k]) for k in eachindex(b.parent) if b.parent[k] > 0)
end
end
