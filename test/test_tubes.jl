# Tube meshes into inlet trees (Vasculature.XCAT.read_mesh_blocks, tube_centerline, tube_forest, graft, drain) and label
# volumes into growth inputs (label_mask, label_blockers): synthetic tubes of known axis and radius, joined into binary trees.
module TubesTest
using Test
using Vasculature: VoxelGrid
using Vasculature.XCAT: Centerline, BaseTree, nnodes, read_mesh_blocks, tube_centerline, tube_forest, graft, drain, label_mask, label_blockers

# the vertices of a tube of radius r from a to b: rings every 0.5 mm, 24 points around
function tube(a, b, r)
    a = collect(Float64, a); b = collect(Float64, b); ax = (b .- a) ./ sqrt(sum(abs2, b .- a))
    e = abs(ax[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]; u = e .- (e' * ax) .* ax; u ./= sqrt(sum(abs2, u))
    v = [ax[2] * u[3] - ax[3] * u[2], ax[3] * u[1] - ax[1] * u[3], ax[1] * u[2] - ax[2] * u[1]]
    L = sqrt(sum(abs2, b .- a))
    return reduce(hcat, [a .+ s .* ax .+ r .* (cos(θ) .* u .+ sin(θ) .* v) for s in 0:0.5:L for θ in range(0, 2π; length = 25)[1:24]])
end
binary(t) = all(count(==(k), t.parent) <= 2 for k in 1:nnodes(t)) && count(==(0), t.parent) == 1
parents_first(t) = all(t.parent[k] < k for k in 2:nnodes(t))

@testset "tube meshes into inlet trees" begin
    # a tube's centreline: along its axis, one node per step, its radius
    c = tube_centerline(tube((0, 0, 0), (40, 0, 0), 2.0); name = "t", step = 2.0)
    @test c.name == "t" && size(c.X, 2) == 21 && all(abs.(c.X[2:3, :]) .< 1.0e-9) && all(abs.(c.r .- 2.0) .< 0.01)
    @test minimum(c.X[1, :]) < 1.5 && maximum(c.X[1, :]) > 38.5

    # a trunk tube with a side tube near its middle and a far tube: two trees, the side joined to the trunk
    trunk = tube_centerline(tube((0, 0, 0), (60, 0, 0), 2.0); name = "a", step = 2.0)
    side = tube_centerline(tube((30, 3, 0), (30, 25, 0), 1.0); name = "a", step = 2.0)
    far = tube_centerline(tube((0, 100, 0), (30, 100, 0), 1.5); name = "a", step = 2.0)
    F = tube_forest([side, far, trunk]; toward = (0, 0, 0), min_nodes = 3)
    @test length(F) == 2 && all(binary, F) && all(parents_first, F)
    big = F[argmax(nnodes.(F))]
    # rooted at the end toward the hub; the side tube hangs from a new node splitting the trunk segment it meets
    @test nnodes(big) == size(trunk.X, 2) + size(side.X, 2) + 1 && (big.X[:, 1] ≈ trunk.X[:, 1] || big.X[:, 1] ≈ trunk.X[:, end]) && big.X[1, 1] < 1.5
    other = F[argmin(nnodes.(F))]; @test other.X[1, 1] < other.X[1, end]                 # its root: the end nearer the hub
    # the side tube hangs from the trunk node nearest its near end
    first_side = findfirst(k -> big.X[2, k] > 2.0, 1:nnodes(big)); joint = big.X[:, big.parent[first_side]]
    @test abs(joint[1] - 30) <= 2.0 && abs(joint[2]) < 1.0e-9
    @test length(tube_forest([side, far, trunk]; toward = (0, 0, 0), min_nodes = 100)) == 0   # stray trees dropped

    # grafted onto a trunk running from its valve (far from the trees) to its distal end beside them
    main = Centerline("pa", [range(-60.0, -5.0; length = 12)'; zeros(1, 12); zeros(1, 12)], [1.0; fill(8.0, 10); 1.0])
    G = graft(main, F; name = "pa")
    @test binary(G) && parents_first(G) && nnodes(G) == 10 + sum(nnodes.(F))   # the two collapsed end rings dropped
    @test G.X[1, 1] < -50 && G.r[2] == 8.0                                     # from the valve; the root edge at the median ring
    @test_throws ArgumentError graft(main, BaseTree[])

    # drained to a chamber: the largest tree takes the others, the root is the nearest surface point
    surface = [fill(-20.0, 1, 5); range(-2.0, 2.0; length = 5)'; zeros(1, 5)]
    D = drain(F, surface; root_radius = 6.0, name = "v")
    @test binary(D) && parents_first(D) && nnodes(D) == 1 + sum(nnodes.(F))
    @test D.X[:, 1] == surface[:, 3] && D.r[1] == 6.0 && D.r[2] == 6.0 && D.parent[2] == 1

    # the .raw format: a name heads each block; a name repeated is another tube; offsets move every point
    path = tempname()
    open(path, "w") do io
        println(io, "arts"); println(io, "0 0 0  1 0 0  0 1 0"); println(io, "veins"); println(io, "0 0 0  0 0 1  1 1 1")
        println(io, "arts"); println(io, "2 2 2  3 2 2  2 3 2"); println(io, "1 1 1  2 2 2  3 3 3")
    end
    B = read_mesh_blocks(path; offset = (10, 0, 0))
    @test length(B["arts"]) == 2 && length(B["veins"]) == 1 && size(B["arts"][2][1]) == (3, 6) && size(B["arts"][2][2]) == (3, 2)
    @test B["arts"][1][1][:, 2] == [11.0, 0.0, 0.0] && B["arts"][1][2] == reshape(1:3, 3, 1)
    @test collect(keys(read_mesh_blocks(path; names = ("veins",)))) == ["veins"]
end

@testset "label volumes into growth inputs" begin
    # 0.2 mm labels: label 6 fills x < 2 mm, label 7 a single voxel at the far corner
    L = zeros(UInt8, 20, 10, 10); L[1:10, :, :] .= 6; L[20, 10, 10] = 7
    g = VoxelGrid((0.4, 0.4, 0.4), 1.0, (4, 2, 2))                           # 1 mm cells, centres at 0.4, 1.4, ... mm
    M = label_mask(L, (0.0, 0.0, 0.0), 0.2, g, (6,))
    @test M[:, 1, 1] == [true, true, false, false] && all(M[1:2, :, :]) && !any(M[3:4, :, :])
    @test label_mask(L, (0.0, 0.0, 0.0), 0.2, g, (6,); fraction = 0.0) |> all  # every cell holds at least none
    B, O, cell = label_blockers(L, (0.0, 0.0, 0.0), 0.2, g, (7,))
    @test cell ≈ 0.4 && count(!=(0), B) == 1 && O[1] ≈ 0.1
    I = Tuple(findfirst(!=(0), B)); @test all(abs.(O .+ (I .- 1) .* cell .- (3.8, 1.8, 1.8)) .<= cell / 2 + 1.0e-9)   # the cell holding (3.8, 1.8, 1.8)
end
end
