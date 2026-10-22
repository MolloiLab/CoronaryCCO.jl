# Penetrating stubs of an imaged tree (Vasculature.XCAT.penetrating_stubs): sites about one spacing apart however finely
# the vessel is sampled, none where the tissue is out of reach or the path is blocked, an imaged side branch counts as a
# stub that is already there, and the tree stays binary.
module StubsTest
using Test
using Vasculature: VoxelGrid
using Vasculature.XCAT: BaseTree, penetrating_stubs, nnodes
@testset "penetrating stubs" begin
    g = VoxelGrid((0.0, 0.0, 0.0), 1.0, (60, 30, 20))              # 1 mm voxels, centres at 0..59, 0..29, 0..19 mm
    wall = falses(g.dims); wall[:, 9:20, :] .= true                # tissue at y = 8..19 mm
    # a straight vessel along x at y = 5, z = 10, sampled every `dx` mm over 50 mm
    vessel(dx; y = 5.0) = (n = round(Int, 50 / dx) + 1; BaseTree("v", collect(0:(n - 1)), [range(5.0, 55.0; length = n)'; fill(y, 1, n); fill(10.0, 1, n)], fill(1.0, n), fill("v", n)))
    stubsites(s, t) = [s.parent[k] for k in (nnodes(t) + 1):nnodes(s)]
    for dx in (0.25, 0.5, 1.0)
        t = vessel(dx); s, rep = penetrating_stubs(t, g, wall)
        x = [t.X[1, k] for k in stubsites(s, t)]
        @test rep.stubs == rep.sites && rep.no_tissue == rep.outside == rep.blocked == 0
        @test all(2.8 - dx / 2 - 1.0e-9 .<= diff(x) .<= 2.8 + dx / 2 + 1.0e-9)    # one spacing apart, to the sampling
        @test all(wall[(round.(Int, s.X[:, k]) .+ 1)...] for k in (nnodes(t) + 1):nnodes(s))   # tips inside the tissue
        @test s.parent[1:nnodes(t)] == t.parent && s.X[:, 1:nnodes(t)] == t.X       # the imaged tree is kept as it was
    end
    # nodes in any order (children before parents, root last): the same sites
    t = vessel(0.5); n = nnodes(t); perm = [n:-1:2; 1]; inv = invperm(perm)       # new node i is old node perm[i]
    tp = BaseTree("v", [t.parent[o] == 0 ? 0 : inv[t.parent[o]] for o in perm], t.X[:, perm], t.r[perm], t.vessel[perm])
    s1, _ = penetrating_stubs(t, g, wall); s2, _ = penetrating_stubs(tp, g, wall)
    @test sort([t.X[1, k] for k in stubsites(s1, t)]) == sort([tp.X[1, k] for k in stubsites(s2, tp)])
    # an organ twice the heart's size: every length doubles (spacing 5.6 mm)
    t = vessel(0.5); s, _ = penetrating_stubs(t, g, wall; scale = 2)
    @test all(5.6 - 0.25 - 1.0e-9 .<= diff([t.X[1, k] for k in stubsites(s, t)]) .<= 5.6 + 0.25 + 1.0e-9)
    # sampled coarser than the spacing: every node is a site
    t = vessel(50 / 18); s, rep = penetrating_stubs(t, g, wall)
    @test rep.sites == nnodes(t) - 1 && rep.stubs == rep.sites
    # out of reach: the tissue is 7 mm away
    s, rep = penetrating_stubs(vessel(0.5; y = 1.0), g, wall; spacing = 2.8, reach = 6.0)
    @test rep.stubs == 0 && rep.no_tissue == rep.sites > 0 && nnodes(s) == nnodes(vessel(0.5))
    # an imaged side branch at x = 30 counts as a stub: the next site is one spacing past it
    t = vessel(0.5); k30 = findfirst(k -> t.X[1, k] == 30.0, 1:nnodes(t))
    tb = BaseTree("v", [t.parent; k30; nnodes(t) + 1], hcat(t.X, [30.0, 12.0, 10.0], [30.0, 18.0, 10.0]), [t.r; 0.5; 0.5], [t.vessel; "b"; "b"])
    s, _ = penetrating_stubs(tb, g, wall; spacing = 2.8, reach = 6.0)
    x = sort([tb.X[1, k] for k in stubsites(s, tb) if tb.vessel[k] == "v"])
    @test 30.0 ∉ x && minimum(filter(>(30.0), x)) ≈ 33.0                 # a node with two children takes none
    @test all(count(==(k), s.parent) <= 2 for k in 1:nnodes(s))         # binary
    # blockers: a fine slab at y = 6.5..7.5 between vessel and tissue
    B = zeros(UInt8, 240, 120, 80); B[:, 27:30, :] .= 1               # 0.25 mm cells, cell (1, 1, 1) centred at 0 mm
    s, rep = penetrating_stubs(vessel(0.5), g, wall; spacing = 2.8, reach = 6.0, blockers = (B, (0.0, 0.0, 0.0), 0.25))
    @test rep.stubs == 0 && rep.blocked == rep.sites > 0
    # forbid: the tips' voxels (1.5 mm past y = 8) are forbidden
    fb = falses(g.dims); fb[:, 10:11, :] .= true
    _, rep = penetrating_stubs(vessel(0.5), g, wall; spacing = 2.8, reach = 6.0, forbid = fb)
    @test rep.stubs == 0 && rep.outside == rep.sites
    # two tissues within reach on either side: successive stubs alternate
    below = falses(g.dims); below[:, 1:3, :] .= true                  # y = 0..2 mm
    t = vessel(0.5); s, rep = penetrating_stubs(t, g, [wall, below]; spacing = 2.8, reach = 6.0)
    side = [s.X[2, k] > 5.0 for k in (nnodes(t) + 1):nnodes(s)]
    @test rep.stubs == rep.sites && all(side[1:2:end]) && !any(side[2:2:end])
    # a vessel running through the tissue: stubs leave sideways (perpendicular, r + depth long), in turning directions
    solid = trues(g.dims)
    t = vessel(0.5; y = 15.0); s, rep = penetrating_stubs(t, g, solid; spacing = 2.8, reach = 6.0)
    st = (nnodes(t) + 1):nnodes(s); D = [s.X[:, k] .- s.X[:, s.parent[k]] for k in st]
    @test rep.stubs == rep.sites == rep.sideways > 0
    @test all(abs(d[1]) < 1.0e-9 && isapprox(sqrt(sum(abs2, d)), 1.0 + 1.5) for d in D)
    @test length(unique(round.(atan.(getindex.(D, 3), getindex.(D, 2)); digits = 3))) == length(D)
    @test_throws DimensionMismatch penetrating_stubs(t, g, falses(3, 3, 3); spacing = 2.8, reach = 6.0)
end
end
