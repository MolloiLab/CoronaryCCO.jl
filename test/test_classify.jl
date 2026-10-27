# The phantom classifiers' XCAT helpers (Vasculature.XCAT: tissue_numbers, nearest_tissue_voxel, nearest_segment,
# closing, great_vessel_transit and its lookups) on small synthetic volumes with known answers.
module ClassifyTest
using Test
using Vasculature: VOX_POOL
using Vasculature.XCAT: tissue_numbers, nearest_tissue_voxel, imaged_segments, nearest_segment, segment_index, closing, great_vessel_transit,
    vessel_volume_at, vessel_tree_at, great_vessel_pool

@testset "phantom classifier helpers" begin
    m = falses(5, 5, 5); m[2, 2, 2] = true; m[4, 4, 4] = true
    n = tissue_numbers(m)
    @test n[2, 2, 2] == 1 && n[4, 4, 4] == 2 && count(!=(0), n) == 2
    @test nearest_tissue_voxel(n, (2, 2, 2)) == 1 && nearest_tissue_voxel(n, (5, 5, 5)) == 2 && nearest_tissue_voxel(n, (2, 2, 4)) == 1
    @test nearest_tissue_voxel(tissue_numbers(falses(5, 5, 5)), (3, 3, 3)) == 0

    tree = (isbase = [true, false, true], a = Float32[0 0 10; 0 0 0; 0 0 0], b = Float32[10 0 10; 0 5 10; 0 0 0], r = Float32[1, 1, 2])
    segs = imaged_segments(tree)
    @test segs[1] == [1, 3]
    j, d = nearest_segment(segs, (5.0, 3.0, 0.0)); @test j == 1 && d ≈ 2
    j, d = nearest_segment(segs, (10.0, 5.0, 0.5)); @test j == 3 && d < 0                       # inside segment 3's lumen
    # the bucket index gives the same answers near the lumens, and none beyond its reach
    ix = segment_index(segs; reach = 3.0, cell = 2.0)
    for p in ((5.0, 3.0, 0.0), (10.0, 5.0, 0.5), (2.0, -1.5, 1.0), (11.0, 9.0, -2.0))
        @test nearest_segment(ix, p)[1] == nearest_segment(segs, p)[1] && nearest_segment(ix, p)[2] ≈ nearest_segment(segs, p)[2]
    end
    @test nearest_segment(ix, (5.0, 30.0, 0.0)) == (0, Inf32) && nearest_segment(ix, (-500.0, 0.0, 0.0)) == (0, Inf32)

    # closing fills a hole narrower than the ball and keeps the outline
    box = falses(31, 31, 31); box[8:22, 8:22, 8:22] .= true; box[15:16, 15:16, 8:22] .= false
    c = closing(box, 3)
    @test all(c[15:16, 15:16, 10:20]) && all(c[box]) && count(c) == count(c[8:22, 8:22, 8:22])   # nothing outside the cube

    # a heart in a line: RV blood (20), pulmonary artery (25) running +x from it, LV blood (19), aorta (28) running +x then
    # turning down at its top
    L = zeros(UInt8, 30, 9, 30)
    L[1, 2, 5] = 20; L[2:20, 2, 5] .= 25
    L[1, 6, 20] = 19; L[2:10, 6, 20] .= 28; L[10, 6, 21:25] .= 28; L[11:15, 6, 25] .= 28; L[15, 6, 10:24] .= 28
    L[1, 8, 2] = 22; L[2:4, 8, 2] .= 24 # an RA and its vein (the IVC/SVC seeds) and a liver voxel to name the sides
    L[25, 8, 2] = 10
    T = great_vessel_transit(L, (0.0, 0.0, 0.0), 1.0)
    @test T.tree[5, 2, 5] == 0x02 && T.tree[5, 6, 20] == 0x01
    @test T.volume[3, 2, 5] < T.volume[10, 2, 5] < T.volume[20, 2, 5]                            # blood nearer the heart grows along it
    @test T.volume[10, 2, 5] ≈ 8 / 1000                                                          # the 8 cells of 1 mm³ before it
    @test T.descending[15, 6, 12] && !T.descending[5, 6, 20]
    @test vessel_volume_at(T, (9.5, 1.0, 4.0), 0x02) ≈ (T.volume[10, 2, 5] + T.volume[11, 2, 5]) / 2   # halfway between two cells
    @test vessel_tree_at(T, (4.0, 1.0, 4.0), (0x01, 0x02)) == 0x02
    k, pool, ml = great_vessel_pool(T, (9.0, 1.0, 4.0), 0x02); @test k == VOX_POOL && pool == 6 && ml ≈ T.volume[10, 2, 5]
    k, pool, ml = great_vessel_pool(T, (14.0, 5.0, 11.0), 0x01); @test pool == 27 && ml ≈ T.volume[15, 6, 12] - T.apex_ml
    @test great_vessel_pool(T, (4.0, 5.0, 19.0), 0x01)[2] == 5
end
end
