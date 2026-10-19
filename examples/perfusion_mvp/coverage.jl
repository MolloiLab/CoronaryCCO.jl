# SPDX-License-Identifier: GPL-3.0-or-later
# Terminal coverage of the myocardium: gaps (distance from each voxel to its nearest terminal)
# and evenness of terminal density.   julia --project=. coverage.jl <tree.bin> <xcat grow dir>
using Printf
import Vasculature: Flow, Territory; using .Flow, .Territory
tree = load_vtree(ARGS[1]); W = ARGS[2]
gl = split.(readlines(joinpath(W, "grid.txt"))); lower = parse.(Int, gl[1]); gd = Tuple(parse.(Int, gl[2]))
lab = Array{UInt8}(undef, gd); read!(joinpath(W, "heart_labels_1mm.u8"), lab)
vidx = findall((lab .== 0x01) .| (lab .== 0x02)); q = Matrix{Float32}(undef, 3, length(vidx))
for (v, I) in enumerate(vidx), d in 1:3
    q[d, v] = lower[d] + I[d] - 1
end
nch = zeros(Int, length(tree.r)); foreach(p -> p > 0 && (nch[p] += 1), tree.parent); terms = findall(==(0), nch)
own = nearest_terminal(tree.b[:, terms], q)
gapd = [sqrt(sum(abs2, q[:, v] .- tree.b[:, terms[own[v]]])) for v in 1:size(q, 2)]; sg = sort(gapd)
dens = Dict{NTuple{3, Int}, Int}(); for j in terms
    k = Tuple(floor.(Int, tree.b[:, j] ./ 4)); dens[k] = get(dens, k, 0) + 1
end
vcnt = Dict{NTuple{3, Int}, Int}(); for v in 1:size(q, 2)
    k = Tuple(floor.(Int, q[:, v] ./ 4)); vcnt[k] = get(vcnt, k, 0) + 1
end
dv = [get(dens, k, 0) / c for (k, c) in vcnt if c >= 48]; m = sum(dv) / length(dv)
@printf(
    "%s: %d terminals; nearest-terminal distance p99 %.2f mm, max %.2f mm; voxels > 1.5 mm: %.3f%%; empty 4 mm cubes %d of %d; density CV %.2f\n",
    basename(ARGS[1]), length(terms), sg[ceil(Int, 0.99end)], sg[end], 100count(>(1.5), gapd) / length(gapd), count(==(0), [get(dens, k, 0) for (k, c) in vcnt if c >= 48]), length(dv), sqrt(sum((dv .- m) .^ 2) / length(dv)) / m
)
