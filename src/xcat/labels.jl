# SPDX-License-Identifier: GPL-3.0-or-later
# Growth inputs from a voxel label volume: an organ's domain on the growth grid, and the fine blocker mask no vessel may
# cross. Any organ, any labels: the caller names the label values.

"""
    label_mask(L, origin, h, g::VoxelGrid, labels; fraction = 0.5) -> Array{Bool, 3}

Cells of the (coarser) grid `g` in which at least `fraction` of the label volume's voxels carry one of `labels`. `L` is a
label volume with voxel (1, 1, 1) centred at `origin` (mm) and spacing `h` (mm); a cell takes the voxels whose centres lie
inside it.
"""
function label_mask(L::AbstractArray{<:Integer, 3}, origin, h::Real, g::VoxelGrid, labels; fraction::Real = 0.5)
    want = falses(Int(typemax(eltype(L))) + 1)
    foreach(l -> want[l + 1] = true, labels)
    M = zeros(Bool, g.dims)
    AK.foreachindex(1:g.dims[3]; min_elems = 1) do k
        for j in 1:g.dims[2], i in 1:g.dims[1]
            c = center(g, i, j, k)
            r = ntuple(d -> (ceil(Int, (c[d] - g.h / 2 - origin[d]) / h - 1.0e-9) + 1):(ceil(Int, (c[d] + g.h / 2 - origin[d]) / h - 1.0e-9)), 3)
            n = 0; m = 0
            for kk in r[3], jj in r[2], ii in r[1]
                checkbounds(Bool, L, ii, jj, kk) || continue
                n += 1; m += want[L[ii, jj, kk] + 1]
            end
            n > 0 && m >= fraction * n && (M[i, j, k] = true)
        end
    end
    return M
end

"""
    label_blockers(L, origin, h, g::VoxelGrid, labels; factor = 2) -> (mask, origin_mm, cell_mm)

The blocker mask over the box of grid `g`, in the form [`grow`](@ref) and [`penetrating_stubs`](@ref) take: cells of
`factor` × `factor` × `factor` label voxels (0.4 mm for XCAT's 0.2 mm), a cell blocked (1) when any of its voxels carries
one of `labels` (chambers, great vessels, airways). `L`, `origin` and `h` are as for [`label_mask`](@ref).
"""
function label_blockers(L::AbstractArray{<:Integer, 3}, origin, h::Real, g::VoxelGrid, labels; factor::Integer = 2)
    want = falses(Int(typemax(eltype(L))) + 1)
    foreach(l -> want[l + 1] = true, labels)
    first_voxel = [round(Int, (g.origin[d] - g.h / 2 - origin[d]) / h) + 1 for d in 1:3]
    cell = factor * h
    cells = [ceil(Int, (g.dims[d] * g.h + g.h) / cell) for d in 1:3]
    B = zeros(UInt8, cells...)
    AK.foreachindex(1:cells[3]; min_elems = 1) do c
        for b in 1:cells[2], a in 1:cells[1], dk in 0:(factor - 1), dj in 0:(factor - 1), di in 0:(factor - 1)
            I = (first_voxel[1] + factor * (a - 1) + di, first_voxel[2] + factor * (b - 1) + dj, first_voxel[3] + factor * (c - 1) + dk)
            checkbounds(Bool, L, I...) || continue
            want[L[I...] + 1] && (B[a, b, c] = 0x01)
        end
    end
    O = ntuple(d -> origin[d] + h * (first_voxel[d] - 1) + cell / 2 - h / 2, 3)
    return (B, O, cell)
end
