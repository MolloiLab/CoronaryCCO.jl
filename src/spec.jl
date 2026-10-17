# SPDX-License-Identifier: GPL-3.0-or-later
# The one growth entry point: an organ-generic specification and `grow`.

"""
    GrowthSpec(; nterm, tissue_ml, wall_clearance = 0.5, gamma = 2.33, root_caps = nothing, root_fixed = false,
               symmetry_floor = Tuple{Float64, Float64}[], seed = 1, T = Float32, device = identity)

Everything that differs between organs when growing a vascular tree; everything else is the same for every organ (see
[`grow`](@ref)).

- `nterm` or `flow_ml_min` (give one): the number of terminals, or the organ's blood flow (mL/min), from which the count
  follows as `flow_ml_min / terminal_flow_ml_min`: every organ is then resolved down to terminals of the same size. For an
  organ that perfuses its own tissue the flow is its O₂ demand by Fick, VO₂ / (arterial O₂ content × O₂ extraction); for
  the lungs it is the cardiac output, since the pulmonary arteries carry all of it to pick up O₂.
- `terminal_flow_ml_min`: the flow of one terminal. Default: the reference heart's resting flow over its terminal count,
  135.6 mL/min / 957,045 (≈ 15.7 µm terminals), the resolution at which the grown coronary tree matches Kassab's
  morphometry.
- `tissue_ml`: the organ's tissue volume (mL).
- `reference_ml`: the tissue volume HCT-CCO's size-dependent constants are used unchanged at (its longest connection,
  territory reach, convergence tolerance, centreline step). They were set for adult human hearts; the default, 148.1 mL,
  is the adult myocardium this package is validated on (XCAT vmale_50). Every other organ gets
  `size_scale = (tissue_ml / reference_ml)^(1/3)`: a larger organ gets proportionally longer reaches. See
  `CCOParams.size_scale`.
- `wall_clearance`: smallest distance (mm) from a terminal to the organ's boundary.
- `gamma`: bifurcation exponent r₀^γ = r₁^γ + r₂^γ; 7/3 from coronary structure-function scaling (Huo & Kassab 2012).
- `root_caps`: root radius per inlet tree (mm), the imaged inlet radii; `nothing` leaves roots to the pressure law.
- `root_fixed`: `true` sets the forest's radius scale so the binding root sits exactly at its cap after growth.
- `symmetry_floor`: measured minimum daughter symmetry by mother diameter, `(µm, ratio)` knots; empty for none.
- `seed`, `T`, `device`: random stream, number type, array device (`MtlArray`, `CuArray`, ... or `identity` for CPU).

Recommended values (XCAT vmale_50; adapt for other subjects):

| | heart arteries | lung arteries | lung veins |
|---|---|---|---|
| `nterm` | 957_045 | 200_000 | 156_000 |
| `tissue_ml` | 148.1 (myocardium) | 2813 (both lungs) | 2813 |
| `wall_clearance` | 0.5 | 0.5 | 0.5 |
| `root_caps` | imaged LCA and RCA ostial radii | imaged main PA radius | imaged vein ostia |
| `root_fixed` | `false` | `true` | `true` |
| `symmetry_floor` | Schwarz et al. 2020 human coronary knots | none | none |
"""
Base.@kwdef struct GrowthSpec
    nterm::Union{Nothing, Int} = nothing
    flow_ml_min::Union{Nothing, Float64} = nothing
    terminal_flow_ml_min::Float64 = 135.6 / 957_045
    tissue_ml::Float64
    reference_ml::Float64 = 148.1
    wall_clearance::Float64 = 0.5
    gamma::Float64 = 2.33
    root_caps::Union{Nothing, Vector{Float64}} = nothing
    root_fixed::Bool = false
    symmetry_floor::Vector{Tuple{Float64, Float64}} = Tuple{Float64, Float64}[]
    seed::Int = 1
    T::Type = Float32
    device::Any = identity
end

"""
    size_scale(spec::GrowthSpec) -> Float32

The organ's length against the reference, `(tissue_ml / reference_ml)^(1/3)`; exactly 1 for the reference heart.
"""
size_scale(s::GrowthSpec) = Float32(cbrt(s.tissue_ml / s.reference_ml))

"""
    terminal_count(spec::GrowthSpec) -> Int

`spec.nterm` when given; otherwise the organ's flow over the flow of one terminal, `flow_ml_min / terminal_flow_ml_min`.
"""
function terminal_count(s::GrowthSpec)
    (s.nterm === nothing) == (s.flow_ml_min === nothing) && throw(ArgumentError("GrowthSpec: give exactly one of nterm and flow_ml_min"))
    return s.nterm !== nothing ? s.nterm : round(Int, s.flow_ml_min / s.terminal_flow_ml_min)
end

"""
    grow(spec::GrowthSpec, domain, inlets; root = nothing, blockers = nothing, territories = nothing,
         base_groups = nothing) -> Engine

Grow a vascular forest from the imaged `inlets` (GXL files, one tree each: two coronary roots, a main pulmonary artery,
several cerebral inlets, ...) into the tissue mask `domain` (a `.vol` file, 1 mm voxels). The algorithm is the same for
every organ: HCT-CCO's staged growth (new vessels attach where its stage schedule allows), sizes tested at their final
radius, no daughter turning back past 90°, terminals sampled continuously, imaged centrelines kept in place, Euclidean
wall distance, full tree volume, batched on the chosen device. Only `spec` differs between organs.

- `root`: a voxel inside the domain to start the domain's distance field from (optional).
- `blockers`: `(mask, origin_mm, cell_mm)`, a fine mask no new vessel may cross (heart chambers and great vessels, lung
  airways, brain ventricles, ...).
- `territories`, `base_groups`: optional per-voxel tree-group bits and the group of each imaged segment, so each imaged
  vessel keeps the tissue nearest to it (see [`set_territories!`](@ref)).
"""
function grow(spec::GrowthSpec, domain, inlets; root = nothing, blockers = nothing, territories = nothing, base_groups = nothing)
    img, _ = read_vol(domain); V = count(>(MASK_THRESHOLD), img)
    # perfusion area per voxel as upstream's sample run sets it (aperf 10000 over its 255,000-voxel case0 domain): HCT-CCO's
    # aperf is an area in voxel units, so it scales with the domain's voxel count
    caps = spec.root_caps === nothing ? 1.588 : spec.root_caps
    return generate(
        domain, inlets; nterm = terminal_count(spec), aperf = 10000 * V / 255000, min_dist = spec.wall_clearance, root = root,
        seed = spec.seed, T = spec.T, gamma = spec.gamma, volume = :full, metric = :euclidean, sampling = :continuous,
        fix_base = true, collision = :final, max_branch_angle = 90, root_radius_max = caps, root_fixed = spec.root_fixed,
        device = spec.device, territories = territories, base_groups = base_groups, forbidden = blockers,
        symmetry_floor = spec.symmetry_floor, size_scale = size_scale(spec)
    )
end

"""
    write_vol(path, mask::AbstractArray{Bool, 3}, lower::NTuple{3, <:Integer}) -> path

Write a Bool mask as a DGtal `.vol` (Version 2, 1 mm voxels): voxel `mask[i, j, k]` is the point `lower .+ (i, j, k) .- 1`.
"""
function write_vol(path::AbstractString, mask::AbstractArray{Bool, 3}, lower::NTuple{3, <:Integer})
    s = size(mask); c = lower .+ (s .- 1) .÷ 2
    open(path, "w") do io
        print(io, "Center-X: $(c[1])\nCenter-Y: $(c[2])\nCenter-Z: $(c[3])\nX: $(s[1])\nY: $(s[2])\nZ: $(s[3])\n")
        print(io, "Voxel-Size: 1\nAlpha-Color: 0\nVoxel-Endian: 0\nInt-Endian: 0123\nVersion: 2\n.\n")
        write(io, UInt8[m ? 0xff : 0x00 for m in mask])
    end
    return path
end

"""
    mesh_domain(path, vertices, triangles; margin = 3) -> path

The growth domain of an organ given as a closed surface mesh (`vertices` 3 × n in mm, `triangles` 3 × m, 1-based):
the voxels of a 1 mm grid whose centres lie inside the mesh, written to `path` as a `.vol` for [`grow`](@ref).
"""
function mesh_domain(path::AbstractString, X::AbstractMatrix, T::AbstractMatrix{<:Integer}; margin::Real = 3)
    g0 = grid_around(X, 1.0; margin)
    lo = round.(Int, g0.origin); g = VoxelGrid(Float64.(lo), 1.0, g0.dims)
    return write_vol(path, voxelize(X, T, g), lo)
end
