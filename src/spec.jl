# SPDX-License-Identifier: GPL-3.0-or-later
# The one growth entry point: an organ-generic specification and `grow`.

"""
    GrowthSpec(; nterm, tissue_ml, wall_clearance = 0.5, gamma = 2.33, root_caps = nothing, root_fixed = false,
               symmetry_floor = Tuple{Float64, Float64}[], seed = 1, T = Float32, device = identity)

Everything that differs between organs when growing a vascular tree; everything else is the same for every organ (see
[`grow`](@ref)).

- `nterm`: terminals to grow. Each feeds an equal share of the tissue; more terminals resolve the tree further.
- `tissue_ml`: the organ's tissue volume (mL). HCT-CCO's size-dependent constants were set on adult human hearts, so they
  scale with the organ's length against that reference, `size_scale = (tissue_ml / 148.1)^(1/3)` (148.1 mL: XCAT
  vmale_50 myocardium): the heart keeps every constant, a larger organ gets longer reaches. See `CCOParams.size_scale`.
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
    nterm::Int
    tissue_ml::Float64
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

The organ's length against the reference heart, `(tissue_ml / 148.1)^(1/3)`; exactly 1 for the reference heart.
"""
size_scale(s::GrowthSpec) = Float32(cbrt(s.tissue_ml / 148.1))

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
    caps = spec.root_caps === nothing ? 1.588 : spec.root_caps
    return generate(
        domain, inlets; nterm = spec.nterm, aperf = 10000 * V / 255000, min_dist = spec.wall_clearance, root = root,
        seed = spec.seed, T = spec.T, gamma = spec.gamma, volume = :full, metric = :euclidean, sampling = :continuous,
        fix_base = true, collision = :final, max_branch_angle = 90, root_radius_max = caps, root_fixed = spec.root_fixed,
        device = spec.device, territories = territories, base_groups = base_groups, forbidden = blockers,
        symmetry_floor = spec.symmetry_floor, size_scale = size_scale(spec)
    )
end
