# SPDX-License-Identifier: GPL-3.0-or-later
# The one growth entry point: an organ-generic specification and `grow`.

"""
    GrowthSpec(; nterm, flow_ml_min, terminal_flow_ml_min, tissue_ml, reference_ml, wall_clearance, gamma, root_caps,
               root_fixed, split, symmetry_floor, seed, T, device)
    GrowthSpec(spec::GrowthSpec; kwargs...)

Everything that differs between organs when growing a vascular tree; everything else is the same for every organ (see
[`grow`](@ref)). Every field is a keyword and every one is required: an organ is an instance (`HEART_GROWTH`,
`HEART_VEIN_GROWTH`, `LUNG_GROWTH`, `LUNG_VEIN_GROWTH`, organs.jl), and a run varies one with the second form, e.g.
`GrowthSpec(LUNG_GROWTH; seed = 2, device = dev)` (`dev` e.g. `MtlArray`).

- `nterm` or `flow_ml_min` (one of them, the other `nothing`): the number of terminals, or the organ's blood flow
  (mL/min), from which the count follows as `flow_ml_min / terminal_flow_ml_min`.
- `terminal_flow_ml_min`: the flow of one terminal, used with `flow_ml_min`.
- `tissue_ml`: the organ's tissue volume (mL).
- `reference_ml`: the tissue volume HCT-CCO's size-dependent constants are used unchanged at (its longest connection,
  territory reach, convergence tolerance, centreline step); they were set for adult human hearts (148.1 mL, the XCAT
  vmale_50 myocardium). Every organ gets `size_scale = (tissue_ml / reference_ml)^(1/3)`. See `CCOParams.size_scale`.
- `wall_clearance`: smallest distance (mm) from a terminal to the organ's boundary.
- `gamma`: bifurcation exponent r₀^γ = r₁^γ + r₂^γ.
- `root_caps`: root radius per inlet tree (mm), the imaged inlet radii; `nothing` leaves roots to the pressure law.
- `root_fixed`: `true` sets the forest's radius scale so the binding root sits exactly at its cap after growth.
- `split`: how a forest's inlet trees share the terminals. `:free`: by the growth's own competition. `:bifurcation`: the
  inlets are the daughters of one origin (the aortic root for the two coronaries), so they split the flow by the same
  bifurcation law as every other branch, share ∝ `root_caps[k]^gamma`: a tree takes terminals only while its share of the
  forest's flow is at or below that, so the trees grow in step and each claims the tissue it reaches. A single tree, or no
  `root_caps`, has nothing to split.
- `symmetry_floor`: measured minimum daughter symmetry by mother diameter, `(µm, ratio)` knots; empty for none.
- `seed`, `T`, `device`: random stream, number type, array device (`MtlArray`, `CuArray`, ... or `identity` for CPU).
"""
struct GrowthSpec
    nterm::Union{Nothing, Int}
    flow_ml_min::Union{Nothing, Float64}
    terminal_flow_ml_min::Float64
    tissue_ml::Float64
    reference_ml::Float64
    wall_clearance::Float64
    gamma::Float64
    root_caps::Union{Nothing, Vector{Float64}}
    root_fixed::Bool
    split::Symbol
    symmetry_floor::Vector{Tuple{Float64, Float64}}
    seed::Int
    T::Type
    device::Any
end
function GrowthSpec(; kwargs...)
    bad = setdiff(keys(kwargs), fieldnames(GrowthSpec))
    isempty(bad) || throw(ArgumentError("unknown GrowthSpec fields $(collect(bad))"))
    miss = setdiff(fieldnames(GrowthSpec), keys(kwargs))
    isempty(miss) || throw(ArgumentError("GrowthSpec: no value for $(collect(miss))"))
    s = GrowthSpec((kwargs[f] for f in fieldnames(GrowthSpec))...)
    s.split in (:free, :bifurcation) || throw(ArgumentError("GrowthSpec: split must be :free or :bifurcation"))
    return s
end
GrowthSpec(s::GrowthSpec; kwargs...) = GrowthSpec(; merge((; (f => getfield(s, f) for f in fieldnames(GrowthSpec))...), values(kwargs))...)

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
    shares = spec.split === :bifurcation && spec.root_caps !== nothing ? (w = spec.root_caps .^ spec.gamma; w ./ sum(w)) : nothing
    return generate(
        domain, inlets; nterm = terminal_count(spec), aperf = 10000 * V / 255000, min_dist = spec.wall_clearance, root = root,
        seed = spec.seed, T = spec.T, gamma = spec.gamma, volume = :full, metric = :euclidean, sampling = :continuous,
        fix_base = true, collision = :final, max_branch_angle = 90, root_radius_max = caps, root_fixed = spec.root_fixed,
        root_shares = shares,
        device = spec.device, territories = territories, base_groups = base_groups, forbidden = blockers,
        symmetry_floor = spec.symmetry_floor, size_scale = size_scale(spec)
    )
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
