### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9cff4f49-4723-43e1-9290-be0f428e0352
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    # a writable font cache for CairoMakie (a notebook worker may have no HOME)
    ENV["XDG_CACHE_HOME"] = get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))
end;

# ╔═╡ 6602ea59-d68c-4c71-8be6-e7c15ddab79e
using PlutoUI: Slider, Select, TableOfContents

# ╔═╡ a4268901-4cf8-4181-a7ff-52aa62e27d6c
using Printf: @sprintf

# ╔═╡ 92dde8bc-c022-4a93-bf00-b6e6a3342fea
using Serialization: serialize, deserialize

# ╔═╡ 0e0a5727-5a5d-4d7f-b27f-018760e1ba8b
using Mmap: mmap

# ╔═╡ fbd087b6-3b8c-4d0a-804c-8b0982c90812
using Statistics: median

# ╔═╡ bfa4faa8-3c77-4464-ae5d-87856a7f961b
md"""
# Cardiac Pipeline

**From the raw XCAT heart to a perfused CT phantom, and the same heart with a 90% LAD stenosis.**

| Stage | What happens |
|---|---|
| **1. Load** | XCAT's 0.2 mm label volume and NURBS surfaces |
| **2. Prepare** | the imaged coronary trees and the myocardium to grow into |
| **3. Grow** | ~1 million arterial terminals and a venous tree |
| **4. Flow** | blood flow from pressures alone, at rest and adenosine stress |
| **5. Transport** | a clinical contrast injection through the whole body |
| **6. Compartments** | the regions ground truth is given for, and their flow |
| **7. Voxelize and Save** | the 0.2 mm perfused phantom |
| **8. Inspect** | the saved phantom: labels, compartments, idealized images, ground truth |
| **9. Stenosis** | a 90% proximal LAD lesion with plaque, as its own phantom |
| **10. Next Steps** | scanning the phantom with BasisSimulator.jl |

Every heavy step is cached in `CACHE_DIR`, so a rerun takes seconds. The XCAT phantom is licensed: its files and anything made from them stay on lab storage.
"""

# ╔═╡ 45b99d76-06b3-436a-9149-0bf0fbaa3a45
md"""
## Setup
"""

# ╔═╡ 99ec1d1a-81b2-4f1b-b9d9-602949ccee63
import Vasculature as V

# ╔═╡ 06f98b38-fbb4-48c9-83c3-a327c6d74052
# ╠═╡ show_logs = false
import CairoMakie as Mke

# ╔═╡ 29a716e1-b776-4949-b130-452855d6708a
import GPUSelect

# ╔═╡ a120dd81-758c-4729-af09-19ff607347f2
TableOfContents()

# ╔═╡ db87cf42-7622-4efc-91e5-3f45da374c65
begin
    # licensed: these files and anything made from them stay on lab storage
    const XCAT_DIR = get(
        ENV,
        "XCAT_DIR",
        joinpath(
            homedir(),
            "Documents",
            "dev",
            "MolloiLab",
            "xcat-data",
            "vmale_50_heart_1600_1400_867",
        ),
    )
    # the full chest, arms included, at the same 0.2 mm: XCAT's heart block inside it and its surface meshes voxelized
    # around it (built once by xcat-data/work/phantom/build_chest.jl), so a scanner sees the whole body
    const XCAT_CHEST_DIR = get(
        ENV,
        "XCAT_CHEST_DIR",
        joinpath(
            homedir(),
            "Documents",
            "dev",
            "MolloiLab",
            "xcat-data",
            "work",
            "phantom",
        ),
    )
    const CACHE_DIR = joinpath(homedir(), "xcat-perfused", "cardiac_pipeline_cache")
    const PHANTOM_DIR = joinpath(homedir(), "xcat-perfused", "cardiac_pipeline")

    mkpath(CACHE_DIR)
    mkpath(PHANTOM_DIR)
end;

# ╔═╡ a1035ac8-ddb9-42a2-b0c2-e86eec8f1769
# ╠═╡ show_logs = false
"Arrays to this machine's GPU (CUDA, Metal, ROCm, oneAPI), or the CPU if there is none."
const DEVICE = GPUSelect.Storage()

# ╔═╡ 72d7878f-c03f-4b24-b5fa-0ef529dbcf68
"Run `compute` once and keep its result in `CACHE_DIR/name.jls`; afterwards load it."
function cached(compute, name)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) && return deserialize(path)

    seconds = @elapsed result = compute()
    serialize(path, result)
    @info "$name: $(round(seconds; digits = 1)) s"
    return result
end

# ╔═╡ ed124ecb-5d54-4641-b775-0f976fe0ec74
md"""
## 1. Load the Raw XCAT Chest

XCAT vmale\_50's whole chest, arms included: one `UInt8` label per 0.2 mm voxel (XCAT's heart block inside it), and the NURBS surfaces of every heart structure, the coronary arteries and veins among them. The phantom is saved on this full chest so a CT scan sees the whole body.
"""

# ╔═╡ 7066ffd5-cdf7-4722-968f-d48ffbc2fd0e
begin
    # the chest grid: origin (mm, the centre of voxel (1, 1, 1)), voxels, voxel size (mm)
    const CHEST_GRID = readlines(joinpath(XCAT_CHEST_DIR, "chest_grid.txt"))
    const XCAT_ORIGIN = Tuple(parse.(Float64, split(CHEST_GRID[1])))
    const XCAT_DIMS = Tuple(parse.(Int, split(CHEST_GRID[2])))
    const XCAT_VOXEL_MM = parse(Float64, CHEST_GRID[3])

    xcat_labels = mmap(
        open(joinpath(XCAT_CHEST_DIR, "chest_labels_0p2.u8")),
        Array{UInt8, 3},
        XCAT_DIMS,
    )
    xcat_surfaces = V.XCAT.read_nrb(joinpath(XCAT_DIR, "vmale_50_heart.nrb"))
end;

# ╔═╡ b5b3a68d-bde2-42e2-8df7-d09d41283847
begin
    "XCAT label groups for pictures: (name, colour, XCAT labels)."
    const LABEL_GROUPS = [
        ("myocardium", :darkorange, (15, 16)),
        ("chamber blood", :mediumpurple, (19, 20, 21, 22)),
        ("great arteries", :firebrick, (25, 28)),
        ("great veins", :royalblue, (24,)),
        ("coronary arteries", :red, (26,)),
        ("coronary veins", :deepskyblue, (27,)),
        ("lung", :gray25, (6,)),
        ("fat", :khaki, (1,)),
        ("muscle", :rosybrown, (4,)),
        ("bone", :ivory, (5, 8, 9)),
        ("other", :gray60, ()),
    ]

    "The picture group of every XCAT label (0 = air)."
    const GROUP_OF_LABEL = let group = fill(UInt8(length(LABEL_GROUPS)), 256)
        group[1] = 0
        for (k, (_, _, labels)) in enumerate(LABEL_GROUPS), label in labels
            group[label + 1] = k
        end
        group
    end

    label_group_colormap() = Mke.cgrad(
        [:black; [colour for (_, colour, _) in LABEL_GROUPS]],
        length(LABEL_GROUPS) + 1;
        categorical = true,
    )
    xcat_slice(z) = clamp(
        round(Int, (z - XCAT_ORIGIN[3]) / XCAT_VOXEL_MM) + 1,
        1,
        XCAT_DIMS[3],
    )
    xcat_x_mm = XCAT_ORIGIN[1] .+ (0:(XCAT_DIMS[1] - 1)) .* XCAT_VOXEL_MM
    xcat_y_mm = XCAT_ORIGIN[2] .+ (0:(XCAT_DIMS[2] - 1)) .* XCAT_VOXEL_MM
end;

# ╔═╡ 7e5c835a-c16c-4c0f-81ae-5410e13fcad8
md"""
z (mm): $(@bind z_xcat Slider(round.(XCAT_ORIGIN[3] .+ (0:5:(XCAT_DIMS[3] - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 255.0, show_value = true))
"""

# ╔═╡ 460614b2-2f7d-4af4-94ad-8051d7c8e70b
begin
    # shades for vessels by radius, pale for the thinnest shown and deep for the largest:
    # arteries red, veins blue
    ARTERY_COLOURS = Mke.cgrad(["#f4a582", "#d6604d", "#b2182b", "#67001f"])
    VEIN_COLOURS = Mke.cgrad(["#92c5de", "#4393c3", "#2166ac", "#053061"])
end;

# ╔═╡ 24dbb378-3f1d-4232-b123-ecb319c51c48
"""
A 3-D axis for the heart in true proportions, seen from the front: head up, the patient's left on the right. A gentle
perspective, no box or grid.
"""
function vessel_axis(position; title)
    ax = Mke.Axis3(
        position;
        title = title,
        aspect = :data,
        azimuth = 1.5π,
        elevation = 0.08π,
        perspectiveness = 0.3,
        viewmode = :fit,
        protrusions = 0,
    )
    Mke.hidedecorations!(ax)
    Mke.hidespines!(ax)
    return ax
end

# ╔═╡ 94428579-cfcc-4c6d-a766-a1a3d36fede9
let k = xcat_slice(z_xcat)
    every = 2
    groups = [
        GROUP_OF_LABEL[xcat_labels[i, j, k] + 1]
            for i in 1:every:XCAT_DIMS[1], j in 1:every:XCAT_DIMS[2]
    ]

    fig = Mke.Figure(size = (1200, 720))
    ax = Mke.Axis(
        fig[1, 1];
        title = "XCAT labels · z = $(z_xcat) mm",
        aspect = Mke.DataAspect(),
        yreversed = true,
    )
    Mke.heatmap!(
        ax,
        xcat_x_mm[1:every:end],
        xcat_y_mm[1:every:end],
        groups;
        colormap = label_group_colormap(),
        colorrange = (-0.5, length(LABEL_GROUPS) + 0.5),
    )
    Mke.hidedecorations!(ax)

    ax3 = vessel_axis(
        fig[1, 2];
        title = "Imaged coronary vessels (NURBS centrelines), seen from the front",
    )
    vessels = [
        ("lad1", :red),
        ("lad2", :red),
        ("lad3", :red),
        ("lcx", :red),
        ("rca1", :red),
        ("rca2", :red),
        ("cor_vein1", :deepskyblue),
        ("cor_vein2", :deepskyblue),
    ]
    for (name, colour) in vessels
        line = V.XCAT.centerline(xcat_surfaces["dias_$name"]; step = 2.0)
        Mke.lines!(
            ax3,
            line.X[1, :],
            line.X[2, :],
            line.X[3, :];
            color = colour,
            linewidth = 3,
        )
    end

    Mke.Legend(
        fig[2, 1:2],
        [Mke.PolyElement(color = colour) for (_, colour, _) in LABEL_GROUPS],
        [name for (name, _, _) in LABEL_GROUPS];
        orientation = :horizontal,
        nbanks = 2,
        framevisible = false,
    )
    fig
end

# ╔═╡ 1b97b382-e826-4989-a7d7-90fd801bcfb6
md"""
## 2. Prepare the Growth Inputs

Growth runs on a 1 mm grid around the pericardium, into the LV and RV myocardium.

- **Inlet trees.** The imaged arteries become two trees (LCA, RCA), each with a **penetrating stub** about every 2.8 mm where it runs beside the wall, so growth starts along the whole vessel. The coronary veins become one tree draining into the coronary sinus.
- **Pass-through voxels.** The imaged lumens and short bridges to the wall: vessels may cross them, terminals may not.
- **Blockers.** Blood pools at 0.4 mm that no new vessel may cross.
"""

# ╔═╡ 985bdb9d-847b-446c-9e54-b9f9fe0be043
growth_grid = let (points, _) = V.XCAT.tessellate(
        xcat_surfaces["dias_pericardium"];
        nu = 60,
        nv = 60,
    )
    lower = floor.(Int, vec(minimum(points; dims = 2))) .- 2
    upper = ceil.(Int, vec(maximum(points; dims = 2))) .+ 2
    V.VoxelGrid(Tuple(Float64.(lower)), 1.0, Tuple(upper .- lower .+ 1))
end

# ╔═╡ c618deff-991a-43f2-a4ca-fbc0b65cf28e
heart_masks = cached("masks") do
    pericardium = V.voxelize(
        V.XCAT.tessellate(
            xcat_surfaces["dias_pericardium"];
            nu = 160,
            nv = 160,
        )...,
        growth_grid,
    )
    V.XCAT.label_masks(
        xcat_labels,
        XCAT_ORIGIN,
        XCAT_VOXEL_MM,
        growth_grid;
        roi = pericardium,
    )
end;

# ╔═╡ 0212424e-a0bb-4820-b933-f6f6a870c5b9
begin
    myocardium = heart_masks.lvmyo .| heart_masks.rvmyo
    forbid = (heart_masks.pool .& .!myocardium) .| heart_masks.lvcav .| heart_masks.rvcav
end;

# ╔═╡ 98c9c13d-a2a5-4316-84bd-715dbd077b88
"""
The myocardium's surface: each surface voxel's centre (mm) and its number among the myocardial voxels, so per-voxel
values (such as MBF) can be painted on it.
"""
MYOCARDIUM_SURFACE = let
    neighbours = (
        CartesianIndex(1, 0, 0),
        CartesianIndex(-1, 0, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 0, 1),
        CartesianIndex(0, 0, -1),
    )
    on_surface(I) = any(
        !checkbounds(Bool, myocardium, I + d) || !myocardium[I + d] for d in neighbours
    )
    voxels = findall(myocardium)
    surface = [n for (n, I) in enumerate(voxels) if on_surface(I)]
    (
        points = [Mke.Point3f(Tuple(voxels[n]) .- 1 .+ growth_grid.origin) for n in surface],
        voxels = surface,
    )
end;

# ╔═╡ 27920760-9210-402b-80d5-d13054040d6b
"Draw every other point of the myocardium's surface as a faint shell around the vessels."
function draw_myocardium_shell!(ax)
    return Mke.scatter!(
        ax,
        MYOCARDIUM_SURFACE.points[1:2:end];
        color = (:gray50, 0.07),
        markersize = 5,
        strokewidth = 0,
    )
end

# ╔═╡ ddfd8b55-868d-4c32-9d96-705b039efd06
"""
Draw a grown tree's segments of radius ≥ `min_radius` (mm): width and shade grow with radius, and the largest vessels
are drawn last so they sit on top.
"""
function draw_vessels!(ax, tree; colormap, min_radius)
    shown = filter(j -> tree.r[j] >= min_radius, sortperm(tree.r))
    segments = [
        Mke.Point3f(tree.a[:, j]) => Mke.Point3f(tree.b[:, j])
            for j in shown
    ]
    radius = tree.r[shown]
    # CairoMakie wants a colour and a width at each end of every segment
    return Mke.linesegments!(
        ax,
        segments;
        color = repeat(log10.(radius); inner = 2),
        colormap = colormap,
        colorrange = (log10(min_radius), log10(maximum(tree.r))),
        linewidth = repeat(1.5 .* sqrt.(radius ./ min_radius); inner = 2),
    )
end

# ╔═╡ fcca3b04-13a7-4274-b1fd-186837613b97
"Blood (chambers and great vessels) at 0.4 mm over the growth grid's box, as `(mask, origin_mm, cell_mm)`: the blockers."
BLOOD_BLOCKERS = cached("blood") do
    first_voxel = [
        round(Int, (growth_grid.origin[d] - 0.5 - XCAT_ORIGIN[d]) / XCAT_VOXEL_MM) + 1
            for d in 1:3
    ]
    cells = [ceil(Int, (growth_grid.dims[d] + 1) / 0.4) for d in 1:3]
    is_blood = falses(256)
    for label in (19, 20, 21, 22, 24, 25, 28)
        is_blood[label + 1] = true
    end

    blood = zeros(UInt8, cells...)
    Threads.@threads for c in 1:cells[3]
        for b in 1:cells[2], a in 1:cells[1], dk in 0:1, dj in 0:1, di in 0:1
            i = first_voxel[1] + 2(a - 1) + di
            j = first_voxel[2] + 2(b - 1) + dj
            k = first_voxel[3] + 2(c - 1) + dk
            (
                1 <= i <= XCAT_DIMS[1] &&
                    1 <= j <= XCAT_DIMS[2] &&
                    1 <= k <= XCAT_DIMS[3]
            ) || continue
            is_blood[xcat_labels[i, j, k] + 1] && (blood[a, b, c] = 0x01)
        end
    end
    origin = ntuple(
        d -> XCAT_ORIGIN[d] + XCAT_VOXEL_MM * (first_voxel[d] - 1) + XCAT_VOXEL_MM / 2,
        3,
    )
    (blood, origin, 0.4)
end;

# ╔═╡ 4af8c280-7734-4920-9683-73ff73ac3ed0
artery_inputs = cached("artery_inputs") do
    names = ["dias_lad1", "dias_lad2", "dias_lad3", "dias_lcx", "dias_rca1", "dias_rca2"]
    trees = V.XCAT.base_forest(
        [V.XCAT.centerline(xcat_surfaces[name]; step = 2.8) for name in names];
        roots = ["dias_lad1", "dias_rca1"],
    )
    stubbed = [
        V.XCAT.penetrating_stubs(
            tree,
            growth_grid,
            [heart_masks.lvmyo, heart_masks.rvmyo];
            domain = myocardium,
            forbid,
            blockers = BLOOD_BLOCKERS,
        )
            for tree in trees
    ]
    lumens = [
        V.XCAT.vessel_tube(tree, growth_grid) .|
            V.XCAT.vessel_bridges(tree, myocardium, growth_grid)
            for tree in trees
    ]
    (
        trees = first.(stubbed),
        stubs = last.(stubbed),
        pass = reduce(.|, lumens) .& .!forbid .& .!myocardium,
    )
end;

# ╔═╡ 28d82d8b-1f94-462f-8b3e-04cb56e127dc
vein_inputs = cached("vein_inputs") do
    vein1 = V.XCAT.centerline(xcat_surfaces["dias_cor_vein1"]; step = 2.8)
    vein2 = V.XCAT.centerline(xcat_surfaces["dias_cor_vein2"]; step = 2.8)

    # the coronary sinus ostium: the point of cor_vein1 nearest the right atrium, joined to it
    # by a 6 mm connector
    atrium, _ = V.XCAT.tessellate(
        xcat_surfaces["dias_ra"];
        nu = 120,
        nv = 120,
    )
    nearest_atrium(p) = begin
        d = vec(sqrt.(sum(abs2, atrium .- p; dims = 1)))
        k = argmin(d)
        (d[k], atrium[:, k])
    end
    ostium = argmin([nearest_atrium(vein1.X[:, k])[1] for k in 1:size(vein1.X, 2)])
    gap, atrium_point = nearest_atrium(vein1.X[:, ostium])
    start =
        vein1.X[:, ostium] .-
        6.0 .* (vein1.X[:, ostium] .- atrium_point) ./ max(gap, 1.0e-6)

    lines = [
        V.XCAT.Centerline(
            "cs_ostium",
            hcat(start, vein1.X[:, ostium]),
            fill(vein1.r[ostium], 2),
        ),
        V.XCAT.Centerline("cor_vein1a", vein1.X[:, ostium:-1:1], vein1.r[ostium:-1:1]),
        V.XCAT.Centerline("cor_vein1b", vein1.X[:, ostium:end], vein1.r[ostium:end]),
        vein2,
    ]
    base = only(
        V.XCAT.base_forest(
            lines;
            roots = ["cs_ostium"],
            attach_tol = 4.0,
        ),
    )
    # tributaries join the imaged veins along their length, as branches join the imaged arteries (imaged vessels are fixed
    # in growth, so without stubs new vessels could join only at their tips)
    tree, stubs = V.XCAT.penetrating_stubs(
        base,
        growth_grid,
        [heart_masks.lvmyo, heart_masks.rvmyo];
        domain = myocardium,
        forbid,
        blockers = BLOOD_BLOCKERS,
    )
    lumen =
        V.XCAT.vessel_tube(tree, growth_grid) .|
        V.XCAT.vessel_bridges(tree, myocardium, growth_grid)
    (tree = tree, stubs = stubs, pass = lumen .& .!forbid .& .!myocardium)
end;

# ╔═╡ 60afcebc-cc3d-482d-a811-3e1544584240
md"""
z (mm): $(@bind z_domain Slider(growth_grid.origin[3]:2:(growth_grid.origin[3] + growth_grid.dims[3] - 1); default = 255, show_value = true))
"""

# ╔═╡ ed8418d2-529c-4fad-8632-4d2ba5ba16e1
let k = clamp(round(Int, z_domain - growth_grid.origin[3]) + 1, 1, growth_grid.dims[3])
    function kind(i, j)
        heart_masks.lvmyo[i, j, k] && return 1
        heart_masks.rvmyo[i, j, k] && return 2
        artery_inputs.pass[i, j, k] && return 3
        vein_inputs.pass[i, j, k] && return 4
        heart_masks.pool[i, j, k] && return 5
        return 0
    end
    kinds = [kind(i, j) for i in 1:growth_grid.dims[1], j in 1:growth_grid.dims[2]]
    xs = growth_grid.origin[1] .+ (0:(growth_grid.dims[1] - 1))
    ys = growth_grid.origin[2] .+ (0:(growth_grid.dims[2] - 1))
    colours = [:darkorange, :goldenrod, :red, :deepskyblue, :mediumpurple]

    fig = Mke.Figure(size = (1200, 640))
    ax = Mke.Axis(
        fig[1, 1];
        title = "Growth domain · z = $(z_domain) mm",
        aspect = Mke.DataAspect(),
        yreversed = true,
    )
    Mke.heatmap!(
        ax,
        xs,
        ys,
        kinds;
        colormap = Mke.cgrad([:black; colours], 6; categorical = true),
        colorrange = (-0.5, 5.5),
    )
    Mke.hidedecorations!(ax)

    stub_counts = join(
        [
            "$name $(r.stubs) of $(r.sites)"
                for (name, r) in zip(["LCA", "RCA", "veins"], [artery_inputs.stubs; vein_inputs.stubs])
        ],
        ", ",
    )
    ax3 = vessel_axis(
        fig[1, 2];
        title = "Inlet trees: imaged (thick) and stubs (thin; $stub_counts)",
    )
    draw_myocardium_shell!(ax3)
    trees = [
        (artery_inputs.trees[1], :red),
        (artery_inputs.trees[2], :orangered),
        (vein_inputs.tree, :deepskyblue),
    ]
    for (tree, colour) in trees, k in 2:length(tree.parent)
        a = tree.X[:, tree.parent[k]]
        b = tree.X[:, k]
        stub = tree.r[k] == 0.1
        Mke.lines!(
            ax3,
            [a[1], b[1]],
            [a[2], b[2]],
            [a[3], b[3]];
            color = stub ? :black : colour,
            linewidth = stub ? 1 : 3,
        )
    end

    Mke.Legend(
        fig[2, 1:2],
        [Mke.PolyElement(color = c) for c in colours],
        [
            "LV myocardium",
            "RV myocardium",
            "artery pass-through",
            "vein pass-through",
            "blood (blocked)",
        ];
        orientation = :horizontal,
        framevisible = false,
    )
    fig
end

# ╔═╡ b7fd02d1-47f4-4252-aca5-1f4a3b0cdd75
md"""
## 3. Grow

HCT-CCO's staged constrained constructive optimization, on the GPU: `HEART_GROWTH` grows 957,045 arterial terminals from the two inlet trees, `HEART_VEIN_GROWTH` 100,000 venules into the coronary sinus.

The two coronaries are daughters of the aortic root, so they split its flow by the same bifurcation law as every other branch: share ∝ (imaged ostium radius)^γ, from XCAT's left main (2.54 mm) and proximal RCA (1.90 mm), so the RCA carries 34%. A tree grows only while it is at or below its share, and only where it is the near artery, so each claims its own territory: the RCA feeds the inferior LV and most of the RV.
"""

# ╔═╡ ec898038-0e7a-4f98-9737-fd711134217c
"Grow one forest: the domain and its pass-through voxels as a `.vol` file, the inlet trees as GXL files."
function grow_forest(spec, trees, pass, name)
    dir = mktempdir()
    lower = Tuple(round.(Int, growth_grid.origin))
    allow = zeros(UInt16, growth_grid.dims)
    allow[pass] .= V.NOTERM # vessels may cross these voxels, no terminal may end in them

    domain = V.write_vol(joinpath(dir, "$name.vol"), myocardium .| pass, lower)
    inlets = [
        V.XCAT.write_base_gxl(joinpath(dir, "$(name)_$k.xml"), tree)
            for (k, tree) in enumerate(trees)
    ]
    engine = V.grow(
        V.GrowthSpec(spec; device = DEVICE),
        domain,
        inlets;
        blockers = BLOOD_BLOCKERS,
        territories = allow,
    )
    return V.vtree(engine)
end

# ╔═╡ 1e11cee2-aa82-45e0-8965-2c09d62c44ae
arteries = cached("arteries") do
    grow_forest(V.HEART_GROWTH, artery_inputs.trees, artery_inputs.pass, "arteries")
end;

# ╔═╡ 13aff1b3-5b21-4270-8eaf-2f2ce16376a0
veins = cached("veins") do
    grow_forest(V.HEART_VEIN_GROWTH, [vein_inputs.tree], vein_inputs.pass, "veins")
end;

# ╔═╡ 9fc11643-e2a0-4e3f-88cb-17f89c5c90e3
md"""
show vessels of radius ≥ $(@bind min_radius_shown Select([0.3 => "0.3 mm", 0.15 => "0.15 mm", 0.08 => "0.08 mm"]; default = 0.15))
"""

# ╔═╡ 5da9b3a3-c69d-424b-9dd0-b4ef94fa0398
let fig = Mke.Figure(size = (1200, 560))
    ax = vessel_axis(
        fig[1, 1];
        title = "Grown arteries (red) and veins (blue), radius ≥ $(min_radius_shown) mm",
    )
    draw_myocardium_shell!(ax)
    draw_vessels!(
        ax,
        veins;
        colormap = VEIN_COLOURS,
        min_radius = min_radius_shown,
    )
    draw_vessels!(
        ax,
        arteries;
        colormap = ARTERY_COLOURS,
        min_radius = min_radius_shown,
    )

    terminal_radius_um = Dict(
        name => begin
            has_child = falses(length(tree.r))
            foreach(p -> p > 0 && (has_child[p] = true), tree.parent)
            1000 .* tree.r[.!has_child]
        end
            for (name, tree) in (("arteries", arteries), ("veins", veins))
    )
    ax2 = Mke.Axis(
        fig[1, 2];
        title = @sprintf(
            "Terminal radius (median %.1f µm arteries, %.1f µm veins)",
            median(terminal_radius_um["arteries"]),
            median(terminal_radius_um["veins"]),
        ),
        xlabel = "µm",
        ylabel = "terminals",
        xscale = log10,
    )
    for (name, colour) in (("arteries", :red3), ("veins", :royalblue))
        Mke.hist!(
            ax2,
            terminal_radius_um[name];
            bins = exp10.(range(0, 2.5; length = 50)),
            color = (colour, 0.5),
            label = name,
        )
    end
    Mke.axislegend(ax2)
    Mke.colsize!(fig.layout, 1, Mke.Relative(0.7))
    fig
end

# ╔═╡ bd4149bd-1699-4f34-845d-5e7272c46571
md"""
## 4. Flow

No flow is prescribed: aortic and coronary-sinus pressures drive blood through every vessel, capillary bed and vein, and ~1 million arterioles set their own tone from their O₂ balance. `flow` solves maximal dilation (stress) and rest. At maximal dilation the imaged epicardial arteries are also 13% wider, as adenosine widens normal human coronaries. Each 1 mm myocardial voxel gets its flow from the terminals that feed it.
"""

# ╔═╡ 4cb380c8-97b1-4fc4-b374-2832ed3f991c
"Centres (mm) of the 1 mm myocardial voxels: the tissue the flow is solved for."
tissue_centres = let voxels = findall(myocardium)
    Float32[growth_grid.origin[d] + voxels[k][d] - 1 for d in 1:3, k in eachindex(voxels)]
end

# ╔═╡ 1279319b-88cc-4214-adae-f8fac1701c52
flows = cached("flow") do
    V.flow(
        V.HEART_FLOW,
        arteries,
        veins,
        tissue_centres;
        to = DEVICE,
    )
end;

# ╔═╡ e442b143-f7da-4e3b-9151-6f81358788c5
"Blood flow (mL/min/g) of every 1 mm myocardial voxel in circuit `circuit`."
function voxel_mbf(circuit)
    terminal_flow = reshape(Float32[circuit.q[j] for j in circuit.terminals], 1, :)
    voxels = collect(Int32(1):Int32(size(tissue_centres, 2)))
    # mL/s per mL of tissue
    per_ml = V.OrganCT.tissue_iodine(terminal_flow, flows.pairs, flows.terr, voxels)
    return vec(per_ml) .* 60 ./ V.HEART_FLOW.rho
end

# ╔═╡ 52e13f82-7260-41db-8ebf-32149123ac4e
"A 1 mm map (x, y) of `values` (one per myocardial voxel) at z (mm); NaN elsewhere."
function myocardium_slice(values, z)
    k = clamp(round(Int, z - growth_grid.origin[3]) + 1, 1, growth_grid.dims[3])
    slice = fill(NaN32, growth_grid.dims[1], growth_grid.dims[2])
    for (n, voxel) in enumerate(findall(myocardium))
        voxel[3] == k && (slice[voxel[1], voxel[2]] = values[n])
    end
    return slice
end

# ╔═╡ ecf9a259-972e-4b13-a1f2-26cab5610bcf
md"""
z (mm): $(@bind z_flow Slider(growth_grid.origin[3]:2:(growth_grid.origin[3] + growth_grid.dims[3] - 1); default = 255, show_value = true))
"""

# ╔═╡ 68d1c5cd-fa2f-4961-9636-f6731ad142d6
let fig = Mke.Figure(size = (1200, 520))
    xs = growth_grid.origin[1] .+ (0:(growth_grid.dims[1] - 1))
    ys = growth_grid.origin[2] .+ (0:(growth_grid.dims[2] - 1))
    for (column, (state, circuit)) in enumerate(
            (("Rest", flows.rest), ("Stress", flows.stress)),
        )
        ax = Mke.Axis(
            fig[1, column];
            title = @sprintf("%s · organ MBF %.2f mL/min/g", state, circuit.resist.mbf),
            aspect = Mke.DataAspect(),
            yreversed = true,
        )
        map = Mke.heatmap!(
            ax,
            xs,
            ys,
            myocardium_slice(voxel_mbf(circuit), z_flow);
            colormap = :turbo,
            colorrange = (0, 6),
            nan_color = :black,
        )
        Mke.hidedecorations!(ax)
        column == 2 && Mke.Colorbar(
            fig[1, 3],
            map;
            label = "MBF (mL/min/g)",
        )
    end
    reserve = flows.stress.resist.mbf / flows.rest.resist.mbf
    Mke.Label(
        fig[0, :],
        @sprintf("Coronary flow reserve %.2f", reserve);
        fontsize = 20,
        font = :bold,
    )
    fig
end

# ╔═╡ e03d6af1-41ef-48c6-b867-4536d946cec4
md"""
## 5. Transport

The lab's current injection (`INJECTION`, below) into the right arm, carried through the whole-body circulation, the lungs and the left heart into the coronary tree, with recirculation. Adenosine stress raises cardiac output 1.71×, so the same injection is diluted into more blood and the arterial peak is lower.
"""

# ╔═╡ 0c4d6e8f-1a3b-4c5d-9e7f-2b4c6d8e0f13
"""
The injection: the lab's current protocol (Hadjiabdolhamid et al. 2024, Shojazadeh et al. 2026), 0.5 mL/kg of contrast at
5 mL/s, then 0.5 mL/kg of a 30:70 contrast:saline chaser at the same rate, for the body model's 70 kg: 35 mL of iohexol
350, then 35 mL at 105 mg/mL. A short bolus keeps the first pass inside the scan window: at stress, iodine starts leaving
the heart about 2.5 s after it arrives.
"""
const INJECTION = let kg = 70, conc = 350
    V.Body.Injection(; phases = [(0.5 * kg, 5, conc), (0.5 * kg, 5, 0.3 * conc)])
end

# ╔═╡ 3efb6e3b-8c83-43e6-8e17-7152d42f6c2d
"Contrast transport through the body and the coronary circuit `circuit`, at rest or adenosine stress."
function transport_state(state, circuit)
    heart = state == "stress" ? V.Body.adenosine() : (co = 81.6f0, beat = 1.0f0)
    spec = V.TransportSpec(
        V.HEART_TRANSPORT;
        phys = circuit.phys,
        tissue_ml = Float32(size(tissue_centres, 2) / 1000),
        heart...,
    )
    return (
        state,
        circuit,
        V.transport(
            spec,
            circuit,
            INJECTION;
            times = V.HEART_PHANTOM.times,
        ),
    )
end

# ╔═╡ 213ad6d9-995f-4b5c-9d87-549824095172
"""
Contrast transport at rest and stress through `circuits`, solved once and kept in `CACHE_DIR/name.jls`; returns that
file. The full results (every vessel's iodine at every frame) take several GB, so the notebook does not hold them:
`save_phantom` reads them while it writes, and the figures use `transport_curves`.
"""
function transport_cache(name, circuits)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) || serialize(
        path,
        [
            transport_state("rest", circuits.rest),
            transport_state("stress", circuits.stress),
        ],
    )
    return path
end

# ╔═╡ c2d2cdfc-a6dd-4892-be8d-9833d47b56b3
TRANSPORT = transport_cache("transport", flows);

# ╔═╡ 5f313e94-82d9-4335-9fae-0d8fcbec0fe8
"Each state's LV blood and whole-myocardium iodine (mg/mL) over time, from a transport file."
function transport_curves(path)
    tissue_ml = size(tissue_centres, 2) / 1000
    return map(deserialize(path)) do (state, _, result)
        held = vec(
            sum(result.tissue; dims = 2) .+
                sum(result.artery; dims = 2) .+
                sum(result.vein; dims = 2),
        )
        (
            state = state,
            blood_times = result.body.t,
            lv = result.body.mean[:lv],
            myocardium = held ./ tissue_ml,
        )
    end
end

# ╔═╡ 3ddafaeb-6b8c-4dfc-a037-20905115a1e7
TRANSPORT_CURVES = transport_curves(TRANSPORT);

# ╔═╡ f5b9f215-e830-4b40-bc5e-d5c3670a6e5d
let fig = Mke.Figure(size = (1200, 420))
    times = V.HEART_PHANTOM.times
    for (column, curve) in enumerate(TRANSPORT_CURVES)
        ax = Mke.Axis(
            fig[1, column];
            title = uppercasefirst(curve.state),
            xlabel = "time after injection (s)",
            ylabel = "iodine (mg/mL)",
        )
        Mke.lines!(
            ax,
            curve.blood_times,
            curve.lv;
            color = :firebrick,
            linewidth = 3,
            label = "LV blood",
        )
        Mke.lines!(
            ax,
            times,
            curve.myocardium;
            color = :darkorange,
            linewidth = 3,
            label = "myocardium (all its iodine)",
        )
        Mke.xlims!(ax, 0, 40)
        Mke.axislegend(ax)
    end
    fig
end

# ╔═╡ cf6401e0-6a84-4016-b458-7a17408173df
md"""
## 6. Compartments

A **compartment** is the region ground truth is given for: a region of the heart with everything in it, tissue, and every artery, vein and lumen that lies in it or overlaps it. Here: the LV and RV walls, each split into the LAD, LCx and RCA territories (the left main counts as LAD), and the **papillary muscles** as their own compartment, left out of the walls.

**Edit the compartments here**, e.g. `split_walls = false` and `split_territories = false` give one compartment, the whole myocardium.

The rule below is applied once, when the phantom is voxelized (section 7), and every compartment number and picture after that comes from the saved phantom (section 8), so there is one source of truth.
"""

# ╔═╡ 105100e4-ca4a-4098-ad30-42f0cb2ad7f0
begin
    split_walls = true # LV and RV apart
    split_territories = true # LAD, LCx and RCA apart
end;

# ╔═╡ 1488815e-8708-46f2-ade5-31b7438b0a41
begin
    wall_names = split_walls ? ["LV", "RV"] : ["Myocardium"]
    territory_names = split_territories ? ["LAD", "LCx", "RCA"] : String[]
    wall_rois =
        isempty(territory_names) ? wall_names :
        ["$w · $t" for w in wall_names for t in territory_names]
    ROI_NAMES = [wall_rois; "Papillary muscles"]

    roi_index(name) = findfirst(==(name), ROI_NAMES)
    ROI_GROUPS = Pair{String, Vector{Int}}[]
    length(wall_rois) > 1 && push!(ROI_GROUPS, "heart" => collect(eachindex(wall_rois)))
    if length(wall_names) > 1 && !isempty(territory_names)
        append!(
            ROI_GROUPS,
            [w => [roi_index("$w · $t") for t in territory_names] for w in wall_names],
        )
        append!(
            ROI_GROUPS,
            [t => [roi_index("$w · $t") for w in wall_names] for t in territory_names],
        )
    end
    ROI_NAMES, ROI_GROUPS
end

# ╔═╡ 2235d2e6-01bc-43ac-b6a4-1a8b5ece93c8
"""
Each arterial terminal's branch (`of`: 1 LAD, 2 LCx, 3 RCA) and the LAD's first segment (`lad`). The LCA splits at its
first imaged bifurcation where both branches feed ≥ 5% of it; the more anterior one is the LAD. Terminals of the left main
before it count as LAD.
"""
coronary_branch = let tree = arteries, n = length(tree.r)
    children = [Int[] for _ in 1:n]
    for j in 1:n
        tree.parent[j] > 0 && push!(children[tree.parent[j]], j)
    end
    root_of(j) = begin
        k = j
        while tree.parent[k] != 0
            k = tree.parent[k]
        end
        k
    end
    depth(j) = begin
        k = j
        d = 0
        while tree.parent[k] != 0
            k = tree.parent[k]
            d += 1
        end
        d
    end

    # terminals below each segment, and their summed y (front to back)
    terminals_below = Float64.(length.(children) .== 0)
    y_below = terminals_below .* tree.b[2, :]
    for j in sortperm(depth.(1:n); rev = true)
        p = tree.parent[j]
        p > 0 || continue
        terminals_below[p] += terminals_below[j]
        y_below[p] += y_below[j]
    end

    roots = unique(root_of.(1:n))
    lca = roots[argmax(terminals_below[roots])]
    j = lca
    lad = lcx = 0
    while lad == 0
        big = [
            c
                for c in children[j]
                if tree.isbase[c] && terminals_below[c] > 0.05 * terminals_below[lca]
        ]
        if length(big) >= 2
            front_to_back = sortperm([y_below[c] / terminals_below[c] for c in big])
            lad, lcx = big[front_to_back[1]], big[front_to_back[2]]
        else
            imaged = [c for c in children[j] if tree.isbase[c]]
            j = imaged[argmax(terminals_below[imaged])]
        end
    end

    below(j, a) = begin
        k = j
        while k != 0
            k == a && return true
            k = tree.parent[k]
        end
        false
    end
    branch_of(j) = root_of(j) != lca ? 3 : below(j, lcx) ? 2 : 1
    (of = branch_of.(flows.rest.terminals), lad = lad)
end;

# ╔═╡ a162fdf2-c33e-4724-8874-b0ab23c96cae
"Each 1 mm voxel's territory (1–3), carried 12 mm out from the myocardium so the lumens on its surface get one too."
territory_map = let weights = zeros(Float32, 3, size(tissue_centres, 2)),
        map = zeros(Int8, growth_grid.dims)
    for (voxel, terminal, w) in flows.pairs
        weights[coronary_branch.of[terminal], voxel] += w
    end
    voxels = findall(myocardium)
    for (n, voxel) in enumerate(voxels)
        map[voxel] = argmax(@view weights[:, n])
    end

    steps = (
        CartesianIndex(1, 0, 0),
        CartesianIndex(-1, 0, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 0, 1),
        CartesianIndex(0, 0, -1),
    )
    front = voxels
    for _ in 1:12
        next = CartesianIndex{3}[]
        for voxel in front, step in steps
            neighbour = voxel + step
            checkbounds(Bool, map, neighbour) && map[neighbour] == 0 || continue
            map[neighbour] = map[voxel]
            push!(next, neighbour)
        end
        front = next
    end
    map
end;

# ╔═╡ 1cb4005e-52ef-4b66-b4ea-02d7e982fc0c
"""
The coronary vessels a wall wraps, as one segmenting a patient's myocardium would include them: from each vessel voxel
(lumen or vessel wall) of `box`, 26 rays run through the vessel to the first voxel that is not vessel. A vessel's cross-section
(its connected voxels in one axial slice) joins the wall whose myocardium more than half of its rays end in, on average;
surface vessels, wrapped on one side only, stay out. Returns 0 (neither), 1 (LV) or 2 (RV) for every voxel of `box`.
"""
function enclosed_vessels(box; reach_mm = 8.0)
    vessel(label) = label == 0x1a || label == 0x1b || label == 0x1e
    directions = [
        (i, j, k) ./ sqrt(i^2 + j^2 + k^2)
            for i in -1:1, j in -1:1, k in -1:1 if (i, j, k) != (0, 0, 0)
    ]
    steps = round(Int, reach_mm / XCAT_VOXEL_MM)
    sides = (
        CartesianIndex(1, 0),
        CartesianIndex(-1, 0),
        CartesianIndex(0, 1),
        CartesianIndex(0, -1),
    )
    wall = zeros(UInt8, size(box))
    Threads.@threads for k in axes(box, 3)
        # the share of each vessel voxel's rays that end in each wall
        lv_share = zeros(Float32, size(box, 1), size(box, 2))
        rv_share = zeros(Float32, size(box, 1), size(box, 2))
        for j in axes(box, 2), i in axes(box, 1)
            vessel(box[i, j, k]) || continue
            for d in directions
                for s in 1:steps
                    q = round.(Int, (i, j, k) .+ s .* d)
                    checkbounds(Bool, box, q...) || break
                    label = box[q...]
                    vessel(label) && continue
                    label == 0x0f && (lv_share[i, j] += 1 / length(directions))
                    label == 0x10 && (rv_share[i, j] += 1 / length(directions))
                    break
                end
            end
        end
        # one decision per cross-section: flood each connected patch of vessel in this slice
        seen = falses(size(box, 1), size(box, 2))
        patch = CartesianIndex{2}[]
        for j in axes(box, 2), i in axes(box, 1)
            (vessel(box[i, j, k]) && !seen[i, j]) || continue
            empty!(patch)
            stack = [CartesianIndex(i, j)]
            seen[i, j] = true
            while !isempty(stack)
                c = pop!(stack)
                push!(patch, c)
                for side in sides
                    n = c + side
                    checkbounds(Bool, seen, n) || continue
                    (seen[n] || !vessel(box[n, k])) && continue
                    seen[n] = true
                    push!(stack, n)
                end
            end
            lv = sum(lv_share[c] for c in patch) / length(patch)
            rv = sum(rv_share[c] for c in patch) / length(patch)
            lv + rv > 0.5 || continue
            for c in patch
                wall[c, k] = lv >= rv ? 0x01 : 0x02
            end
        end
    end
    return wall
end

# ╔═╡ a648cc82-ba31-49ce-8c88-b670dbf3a513
"""
Where the LV and RV walls would be, at 0.2 mm: each wall closed by 3 mm, so a vessel lumen in or overlapping it is in it,
and the vessels each wall wraps (`enclosed_vessels`). And the papillary muscles: the LV myocardium that a 10 mm closing of
the LV cavity fills (thin muscle surrounded by blood).
"""
WALL_ENVELOPES = cached("wall_envelopes_enclosed") do
    walls = findall(label -> label == 0x0f || label == 0x10, xcat_labels)
    margin = ceil(Int, 12 / XCAT_VOXEL_MM)
    lower = max.(Tuple(minimum(walls)) .- margin, 1)
    upper = min.(Tuple(maximum(walls)) .+ margin, XCAT_DIMS)
    box = xcat_labels[CartesianIndices(UnitRange.(lower, upper))]
    (
        # the box's first voxel centre in mm, so the masks stay right whatever grid the labels come on
        origin_mm = XCAT_ORIGIN .+ (lower .- 1) .* XCAT_VOXEL_MM,
        lv = V.XCAT.closing(box .== 0x0f, 3 / XCAT_VOXEL_MM),
        rv = V.XCAT.closing(box .== 0x10, 3 / XCAT_VOXEL_MM),
        papillary = V.XCAT.closing(box .== 0x13, 10 / XCAT_VOXEL_MM) .& (box .== 0x0f),
        enclosed = enclosed_vessels(box),
    )
end;

# ╔═╡ ba8d75ec-ca1a-4952-970b-a31cdd5be016
"""
    make_roi(labels, envelopes, territories, grid, n_rois; split_walls, split_territories) -> roi

The compartment rule `roi(p)` for the 0.2 mm voxel at `p` (0 = none): the papillary muscles (the last compartment); else its
wall (XCAT's label, or for a vessel lumen or vessel wall the wall envelope it lies in) and its 1 mm voxel's territory.
"""
function make_roi(
        labels,
        envelopes,
        territories,
        grid,
        n_rois;
        split_walls,
        split_territories,
    )
    lv = envelopes.lv
    rv = envelopes.rv
    papillary = envelopes.papillary
    return function roi(p)
        voxel = round.(Int, (p .- XCAT_ORIGIN) ./ XCAT_VOXEL_MM) .+ 1
        all(1 .<= voxel .<= XCAT_DIMS) || return 0
        label = labels[voxel...]
        k = round.(Int, (p .- envelopes.origin_mm) ./ XCAT_VOXEL_MM) .+ 1
        inside = all(1 .<= k .<= size(lv))
        inside && papillary[k...] && return n_rois

        wall = label == 0x0f ? 1 : label == 0x10 ? 2 : 0
        # coronary artery, vein or vessel wall: in a wall's envelope, or wrapped by a wall
        if wall == 0 && inside && label in (0x1a, 0x1b, 0x1e)
            wall = lv[k...] ? 1 : rv[k...] ? 2 : Int(envelopes.enclosed[k...])
        end
        wall == 0 && return 0

        g = round.(Int, p .- grid.origin) .+ 1
        territory = all(1 .<= g .<= size(territories)) ? Int(territories[g...]) : 0
        territory == 0 && return 0

        wall = split_walls ? wall : 1
        return split_territories ? 3 * (wall - 1) + territory : wall
    end
end

# ╔═╡ 48514abd-afbb-44c0-8f3d-f3293da05204
roi = make_roi(
    xcat_labels,
    WALL_ENVELOPES,
    territory_map,
    growth_grid,
    length(ROI_NAMES);
    split_walls,
    split_territories,
)

# ╔═╡ a19d475f-edfe-41c7-ad61-2ea8bdf170bd
"""
Categorical colours of the compartments: hue for the territory (LAD red, LCx green, RCA blue), dark for the LV and light
for the RV, gold for the papillary muscles.
"""
function roi_colormap()
    lv = ["#c62828", "#2e7d32", "#1565c0"]
    rv = ["#ef9a9a", "#a5d6a7", "#90caf9"]
    walls =
        split_territories ? (split_walls ? [lv; rv] : lv) :
        (split_walls ? ["#c62828", "#90caf9"] : ["#c62828"])
    return Mke.cgrad(
        [:black; walls; "#ffc400"],
        length(ROI_NAMES) + 1;
        categorical = true,
    )
end

# ╔═╡ 44638d83-bacc-489c-bbc1-f28ce8ff43ac
md"""
## 7. Voxelize and Save

The 0.2 mm phantom of the whole XCAT heart, one folder per scenario. Every voxel gets a **label** (its material and iodine curve) and a **compartment**. Every label lies in one compartment, so a compartment's iodine and ground truth are sums over its labels, and painting conserves iodine exactly: the build fails if the heart's iodine ledger does not balance.

| File | What it holds |
|---|---|
| `labels.u16.zst` | every voxel's label |
| `rois.u8.zst` | every voxel's compartment |
| `mbf_rest.f32.zst`, `mbf_stress.f32.zst` | every voxel's ground-truth MBF (mL/min/g) |
| `label_rest.csv`, `label_stress.csv` | each label's compartment, material, ground-truth MBF and iodine every second |
| `ledger_rest.csv`, `ledger_stress.csv` | the heart's iodine every second: delivered, held, out through the veins, painted |
| `phantom.toml` | the grid, the compartment names and groups, the tracking point |
"""

# ╔═╡ ec88e71f-17a6-42ab-ae1d-131b91a396df
begin
    "1 mm myocardial voxel number of each growth-grid voxel (0: not myocardium)."
    const TISSUE_VOXEL_NUMBER = V.XCAT.tissue_numbers(myocardium)

    # the imaged segments of the grown trees, for naming XCAT's drawn coronary lumens
    const IMAGED_ARTERY_SEGMENTS = V.XCAT.imaged_segments(arteries)
    const IMAGED_VEIN_SEGMENTS = V.XCAT.imaged_segments(veins)
end;

# ╔═╡ 5b0c2f7e-3d41-4a8e-9b6f-2c7e1d9a4f10
"""
Which great vessel each vessel voxel belongs to, and where it lies along it, on a 1 mm grid: its tree and the blood (mL)
between it and the heart. The trees grow by connection through XCAT's vessel labels from where they meet the heart: the
aorta from the LV, the pulmonary arteries from the RV, the venae cavae and brachiocephalic veins from the RA, the pulmonary
veins from the LA (`V.XCAT.great_vessel_transit`). Each state's transit is volume over flow (`V.pool_delay`).

The aorta is timed by the body model's volume, not XCAT's. XCAT vmale_50's aorta is about three standard deviations
narrower than adult men's (21–22 mm ascending and 15 mm descending, against 33 ± 4 and 24 ± 3 mm; Wolak et al. 2008), so its
lumen holds 26.7 mL from the valve to the arch's top. The whole-body model that makes the arterial curves holds 60 mL there
(aortic root 25 mL, arch 35 mL). Each aortic voxel keeps its place along XCAT's aorta, and the blood upstream of it is
scaled by that ratio (`V.XCAT.scale_tree`), so the arterial curves and their timing come from one body.
"""
GREAT_VESSEL_TRANSIT = let
    xcat = cached("great_vessel_transit") do
        V.XCAT.great_vessel_transit(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM)
    end
    body = V.circulation(V.HEART_TRANSPORT)
    volume(segment) = body.segs[findfirst(s -> s.name === segment, body.segs)].v
    V.XCAT.scale_tree(xcat, 0x01, (volume(:aortic_root) + volume(:arch)) / xcat.apex_ml)
end;

# ╔═╡ 2f6d8c41-7a3e-4b9c-8e15-d4a7c0b93e62
"The heart's phantom spec: the body pools, with the descending aorta timed from this subject's arch."
PHANTOM_SPEC = let pools = copy(V.HEART_PHANTOM.pools)
    i = findfirst(p -> p.name == "descending aorta", pools)
    d = pools[i]
    pools[i] = V.Pool(
        d.name,
        d.material,
        d.segment,
        d.what;
        through = d.through,
        lead_ml = GREAT_VESSEL_TRANSIT.apex_ml,
        lead_through = d.lead_through,
    )
    V.PhantomSpec(V.HEART_PHANTOM; pools)
end;

# ╔═╡ 0fbff6a0-6c59-4259-9974-b51529affa05
"""
    make_classifier(labels, grid, numbers, artery_segments, vein_segments, transit) -> classify

The classifier `classify(p)`: what is at the 0.2 mm voxel centre `p`. Myocardium is its 1 mm voxel, a coronary lumen its
imaged segment, a great-vessel voxel its vessel's pool and its place along the vessel (`transit`); every other XCAT
label is decided by the label alone (a body pool, `V.XCAT_POOL_OF_LABEL`, or a static material). Everything it reads is
captured with its type, so it compiles to fast code.
"""
function make_classifier(
        labels,
        grid,
        numbers,
        artery_segments,
        vein_segments,
        transit,
    )
    # by the label alone: 0 needs the voxel's geometry
    by_label = fill((V.VOX_POOL, Int32(24), 0.0f0), 256)
    needs_geometry = falses(256)
    for label in 0x00:0xff
        if label in (0x0f, 0x10, 0x1a, 0x1b, 0x18, 0x19, 0x1c, 0x7c, 0x7d)
            needs_geometry[label + 1] = true
        elseif haskey(V.XCAT_POOL_OF_LABEL, label)
            by_label[label + 1] = (V.VOX_POOL, Int32(V.XCAT_POOL_OF_LABEL[label]), 0.0f0)
        elseif haskey(V.XCAT_STATIC, Int(label))
            by_label[label + 1] = (V.VOX_STATIC, Int32(label), 0.0f0)
        end
    end
    vessel_pool(p, t) = V.XCAT.great_vessel_pool(transit, p, t)

    return function classify(p)
        voxel = round.(Int, (p .- XCAT_ORIGIN) ./ XCAT_VOXEL_MM) .+ 1
        all(1 .<= voxel .<= XCAT_DIMS) || return (V.VOX_STATIC, Int32(0), 0.0f0)
        label = labels[voxel...]
        needs_geometry[label + 1] || return by_label[label + 1]

        if label == 0x0f || label == 0x10
            number = V.XCAT.nearest_tissue_voxel(numbers, round.(Int, p .- grid.origin) .+ 1)
            return number > 0 ?
                (V.VOX_TISSUE, number, 0.0f0) :
                (V.VOX_POOL, Int32(26), 0.0f0)
        end
        label == 0x1a && return (
            V.VOX_VESSEL,
            Int32(first(V.XCAT.nearest_segment(artery_segments, p))),
            0.0f0,
        )
        label == 0x1b && return (
            V.VOX_VESSEL,
            -Int32(first(V.XCAT.nearest_segment(vein_segments, p))),
            0.0f0,
        )
        label == 0x1c && return vessel_pool(p, 0x01)
        label == 0x19 && return vessel_pool(p, V.XCAT.vessel_tree_at(transit, p, (0x01, 0x02)))
        label == 0x7d && return vessel_pool(p, 0x02)
        label == 0x7c && return vessel_pool(p, 0x05)
        return vessel_pool(p, V.XCAT.vessel_tree_at(transit, p, (0x03, 0x04, 0x05, 0x06, 0x07)))
    end
end

# ╔═╡ 23550e90-a7a1-4923-b7ca-f0910276f80a
classify = make_classifier(
    xcat_labels,
    growth_grid,
    TISSUE_VOXEL_NUMBER,
    IMAGED_ARTERY_SEGMENTS,
    IMAGED_VEIN_SEGMENTS,
    GREAT_VESSEL_TRANSIT,
)

# ╔═╡ b7f0823b-85d9-4ed5-8096-9818ec1877f8
"""
The bolus-tracking points: the descending aorta at the myocardium's mid level, and the ascending aorta 15 mm above where it
leaves the heart (clear of the root's sinuses), the centre of its cross-section there.
"""
TRACKING_POINTS = let voxels = findall(myocardium)
    z_mid =
        (minimum(v[3] for v in voxels) + maximum(v[3] for v in voxels)) / 2 +
        growth_grid.origin[3] - 1
    k = xcat_slice(z_mid)
    aorta = [
        (i, j)
            for j in 1:XCAT_DIMS[2], i in 1:XCAT_DIMS[1]
            if xcat_labels[i, j, k] == 0x1c
    ]
    back = maximum(last, aorta)
    descending = [
        (i, j)
            for (i, j) in aorta
            if j > back - 15 / XCAT_VOXEL_MM # within 15 mm of its back
    ]
    descending_point = (
        XCAT_ORIGIN[1] + (sum(first, descending) / length(descending) - 1) * XCAT_VOXEL_MM,
        XCAT_ORIGIN[2] + (sum(last, descending) / length(descending) - 1) * XCAT_VOXEL_MM,
        XCAT_ORIGIN[3] + (k - 1) * XCAT_VOXEL_MM,
    )

    # the ascending aorta: the aortic tree's cells before the arch's top, 15 mm above the lowest of them
    transit = GREAT_VESSEL_TRANSIT
    step = round(Int, transit.h / XCAT_VOXEL_MM)
    ascending = [
        c
            for c in CartesianIndices(transit.tree)
            if transit.tree[c] == 0x01 && !transit.descending[c] && !isnan(transit.volume[c]) &&
            transit.volume[c] < transit.apex_ml && xcat_labels[((Tuple(c) .- 1) .* step .+ 1)...] == 0x1c
    ]
    z_low = minimum(c[3] for c in ascending)
    level = [c for c in ascending if c[3] == z_low + round(Int, 15 / transit.h)]
    ascending_point = Tuple(
        transit.origin_mm[d] + (sum(c[d] for c in level) / length(level) - 1) * transit.h for d in 1:3
    )
    ["descending aorta" => descending_point, "ascending aorta" => ascending_point]
end

# ╔═╡ 2fe6d117-a55a-412c-bf8b-10cd9444056e
"""
Write one scenario's phantom into `PHANTOM_DIR/name`, once, from its transport file (read here, and let go once the
phantom is written).
"""
function save_phantom(
        name,
        transport_file;
        spec = PHANTOM_SPEC,
        classify = classify,
    )
    dir = joinpath(PHANTOM_DIR, name)
    isfile(joinpath(dir, "phantom.toml")) && return dir

    scans = deserialize(transport_file)
    seconds = @elapsed V.phantom(
        spec,
        scans,
        classify,
        V.VoxelGrid(XCAT_ORIGIN, XCAT_VOXEL_MM, XCAT_DIMS);
        pairs = flows.pairs,
        terr = flows.terr,
        nvox = size(tissue_centres, 2),
        out = dir,
        roi = roi,
        roi_names = ROI_NAMES,
        roi_groups = ROI_GROUPS,
        tracking = TRACKING_POINTS,
        facts = ["organ" => "heart (XCAT vmale_50, 0.2 mm)", "model" => "Vasculature.jl"],
        verbose = true,
    )
    @info "$name: $(round(seconds; digits = 1)) s"
    scans = nothing
    GC.gc()
    return dir
end

# ╔═╡ 487f30c5-950f-4177-b5b8-4c50cf4ca156
HEALTHY_DIR = save_phantom("healthy", TRANSPORT);

# ╔═╡ ef9926c5-a216-4323-abe3-80796ebab79a
md"""
## 8. Inspect the Saved Phantom

Only the heart's slab is read back. Labels and compartments are categorical; the idealized 70 keV image (the phantom as a perfect scanner would see it) is gray; iodine and ground-truth MBF use `:turbo`.

Ground-truth MBF is shown per compartment, the number a CT measurement of that region is compared with: the flow the compartment's capillary beds receive over the mass of everything in it, tissue and the vessels in it alike. It covers the whole compartment, vessels included, as a CT region does.
"""

# ╔═╡ 95d13834-2f64-4318-a83b-7a447e8f9e5b
"""
The heart's box of a saved phantom, the only part read back: the growth grid's (the pericardium's) box and the
descending aorta's tracking point, with 15 mm to spare. Returns the bundle (labels and
compartments) and its first z (mm). Ground-truth MBF is each compartment's own (its beds' flow over the mass of
everything in it, tissue and vessels), from the bundle's tables, so the `mbf_<state>.f32.zst` volumes are not read back.
"""
function inspect(dir)
    tracking = V.DynamicPhantom.load(dir; volume = false).meta["tracking"]
    aorta_mm = tracking["descending aorta"]["point_mm"]
    heart_lower_mm = collect(growth_grid.origin)
    heart_upper_mm = heart_lower_mm .+ collect(growth_grid.dims) .- 1
    lower_mm = [min.(heart_lower_mm[1:2], aorta_mm[1:2]); heart_lower_mm[3]] .- 15
    upper_mm = [max.(heart_upper_mm[1:2], aorta_mm[1:2]); heart_upper_mm[3]] .+ 15
    function voxel(mm, d)
        number = round(Int, (mm - XCAT_ORIGIN[d]) / XCAT_VOXEL_MM) + 1
        return clamp(number, 1, XCAT_DIMS[d])
    end
    bundle = V.DynamicPhantom.load(
        dir;
        rois = true,
        zrange = voxel(lower_mm[3], 3):voxel(upper_mm[3], 3),
        xrange = voxel(lower_mm[1], 1):voxel(upper_mm[1], 1),
        yrange = voxel(lower_mm[2], 2):voxel(upper_mm[2], 2),
    )
    return (bundle = bundle, z0 = bundle.origin_mm[3])
end

# ╔═╡ 2b074062-ec5b-42d2-9b94-95329585ea3d
healthy = inspect(HEALTHY_DIR);

# ╔═╡ 3fc08b34-bd48-4687-98bb-efdfb79bd05e
"""
One phantom slice in six panels, each with its own scale underneath: labels by kind, compartments and ground-truth MBF in
`state`; the idealized 70 keV image and iodine at time `t`; ground-truth MBF in the other state.
"""
function slice_figure(view, state, t, z)
    bundle = view.bundle
    k = clamp(round(Int, (z - view.z0) / XCAT_VOXEL_MM) + 1, 1, size(bundle.labels, 3))
    is, js = axes(bundle.labels, 1), axes(bundle.labels, 2)
    xs = bundle.origin_mm[1] .+ (is .- 1) .* XCAT_VOXEL_MM
    ys = bundle.origin_mm[2] .+ (js .- 1) .* XCAT_VOXEL_MM
    labels = bundle.labels[is, js, k]
    other = state == "rest" ? "stress" : "rest"

    # per-label lookups: kind, idealized HU and iodine at time t
    kind_colours = [:gray30, :gray65, :mediumpurple, :darkorange, :red, :deepskyblue]
    kind_names = [
        "static (air, cartilage)",
        "body tissue",
        "blood",
        "myocardium",
        "artery",
        "vein",
    ]
    function kind_code(row)
        row.kind == "static" && return 1
        row.kind == "pool" && return occursin("blood", row.base) ? 3 : 2
        return Dict("tissue" => 4, "artery" => 5, "vein" => 6)[row.kind]
    end
    kind_of = zeros(UInt8, 65536)
    for (label, row) in bundle.table
        kind_of[label + 1] = label == 0 ? 0 : kind_code(row)
    end
    # each compartment's ground truth over all of it, tissue and vessels: its beds' flow over the mass of everything in it
    mbf = Dict(st => zeros(Float32, length(ROI_NAMES) + 1) for st in ("rest", "stress"))
    for st in ("rest", "stress"), (r, name) in enumerate(ROI_NAMES)
        mbf[st][r + 1] = bundle.roi_table[(st, name)].mbf
    end
    rois = bundle.rois[is, js, k]
    hu = fill(-1000.0f0, 65536)
    for (label, value) in V.DynamicPhantom.hu(bundle, state, t)
        hu[label + 1] = value
    end
    iodine = zeros(Float32, 65536)
    for (label, value) in V.DynamicPhantom.iodine(bundle, state, t)
        iodine[label + 1] = value
    end

    fig = Mke.Figure(; size = (1500, 1120))
    function panel(
            row,
            column,
            title,
            image;
            colormap,
            colorrange,
            label = nothing,
            legend = nothing,
        )
        cell = fig[row, column] = Mke.GridLayout()
        ax = Mke.Axis(
            cell[1, 1];
            title,
            aspect = Mke.DataAspect(),
            yreversed = true,
        )
        map = Mke.heatmap!(
            ax,
            xs,
            ys,
            image;
            colormap,
            colorrange,
        )
        Mke.hidedecorations!(ax)
        if legend === nothing
            Mke.Colorbar(
                cell[2, 1],
                map;
                label,
                vertical = false,
                flipaxis = false,
                width = Mke.Relative(0.9),
                tellwidth = false,
            )
        else
            colours, names = legend
            Mke.Legend(
                cell[2, 1],
                [Mke.PolyElement(; color = c) for c in colours],
                names;
                orientation = :horizontal,
                nbanks = 3,
                framevisible = false,
                tellwidth = false,
            )
        end
        return ax
    end

    panel(
        1,
        1,
        "Labels by kind",
        [kind_of[l + 1] for l in labels];
        colormap = Mke.cgrad([:black; kind_colours], 7; categorical = true),
        colorrange = (-0.5, 6.5),
        legend = (kind_colours, kind_names),
    )
    panel(
        1,
        2,
        "Compartments",
        rois;
        colormap = roi_colormap(),
        colorrange = (-0.5, length(ROI_NAMES) + 0.5),
        legend = (roi_colormap().colors[2:end], ROI_NAMES),
    )
    panel(
        1,
        3,
        "Ground-truth MBF by compartment, $state",
        [mbf[state][r + 1] for r in rois];
        colormap = :turbo,
        colorrange = (0, 6),
        label = "mL/min/g",
    )
    panel(
        2,
        1,
        "Idealized 70 keV, $state, t = $t s",
        [hu[l + 1] for l in labels];
        colormap = :grays,
        colorrange = (-200, 600),
        label = "HU",
    )
    panel(
        2,
        2,
        "Iodine, $state, t = $t s",
        [iodine[l + 1] for l in labels];
        colormap = :turbo,
        colorrange = (0, 8),
        label = "mg/mL",
    )
    panel(
        2,
        3,
        "Ground-truth MBF by compartment, $other",
        [mbf[other][r + 1] for r in rois];
        colormap = :turbo,
        colorrange = (0, 6),
        label = "mL/min/g",
    )
    Mke.Label(
        fig[0, :],
        "z = $(round(z; digits = 1)) mm";
        fontsize = 20,
        font = :bold,
    )
    # three equal columns: no legend or colour bar widens its own
    for column in 1:3
        Mke.colsize!(fig.layout, column, Mke.Relative(1 / 3))
    end
    return fig
end

# ╔═╡ 53be3442-4eb3-4547-95b6-5fc32f670646
md"""
* state $(@bind healthy_state Select(["rest", "stress"]))
* time $(@bind healthy_time Slider(0:1:40; default = 20, show_value = true)) s
* z $(@bind healthy_z Slider(round.(healthy.z0 .+ (0:5:(size(healthy.bundle.labels, 3) - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 259.0, show_value = true)) mm
"""

# ╔═╡ e44cab6e-9bc9-44e7-a481-08938d5fac32
slice_figure(healthy, healthy_state, healthy_time, healthy_z)

# ╔═╡ 3c468348-759d-463a-bd1c-9ded251211cc
"Each compartment's iodine, with the tracking arrival and peak, and the heart's iodine ledger at stress."
function curves_figure(view)
    bundle = view.bundle
    fig = Mke.Figure(size = (1500, 540))
    for (column, state) in enumerate(("rest", "stress"))
        ax = Mke.Axis(
            fig[1, column];
            title = "Compartments, $state",
            xlabel = "time after injection (s)",
            ylabel = "iodine (mg)",
        )
        for name in ROI_NAMES
            row = bundle.roi_table[(state, name)]
            Mke.lines!(
                ax,
                bundle.times,
                row.mass;
                linewidth = 2,
                label = @sprintf("%s · %.2f mL/min/g", name, row.mbf),
            )
        end
        timing = V.DynamicPhantom.scan_times(bundle; state, roi = "roi:descending aorta")
        Mke.vlines!(
            ax,
            [timing.arrival, timing.peak_time];
            color = :gray40,
            linestyle = :dash,
        )
        Mke.xlims!(ax, 0, 40)
        column == 1 && Mke.axislegend(
            ax;
            position = :lt,
            labelsize = 11,
        )
    end

    ledger = bundle.ledger["stress"]
    ax = Mke.Axis(
        fig[1, 3];
        title = "Iodine ledger, stress",
        xlabel = "time after injection (s)",
        ylabel = "iodine (mg)",
    )
    entries = [
        ("delivered_mg", :black, "delivered", :solid),
        ("held_mg", :darkorange, "held", :solid),
        ("outflow_mg", :royalblue, "out through the veins", :solid),
        ("painted_mg", :red, "painted", :dash),
    ]
    for (column, colour, name, style) in entries
        Mke.lines!(
            ax,
            ledger["t_s"],
            ledger[column];
            color = colour,
            linewidth = 2,
            linestyle = style,
            label = name,
        )
    end
    Mke.xlims!(ax, 0, 40)
    Mke.axislegend(ax; position = :lt)
    return fig
end

# ╔═╡ ab9d2337-2a36-4416-97ad-e9ad2c753d4b
"""
The share of the delivered iodine that has already left the heart through the coronary sinus, 1, 2 and 3 s after
arrival, 1 s before the arterial peak and at the peak.
"""
function outflow_check(view)
    bundle = view.bundle
    function share_out(ledger, t)
        ts = ledger["t_s"]
        k = clamp(searchsortedlast(ts, t), 1, length(ts) - 1)
        w = (t - ts[k]) / (ts[k + 1] - ts[k])
        out = (1 - w) * ledger["outflow_mg"][k] + w * ledger["outflow_mg"][k + 1]
        delivered = (1 - w) * ledger["delivered_mg"][k] + w * ledger["delivered_mg"][k + 1]
        return out / max(delivered, eps())
    end
    rows = map(("rest", "stress")) do state
        ledger = bundle.ledger[state]
        timing = V.DynamicPhantom.scan_times(bundle; state, roi = "roi:descending aorta")
        times = [timing.arrival .+ [1, 2, 3]; timing.peak_time - 1; timing.peak_time]
        shares = [@sprintf("%.0f%%", 100 * share_out(ledger, t)) for t in times]
        @sprintf(
            "| **%s** (arrival at %.1f s, peak at %.1f s) | %s |",
            state,
            timing.arrival,
            timing.peak_time,
            join(shares, " | "),
        )
    end
    text = """
    !!! note "Iodine already out of the heart"
        | | arrival + 1 s | arrival + 2 s | arrival + 3 s | 1 s before the peak | arterial peak |
        |---|---|---|---|---|---|
    $(join("    " .* rows, "\n"))

        A measurement that counts the iodine the heart has taken up misses this share.
    """
    return Markdown.parse(text)
end

# ╔═╡ fa0357e5-2d1e-4354-a12d-7e1c8c308c3a
curves_figure(healthy)

# ╔═╡ 173bb953-3b21-4c4a-aa00-9c838e0b21e1
outflow_check(healthy)


# ╔═╡ 948065f2-a977-4980-9eba-dec19b754a5a
md"""
## 9. A 90% LAD Stenosis

The same heart with a lesion over the first 10 mm of the proximal LAD, graded on the imaged lumen. Flow is solved again with the healthy tone constant, so the arterioles downstream relax on their own. The lesion is drawn with plaque (fibrous tissue, a lipid core, calcification) in place of the lumen it takes, and saved as its own scenario.
"""

# ╔═╡ 872d6939-4654-4b12-82be-58aed8a9f85a
const DIAMETER_STENOSIS = 0.9

# ╔═╡ 168b517e-ea3e-4933-b0b3-1f85f31a8512
"The lesion along the LAD's first 10 mm of imaged path: its segments, their imaged radii and lengths, and `circuit`'s list."
LESION = let tree = arteries, n = length(tree.r)
    children = [Int[] for _ in 1:n]
    for j in 1:n
        tree.parent[j] > 0 && push!(children[tree.parent[j]], j)
    end
    segment_length(j) = sqrt(sum(abs2, tree.b[:, j] .- tree.a[:, j]))

    path = Int[]
    covered = 0.0f0
    j = coronary_branch.lad
    while covered < 10.0f0
        push!(path, j)
        covered += segment_length(j)
        imaged = [c for c in children[j] if tree.isbase[c]]
        isempty(imaged) && break
        j = imaged[argmax([segment_length(c) * tree.r[c]^2 for c in imaged])]
    end

    # the imaged lumen radius at each segment: XCAT's LAD centreline radius nearest to it
    lad = V.XCAT.centerline(xcat_surfaces["dias_lad1"]; step = 0.5)
    imaged_radius(j) = lad.r[
        argmin(vec(sum(abs2, lad.X .- (tree.a[:, j] .+ tree.b[:, j]) ./ 2; dims = 1))),
    ]
    radii = Float32[imaged_radius(j) for j in path]
    lengths = Float32[segment_length(j) for j in path]
    area_fraction(j, r) = min(1.0f0, (r * (1 - DIAMETER_STENOSIS) / tree.r[j])^2)
    (
        path = path,
        radii = radii,
        stenosis = [(j, area_fraction(j, r), l) for (j, r, l) in zip(path, radii, lengths)],
    )
end

# ╔═╡ cbf27751-a278-460d-9b99-5db1be2e9bf1
flows_lad90 = cached("flow_lad90") do
    stress = V.Flow.circuit(
        arteries,
        flows.terr;
        state = :stress,
        veins,
        vleaf = flows.vleaf,
        reference = flows.stress,
        kappa = flows.stress.resist.kappa,
        stenosis = LESION.stenosis,
        phys = flows.stress.phys,
        to = DEVICE,
    )
    rest = V.Flow.circuit(
        arteries,
        flows.terr;
        state = :rest,
        veins,
        vleaf = flows.vleaf,
        reference = flows.rest,
        kappa = flows.stress.resist.kappa,
        stenosis = LESION.stenosis,
        phys = flows.rest.phys,
        demand = flows.demand,
        to = DEVICE,
    )
    (rest = rest, stress = stress)
end;

# ╔═╡ aed22c26-9c9b-4fb7-b5c7-b63228c511b7
TRANSPORT_LAD90 = transport_cache("transport_lad90", flows_lad90);

# ╔═╡ 311783a3-e1b0-4ccc-b43f-ebb4ce97c47d
let fig = Mke.Figure(size = (1500, 560))
    ax = vessel_axis(
        fig[1, 1];
        title = "Arteries ≥ 0.15 mm, the 90% LAD lesion in gold",
    )
    draw_myocardium_shell!(ax)
    draw_vessels!(
        ax,
        arteries;
        colormap = ARTERY_COLOURS,
        min_radius = 0.15,
    )
    lesion = [
        Mke.Point3f(arteries.a[:, j]) => Mke.Point3f(arteries.b[:, j])
            for j in LESION.path
    ]
    Mke.linesegments!(
        ax,
        lesion;
        color = :gold,
        linewidth = 9,
    )

    panels = (
        ("Healthy · stress MBF", flows.stress, 2),
        ("90% LAD · stress MBF", flows_lad90.stress, 3),
    )
    for (title, circuit, column) in panels
        ax = vessel_axis(fig[1, column]; title = title)
        surface = Mke.scatter!(
            ax,
            MYOCARDIUM_SURFACE.points;
            color = voxel_mbf(circuit)[MYOCARDIUM_SURFACE.voxels],
            colormap = :turbo,
            colorrange = (0, 6),
            markersize = 4,
            strokewidth = 0,
        )
        column == 3 && Mke.Colorbar(fig[1, 4], surface; label = "MBF (mL/min/g)")
    end
    fig
end

# ╔═╡ 6ba85ada-0a3f-4d20-a5c7-555395dadb37
begin
    "Plaque materials, as static labels: fibrous tissue (collagen), a lipid core, calcification (hydroxyapatite)."
    const PLAQUE = Dict(40 => "basis_collagen", 41 => "basis_lipid", 42 => "hydroxyapatite")
    const STENOSIS_SPEC = V.PhantomSpec(
        PHANTOM_SPEC;
        materials = merge(V.XCAT_STATIC, PLAQUE),
    )

    "Two unit vectors perpendicular to `axis` and to each other."
    function perpendiculars(axis)
        a = axis ./ sqrt(sum(abs2, axis))
        e = abs(a[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        u = e .- sum(e .* a) .* a
        u ./= sqrt(sum(abs2, u))
        v = [
            a[2] * u[3] - a[3] * u[2],
            a[3] * u[1] - a[1] * u[3],
            a[1] * u[2] - a[2] * u[1],
        ]
        return u, v
    end

    """
    The plaque at `p` inside XCAT's LAD lumen along the lesion: calcium in a 60° arc at the outer rim, a lipid core behind the
    cap on the opposite side, fibrous tissue elsewhere; `nothing` where `p` is lumen or outside the lesion.
    """
    function plaque(p)
        for (j, r) in zip(LESION.path, LESION.radii)
            a = arteries.a[:, j]
            b = arteries.b[:, j]
            axis = b .- a
            s = clamp(sum((p .- a) .* axis) / sum(abs2, axis), 0, 1)
            offset = p .- (a .+ s .* axis)
            distance = sqrt(sum(abs2, offset))
            (0 < s < 1 && distance <= r + 0.2) || continue

            lumen = r * (1 - DIAMETER_STENOSIS)
            distance <= lumen && return nothing
            u, v = perpendiculars(axis)
            angle = atan(sum(offset .* v), sum(offset .* u))
            # 0 at the lumen, 1 at the vessel wall
            depth = (distance - lumen) / max(r - lumen, 1.0e-6)
            abs(angle) < π / 6 && depth > 0.6 && return 42
            abs(abs(angle) - π) < π / 3 && 0.25 < depth < 0.75 && return 41
            return 40
        end
        return nothing
    end

    "The box (mm) around the lesion's segments and their imaged radius: plaque is only looked for inside it."
    const LESION_BOX = let
        ends = reduce(hcat, [hcat(arteries.a[:, j], arteries.b[:, j]) for j in LESION.path])
        reach = maximum(LESION.radii) + 0.2
        (vec(minimum(ends; dims = 2)) .- reach, vec(maximum(ends; dims = 2)) .+ reach)
    end

    """
    The healthy classifier `healthy`, with plaque where the lesion narrows XCAT's LAD lumen. A factory, like the healthy
    one: the returned closure holds `healthy` and the lesion's box as concrete values, so the billions of voxel calls
    are compiled for them rather than looked up as untyped notebook globals each time.
    """
    function make_plaque_classifier(healthy, box)
        lower, upper = box
        return function (p)
            what = healthy(p)
            if what[1] == V.VOX_VESSEL && what[2] > 0 && all(lower .<= p .<= upper)
                material = plaque(collect(Float64, p))
                material !== nothing && return (V.VOX_STATIC, Int32(material), 0.0f0)
            end
            return what
        end
    end

    classify_with_plaque = make_plaque_classifier(classify, LESION_BOX)
end;

# ╔═╡ 0bf1b83a-bda5-4802-b067-b55197d8784b
STENOSIS_DIR = save_phantom(
    "lad_stenosis_90",
    TRANSPORT_LAD90;
    spec = STENOSIS_SPEC,
    classify = classify_with_plaque,
);

# ╔═╡ c2544222-28dd-4765-9529-ab84a8716d5a
stenosis = inspect(STENOSIS_DIR);

# ╔═╡ b458071e-51db-4dcb-8e46-2f751e9ac1d0
md"""
* state $(@bind stenosis_state Select(["stress", "rest"]))
* time $(@bind stenosis_time Slider(0:1:40; default = 20, show_value = true)) s
* z $(@bind stenosis_z Slider(round.(stenosis.z0 .+ (0:5:(size(stenosis.bundle.labels, 3) - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 280.0, show_value = true)) mm
"""

# ╔═╡ 64d1941b-53a3-4d2d-b16f-7ec998c569a4
slice_figure(stenosis, stenosis_state, stenosis_time, stenosis_z)

# ╔═╡ 09e26b54-5bec-4e4f-8e2c-9f5618bde86c
curves_figure(stenosis)

# ╔═╡ 48ee3e74-16d3-496e-9b9d-1f944cfff23e
outflow_check(stenosis)

# ╔═╡ a82aaa5f-af8d-4c94-9173-936dd0b36c0e
let fig = Mke.Figure(size = (1200, 480))
    for (column, state) in enumerate(("rest", "stress"))
        ax = Mke.Axis(
            fig[1, column];
            title = "Ground-truth MBF, $state",
            ylabel = "mL/min/g",
            xticks = (eachindex(ROI_NAMES), ROI_NAMES),
            xticklabelrotation = π / 5,
        )
        healthy_mbf = [healthy.bundle.roi_table[(state, name)].mbf for name in ROI_NAMES]
        stenosis_mbf = [stenosis.bundle.roi_table[(state, name)].mbf for name in ROI_NAMES]
        Mke.barplot!(
            ax,
            eachindex(ROI_NAMES) .- 0.2,
            healthy_mbf;
            width = 0.38,
            color = :gray60,
            label = "healthy",
        )
        Mke.barplot!(
            ax,
            eachindex(ROI_NAMES) .+ 0.2,
            stenosis_mbf;
            width = 0.38,
            color = :firebrick,
            label = "90% LAD",
        )
        column == 1 && Mke.axislegend(ax; position = :lt)
    end
    fig
end

# ╔═╡ 40fdee14-9c29-464e-a674-04dc149ac5e9
md"""
## 10. Next Steps: Scanning a Phantom with BasisSimulator.jl

A standalone script: the saved phantom scanned at any time `t` on the validated GE Revolution Apex Elite model (BasisSimulator.jl 0.15+), here at the bolus arrival and at the arterial peak. Ground truth for any compartment comes straight from the phantom.

```julia
import Vasculature as V
import BasisSimulator as BS
import GPUSelect

const GPU = GPUSelect.Storage()
const DIR = expanduser("~/xcat-perfused/cardiac_pipeline/healthy") # or lad_stenosis_90
const STATE = "rest" # or "stress"

# ---- the GE Revolution Apex Elite model (MolloiLab/basis-verification, validated against
# the clinical scanner)
const GE_MODEL = (
    source_to_isocenter = 625.6,
    source_to_detector = 1100.0,
    detector_rows = 256,
    detector_cols = 834,
    detector_row_size = 0.625,
    detector_col_size = 0.6,
    detector_col_offset = 0.25,
    detector_shape = :arc,
    target_angle = 10.0,
    flat_filter_material = :aluminum,
    flat_filter_thickness = 2.5,
    bowtie_filter = :ge_revolution_large,
    detector_material = :lumex,
    detector_depth = 3.0,
    fill_factor_row = 0.9,
    fill_factor_col = 0.9,
    electronic_noise = 1500,
    detection_gain = 10.0,
)
# the scanner's large spot at high tube power
focal_spot(mA) = mA >= 250 ? (1.2, 1.6) : (0.7, 1.0)
scanner(mA) = BS.EICTScanner(;
    GE_MODEL...,
    focal_spot_width = focal_spot(mA)[1],
    focal_spot_length = focal_spot(mA)[2],
)
const SIM = BS.SimOptions(
    projector = :dd_fast,
    detector_efficiency_mode = :mc_lut,
    view_samples = 5,
)
with_seed(seed) = BS.SimOptions(;
    (f => getfield(SIM, f) for f in fieldnames(BS.SimOptions) if f !== :seed)...,
    seed,
)

# the STANDARD kernel: an apodized ramp fitted to the clinical one, on a 22 cm cardiac field
# (512², 0.43 mm)
const FOV_CM, SCAN_FOV_CM = 22.0, 30.0
const KNOTS = Tuple(range(0.0, 1.0; length = 11))
const WINDOW = (fc = 0.7 * (10FOV_CM / 512) / (350 / 512), p = 1.5)
const STANDARD = BS.CustomFilter(
    KNOTS,
    Tuple(round(exp(-(x / WINDOW.fc)^WINDOW.p); digits = 5) for x in KNOTS),
)

# ---- one axial rotation of the phantom slab at time t: project, then water BHC, FDK and HU
# on the cardiac field
function scan(bundle, t; mA, collimation_mm, slice_mm, seed, field_centre_mm)
    # ≤ 64 materials, every compartment's iodine kept
    P0, _ = V.DynamicPhantom.bs_phantom(
        bundle,
        BS;
        state = STATE,
        t,
    )
    phantom = BS.Phantom(
        GPU(P0.mask),
        Dict(i - 1 => m for (i, m) in enumerate(P0.materials)),
        P0.voxel_size,
    )
    protocol = BS.CTProtocol(;
        kVp = 100,
        mA,
        views = 984,
        rotation_time = 0.35,
        collimation_mm,
        additional_filters = [("Al", 4.5)],
    )
    sim = with_seed(seed)
    recon = BS.ReconOptions(;
        matrix_size = (512, 512, round(Int, collimation_mm / slice_mm)),
        fov_cm = SCAN_FOV_CM,
        z_cm = collimation_mm / 10,
    )
    ws = BS.create_eict_workspace(scanner(mA), protocol, sim, recon, phantom)
    dose = BS.simulate!(ws, phantom, protocol, sim).dose
    bhc = BS.calibrate_bhc_water(
        sim,
        protocol;
        scanner = scanner(mA),
        geom = ws.geom,
    )
    sinogram = GPU(BS.apply_bhc_water(GPU(ws.sinogram), bhc))
    # the cardiac field: every source and detector position moved by −(field centre −
    # isocentre), exact for FDK
    g = ws.geom
    d = [field_centre_mm[1] / 10, field_centre_mm[2] / 10, 0.0]
    geom = BS.CTGeometry(
        g.SAD,
        g.SDD,
        g.n_angles,
        g.n_rows,
        g.n_cols,
        g.pixel_size,
        g.pixel_row_size,
        g.angles,
        g.source_positions .- d,
        g.detector_centers .- d,
        g.detector_u,
        g.detector_v,
        (FOV_CM, FOV_CM, g.fov[3]),
        g.pitch,
        g.table_feed,
        g.detector_shape,
        g.column_offset,
    )
    wr = BS.create_fdk_recon_workspace(
        sinogram,
        geom,
        recon.matrix_size;
        filter = STANDARD,
    )
    hu = Float32.(
        BS.to_hounsfield(
            Array(BS.reconstruct!(wr, sinogram, geom));
            μ_water = bhc.μ_water_ref,
        ),
    )
    BS.release_backend!(ws)
    BS.release_backend!((wr, sinogram))
    BS.release_backend!(phantom.mask)
    return hu, dose
end

# ---- the protocol. The phantom is centred on the isocentre; positions below are in the
# phantom's mm frame minus its centre.
meta = V.DynamicPhantom.load(DIR; volume = false)
# the tracking point (x, y, z mm)
track = meta.meta["tracking"]["descending aorta"]["point_mm"]
# arrival (+80 HU at the point) and the arterial peak
timing = V.DynamicPhantom.scan_times(meta; state = STATE, roi = "roi:descending aorta")
# two scans: 100 mm over the heart, 0.625 mm slices, 500 mA, at the arrival and at the peak
grid = meta.meta["grid"]
# the phantom's centre (x, y), on the isocentre
centre = grid["origin_mm"][1:2] .+ (grid["dims"][1:2] .- 1) .* grid["voxel_mm"] ./ 2
slices(zlo, zhi) = (
    (round(Int, (zlo - grid["origin_mm"][3]) / grid["voxel_mm"]) + 1):
    round(Int, (zhi - grid["origin_mm"][3]) / grid["voxel_mm"])
)
const HEART_Z = (207.4, 307.4) # mm: the 100 mm beam over the heart
heart = V.DynamicPhantom.load(
    DIR;
    zrange = max(1, first(slices(HEART_Z[1] - 15, 0))):last(slices(0, HEART_Z[2] + 15)),
)
# the 22 cm field on the heart and the descending aorta
field_centre_mm = (151.6, 227.8) .- Tuple(centre)
early, early_dose = scan(
    heart,
    timing.arrival;
    mA = 500,
    collimation_mm = 100.0,
    slice_mm = 0.625,
    seed = 7,
    field_centre_mm,
)
peak, peak_dose = scan(
    heart,
    timing.peak_time;
    mA = 500,
    collimation_mm = 100.0,
    slice_mm = 0.625,
    seed = 8,
    field_centre_mm,
)

# ---- ground truth: any compartment or group is a sum over its labels
# .mbf (mL/min/g), .mass (mg of iodine every second), .mass_g, .voxels
truth = heart.roi_table[(STATE, "LV")]
```

`V.DynamicPhantom.hu(bundle, state, t)` gives each label's idealized 70 keV HU, the scanner-free reference.
"""

# ╔═╡ Cell order:
# ╟─bfa4faa8-3c77-4464-ae5d-87856a7f961b
# ╟─45b99d76-06b3-436a-9149-0bf0fbaa3a45
# ╠═9cff4f49-4723-43e1-9290-be0f428e0352
# ╠═99ec1d1a-81b2-4f1b-b9d9-602949ccee63
# ╠═06f98b38-fbb4-48c9-83c3-a327c6d74052
# ╠═29a716e1-b776-4949-b130-452855d6708a
# ╠═6602ea59-d68c-4c71-8be6-e7c15ddab79e
# ╠═a4268901-4cf8-4181-a7ff-52aa62e27d6c
# ╠═92dde8bc-c022-4a93-bf00-b6e6a3342fea
# ╠═0e0a5727-5a5d-4d7f-b27f-018760e1ba8b
# ╠═fbd087b6-3b8c-4d0a-804c-8b0982c90812
# ╠═a120dd81-758c-4729-af09-19ff607347f2
# ╠═db87cf42-7622-4efc-91e5-3f45da374c65
# ╠═a1035ac8-ddb9-42a2-b0c2-e86eec8f1769
# ╠═72d7878f-c03f-4b24-b5fa-0ef529dbcf68
# ╟─ed124ecb-5d54-4641-b775-0f976fe0ec74
# ╠═7066ffd5-cdf7-4722-968f-d48ffbc2fd0e
# ╠═b5b3a68d-bde2-42e2-8df7-d09d41283847
# ╟─7e5c835a-c16c-4c0f-81ae-5410e13fcad8
# ╠═460614b2-2f7d-4af4-94ad-8051d7c8e70b
# ╠═24dbb378-3f1d-4232-b123-ecb319c51c48
# ╟─94428579-cfcc-4c6d-a766-a1a3d36fede9
# ╟─1b97b382-e826-4989-a7d7-90fd801bcfb6
# ╠═985bdb9d-847b-446c-9e54-b9f9fe0be043
# ╠═c618deff-991a-43f2-a4ca-fbc0b65cf28e
# ╠═0212424e-a0bb-4820-b933-f6f6a870c5b9
# ╠═98c9c13d-a2a5-4316-84bd-715dbd077b88
# ╠═27920760-9210-402b-80d5-d13054040d6b
# ╠═ddfd8b55-868d-4c32-9d96-705b039efd06
# ╠═fcca3b04-13a7-4274-b1fd-186837613b97
# ╠═4af8c280-7734-4920-9683-73ff73ac3ed0
# ╠═28d82d8b-1f94-462f-8b3e-04cb56e127dc
# ╟─60afcebc-cc3d-482d-a811-3e1544584240
# ╟─ed8418d2-529c-4fad-8632-4d2ba5ba16e1
# ╟─b7fd02d1-47f4-4252-aca5-1f4a3b0cdd75
# ╠═ec898038-0e7a-4f98-9737-fd711134217c
# ╠═1e11cee2-aa82-45e0-8965-2c09d62c44ae
# ╠═13aff1b3-5b21-4270-8eaf-2f2ce16376a0
# ╟─9fc11643-e2a0-4e3f-88cb-17f89c5c90e3
# ╟─5da9b3a3-c69d-424b-9dd0-b4ef94fa0398
# ╟─bd4149bd-1699-4f34-845d-5e7272c46571
# ╠═4cb380c8-97b1-4fc4-b374-2832ed3f991c
# ╠═1279319b-88cc-4214-adae-f8fac1701c52
# ╠═e442b143-f7da-4e3b-9151-6f81358788c5
# ╠═52e13f82-7260-41db-8ebf-32149123ac4e
# ╟─ecf9a259-972e-4b13-a1f2-26cab5610bcf
# ╟─68d1c5cd-fa2f-4961-9636-f6731ad142d6
# ╟─e03d6af1-41ef-48c6-b867-4536d946cec4
# ╠═0c4d6e8f-1a3b-4c5d-9e7f-2b4c6d8e0f13
# ╠═3efb6e3b-8c83-43e6-8e17-7152d42f6c2d
# ╠═213ad6d9-995f-4b5c-9d87-549824095172
# ╠═c2d2cdfc-a6dd-4892-be8d-9833d47b56b3
# ╠═5f313e94-82d9-4335-9fae-0d8fcbec0fe8
# ╠═3ddafaeb-6b8c-4dfc-a037-20905115a1e7
# ╟─f5b9f215-e830-4b40-bc5e-d5c3670a6e5d
# ╟─cf6401e0-6a84-4016-b458-7a17408173df
# ╠═105100e4-ca4a-4098-ad30-42f0cb2ad7f0
# ╠═1488815e-8708-46f2-ade5-31b7438b0a41
# ╠═2235d2e6-01bc-43ac-b6a4-1a8b5ece93c8
# ╠═a162fdf2-c33e-4724-8874-b0ab23c96cae
# ╠═1cb4005e-52ef-4b66-b4ea-02d7e982fc0c
# ╠═a648cc82-ba31-49ce-8c88-b670dbf3a513
# ╠═ba8d75ec-ca1a-4952-970b-a31cdd5be016
# ╠═48514abd-afbb-44c0-8f3d-f3293da05204
# ╠═a19d475f-edfe-41c7-ad61-2ea8bdf170bd
# ╟─44638d83-bacc-489c-bbc1-f28ce8ff43ac
# ╠═ec88e71f-17a6-42ab-ae1d-131b91a396df
# ╠═5b0c2f7e-3d41-4a8e-9b6f-2c7e1d9a4f10
# ╠═2f6d8c41-7a3e-4b9c-8e15-d4a7c0b93e62
# ╠═0fbff6a0-6c59-4259-9974-b51529affa05
# ╠═23550e90-a7a1-4923-b7ca-f0910276f80a
# ╠═b7f0823b-85d9-4ed5-8096-9818ec1877f8
# ╠═2fe6d117-a55a-412c-bf8b-10cd9444056e
# ╠═487f30c5-950f-4177-b5b8-4c50cf4ca156
# ╟─ef9926c5-a216-4323-abe3-80796ebab79a
# ╠═95d13834-2f64-4318-a83b-7a447e8f9e5b
# ╠═2b074062-ec5b-42d2-9b94-95329585ea3d
# ╠═3fc08b34-bd48-4687-98bb-efdfb79bd05e
# ╟─53be3442-4eb3-4547-95b6-5fc32f670646
# ╟─e44cab6e-9bc9-44e7-a481-08938d5fac32
# ╠═3c468348-759d-463a-bd1c-9ded251211cc
# ╠═ab9d2337-2a36-4416-97ad-e9ad2c753d4b
# ╟─fa0357e5-2d1e-4354-a12d-7e1c8c308c3a
# ╟─173bb953-3b21-4c4a-aa00-9c838e0b21e1
# ╟─948065f2-a977-4980-9eba-dec19b754a5a
# ╠═872d6939-4654-4b12-82be-58aed8a9f85a
# ╠═168b517e-ea3e-4933-b0b3-1f85f31a8512
# ╠═cbf27751-a278-460d-9b99-5db1be2e9bf1
# ╠═aed22c26-9c9b-4fb7-b5c7-b63228c511b7
# ╟─311783a3-e1b0-4ccc-b43f-ebb4ce97c47d
# ╠═6ba85ada-0a3f-4d20-a5c7-555395dadb37
# ╠═0bf1b83a-bda5-4802-b067-b55197d8784b
# ╠═c2544222-28dd-4765-9529-ab84a8716d5a
# ╟─b458071e-51db-4dcb-8e46-2f751e9ac1d0
# ╟─64d1941b-53a3-4d2d-b16f-7ec998c569a4
# ╟─09e26b54-5bec-4e4f-8e2c-9f5618bde86c
# ╟─48ee3e74-16d3-496e-9b9d-1f944cfff23e
# ╟─a82aaa5f-af8d-4c94-9173-936dd0b36c0e
# ╟─40fdee14-9c29-464e-a674-04dc149ac5e9
