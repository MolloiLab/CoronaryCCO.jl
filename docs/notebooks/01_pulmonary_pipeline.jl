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

# ╔═╡ 53ef55e4-f46a-45b5-8f84-ecc59711a13d
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    # a writable font cache for CairoMakie (a notebook worker may have no HOME)
    ENV["XDG_CACHE_HOME"] = get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))
end;

# ╔═╡ ef4704c9-3438-48e1-9f42-2f97b9f28acc
using PlutoUI: Slider, Select, TableOfContents

# ╔═╡ 4152495c-72ba-4321-8ee2-4e4754ff2797
using Printf: @sprintf

# ╔═╡ a9210218-f94e-45a4-bde7-c4f548a4eae6
using Serialization: serialize, deserialize

# ╔═╡ 07c920f4-5c4a-4653-af94-284f616f8256
using Mmap: mmap

# ╔═╡ cf8017b9-c58c-465b-8945-a13acd593609
using Statistics: median, quantile

# ╔═╡ d0311ce2-5aca-446b-bf46-44ea40dd8921
md"""
# Pulmonary Pipeline

**From the raw XCAT chest to a perfused lung phantom, and the same lungs with a left lower lobe artery occluded.**

| Stage | What happens |
|---|---|
| **1. Load** | XCAT's 0.2 mm label volume, its NURBS heart and its lung and vessel meshes |
| **2. Prepare** | the imaged pulmonary arteries and veins, the airways and lobes, the lungs to grow into |
| **3. Grow** | one arterial terminal and one venule per acinus, 29,912 of each |
| **4. Flow** | blood flow from the pulmonary arterial and left atrial pressures alone |
| **5. Transport** | a clinical contrast injection through the whole body |
| **6. Compartments** | the lobes ground truth is given for, and their flow |
| **7. Voxelize and Save** | the 0.2 mm perfused phantom |
| **8. Inspect** | the saved phantom: labels, compartments, idealized images, ground truth |
| **9. Occlusion** | a left lower lobe artery blocked, as its own phantom |
| **10. Next Steps** | scanning the phantom with BasisSimulator.jl |

Every heavy step is cached in `CACHE_DIR`, so a rerun takes seconds. The XCAT phantom is licensed: its files and anything made from them stay on lab storage.
"""

# ╔═╡ 21cabb67-46bc-4458-b215-d2484b257de6
md"""
## Setup
"""

# ╔═╡ 90ea862c-5c67-49ab-920a-d5acb3244766
import TOML

# ╔═╡ 5fd8df67-aca2-4daa-ab7d-3a3620dff99a
import Vasculature as V

# ╔═╡ 5f1b74ca-efbd-4370-bda2-2a499cf7be7c
# ╠═╡ show_logs = false
import CairoMakie as Mke

# ╔═╡ 58d7db22-8641-41c2-aed3-1a23e3bd6883
import GPUSelect

# ╔═╡ 2e61946c-82b1-402f-8ef7-e26f22dbd5b2
TableOfContents()

# ╔═╡ bdf025f0-bc0b-436f-9581-23318ebe26fa
begin
    # licensed: these files and anything made from them stay on lab storage
    const XCAT_DIR = get(
        ENV,
        "XCAT_DIR",
        joinpath(homedir(), "Documents", "dev", "MolloiLab", "xcat-data", "vmale_50_heart_1600_1400_867"),
    )
    # the full chest, arms included, at the same 0.2 mm: XCAT's heart block inside it and its surface meshes voxelized
    # around it (built once by xcat-data/work/phantom/build_chest.jl), so a scanner sees the whole body
    const XCAT_CHEST_DIR = get(
        ENV,
        "XCAT_CHEST_DIR",
        joinpath(homedir(), "Documents", "dev", "MolloiLab", "xcat-data", "work", "phantom"),
    )
    const CACHE_DIR = joinpath(homedir(), "xcat-perfused", "pulmonary_pipeline_cache")
    const PHANTOM_DIR = joinpath(homedir(), "xcat-perfused", "pulmonary_pipeline")

    mkpath(CACHE_DIR)
    mkpath(PHANTOM_DIR)
end;

# ╔═╡ 7d6dc0e9-93bb-4e12-96f9-21efbcbf3d38
# ╠═╡ show_logs = false
"Arrays to this machine's GPU (CUDA, Metal, ROCm, oneAPI), or the CPU if there is none."
const DEVICE = GPUSelect.Storage()

# ╔═╡ c69a5e75-f4ef-424f-8876-f26ad3501180
"Run `compute` once and keep its result in `CACHE_DIR/name.jls`; afterwards load it."
function cached(compute, name)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) && return deserialize(path)

    seconds = @elapsed result = compute()
    serialize(path, result)
    @info "$name: $(round(seconds; digits = 1)) s"
    return result
end

# ╔═╡ 440ec9de-6961-4ed1-9bf4-43e0dbfea35c
"The six face neighbours of a voxel."
const SIX_STEPS = (
    CartesianIndex(1, 0, 0),
    CartesianIndex(-1, 0, 0),
    CartesianIndex(0, 1, 0),
    CartesianIndex(0, -1, 0),
    CartesianIndex(0, 0, 1),
    CartesianIndex(0, 0, -1),
)

# ╔═╡ dddaae8d-18cb-4213-949c-f72ea7d8d717
md"""
## 1. Load the Raw XCAT Chest

XCAT vmale\_50's whole chest, arms included: one `UInt8` label per 0.2 mm voxel (XCAT's heart block inside it), the NURBS surfaces of the heart (the main pulmonary artery and the left atrium among them), and XCAT's triangle meshes of the two lungs and of the vessels inside them: the pulmonary arteries and veins, each drawn as a few hundred short tubes. XCAT also labels the airways, their walls and the air inside them. The phantom is saved on the full chest so a CT scan sees the whole body.
"""

# ╔═╡ 26bbaa0c-82be-415f-a934-7e86734b719d
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

    # XCAT's surface meshes lie this far (mm) from its NURBS and label frame: the shift the chest volume was built with
    # (its lungs land on XCAT's own lung labels, section 2)
    const MESH_OFFSET = (131.4, 238.4, 193.8)
    const ARTERY_TUBES = ("right_pulmonary_arts", "left_pulmonary_arts")
    const VEIN_TUBES = ("right_pulmonary_veins", "left_pulmonary_veins")
    lung_meshes = V.XCAT.read_mesh_blocks(
        joinpath(XCAT_DIR, "vmale_50_heart_1_organs.raw");
        names = ("rlung", "llung"),
        offset = MESH_OFFSET,
    )
    vessel_tubes = V.XCAT.read_mesh_blocks(
        joinpath(XCAT_DIR, "vmale_50_heart_1_vessels.raw");
        names = (ARTERY_TUBES..., VEIN_TUBES...),
        offset = MESH_OFFSET,
    )
end;

# ╔═╡ b3aecd39-f506-4538-8fcd-4bdcb22f595c
begin
    "XCAT label groups for pictures: (name, colour, XCAT labels)."
    const LABEL_GROUPS = [
        ("lung", :gray40, (6,)),
        ("airways", :paleturquoise, (7, 11)),
        ("arteries", :red, (25, 28, 125)),
        ("veins", :royalblue, (24, 124)),
        ("heart blood", :mediumpurple, (19, 20, 21, 22)),
        ("myocardium", :darkorange, (15, 16, 17, 18)),
        ("fat", :khaki, (1,)),
        ("muscle", :rosybrown, (4, 13)),
        ("bone", :ivory, (5, 8, 9)),
        ("other", :gray65, ()),
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
    xcat_slice(z) = clamp(round(Int, (z - XCAT_ORIGIN[3]) / XCAT_VOXEL_MM) + 1, 1, XCAT_DIMS[3])
    xcat_x_mm = XCAT_ORIGIN[1] .+ (0:(XCAT_DIMS[1] - 1)) .* XCAT_VOXEL_MM
    xcat_y_mm = XCAT_ORIGIN[2] .+ (0:(XCAT_DIMS[2] - 1)) .* XCAT_VOXEL_MM

    "The XCAT label at `p` (mm); 0 outside the chest."
    function xcat_label(p)
        k = round.(Int, (p .- XCAT_ORIGIN) ./ XCAT_VOXEL_MM) .+ 1
        return all(1 .<= k .<= XCAT_DIMS) ? xcat_labels[k...] : 0x00
    end
end;

# ╔═╡ 2fe11641-05cb-4803-9c8d-a30328bc6235
md"""
z (mm): $(@bind z_xcat Slider(round.(XCAT_ORIGIN[3] .+ (0:5:(XCAT_DIMS[3] - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 290.0, show_value = true))
"""

# ╔═╡ c5f15420-4d45-4054-873f-282a655f9442
begin
    # shades for vessels by radius, pale for the thinnest shown and deep for the largest:
    # arteries red, veins blue
    ARTERY_COLOURS = Mke.cgrad(["#f4a582", "#d6604d", "#b2182b", "#67001f"])
    VEIN_COLOURS = Mke.cgrad(["#92c5de", "#4393c3", "#2166ac", "#053061"])
end;

# ╔═╡ aeb894e0-eea5-416d-a0c1-04f02b797b7d
"""
A 3-D axis for the chest in true proportions, seen from the front: head up, the patient's left on the right. A gentle
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
        viewmode = :fitzoom, # zoomed to the data, not to the box it could rotate in
        protrusions = 0,
    )
    Mke.hidedecorations!(ax)
    Mke.hidespines!(ax)
    return ax
end

# ╔═╡ 44aca054-f174-4b36-aed6-30b5fd4fb6df
let k = xcat_slice(z_xcat)
    every = 2
    groups = [GROUP_OF_LABEL[xcat_labels[i, j, k] + 1] for i in 1:every:XCAT_DIMS[1], j in 1:every:XCAT_DIMS[2]]

    fig = Mke.Figure(size = (1200, 720))
    ax = Mke.Axis(fig[1, 1]; title = "XCAT labels · z = $(z_xcat) mm", aspect = Mke.DataAspect(), yreversed = true)
    Mke.heatmap!(
        ax,
        xcat_x_mm[1:every:end],
        xcat_y_mm[1:every:end],
        groups;
        colormap = label_group_colormap(),
        colorrange = (-0.5, length(LABEL_GROUPS) + 0.5),
    )
    Mke.hidedecorations!(ax)

    ax3 = vessel_axis(fig[1, 2]; title = "XCAT's lung vessel tubes and the main pulmonary artery, seen from the front")
    # every tube of a system in one call, broken by NaN between tubes (one call per tube is slow)
    for (names, colour) in ((ARTERY_TUBES, :red), (VEIN_TUBES, :royalblue))
        points = Mke.Point3f[]
        for name in names, (vertices, _) in vessel_tubes[name]
            line = V.XCAT.tube_centerline(vertices; step = 3.0)
            append!(points, [Mke.Point3f(line.X[:, k]) for k in axes(line.X, 2)])
            push!(points, Mke.Point3f(NaN, NaN, NaN))
        end
        Mke.lines!(ax3, points; color = colour, linewidth = 1.5)
    end
    trunk = V.XCAT.centerline(xcat_surfaces["dias_pul_art"]; step = 3.0)
    Mke.lines!(ax3, trunk.X[1, :], trunk.X[2, :], trunk.X[3, :]; color = :darkred, linewidth = 6)

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

# ╔═╡ 4d68e042-f2da-44b9-88f2-d86ed350d6c7
md"""
## 2. Prepare the Growth Inputs

Growth runs on a 1 mm grid around both lungs, into the lung voxels.

- **Inlet trees.** The main pulmonary artery from its valve, with XCAT's lung artery tubes joined to it, is one arterial tree. Each lung's vein tubes, joined and carried to the left atrium, are one venous tree. Every imaged vessel gets **stubs**, short side branches one vessel segment apart, so growth starts along its whole length, as real pulmonary vessels branch along theirs. The segment length is human pulmonary morphometry's for vessels of the imaged ones' calibre (Huang et al. 1996).
- **Lobes.** XCAT draws no fissures, so the lobes come from the airways, as a lobe is assigned without fissures: the lobar bronchi are found from the airway tree's branching, and each lung voxel takes the lobe of the nearest lobar airway, by the distance through the lung.
- **Blockers.** No new vessel may cross the airways, the heart's chambers, the aorta, the other system's imaged lumens, or a **fissure**: the boundary between two lobes, away from the hilum where the lobar vessels and bronchi enter.
"""

# ╔═╡ f2c25bd3-fd90-4742-9996-067b2f16b04c
growth_grid = let points = hcat(first(only(lung_meshes["rlung"])), first(only(lung_meshes["llung"])))
    lower = floor.(Int, vec(minimum(points; dims = 2))) .- 2
    upper = ceil.(Int, vec(maximum(points; dims = 2))) .+ 2
    V.VoxelGrid(Tuple(Float64.(lower)), 1.0, Tuple(upper .- lower .+ 1))
end

# ╔═╡ abb33992-3cb7-41ff-859f-f55c8d6d9d3b
"The lung: the 1 mm cells that are mostly XCAT's lung (label 6), the voxels the phantom paints as lung."
LUNG = cached("lung") do
    V.XCAT.label_mask(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, growth_grid, (0x06,))
end;

# ╔═╡ 8cbd6bf0-9455-4aa5-8fde-d21cea7b7111
let meshes = V.voxelize(only(lung_meshes["rlung"])..., growth_grid) .| V.voxelize(only(lung_meshes["llung"])..., growth_grid)
    overlap = 2 * count(LUNG .& meshes) / (count(LUNG) + count(meshes))
    md"""
    The lung: $(round(count(LUNG) / 1000; digits = 1)) mL (XCAT's lung meshes: $(round(count(meshes) / 1000; digits = 1)) mL; Dice $(round(overlap; digits = 3)) between the two, so the meshes sit on the labels).
    """
end

# ╔═╡ 27f9d235-e736-4ce8-b7f5-68e20ca83d73
"""
Stub spacing and length, from human pulmonary morphometry (Huang et al. 1996, `test/reference`): the segment length (between
two branchings) of a vessel of the imaged vessels' length-weighted median diameter, and of one two orders smaller, the
typical side branch, as the stub's length.
"""
STUB_LENGTHS = let
    huang = TOML.parsefile(
        joinpath(pkgdir(V), "test", "reference", "huang1996_human_pulmonary_morphometry.toml"),
    )["data"]
    # segment length (mm) at a diameter (mm) `below` orders smaller, log-linear between orders
    function segment_length(table, diameter; below = 0)
        d = log.(table["diameter_mean_mm"])
        l = log.(table["length_mean_mm"])
        x = clamp(log(diameter), d[1], d[end])
        k = clamp(searchsortedlast(d, x), 1, length(d) - 1)
        order = clamp(k + (x - d[k]) / (d[k + 1] - d[k]) - below, 1, length(d))
        j = clamp(floor(Int, order), 1, length(d) - 1)
        w = order - j
        return exp((1 - w) * l[j] + w * l[j + 1])
    end
    function median_diameter(names)
        lines = [V.XCAT.tube_centerline(points) for name in names for (points, _) in vessel_tubes[name]]
        lengths = [sum(sqrt.(sum(abs2, diff(line.X; dims = 2); dims = 1))) for line in lines]
        diameters = [2 * line.r[1] for line in lines]
        order = sortperm(diameters)
        covered = cumsum(lengths[order])
        return diameters[order][findfirst(>=(covered[end] / 2), covered)]
    end
    artery = median_diameter(ARTERY_TUBES)
    vein = median_diameter(VEIN_TUBES)
    (
        artery = (
            diameter = artery,
            spacing = segment_length(huang["arterial_segments"], artery),
            depth = segment_length(huang["arterial_segments"], artery; below = 2),
        ),
        vein = (
            diameter = vein,
            spacing = segment_length(huang["venous_segments"], vein),
            depth = segment_length(huang["venous_segments"], vein; below = 2),
        ),
    )
end

# ╔═╡ c4b8f772-02aa-4cac-9a1b-d4d4540585f1
"""
The arterial inlet tree: XCAT's main pulmonary artery (its NURBS surface, from the valve), with the lung artery tubes joined
into trees and grafted onto it (`V.XCAT.graft`). Centrelines are sampled at the stub spacing.
"""
artery_base = cached("artery_base") do
    step = STUB_LENGTHS.artery.spacing
    trunk = V.XCAT.centerline(xcat_surfaces["dias_pul_art"]; step)
    hub = vec(sum(trunk.X; dims = 2)) ./ size(trunk.X, 2)
    lines = [V.XCAT.tube_centerline(points; name, step) for name in ARTERY_TUBES for (points, _) in vessel_tubes[name]]
    V.XCAT.graft(trunk, V.XCAT.tube_forest(lines; toward = hub); name = "pulmonary_arteries")
end;

# ╔═╡ 03d82855-b9d6-43e5-96a2-f869776c5463
"""
The venous inlet trees, one per lung: its vein tubes joined into trees, the largest taking the others, carried straight to
the nearest blood of the left atrium (XCAT's label, so the veins end where the atrium's blood begins; `V.XCAT.drain`). The
ostium's radius is the lung's two pulmonary veins (Huang 1996 venous order 15, 12.97 mm) as one cube-law equivalent.
"""
vein_bases = cached("vein_bases") do
    step = STUB_LENGTHS.vein.spacing
    points, _ = V.XCAT.tessellate(xcat_surfaces["dias_la"]; nu = 60, nv = 60)
    grid = V.grid_around(points, 1.0; margin = 10)
    atrium = V.XCAT.label_mask(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, grid, (0x15,))
    wall = [c for c in findall(atrium) if any(s -> !checkbounds(Bool, atrium, c + s) || !atrium[c + s], SIX_STEPS)]
    surface = reduce(hcat, [collect(V.center(grid, Tuple(c)...)) for c in wall])
    hub = vec(sum(surface; dims = 2)) ./ size(surface, 2)
    [
        V.XCAT.drain(
            V.XCAT.tube_forest(
                [V.XCAT.tube_centerline(points; name, step) for (points, _) in vessel_tubes[name]];
                toward = hub,
            ),
            surface;
            root_radius = cbrt(2 * (12.97 / 2)^3),
            name,
        ) for name in VEIN_TUBES
    ]
end;

# ╔═╡ 1f2278a1-9326-417c-b472-405539f41c50
"""
The airway tree, from XCAT's airway labels (wall 7 and air 11) on a 0.6 mm grid: its distance from the trachea's top through
the airway, cut into 3 mm shells; each shell's connected pieces are the tree's nodes, each joined to the piece of the shell
before it that it touches. Returns the grid, the airway mask, each node's cells, parent and shell, and each cell's node.
"""
function airway_tree(; cell = 0.6, shell_mm = 3.0)
    grid = V.VoxelGrid(growth_grid.origin, cell, Tuple(ceil.(Int, growth_grid.dims ./ cell)))
    airway = V.XCAT.label_mask(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, grid, (0x07, 0x0b); fraction = 1.0e-9)
    steps = [CartesianIndex(i, j, k) for i in -1:1, j in -1:1, k in -1:1 if (i, j, k) != (0, 0, 0)]
    # keep the largest connected piece
    piece = zeros(Int32, size(airway))
    sizes = Int[]
    for seed in findall(airway)
        piece[seed] == 0 || continue
        push!(sizes, 0)
        stack = [seed]
        piece[seed] = length(sizes)
        while !isempty(stack)
            c = pop!(stack)
            sizes[end] += 1
            for s in steps
                n = c + s
                checkbounds(Bool, airway, n) && airway[n] && piece[n] == 0 || continue
                piece[n] = length(sizes)
                push!(stack, n)
            end
        end
    end
    airway .= piece .== argmax(sizes)
    # distance (cells) from the trachea's top, through the airway
    top = maximum(c[3] for c in findall(airway))
    distance = fill(Int32(-1), size(airway))
    queue = [c for c in findall(airway) if c[3] >= top - 1]
    foreach(c -> distance[c] = 0, queue)
    head = 1
    while head <= length(queue)
        c = queue[head]
        head += 1
        for s in steps
            n = c + s
            checkbounds(Bool, airway, n) && airway[n] && distance[n] < 0 || continue
            distance[n] = distance[c] + 1
            push!(queue, n)
        end
    end
    shell_of(c) = floor(Int, distance[c] * cell / shell_mm)
    # nodes: the connected pieces of each shell
    node = zeros(Int32, size(airway))
    members = Vector{Vector{CartesianIndex{3}}}()
    shell = Int[]
    for seed in queue
        node[seed] == 0 || continue
        push!(members, CartesianIndex{3}[])
        push!(shell, shell_of(seed))
        stack = [seed]
        node[seed] = length(members)
        while !isempty(stack)
            c = pop!(stack)
            push!(members[end], c)
            for s in steps
                n = c + s
                checkbounds(Bool, airway, n) && airway[n] && node[n] == 0 && shell_of(n) == shell[end] || continue
                node[n] = length(members)
                push!(stack, n)
            end
        end
    end
    parent = zeros(Int, length(members))
    for (id, cells) in enumerate(members)
        for c in cells
            k = findfirst(s -> checkbounds(Bool, airway, c + s) && airway[c + s] && shell[node[c + s]] == shell[id] - 1, steps)
            k === nothing || (parent[id] = node[c + steps[k]]; break)
        end
    end
    return (; grid, airway, members, parent, shell)
end

# ╔═╡ 6de4ed07-39b0-440a-b56f-c3ca02ea2d88
AIRWAYS = cached("airways") do
    airway_tree()
end;

# ╔═╡ e79d69d3-719f-429f-9346-4b5448003243
"""
The five lobar bronchi, as nodes of the airway tree. Walking down from the trachea, a split is the first node whose second
largest branch holds at least `share` of the airway below it. The carina's two branches are the main bronchi (the right at
smaller x). The right main bronchus's first split gives off the upper lobe bronchus, the branch that goes up; the bronchus
intermedius then splits into the middle lobe bronchus (the anterior one) and the lower lobe bronchus. The left main bronchus
splits into the upper and lower lobe bronchi.
"""
LOBAR_BRONCHI = let tree = AIRWAYS, n = length(tree.members)
    children = [Int[] for _ in 1:n]
    foreach(i -> tree.parent[i] > 0 && push!(children[tree.parent[i]], i), 1:n)
    # airway cells below each node, and their summed positions
    below = Float64.(length.(tree.members))
    position = zeros(3, n)
    for i in 1:n
        position[:, i] = sum(collect(Tuple(c)) for c in tree.members[i])
    end
    for i in sortperm(tree.shell; rev = true)
        p = tree.parent[i]
        p > 0 || continue
        below[p] += below[i]
        position[:, p] .+= position[:, i]
    end
    centre(i) = collect(tree.grid.origin) .+ (position[:, i] ./ below[i] .- 1) .* tree.grid.h
    function first_split(i, share)
        while true
            branches = sort(children[i]; by = c -> below[c], rev = true)
            isempty(branches) && error("no airway split below node $i")
            length(branches) >= 2 && below[branches[2]] >= share * below[i] && return branches[1:2]
            i = branches[1]
        end
        return
    end
    # a split's two branches ordered by one coordinate of their centres (1 x, 2 y front to back, 3 z up)
    ordered(pair, axis; rev = false) = pair[sortperm([centre(i)[axis] for i in pair]; rev)]
    trachea = argmax([tree.parent[i] == 0 ? below[i] : 0.0 for i in 1:n])
    right, left = ordered(first_split(trachea, 0.25), 1)
    upper, intermedius = ordered(first_split(right, 0.15), 3; rev = true)
    middle, lower = ordered(first_split(intermedius, 0.15), 2)
    left_upper, left_lower = ordered(first_split(left, 0.15), 3; rev = true)
    (
        names = ["right upper", "right middle", "right lower", "left upper", "left lower"],
        nodes = [upper, middle, lower, left_upper, left_lower],
        children = children,
    )
end;

# ╔═╡ 6e84d94d-7e21-494b-8cab-e0cb9119b4de
"""
Each 1 mm lung voxel's lobe (1-5, `LOBAR_BRONCHI.names`; 0 outside the lung): the lobe of the nearest lobar airway, by the
distance through the lung (6-connected, all lobes spreading at once from their own airways).
"""
LOBES = cached("lobes") do
    lobe = zeros(Int8, growth_grid.dims)
    for (k, node) in enumerate(LOBAR_BRONCHI.nodes)
        stack = [node]
        while !isempty(stack)
            i = pop!(stack)
            append!(stack, LOBAR_BRONCHI.children[i])
            for c in AIRWAYS.members[i]
                p = AIRWAYS.grid.origin .+ (Tuple(c) .- 1) .* AIRWAYS.grid.h
                q = round.(Int, (p .- growth_grid.origin) ./ growth_grid.h) .+ 1
                all(1 .<= q .<= growth_grid.dims) && (lobe[q...] = k)
            end
        end
    end
    queue = findall(!=(0), lobe)
    head = 1
    while head <= length(queue)
        c = queue[head]
        head += 1
        for s in SIX_STEPS
            n = c + s
            checkbounds(Bool, lobe, n) && lobe[n] == 0 && LUNG[n] || continue
            lobe[n] = lobe[c]
            push!(queue, n)
        end
    end
    lobe .* LUNG
end;

# ╔═╡ 95aca34c-b6bd-442c-a95c-8a2cee41971e
"""
    fissures(lobes, airways, grid; hilum_mm = 5.0) -> BitArray

The fissures on `grid`: lung voxels of `lobes` beside a voxel of another lobe, except within `hilum_mm` of an airway
(`airways`, as `airway_tree` returns it), where the lobar vessels and bronchi enter the lobes.
"""
function fissures(lobes, airways, grid; hilum_mm = 5.0)
    near_airway = falses(size(lobes))
    r = ceil(Int, hilum_mm)
    for c in findall(airways.airway)
        p = airways.grid.origin .+ (Tuple(c) .- 1) .* airways.grid.h
        q = round.(Int, (p .- grid.origin) ./ grid.h) .+ 1
        for k in -r:r, j in -r:r, i in -r:r
            i^2 + j^2 + k^2 <= hilum_mm^2 || continue
            n = q .+ (i, j, k)
            all(1 .<= n .<= size(lobes)) && (near_airway[n...] = true)
        end
    end
    fissure = falses(size(lobes))
    for c in findall(!=(0), lobes)
        near_airway[c] && continue
        other(s) = checkbounds(Bool, lobes, c + s) && lobes[c + s] != 0 && lobes[c + s] != lobes[c]
        any(other, SIX_STEPS) && (fissure[c] = true)
    end
    return fissure
end

# ╔═╡ 80f43c52-f003-4912-be1f-c485974abad6
"The labels no new vessel may cross, for both trees: the airways (wall 7, air 11), the heart's chambers and the aorta."
const AIRWAYS_AND_HEART = (0x07, 0x0b, 0x13, 0x14, 0x15, 0x16, 0x1c)

# ╔═╡ 5a8c1f3e-2d4b-4e6a-9f7c-1b3d5e7a9c02
"A tree's blockers `(mask, origin, cell)` with every cell whose centre lies in a fissure voxel (on `grid`) added."
function with_fissures((mask, origin, cell), fissure, grid)
    for c in CartesianIndices(mask)
        p = origin .+ (Tuple(c) .- 1) .* cell
        k = round.(Int, (p .- grid.origin) ./ grid.h) .+ 1
        all(1 .<= k .<= size(fissure)) && fissure[k...] && (mask[c] = 0x01)
    end
    return (mask, origin, cell)
end

# ╔═╡ 9def37df-4fe1-4716-bb70-1845a0433a67
"""
The blockers, at 0.4 mm, one set per tree: `AIRWAYS_AND_HEART`, the other system's imaged lumens (veins 24 and 124 for the
arteries, arteries 25 and 125 for the veins) and the fissures.
"""
BLOCKERS = cached("blockers") do
    fissure = fissures(LOBES, AIRWAYS, growth_grid)
    blocked(labels) = with_fissures(
        V.XCAT.label_blockers(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, growth_grid, labels),
        fissure,
        growth_grid,
    )
    (
        artery = blocked((AIRWAYS_AND_HEART..., 0x18, 0x7c)),
        vein = blocked((AIRWAYS_AND_HEART..., 0x19, 0x7d)),
        fissure = fissure,
    )
end;

# ╔═╡ d4b22c25-1df4-4b6f-bf71-188d026fbfec
"The 1 mm cells no terminal or stub may end in: mostly airway, chamber or aortic blood."
FORBID = V.XCAT.label_mask(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, growth_grid, AIRWAYS_AND_HEART);

# ╔═╡ 7fd7dcf9-fe92-4d4b-86e1-0ca90317628b
"""
An inlet tree with its stubs (`V.XCAT.penetrating_stubs` at the morphometric spacing; where a vessel runs through the lung,
as almost all of them do, the stub leaves sideways) and its pass-through voxels: its lumens outside the lung, which new
vessels may cross but no terminal may end in.
"""
function inlet(tree, lengths, blockers)
    stubbed, report = V.XCAT.penetrating_stubs(
        tree,
        growth_grid,
        LUNG;
        spacing = lengths.spacing,
        depth = lengths.depth,
        reach = lengths.spacing,
        domain = LUNG,
        forbid = FORBID,
        blockers,
    )
    return (tree = stubbed, report = report, pass = V.XCAT.vessel_tube(stubbed, growth_grid) .& .!FORBID .& .!LUNG)
end

# ╔═╡ f1bd3790-c7e3-43d4-9998-aa4bdb669a40
artery_inputs = cached("artery_inputs") do
    inlet(artery_base, STUB_LENGTHS.artery, BLOCKERS.artery)
end;

# ╔═╡ da727f7d-e412-46a1-9154-823facbf1f98
vein_inputs = cached("vein_inputs") do
    [inlet(tree, STUB_LENGTHS.vein, BLOCKERS.vein) for tree in vein_bases]
end;

# ╔═╡ a362e6a5-c0e2-412b-b96d-9c13609d4d73
"Lobe colours: right lobes warm, left lobes cool."
const LOBE_COLOURS = ["#e41a1c", "#ff7f00", "#ffd92f", "#377eb8", "#984ea3"]

# ╔═╡ cfdb51b7-71df-4da8-a11d-41b8ee1964fc
md"""
y (mm, front to back): $(@bind y_domain Slider(growth_grid.origin[2]:2:(growth_grid.origin[2] + growth_grid.dims[2] - 1); default = 262, show_value = true))
"""

# ╔═╡ 0fb6d8c3-c1c2-4f6d-96c3-1a47a00458a5
let j = clamp(round(Int, y_domain - growth_grid.origin[2]) + 1, 1, growth_grid.dims[2])
    function kind(i, k)
        artery_inputs.pass[i, j, k] && return 6
        any(v -> v.pass[i, j, k], vein_inputs) && return 7
        BLOCKERS.fissure[i, j, k] && return 8
        LOBES[i, j, k] > 0 && return Int(LOBES[i, j, k])
        FORBID[i, j, k] && return 9
        return 0
    end
    kinds = [kind(i, k) for i in 1:growth_grid.dims[1], k in 1:growth_grid.dims[3]]
    xs = growth_grid.origin[1] .+ (0:(growth_grid.dims[1] - 1))
    zs = growth_grid.origin[3] .+ (0:(growth_grid.dims[3] - 1))
    colours = [LOBE_COLOURS; :red; :deepskyblue; :white; :gray50]

    fig = Mke.Figure(size = (1200, 680))
    ax = Mke.Axis(fig[1, 1]; title = "Growth domain · coronal y = $(y_domain) mm", aspect = Mke.DataAspect())
    Mke.heatmap!(ax, xs, zs, kinds; colormap = Mke.cgrad([:black; colours], 10; categorical = true), colorrange = (-0.5, 9.5))
    Mke.hidedecorations!(ax)

    stubs = sum(r.stubs for r in [artery_inputs.report; [v.report for v in vein_inputs]])
    ax3 = vessel_axis(fig[1, 2]; title = "Inlet trees: imaged (thick), $stubs stubs (thin)")
    # each tree's imaged segments in one call and its stubs in another (one call per segment is far too slow)
    for (tree, colour) in [(artery_inputs.tree, :red); [(v.tree, :deepskyblue) for v in vein_inputs]]
        segment(k) = Mke.Point3f(tree.X[:, tree.parent[k]]) => Mke.Point3f(tree.X[:, k])
        stub = [tree.r[k] == 0.1 for k in eachindex(tree.parent)]
        Mke.linesegments!(ax3, [segment(k) for k in 2:length(tree.parent) if !stub[k]]; color = colour, linewidth = 2)
        Mke.linesegments!(ax3, [segment(k) for k in 2:length(tree.parent) if stub[k]]; color = :black, linewidth = 0.5)
    end
    Mke.Legend(
        fig[2, 1:2],
        [Mke.PolyElement(color = c, strokecolor = :gray40, strokewidth = 1) for c in colours],
        [LOBAR_BRONCHI.names; "artery pass-through"; "vein pass-through"; "fissure"; "airway or heart (blocked)"];
        orientation = :horizontal,
        nbanks = 2,
        framevisible = false,
    )
    fig
end

# ╔═╡ 8960d456-a118-462d-be93-db4b884d85d6
let names = LOBAR_BRONCHI.names, total = count(LUNG)
    rows = join(["| $(names[k]) | $(round(100 * count(==(k), LOBES) / total; digits = 1))% |" for k in 1:5], "\n")
    Markdown.parse(
        """
        | Lobe | share of the lung's volume |
        |---|---|
        $rows

        Fissure voxels: $(count(BLOCKERS.fissure)) ($(round(100 * count(BLOCKERS.fissure) / total; digits = 1))% of the lung). Stubs every $(round(STUB_LENGTHS.artery.spacing; digits = 2)) mm on the arteries (median imaged diameter $(round(STUB_LENGTHS.artery.diameter; digits = 2)) mm; $(artery_inputs.report.stubs) of $(artery_inputs.report.sites) sites took one), every $(round(STUB_LENGTHS.vein.spacing; digits = 2)) mm on the veins ($(round(STUB_LENGTHS.vein.diameter; digits = 2)) mm; $(sum(v.report.stubs for v in vein_inputs)) of $(sum(v.report.sites for v in vein_inputs)) sites).
        """
    )
end

# ╔═╡ 611b7d1b-991d-45eb-987e-c390b6498a04
md"""
## 3. Grow

HCT-CCO's staged constrained constructive optimization, on the GPU: `LUNG_GROWTH` grows one arterial terminal per acinus (29,912, by micro-CT stereology of human lungs) from the main pulmonary artery, and `LUNG_VEIN_GROWTH` one venule per acinus into each lung's pulmonary veins. The vessels inside each acinus are not grown: the flow model carries them as Huang's per-order tables below each terminal.
"""

# ╔═╡ 2bf2e225-251c-44c2-a310-fa151e5b4971
"Grow one forest: the domain and its pass-through voxels as a `.vol` file, the inlet trees as GXL files."
function grow_forest(spec, trees, pass, blockers, name)
    dir = mktempdir()
    lower = Tuple(round.(Int, growth_grid.origin))
    allow = zeros(UInt16, growth_grid.dims)
    allow[pass] .= V.NOTERM # vessels may cross these voxels, no terminal may end in them

    domain = V.write_vol(joinpath(dir, "$name.vol"), LUNG .| pass, lower)
    inlets = [V.XCAT.write_base_gxl(joinpath(dir, "$(name)_$k.xml"), tree) for (k, tree) in enumerate(trees)]
    engine = V.grow(V.GrowthSpec(spec; device = DEVICE), domain, inlets; blockers, territories = allow)
    return V.vtree(engine)
end

# ╔═╡ 05ecf689-b662-41ba-b3fd-2b39e3060317
arteries = cached("arteries") do
    # the preset is this subject's: its tissue is these lung voxels and its root the trunk's median ring
    @assert V.LUNG_GROWTH.tissue_ml ≈ count(LUNG) / 1000 && V.LUNG_GROWTH.root_caps[1] ≈ artery_base.r[2]
    grow_forest(V.LUNG_GROWTH, [artery_inputs.tree], artery_inputs.pass, BLOCKERS.artery, "arteries")
end;

# ╔═╡ f909bf22-0953-4b41-9822-036f7042ccca
veins = cached("veins") do
    @assert all(V.LUNG_VEIN_GROWTH.root_caps .≈ [tree.r[1] for tree in vein_bases])
    grow_forest(
        V.LUNG_VEIN_GROWTH,
        [v.tree for v in vein_inputs],
        reduce(.|, [v.pass for v in vein_inputs]),
        BLOCKERS.vein,
        "veins",
    )
end;

# ╔═╡ 5a45700e-ddd3-40d4-8504-4712118904e9
"""
The lung's surface: each surface voxel's centre (mm) and its number among the lung voxels, so per-voxel values (such as
perfusion) can be painted on it.
"""
LUNG_SURFACE = let
    function surface_of(lung, origin)
        voxels = findall(lung)
        on_surface(I) = any(s -> !checkbounds(Bool, lung, I + s) || !lung[I + s], SIX_STEPS)
        surface = [n for (n, I) in enumerate(voxels) if on_surface(I)]
        return (points = [Mke.Point3f(Tuple(voxels[n]) .- 1 .+ origin) for n in surface], voxels = surface)
    end
    surface_of(LUNG, growth_grid.origin)
end;

# ╔═╡ daa1c63b-b988-4d38-ad53-c22630d11af2
"Draw every fourth point of the lungs' surface as a faint shell around the vessels."
function draw_lung_shell!(ax)
    return Mke.scatter!(ax, LUNG_SURFACE.points[1:4:end]; color = (:gray50, 0.04), markersize = 4, strokewidth = 0)
end

# ╔═╡ d5e4c316-fadc-4dbe-8f51-471929992f45
"""
Draw a grown tree's segments of radius ≥ `min_radius` (mm): width and shade grow with radius, and the largest vessels
are drawn last so they sit on top.
"""
function draw_vessels!(ax, tree; colormap, min_radius)
    shown = filter(j -> tree.r[j] >= min_radius, sortperm(tree.r))
    segments = [Mke.Point3f(tree.a[:, j]) => Mke.Point3f(tree.b[:, j]) for j in shown]
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

# ╔═╡ f75080e9-6836-4f3e-839d-cbc53c637a93
md"""
show vessels of radius ≥ $(@bind min_radius_shown Select([1.0 => "1 mm", 0.5 => "0.5 mm", 0.25 => "0.25 mm"]; default = 0.5))
"""

# ╔═╡ b44fc616-883b-47b7-a8d1-4d60879cc96d
let fig = Mke.Figure(size = (1300, 760), figure_padding = 4)
    ax = vessel_axis(fig[1, 1]; title = "Grown arteries (red) and veins (blue), radius ≥ $(min_radius_shown) mm")
    draw_lung_shell!(ax)
    draw_vessels!(ax, veins; colormap = VEIN_COLOURS, min_radius = min_radius_shown)
    draw_vessels!(ax, arteries; colormap = ARTERY_COLOURS, min_radius = min_radius_shown)

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
            "Terminal radius (median %.0f µm arteries, %.0f µm veins)",
            median(terminal_radius_um["arteries"]),
            median(terminal_radius_um["veins"]),
        ),
        xlabel = "µm",
        ylabel = "terminals",
        xscale = log10,
    )
    for (name, colour) in (("arteries", :red3), ("veins", :royalblue))
        Mke.hist!(ax2, terminal_radius_um[name]; bins = exp10.(range(1, 3; length = 50)), color = (colour, 0.5), label = name)
    end
    Mke.axislegend(ax2)
    Mke.colsize!(fig.layout, 1, Mke.Relative(0.7))
    fig
end

# ╔═╡ 116e3f65-4c0f-4edb-a35f-a6abcf8697a5
md"""
## 4. Flow

No flow is prescribed: the mean pulmonary arterial pressure (17 mmHg) and the left atrial pressure (11 mmHg) drive blood through every vessel, the alveolar capillary sheets and the veins. Normoxic pulmonary vessels carry little tone, so the lung is **passive**: vessels distend with their pressure and the capillary sheets with theirs, and there is one state. Each 1 mm lung voxel gets its flow from the terminals that feed it.

The lung carries the whole cardiac output, so in this phantom the cardiac output **is** the lung's own flow (section 5): the body and the lung carry the same blood.
"""

# ╔═╡ 27a52b36-0fa9-4f86-bcd0-4cf5824daecb
"Centres (mm) of the 1 mm lung voxels: the tissue the flow is solved for."
tissue_centres = let voxels = findall(LUNG), origin = growth_grid.origin
    Float32[origin[d] + voxels[k][d] - 1 for d in 1:3, k in eachindex(voxels)]
end

# ╔═╡ 55a43d63-8f38-4962-b428-9c0f28e814c7
flows = cached("flow") do
    V.flow(V.LUNG_FLOW, arteries, veins, tissue_centres; to = DEVICE)
end;

# ╔═╡ ab72236e-a88a-4eb6-accf-8d67de2ac41f
"Blood flow (mL/min per mL of lung) of every 1 mm lung voxel in circuit `circuit`."
function voxel_perfusion(circuit)
    terminal_flow = reshape(Float32[circuit.q[j] for j in circuit.terminals], 1, :)
    voxels = collect(Int32(1):Int32(size(tissue_centres, 2)))
    # mL/s per mL of tissue
    return vec(V.OrganCT.tissue_iodine(terminal_flow, flows.pairs, flows.terr, voxels)) .* 60
end

# ╔═╡ a40e8ed4-a6f0-4de7-8a9f-c0002a81dab2
"A 1 mm coronal map (x, z) of `values` (one per lung voxel) at y (mm); NaN elsewhere."
function lung_coronal(values, y)
    j = clamp(round(Int, y - growth_grid.origin[2]) + 1, 1, growth_grid.dims[2])
    slice = fill(NaN32, growth_grid.dims[1], growth_grid.dims[3])
    for (n, voxel) in enumerate(findall(LUNG))
        voxel[2] == j && (slice[voxel[1], voxel[3]] = values[n])
    end
    return slice
end

# ╔═╡ 32162101-95fb-4558-a8d7-57590ad05179
md"""
y (mm, front to back): $(@bind y_flow Slider(growth_grid.origin[2]:2:(growth_grid.origin[2] + growth_grid.dims[2] - 1); default = 262, show_value = true))
"""

# ╔═╡ 7ba30cb9-e856-4e5f-a50b-4c14a0c517cd
let fig = Mke.Figure(size = (1200, 560)), C = flows.rest
    xs = growth_grid.origin[1] .+ (0:(growth_grid.dims[1] - 1))
    zs = growth_grid.origin[3] .+ (0:(growth_grid.dims[3] - 1))
    ax = Mke.Axis(
        fig[1, 1];
        title = @sprintf("Lung blood flow · %.2f L/min · coronal y = %d mm", 60 * V.organ_flow(C) / 1000, y_flow),
        aspect = Mke.DataAspect(),
    )
    map = Mke.heatmap!(
        ax,
        xs,
        zs,
        lung_coronal(voxel_perfusion(C), y_flow);
        colormap = :turbo,
        colorrange = (0, 4),
        nan_color = :black,
    )
    Mke.hidedecorations!(ax)
    Mke.Colorbar(fig[1, 2], map; label = "mL/min per mL of lung")
    fig
end

# ╔═╡ 56f6f8b4-669c-4bb8-9e3e-b693c342d67c
let C = flows.rest
    v = C.volumes
    # the trees' and sheets' blood, imaged vessels included; its mean transit is volume over flow (central volume principle)
    blood_ml = v.art_resolved + v.art_unresolved + v.capillary + v.ven_resolved + v.ven_unresolved
    Markdown.parse(
        @sprintf(
            """
            | | model | measured |
            |---|---|---|
            | lung blood flow (L/min) | %.2f | 4.8 ± 0.8 at rest (Bentley 2020) |
            | pressure drop: arteries · capillary sheets · veins | %.2f · %.2f · %.2f | arteries ≈ veins, sheets < 0.16 (Hakim 1982, dog lobes) |
            | pulmonary blood volume (mL) | %.0f | 526 ± 87 (Nelsson 2021), 471 ± 122 (Kanski 2013) |
            | mean pulmonary transit, valve to left atrium (s) | %.1f | 5.7 ± 1.4 (Nelsson 2021), 5.3 (Dock 1961) |

            The flow is about a quarter below the measured mean, within two standard deviations, and the grown veins take more of the pressure drop than the arteries. Both come from the venous tree: it is grown with the arteries' branching exponent from two ostia per lung (no lung venous exponent is measured), which makes its small veins narrower than human casts show. The blood volume is close to the measured, so the lower flow makes the mean transit longer than measured: iodine stays in the lungs longer, and less of it has left by any scan time than in a real lung. These are stated, not tuned: the flow is a prediction from pressures and geometry, and the phantom's ground truth is whatever flow the model gives.
            """,
            60 * V.organ_flow(C) / 1000,
            C.resist.art_share,
            C.resist.cap_share,
            C.resist.ven_share,
            blood_ml,
            blood_ml / V.organ_flow(C),
        ),
    )
end

# ╔═╡ d5da5d67-513f-49d9-a70b-2bfead33c5fd
md"""
## 5. Transport

50 mL of iohexol 350 at 5 mL/s and 40 mL of saline into the right arm, carried through the whole-body circulation, the right heart and the lung's own trees and capillary sheets, then the left heart and the body, with recirculation. The cardiac output is the lung's solved flow (`V.organ_flow`), so the bolus is diluted into the blood the lung really carries.
"""

# ╔═╡ 6817b779-c89b-47ce-b07e-ad48831b78cb
"Contrast transport through the body and the lung circuit `circuit`, in its one (passive) state."
function transport_state(circuit)
    spec = V.TransportSpec(V.LUNG_TRANSPORT; phys = circuit.phys, tissue_ml = Float32(size(tissue_centres, 2) / 1000))
    return ("rest", circuit, V.transport(spec, circuit, V.Body.Injection(); times = V.LUNG_PHANTOM.times))
end

# ╔═╡ 45a3f241-cc85-4202-ad0e-00de0342f5fd
"""
Contrast transport through `circuit`, solved once and kept in `CACHE_DIR/name.jls`; returns that file. The full result
(every vessel's iodine at every frame) is large, so the notebook does not hold it: `save_phantom` reads it while it writes,
and the figures use `transport_curves`.
"""
function transport_cache(name, circuit)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) || serialize(path, [transport_state(circuit)])
    return path
end

# ╔═╡ f6016974-f880-4c69-90fc-37a14b4ff4ce
TRANSPORT = transport_cache("transport", flows.rest);

# ╔═╡ e5b80139-a62d-4497-843b-77d4e1d2a6cf
"The RV and LA blood and the whole lung's iodine (mg/mL) over time, and the cardiac output, from a transport file."
function transport_curves(path)
    tissue_ml = size(tissue_centres, 2) / 1000
    _, circuit, result = only(deserialize(path))
    held = vec(sum(result.tissue; dims = 2) .+ sum(result.artery; dims = 2) .+ sum(result.vein; dims = 2))
    return (
        blood_times = result.body.t,
        rv = result.body.mean[:rv],
        la = result.body.mean[:la],
        lung = held ./ tissue_ml,
        co_l_min = 60 * result.body.circ.co / 1000,
    )
end

# ╔═╡ 5e4a06c3-bed3-4a8e-9550-dd57bb4a33de
TRANSPORT_CURVES = transport_curves(TRANSPORT);

# ╔═╡ 9f8d5635-baf7-43f2-ba01-6d006b3d6a77
let fig = Mke.Figure(size = (900, 420)), curves = TRANSPORT_CURVES
    ax = Mke.Axis(
        fig[1, 1];
        title = @sprintf("Cardiac output %.2f L/min (the lung's flow)", curves.co_l_min),
        xlabel = "time after injection (s)",
        ylabel = "iodine (mg/mL)",
    )
    Mke.lines!(ax, curves.blood_times, curves.rv; color = :royalblue, linewidth = 3, label = "RV blood")
    Mke.lines!(ax, curves.blood_times, curves.la; color = :firebrick, linewidth = 3, label = "LA blood")
    Mke.lines!(ax, V.LUNG_PHANTOM.times, curves.lung; color = :gray30, linewidth = 3, label = "lung (all its iodine)")
    Mke.xlims!(ax, 0, 40)
    Mke.axislegend(ax)
    fig
end

# ╔═╡ b37414bb-2941-4a47-85c3-0c10c246cc0c
md"""
## 6. Compartments

A **compartment** is the region ground truth is given for: a region of the lungs with everything in it, the lung tissue and every artery, vein and lumen that lies in it or overlaps it. Here: the five lobes. The lungs, and the right and left lung, are groups of them.

**Edit the compartments here**: `split_lobes = false` gives one compartment per lung.

The rule below is applied once, when the phantom is voxelized (section 7), and every compartment number and picture after that comes from the saved phantom (section 8), so there is one source of truth.
"""

# ╔═╡ df85dc89-82d9-4d4d-9665-8282155840fb
split_lobes = true # the five lobes apart; false: one compartment per lung

# ╔═╡ 4e169ac1-210c-49d8-8f5c-d7c62084fb18
begin
    ROI_NAMES = split_lobes ? copy(LOBAR_BRONCHI.names) : ["right lung", "left lung"]
    ROI_GROUPS =
        split_lobes ? ["lungs" => [1, 2, 3, 4, 5], "right lung" => [1, 2, 3], "left lung" => [4, 5]] :
        ["lungs" => [1, 2]]
    ROI_NAMES, ROI_GROUPS
end

# ╔═╡ 6b9d2a4f-3e5c-4f7b-8a1d-2c4e6f8b0d13
"`labels` with each nonzero label carried `steps` voxels out into the zeros around it, nearest first."
function carried(labels, steps)
    out = copy(labels)
    front = findall(!=(0), out)
    for _ in 1:steps
        next = CartesianIndex{3}[]
        for c in front, s in SIX_STEPS
            n = c + s
            checkbounds(Bool, out, n) && out[n] == 0 || continue
            out[n] = out[c]
            push!(next, n)
        end
        front = next
    end
    return out
end

# ╔═╡ 31928ce9-9fa0-402a-ba97-cfc37d5898e5
"""
Where the lungs would be with their vessels: XCAT's lung closed by 5 mm on the 1 mm grid, so a vessel lumen inside or
overlapping the lung is in it (the hilar vessels outside the lung are not); and the lobes carried 15 mm out from the lung, so
those lumens get the lobe they lie in.
"""
LUNG_ENVELOPE = cached("lung_envelope") do
    envelope = V.XCAT.closing(V.XCAT.label_mask(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM, growth_grid, (0x06,); fraction = 1.0e-9), 5)
    (inside = envelope, lobe = carried(LOBES, 15))
end;

# ╔═╡ 4b166f97-4217-4a97-95b2-edab3ff8aa46
"""
    make_roi(labels, envelope, grid, n_rois; split_lobes) -> roi

The compartment rule `roi(p)` for the 0.2 mm voxel at `p` (0 = none): lung tissue takes its lobe; a vessel lumen takes the
lobe it lies in when it is inside the lungs' envelope (both on `grid`). (Only tissue and vessel voxels take a compartment.)
Everything it reads is captured with its type, so it compiles to fast code.
"""
function make_roi(labels, envelope, grid, n_rois; split_lobes)
    inside = envelope.inside
    lobe = envelope.lobe
    origin = grid.origin
    return function roi(p)
        voxel = round.(Int, (p .- XCAT_ORIGIN) ./ XCAT_VOXEL_MM) .+ 1
        all(1 .<= voxel .<= XCAT_DIMS) || return 0
        k = round.(Int, p .- origin) .+ 1
        all(1 .<= k .<= size(lobe)) || return 0
        (labels[voxel...] == 0x06 || inside[k...]) || return 0
        l = Int(lobe[k...])
        l == 0 && return 0
        return split_lobes ? l : (l <= 3 ? 1 : 2)
    end
end

# ╔═╡ 08f3aa33-2562-424a-8526-349c90fdd559
roi = make_roi(xcat_labels, LUNG_ENVELOPE, growth_grid, length(ROI_NAMES); split_lobes)

# ╔═╡ b83e74d9-3d5e-4059-b8e6-0e6307b119bf
"Categorical colours of the compartments (`LOBE_COLOURS`, or one per lung)."
roi_colormap() = Mke.cgrad(
    [:black; split_lobes ? LOBE_COLOURS : ["#e41a1c", "#377eb8"]],
    length(ROI_NAMES) + 1;
    categorical = true,
)

# ╔═╡ ac7c95a2-969d-41f6-bf40-b299016a01ea
md"""
## 7. Voxelize and Save

The 0.2 mm phantom of the whole XCAT chest, one folder per scenario. Every voxel gets a **label** (its material and iodine curve) and a **compartment**. Every label lies in one compartment, so a compartment's iodine and ground truth are sums over its labels, and painting conserves iodine exactly: the build fails if the lungs' iodine ledger does not balance.

| File | What it holds |
|---|---|
| `labels.u16.zst` | every voxel's label |
| `rois.u8.zst` | every voxel's compartment |
| `mbf_rest.f32.zst` | every voxel's ground-truth perfusion (mL/min/g) |
| `label_rest.csv` | each label's compartment, material, ground-truth perfusion and iodine every second |
| `ledger_rest.csv` | the lungs' iodine every second: delivered, held, out through the pulmonary veins, painted |
| `phantom.toml` | the grid, the compartment names and groups, the tracking points |
"""

# ╔═╡ e54d8f12-2cb3-4c9b-a1ef-a68e4ad5ec3a
begin
    "1 mm lung voxel number of each growth-grid voxel (0: not lung)."
    const TISSUE_VOXEL_NUMBER = V.XCAT.tissue_numbers(LUNG)

    # the imaged segments of the grown trees, indexed for naming XCAT's drawn pulmonary lumens near them
    const NEAR_MM = 2.0 # a drawn lumen is an imaged segment's when it lies within this of the segment's surface
    # (the arteries' reaches 10 mm: inside the main pulmonary artery's surface every voxel is the trunk's, wherever its
    # grown radius puts the segment's own surface)
    const ARTERY_INDEX = V.XCAT.segment_index(V.XCAT.imaged_segments(arteries); reach = 10.0)
    const VEIN_INDEX = V.XCAT.segment_index(V.XCAT.imaged_segments(veins); reach = NEAR_MM)
end;

# ╔═╡ aeeda182-58d5-45f5-b269-bb9bb53a54d3
"""
Which great vessel each vessel voxel belongs to, and where it lies along it, on a 1 mm grid: its tree and the blood (mL)
between it and the heart (`V.XCAT.great_vessel_transit`). Each state's transit is volume over flow (`V.pool_delay`).
"""
GREAT_VESSEL_TRANSIT = cached("great_vessel_transit") do
    V.XCAT.great_vessel_transit(xcat_labels, XCAT_ORIGIN, XCAT_VOXEL_MM)
end;

# ╔═╡ 75243d5b-b7b8-4b71-95c4-8648072543b4
"""
The main pulmonary artery's lumen, at 0.4 mm: inside XCAT's `dias_pul_art` surface. XCAT's label volume paints this lumen
as the right ventricle's blood (label 20), so its surface, not its label, tells the trunk from the ventricle; its volume
(mL) is the blood between the valve and the pulmonary arteries the label volume does show.
"""
PULMONARY_TRUNK = cached("pulmonary_trunk") do
    points, triangles = V.XCAT.tessellate(xcat_surfaces["dias_pul_art"]; nu = 160, nv = 160)
    grid = V.grid_around(points, 0.4; margin = 2)
    inside = V.voxelize(points, triangles, grid)
    (grid = grid, inside = inside, ml = count(inside) * 0.4^3 / 1000)
end;

# ╔═╡ 7d1e3a52-4b6c-4f0e-9a21-3c5e8f7b9d14
"""
The lungs' phantom spec: the body pools, with the descending aorta timed from this subject's arch and the pulmonary
arteries the label volume shows timed from the valve (the trunk's blood before them).
"""
PHANTOM_SPEC = let pools = copy(V.LUNG_PHANTOM.pools)
    function lead(name, ml)
        i = findfirst(p -> p.name == name, pools)
        d = pools[i]
        return pools[i] = V.Pool(d.name, d.material, d.segment, d.what; through = d.through, lead_ml = ml, lead_through = d.lead_through)
    end
    lead("descending aorta", GREAT_VESSEL_TRANSIT.apex_ml)
    lead("pulmonary arteries", Float32(PULMONARY_TRUNK.ml))
    V.PhantomSpec(V.LUNG_PHANTOM; pools)
end;

# ╔═╡ 9a5e0d2b-6717-4c81-9684-d2a0cb94718e
"""
    make_classifier(labels, grid, numbers, artery_index, vein_index, transit, trunk) -> classify

The classifier `classify(p)`: what is at the 0.2 mm voxel centre `p`. Lung tissue is its 1 mm voxel; a pulmonary vessel
lumen (XCAT's labels 25 and 125 on the pulmonary arterial tree, 24 and 124 on the pulmonary venous tree, and the ventricle's
label inside the main pulmonary artery's surface, `trunk`) is the imaged segment it lies in or beside, else its vessel's pool
and its place along the vessel; every other great-vessel voxel is its
pool; every other XCAT label is decided by the label alone (a body pool, `V.XCAT_POOL_OF_LABEL`, the heart lumped, or a
static material). Everything it reads is captured with its type, so it compiles to fast code.
"""
function make_classifier(labels, grid, numbers, artery_index, vein_index, transit, trunk)
    # by the label alone: 0 needs the voxel's geometry
    by_label = fill((V.VOX_POOL, Int32(24), 0.0f0), 256)
    needs_geometry = falses(256)
    for label in 0x00:0xff
        if label in (0x06, 0x14, 0x18, 0x19, 0x1c, 0x7c, 0x7d)
            needs_geometry[label + 1] = true
        elseif haskey(V.XCAT_POOL_OF_LABEL, label)
            by_label[label + 1] = (V.VOX_POOL, Int32(V.XCAT_POOL_OF_LABEL[label]), 0.0f0)
        elseif haskey(V.XCAT_STATIC, Int(label))
            by_label[label + 1] = (V.VOX_STATIC, Int32(label), 0.0f0)
        end
    end
    vessel_pool(p, t) = V.XCAT.great_vessel_pool(transit, p, t)
    ventricle = (V.VOX_POOL, Int32(V.XCAT_POOL_OF_LABEL[0x14]), 0.0f0)
    # a pulmonary lumen: the imaged segment it lies in or beside, else its pool along the vessel
    function imaged_or_pool(p, index, sign, t)
        j, distance = V.XCAT.nearest_segment(index, p)
        return j > 0 && distance <= NEAR_MM ? (V.VOX_VESSEL, sign * Int32(j), 0.0f0) : vessel_pool(p, t)
    end

    in_trunk(p) = (
        k = round.(Int, (p .- trunk.grid.origin) ./ trunk.grid.h) .+ 1;
        all(1 .<= k .<= size(trunk.inside)) && trunk.inside[k...]
    )

    return function classify(p)
        voxel = round.(Int, (p .- XCAT_ORIGIN) ./ XCAT_VOXEL_MM) .+ 1
        all(1 .<= voxel .<= XCAT_DIMS) || return (V.VOX_STATIC, Int32(0), 0.0f0)
        label = labels[voxel...]
        needs_geometry[label + 1] || return by_label[label + 1]

        # the main pulmonary artery, painted as the ventricle's blood: its imaged trunk segment
        if label == 0x14
            in_trunk(p) || return ventricle
            j, _ = V.XCAT.nearest_segment(artery_index, p)
            return j > 0 ? (V.VOX_VESSEL, Int32(j), 0.0f0) : ventricle
        end
        if label == 0x06
            number = V.XCAT.nearest_tissue_voxel(numbers, round.(Int, p .- grid.origin) .+ 1)
            return number > 0 ? (V.VOX_TISSUE, number, 0.0f0) : (V.VOX_POOL, Int32(26), 0.0f0)
        end
        label == 0x7d && return imaged_or_pool(p, artery_index, Int32(1), 0x02)
        label == 0x7c && return imaged_or_pool(p, vein_index, Int32(-1), 0x05)
        label == 0x1c && return vessel_pool(p, 0x01)
        if label == 0x19
            t = V.XCAT.vessel_tree_at(transit, p, (0x01, 0x02))
            return t == 0x02 ? imaged_or_pool(p, artery_index, Int32(1), t) : vessel_pool(p, t)
        end
        t = V.XCAT.vessel_tree_at(transit, p, (0x03, 0x04, 0x05, 0x06, 0x07))
        return t == 0x05 ? imaged_or_pool(p, vein_index, Int32(-1), t) : vessel_pool(p, t)
    end
end

# ╔═╡ 01e67f8a-e865-40a9-979b-51da0a3e3f5d
classify = make_classifier(
    xcat_labels,
    growth_grid,
    TISSUE_VOXEL_NUMBER,
    ARTERY_INDEX,
    VEIN_INDEX,
    GREAT_VESSEL_TRANSIT,
    PULMONARY_TRUNK,
)

# ╔═╡ 8c0e3b5a-4f6d-4a8c-9b2e-3d5f7a9c1e24
"The point (mm) of label `label`'s voxels deepest inside them (farthest from their surface), on a 1 mm grid."
function deepest_point(labels, label)
    step = round(Int, 1 / XCAT_VOXEL_MM)
    coarse = labels[1:step:end, 1:step:end, 1:step:end]
    blood = coarse .== label
    depth = fill(-1, size(blood))
    queue = [c for c in findall(blood) if any(s -> !checkbounds(Bool, blood, c + s) || !blood[c + s], SIX_STEPS)]
    foreach(c -> depth[c] = 0, queue)
    head = 1
    while head <= length(queue)
        c = queue[head]
        head += 1
        for s in SIX_STEPS
            n = c + s
            checkbounds(Bool, blood, n) && blood[n] && depth[n] < 0 || continue
            depth[n] = depth[c] + 1
            push!(queue, n)
        end
    end
    return Tuple(XCAT_ORIGIN .+ (Tuple(argmax(depth)) .- 1) .* (step * XCAT_VOXEL_MM))
end

# ╔═╡ 8811c92e-4cb3-4c0c-add3-6e7da17660a0
"""
The bolus-tracking points: the main pulmonary artery at the middle of its imaged trunk (the centre of its lumen there), and
the right ventricle at its blood's deepest point (farthest from its wall, on a 1 mm grid).
"""
TRACKING_POINTS = let
    trunk = findall(==("dias_pul_art"), artery_base.vessel)
    artery = Tuple(artery_base.X[:, trunk[(length(trunk) + 1) ÷ 2]])
    ["main pulmonary artery" => artery, "right ventricle" => deepest_point(xcat_labels, 0x14)]
end

# ╔═╡ 15d45948-dcea-4ac6-a552-7c23db9a5988
"""
Write one scenario's phantom into `PHANTOM_DIR/name`, once, from its transport file (read here, and let go once the
phantom is written).
"""
function save_phantom(name, transport_file; spec = PHANTOM_SPEC, classify = classify)
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
        facts = ["organ" => "lungs (XCAT vmale_50, 0.2 mm)", "model" => "Vasculature.jl"],
        verbose = true,
    )
    @info "$name: $(round(seconds; digits = 1)) s"
    scans = nothing
    GC.gc()
    return dir
end

# ╔═╡ e75df31a-d1ec-4690-8b58-f8177b34886c
HEALTHY_DIR = save_phantom("healthy", TRANSPORT);

# ╔═╡ 9a8d5350-20ea-4e6a-820a-0e4d8a3045d1
md"""
## 8. Inspect the Saved Phantom

Only part of the phantom is read back: an 80 mm slab through the hila for the axial views, and one coronal plane through the main pulmonary artery. Labels and compartments are categorical; the idealized 70 keV image (the phantom as a perfect scanner would see it) is gray, in a mediastinal and a lung window; iodine and ground-truth perfusion use `:turbo`.

Ground-truth perfusion is shown per compartment, the number a CT measurement of that region is compared with: the flow the compartment's capillary beds receive over the mass of everything in it, lung tissue and the vessels in it alike.
"""

# ╔═╡ 69f6db9d-70cf-4c42-8ad4-684b798cdb3e
"""
The lungs' box of a saved phantom, 80 mm of it about the main pulmonary artery's tracking point (z), the only part read back.
Returns the bundle (labels and compartments) and its first z (mm). Ground truth is each compartment's own, from the bundle's
tables.
"""
function inspect(dir; half_mm = 40)
    point = V.DynamicPhantom.load(dir; volume = false).meta["tracking"]["main pulmonary artery"]["point_mm"]
    lower = collect(growth_grid.origin) .- 5
    upper = collect(growth_grid.origin) .+ collect(growth_grid.dims) .+ 5
    voxel(mm, d) = clamp(round(Int, (mm - XCAT_ORIGIN[d]) / XCAT_VOXEL_MM) + 1, 1, XCAT_DIMS[d])
    bundle = V.DynamicPhantom.load(
        dir;
        rois = true,
        zrange = voxel(point[3] - half_mm, 3):voxel(point[3] + half_mm, 3),
        xrange = voxel(lower[1], 1):voxel(upper[1], 1),
        yrange = voxel(lower[2], 2):voxel(upper[2], 2),
    )
    return (bundle = bundle, z0 = bundle.origin_mm[3])
end

# ╔═╡ 0c29d1f3-3e97-49c3-9ccf-2ed17a205c4a
healthy = inspect(HEALTHY_DIR);

# ╔═╡ 70f5505d-1073-4b93-82b7-3c7674b657fd
"Per-label lookups of a bundle at time `t`: kind (for pictures), idealized 70 keV HU and iodine (mg/mL)."
function label_lookups(bundle, t)
    kind_of = zeros(UInt8, 65536)
    for (label, row) in bundle.table
        label == 0 && continue
        kind_of[label + 1] =
            row.kind == "static" ? 1 :
            row.kind == "pool" ? (occursin("blood", row.base) ? 3 : 2) :
            Dict("tissue" => 4, "artery" => 5, "vein" => 6)[row.kind]
    end
    hu = fill(-1000.0f0, 65536)
    for (label, value) in V.DynamicPhantom.hu(bundle, "rest", t)
        hu[label + 1] = value
    end
    iodine = zeros(Float32, 65536)
    for (label, value) in V.DynamicPhantom.iodine(bundle, "rest", t)
        iodine[label + 1] = value
    end
    return (; kind_of, hu, iodine)
end

# ╔═╡ 5fe0c8fa-fe29-4114-b0d1-3076eb730ede
# the kinds of label, for pictures
const KIND_COLOURS, KIND_NAMES = (
    [:gray30, :gray65, :mediumpurple, :gray45, :red, :deepskyblue],
    ["static (air, cartilage)", "body tissue", "blood", "lung", "artery", "vein"],
)

# ╔═╡ c238aa74-eace-4853-8b12-80a69ae1cd3a
"""
One heatmap panel with its own scale underneath (a colour bar, or a legend for categories), at `fig[row, column]`; `xs`
and `ys` (mm) span `image`.
"""
function panel!(fig, row, column, title, xs, ys, image; colormap, colorrange, label = nothing, legend = nothing, yreversed = true)
    cell = fig[row, column] = Mke.GridLayout()
    ax = Mke.Axis(cell[1, 1]; title, aspect = Mke.DataAspect(), yreversed)
    map = Mke.heatmap!(ax, xs, ys, image; colormap, colorrange)
    Mke.hidedecorations!(ax)
    if legend === nothing
        Mke.Colorbar(cell[2, 1], map; label, vertical = false, flipaxis = false, width = Mke.Relative(0.9), tellwidth = false)
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

# ╔═╡ d99468dc-2991-44f5-becc-b8edf15f5778
"""
One axial phantom slice in six panels: labels by kind, compartments and ground-truth perfusion; the idealized 70 keV image
at time `t` in a mediastinal and a lung window, and iodine at `t`.
"""
function slice_figure(view, t, z)
    bundle = view.bundle
    k = clamp(round(Int, (z - view.z0) / XCAT_VOXEL_MM) + 1, 1, size(bundle.labels, 3))
    is, js = axes(bundle.labels, 1), axes(bundle.labels, 2)
    xs = bundle.origin_mm[1] .+ (is .- 1) .* XCAT_VOXEL_MM
    ys = bundle.origin_mm[2] .+ (js .- 1) .* XCAT_VOXEL_MM
    labels = bundle.labels[is, js, k]
    rois = bundle.rois[is, js, k]
    look = label_lookups(bundle, t)
    # each compartment's ground truth over all of it: its beds' flow over the mass of everything in it
    truth = [0.0f0; [Float32(bundle.roi_table[("rest", name)].mbf) for name in ROI_NAMES]]

    fig = Mke.Figure(; size = (1500, 1000))
    kinds = [look.kind_of[l + 1] for l in labels]
    panel!(
        fig, 1, 1, "Labels by kind", xs, ys, kinds;
        colormap = Mke.cgrad([:black; KIND_COLOURS], 7; categorical = true), colorrange = (-0.5, 6.5), legend = (KIND_COLOURS, KIND_NAMES)
    )
    panel!(
        fig, 1, 2, "Compartments", xs, ys, rois;
        colormap = roi_colormap(), colorrange = (-0.5, length(ROI_NAMES) + 0.5), legend = (roi_colormap().colors[2:end], ROI_NAMES)
    )
    panel!(
        fig, 1, 3, "Ground-truth perfusion by compartment", xs, ys, [truth[r + 1] for r in rois];
        colormap = :turbo, colorrange = (0, 12), label = "mL/min/g"
    )
    hu = [look.hu[l + 1] for l in labels]
    panel!(fig, 2, 1, "Idealized 70 keV, t = $t s", xs, ys, hu; colormap = :grays, colorrange = (-200, 600), label = "HU")
    panel!(fig, 2, 2, "Idealized 70 keV, lung window, t = $t s", xs, ys, hu; colormap = :grays, colorrange = (-1000, 200), label = "HU")
    panel!(
        fig, 2, 3, "Iodine, t = $t s", xs, ys, [look.iodine[l + 1] for l in labels];
        colormap = :turbo, colorrange = (0, 8), label = "mg/mL"
    )
    Mke.Label(fig[0, :], "z = $(round(z; digits = 1)) mm"; fontsize = 20, font = :bold)
    for column in 1:3
        Mke.colsize!(fig.layout, column, Mke.Relative(1 / 3))
    end
    return fig
end

# ╔═╡ 74ad599e-a3f3-423a-850b-ed84438f8dfe
md"""
* time $(@bind healthy_time Slider(0:1:40; default = 10, show_value = true)) s
* z $(@bind healthy_z Slider(round.(healthy.z0 .+ (0:10:(size(healthy.bundle.labels, 3) - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 290.0, show_value = true)) mm
"""

# ╔═╡ aa982a9c-9ff8-4e46-97d8-65c13a957f3c
slice_figure(healthy, healthy_time, healthy_z)

# ╔═╡ 2abf3ae2-4e6a-4e90-9f2e-df9a68f29a56
"""
One coronal plane of a saved phantom, through the main pulmonary artery's tracking point, at time `t`: compartments,
ground-truth perfusion by compartment, the idealized 70 keV image in a lung window, and iodine. Reads the whole phantom
once for this one plane.
"""
function coronal_figure(dir, t)
    point = V.DynamicPhantom.load(dir; volume = false).meta["tracking"]["main pulmonary artery"]["point_mm"]
    j = clamp(round(Int, (point[2] - XCAT_ORIGIN[2]) / XCAT_VOXEL_MM) + 1, 1, XCAT_DIMS[2])
    voxel(mm, d) = clamp(round(Int, (mm - XCAT_ORIGIN[d]) / XCAT_VOXEL_MM) + 1, 1, XCAT_DIMS[d])
    bundle = V.DynamicPhantom.load(
        dir;
        rois = true,
        yrange = j:j,
        xrange = voxel(growth_grid.origin[1] - 5, 1):voxel(growth_grid.origin[1] + growth_grid.dims[1] + 5, 1),
        zrange = voxel(growth_grid.origin[3] - 5, 3):voxel(growth_grid.origin[3] + growth_grid.dims[3] + 5, 3),
    )
    xs = bundle.origin_mm[1] .+ (axes(bundle.labels, 1) .- 1) .* XCAT_VOXEL_MM
    zs = bundle.origin_mm[3] .+ (axes(bundle.labels, 3) .- 1) .* XCAT_VOXEL_MM
    labels = bundle.labels[:, 1, :]
    rois = bundle.rois[:, 1, :]
    look = label_lookups(bundle, t)
    truth = [0.0f0; [Float32(bundle.roi_table[("rest", name)].mbf) for name in ROI_NAMES]]

    fig = Mke.Figure(; size = (1600, 560))
    panel!(
        fig, 1, 1, "Compartments", xs, zs, rois;
        colormap = roi_colormap(), colorrange = (-0.5, length(ROI_NAMES) + 0.5), legend = (roi_colormap().colors[2:end], ROI_NAMES), yreversed = false
    )
    panel!(
        fig, 1, 2, "Ground-truth perfusion by compartment", xs, zs, [truth[r + 1] for r in rois];
        colormap = :turbo, colorrange = (0, 12), label = "mL/min/g", yreversed = false
    )
    panel!(
        fig, 1, 3, "Idealized 70 keV, lung window, t = $t s", xs, zs, [look.hu[l + 1] for l in labels];
        colormap = :grays, colorrange = (-1000, 200), label = "HU", yreversed = false
    )
    panel!(
        fig, 1, 4, "Iodine, t = $t s", xs, zs, [look.iodine[l + 1] for l in labels];
        colormap = :turbo, colorrange = (0, 8), label = "mg/mL", yreversed = false
    )
    Mke.Label(fig[0, :], "coronal y = $(round(point[2]; digits = 1)) mm"; fontsize = 20, font = :bold)
    for column in 1:4
        Mke.colsize!(fig.layout, column, Mke.Relative(1 / 4))
    end
    return fig
end

# ╔═╡ 8d571f7f-8a17-4191-914e-be250d44f00e
coronal_figure(HEALTHY_DIR, 14)

# ╔═╡ 74e87670-dc56-477e-af66-4667a0502bc6
"Each compartment's iodine, with the tracking arrival and peak, and the lungs' iodine ledger."
function curves_figure(view)
    bundle = view.bundle
    fig = Mke.Figure(size = (1400, 520))
    ax = Mke.Axis(fig[1, 1]; title = "Compartments", xlabel = "time after injection (s)", ylabel = "iodine (mg)")
    for (k, name) in enumerate(ROI_NAMES)
        row = bundle.roi_table[("rest", name)]
        Mke.lines!(
            ax,
            bundle.times,
            row.mass;
            linewidth = 2,
            color = roi_colormap().colors[k + 1],
            label = @sprintf("%s · %.2f mL/min/g", name, row.mbf),
        )
    end
    timing = V.DynamicPhantom.scan_times(bundle; state = "rest", roi = "roi:main pulmonary artery")
    Mke.vlines!(ax, [timing.arrival, timing.peak_time]; color = :gray40, linestyle = :dash)
    Mke.xlims!(ax, 0, 40)
    Mke.axislegend(ax; position = :rt, labelsize = 11)

    ledger = bundle.ledger["rest"]
    ax = Mke.Axis(fig[1, 2]; title = "Iodine ledger", xlabel = "time after injection (s)", ylabel = "iodine (mg)")
    entries = [
        ("delivered_mg", :black, "delivered", :solid),
        ("held_mg", :darkorange, "held", :solid),
        ("outflow_mg", :royalblue, "out through the pulmonary veins", :solid),
        ("painted_mg", :red, "painted", :dash),
    ]
    for (column, colour, name, style) in entries
        Mke.lines!(ax, ledger["t_s"], ledger[column]; color = colour, linewidth = 2, linestyle = style, label = name)
    end
    Mke.xlims!(ax, 0, 40)
    Mke.axislegend(ax; position = :lt)
    return fig
end

# ╔═╡ 8c221ce6-b5dd-4184-9271-da2f65e50f84
"""
The share of the delivered iodine that has already left the lungs through the pulmonary veins, 1, 2 and 3 s after arrival
and at the arterial peak (bolus tracking in the main pulmonary artery).
"""
function outflow_check(view)
    bundle = view.bundle
    ledger = bundle.ledger["rest"]
    function share_out(t)
        ts = ledger["t_s"]
        k = clamp(searchsortedlast(ts, t), 1, length(ts) - 1)
        w = (t - ts[k]) / (ts[k + 1] - ts[k])
        out = (1 - w) * ledger["outflow_mg"][k] + w * ledger["outflow_mg"][k + 1]
        delivered = (1 - w) * ledger["delivered_mg"][k] + w * ledger["delivered_mg"][k + 1]
        return out / max(delivered, eps())
    end
    timing = V.DynamicPhantom.scan_times(bundle; state = "rest", roi = "roi:main pulmonary artery")
    shares = [@sprintf("%.0f%%", 100 * share_out(t)) for t in [timing.arrival .+ [1, 2, 3]; timing.peak_time]]
    text = """
    !!! note "Iodine already out of the lungs"
        | | arrival + 1 s | arrival + 2 s | arrival + 3 s | arterial peak |
        |---|---|---|---|---|
        | arrival at $(round(timing.arrival; digits = 1)) s, peak at $(round(timing.peak_time; digits = 1)) s | $(join(shares, " | ")) |

        A measurement that counts the iodine the lungs have taken up misses this share.
    """
    return Markdown.parse(text)
end

# ╔═╡ a8018905-b9ec-4432-b672-40407e3d54d6
curves_figure(healthy)

# ╔═╡ 1a1c84fa-19f7-4091-92f1-846bfa063dd9
outflow_check(healthy)

# ╔═╡ 6510bc2c-7878-427b-9c90-b3fe8fe0884b
md"""
## 9. A Left Lower Lobe Artery Occluded

The same lungs with one artery of the left lower lobe blocked, as a balloon blocks it in the lab's swine studies: the largest imaged artery whose terminals lie in the left lower lobe and that feeds at most a third of it. The solver needs a finite resistance, so the blocked lumen keeps 0.1% of its area; it passes about 0.3% of its healthy flow. The right ventricle keeps the cardiac output, as it does after an acute embolism, so the pulmonary artery pressure rises until the smaller bed carries it. Saved as its own scenario.
"""

# ╔═╡ 9d1f4c6b-5a7e-4b9d-8c3f-4e6a8b0d2f35
"""
The occluded artery of `tree`: the imaged segment that feeds the most of lobe `lobe` (`lobes` on `grid`) among those whose
terminals lie ≥ 95% in it and that feed at most a third of it; with `circuit`'s narrowing to `area` of its lumen.
"""
function occluded_artery(tree, lobes, grid; lobe, area)
    n = length(tree.r)
    lobe_at(p) = (k = round.(Int, (p .- grid.origin) ./ grid.h) .+ 1; all(1 .<= k .<= grid.dims) ? Int(lobes[k...]) : 0)
    children = [Int[] for _ in 1:n]
    foreach(j -> tree.parent[j] > 0 && push!(children[tree.parent[j]], j), 1:n)
    depth = [
        (
            k = j; d = 0; while tree.parent[k] != 0
                k = tree.parent[k]; d += 1
            end; d
        ) for j in 1:n
    ]
    counts = zeros(Int, 5, n)
    for j in 1:n
        isempty(children[j]) && (l = lobe_at(tree.b[:, j]); l > 0 && (counts[l, j] += 1))
    end
    for j in sortperm(depth; rev = true)
        p = tree.parent[j]
        p > 0 && (counts[:, p] .+= counts[:, j])
    end
    total = sum(counts[lobe, j] for j in 1:n if isempty(children[j]))
    candidates = [
        j for j in 1:n if tree.isbase[j] && sum(counts[:, j]) > 0 &&
            counts[lobe, j] >= 0.95 * sum(counts[:, j]) && counts[lobe, j] <= total / 3
    ]
    j = candidates[argmax([counts[lobe, c] for c in candidates])]
    return (segment = j, share = counts[lobe, j] / total, stenosis = [(j, area, 10.0)])
end

# ╔═╡ 2d09b1bf-4ba6-4e97-98d4-575d2bb8c69f
"""
The occluded artery: of the imaged arterial segments whose terminals lie ≥ 95% in the left lower lobe and that feed at most
a third of it, the one feeding the most; its lumen closed to `area` of its own over its length (at most 10 mm), as `circuit`
takes a narrowing.
"""
OCCLUSION = occluded_artery(arteries, LOBES, growth_grid; lobe = 5, area = 1.0e-3)

# ╔═╡ 72e11792-5070-47b2-87d5-08c51895c710
"""
The occluded lungs' circuit. The right ventricle keeps the cardiac output: the inlet pressure is raised (by the secant rule)
until the lungs carry the healthy flow again.
"""
flows_occluded = cached("flow_occluded") do
    target = V.organ_flow(flows.rest)
    solve(p_in) = V.Flow.circuit(
        arteries,
        flows.terr;
        state = :rest,
        veins,
        vleaf = flows.vleaf,
        reference = flows.rest,
        phys = V.FlowSpec(flows.rest.phys; p_in = Float32(p_in)),
        stenosis = OCCLUSION.stenosis,
        to = DEVICE,
    )
    p0 = Float64(flows.rest.phys.p_in)
    c0 = solve(p0)
    p1 = p0 + 1
    c1 = solve(p1)
    for _ in 1:12
        f0 = V.organ_flow(c0) - target
        f1 = V.organ_flow(c1) - target
        abs(f1) <= 1.0e-4 * target && break
        p0, c0, p1 = p1, c1, p1 - f1 * (p1 - p0) / (f1 - f0)
        c1 = solve(p1)
    end
    (rest = c1, stress = c1, p_in = p1)
end;

# ╔═╡ 73bb4a14-b592-4da3-a025-a7acb725aa58
TRANSPORT_OCCLUDED = transport_cache("transport_occluded", flows_occluded.rest);

# ╔═╡ 41071a00-8b18-4c1c-9e9e-aba7a2cbd257
let fig = Mke.Figure(size = (1500, 560))
    ax = vessel_axis(fig[1, 1]; title = "Arteries ≥ 1 mm, the occluded one in gold")
    draw_lung_shell!(ax)
    draw_vessels!(ax, arteries; colormap = ARTERY_COLOURS, min_radius = 1.0)
    j = OCCLUSION.segment
    Mke.linesegments!(ax, [Mke.Point3f(arteries.a[:, j]) => Mke.Point3f(arteries.b[:, j])]; color = :gold, linewidth = 12)

    panels = (("Healthy", flows.rest, 2), ("Occluded", flows_occluded.rest, 3))
    for (title, circuit, column) in panels
        ax = vessel_axis(fig[1, column]; title = "$title · lung blood flow, seen from the front")
        surface = Mke.scatter!(
            ax,
            LUNG_SURFACE.points[1:2:end];
            color = voxel_perfusion(circuit)[LUNG_SURFACE.voxels[1:2:end]],
            colormap = :turbo,
            colorrange = (0, 4),
            markersize = 3,
            strokewidth = 0,
        )
        column == 3 && Mke.Colorbar(fig[1, 4], surface; label = "mL/min per mL of lung")
    end
    Mke.Label(
        fig[0, :],
        @sprintf(
            "The occluded artery feeds %.0f%% of the left lower lobe · pulmonary artery pressure %.1f → %.1f mmHg",
            100 * OCCLUSION.share,
            flows.rest.phys.p_in,
            flows_occluded.p_in,
        );
        fontsize = 18,
        font = :bold,
    )
    fig
end

# ╔═╡ 71283248-3d0e-490e-97ff-6fed714ea769
OCCLUDED_DIR = save_phantom("lll_occlusion", TRANSPORT_OCCLUDED);

# ╔═╡ 4ac5be16-921c-42c6-8d79-8d6445946ab2
occluded = inspect(OCCLUDED_DIR);

# ╔═╡ b44c76e6-18a2-465c-b293-722318aacc60
md"""
* time $(@bind occluded_time Slider(0:1:40; default = 10, show_value = true)) s
* z $(@bind occluded_z Slider(round.(occluded.z0 .+ (0:10:(size(occluded.bundle.labels, 3) - 1)) .* XCAT_VOXEL_MM; digits = 1); default = 270.0, show_value = true)) mm
"""

# ╔═╡ 5a62b1f2-d47f-400d-a669-a79c6cbcfd57
slice_figure(occluded, occluded_time, occluded_z)

# ╔═╡ ce3898c4-2b0e-48b8-a765-ee5fe8b69b5f
coronal_figure(OCCLUDED_DIR, 14)

# ╔═╡ d5172879-75ed-4750-87be-4920faa3718f
curves_figure(occluded)

# ╔═╡ 2d80a57b-7c8e-4d38-a196-29df385aeb09
outflow_check(occluded)

# ╔═╡ ba7b7bc0-6d9f-4292-9c59-47ca6ba3d66c
let fig = Mke.Figure(size = (1000, 460))
    names = [ROI_NAMES; [g for (g, _) in ROI_GROUPS]]
    ax = Mke.Axis(
        fig[1, 1];
        title = "Ground-truth perfusion",
        ylabel = "mL/min/g",
        xticks = (eachindex(names), names),
        xticklabelrotation = π / 5,
    )
    healthy_truth = [healthy.bundle.roi_table[("rest", name)].mbf for name in names]
    occluded_truth = [occluded.bundle.roi_table[("rest", name)].mbf for name in names]
    Mke.barplot!(ax, eachindex(names) .- 0.2, healthy_truth; width = 0.38, color = :gray60, label = "healthy")
    Mke.barplot!(ax, eachindex(names) .+ 0.2, occluded_truth; width = 0.38, color = :firebrick, label = "left lower lobe artery occluded")
    Mke.axislegend(ax; position = :lt)
    fig
end

# ╔═╡ 30746d3d-3c69-42aa-a58d-743161d21444
md"""
## 10. Next Steps: Scanning a Phantom with BasisSimulator.jl

A standalone script: the saved phantom scanned at any time `t` on the validated GE Revolution Apex Elite model (BasisSimulator.jl 0.15+), here at the bolus arrival in the main pulmonary artery and at its peak. The detector covers 160 mm, so one rotation centred on the hila covers most of the lungs (their apices and bases lie outside it); a 160 mm slab of the full chest takes about 17 GB while it is scanned. Ground truth for any compartment comes straight from the phantom.

```julia
import Vasculature as V
import BasisSimulator as BS
import GPUSelect

const GPU = GPUSelect.Storage()
const DIR = expanduser("~/xcat-perfused/pulmonary_pipeline/healthy") # or lll_occlusion
const STATE = "rest" # the lungs have one state

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

# the STANDARD kernel: an apodized ramp fitted to the clinical one, on a 40 cm field over both
# lungs (512², 0.78 mm)
const FOV_CM, SCAN_FOV_CM = 40.0, 50.0
const KNOTS = Tuple(range(0.0, 1.0; length = 11))
const WINDOW = (fc = 0.7 * (10FOV_CM / 512) / (350 / 512), p = 1.5)
const STANDARD = BS.CustomFilter(
    KNOTS,
    Tuple(round(exp(-(x / WINDOW.fc)^WINDOW.p); digits = 5) for x in KNOTS),
)

# ---- one axial rotation of the phantom slab at time t: project, then water BHC, FDK and HU
# on the lung field
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
    # the lung field: every source and detector position moved by −(field centre −
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
# the tracking point (x, y, z mm): the main pulmonary artery
track = meta.meta["tracking"]["main pulmonary artery"]["point_mm"]
# arrival (+80 HU at the point) and the arterial peak
timing = V.DynamicPhantom.scan_times(meta; state = STATE, roi = "roi:main pulmonary artery")
# two scans: 160 mm about the tracking point, 0.625 mm slices, 300 mA, at the arrival and at the peak
grid = meta.meta["grid"]
# the phantom's centre (x, y), on the isocentre
centre = grid["origin_mm"][1:2] .+ (grid["dims"][1:2] .- 1) .* grid["voxel_mm"] ./ 2
slices(zlo, zhi) = (
    (round(Int, (zlo - grid["origin_mm"][3]) / grid["voxel_mm"]) + 1):
    round(Int, (zhi - grid["origin_mm"][3]) / grid["voxel_mm"])
)
const LUNG_Z = (track[3] - 80, track[3] + 80) # mm: the 160 mm beam about the hila
lungs = V.DynamicPhantom.load(
    DIR;
    zrange = max(1, first(slices(LUNG_Z[1] - 15, 0))):last(slices(0, LUNG_Z[2] + 15)),
)
# the 40 cm field on both lungs (x, y mm of their middle)
field_centre_mm = (136.0, 244.0) .- Tuple(centre)
early, early_dose = scan(
    lungs,
    timing.arrival;
    mA = 300,
    collimation_mm = 160.0,
    slice_mm = 0.625,
    seed = 7,
    field_centre_mm,
)
peak, peak_dose = scan(
    lungs,
    timing.peak_time;
    mA = 300,
    collimation_mm = 160.0,
    slice_mm = 0.625,
    seed = 8,
    field_centre_mm,
)

# ---- ground truth: any compartment or group is a sum over its labels
# .mbf (mL/min/g), .mass (mg of iodine every second), .mass_g, .voxels
truth = lungs.roi_table[(STATE, "lungs")]
```

`V.DynamicPhantom.hu(bundle, state, t)` gives each label's idealized 70 keV HU, the scanner-free reference.
"""

# ╔═╡ Cell order:
# ╟─d0311ce2-5aca-446b-bf46-44ea40dd8921
# ╟─21cabb67-46bc-4458-b215-d2484b257de6
# ╠═53ef55e4-f46a-45b5-8f84-ecc59711a13d
# ╠═ef4704c9-3438-48e1-9f42-2f97b9f28acc
# ╠═4152495c-72ba-4321-8ee2-4e4754ff2797
# ╠═a9210218-f94e-45a4-bde7-c4f548a4eae6
# ╠═07c920f4-5c4a-4653-af94-284f616f8256
# ╠═cf8017b9-c58c-465b-8945-a13acd593609
# ╠═90ea862c-5c67-49ab-920a-d5acb3244766
# ╠═5fd8df67-aca2-4daa-ab7d-3a3620dff99a
# ╠═5f1b74ca-efbd-4370-bda2-2a499cf7be7c
# ╠═58d7db22-8641-41c2-aed3-1a23e3bd6883
# ╠═2e61946c-82b1-402f-8ef7-e26f22dbd5b2
# ╠═bdf025f0-bc0b-436f-9581-23318ebe26fa
# ╠═7d6dc0e9-93bb-4e12-96f9-21efbcbf3d38
# ╠═c69a5e75-f4ef-424f-8876-f26ad3501180
# ╠═440ec9de-6961-4ed1-9bf4-43e0dbfea35c
# ╟─dddaae8d-18cb-4213-949c-f72ea7d8d717
# ╠═26bbaa0c-82be-415f-a934-7e86734b719d
# ╠═b3aecd39-f506-4538-8fcd-4bdcb22f595c
# ╟─2fe11641-05cb-4803-9c8d-a30328bc6235
# ╠═c5f15420-4d45-4054-873f-282a655f9442
# ╠═aeb894e0-eea5-416d-a0c1-04f02b797b7d
# ╟─44aca054-f174-4b36-aed6-30b5fd4fb6df
# ╟─4d68e042-f2da-44b9-88f2-d86ed350d6c7
# ╠═f2c25bd3-fd90-4742-9996-067b2f16b04c
# ╠═abb33992-3cb7-41ff-859f-f55c8d6d9d3b
# ╟─8cbd6bf0-9455-4aa5-8fde-d21cea7b7111
# ╠═27f9d235-e736-4ce8-b7f5-68e20ca83d73
# ╠═c4b8f772-02aa-4cac-9a1b-d4d4540585f1
# ╠═03d82855-b9d6-43e5-96a2-f869776c5463
# ╠═1f2278a1-9326-417c-b472-405539f41c50
# ╠═6de4ed07-39b0-440a-b56f-c3ca02ea2d88
# ╠═e79d69d3-719f-429f-9346-4b5448003243
# ╠═6e84d94d-7e21-494b-8cab-e0cb9119b4de
# ╠═95aca34c-b6bd-442c-a95c-8a2cee41971e
# ╠═80f43c52-f003-4912-be1f-c485974abad6
# ╠═9def37df-4fe1-4716-bb70-1845a0433a67
# ╠═5a8c1f3e-2d4b-4e6a-9f7c-1b3d5e7a9c02
# ╠═d4b22c25-1df4-4b6f-bf71-188d026fbfec
# ╠═7fd7dcf9-fe92-4d4b-86e1-0ca90317628b
# ╠═f1bd3790-c7e3-43d4-9998-aa4bdb669a40
# ╠═da727f7d-e412-46a1-9154-823facbf1f98
# ╠═a362e6a5-c0e2-412b-b96d-9c13609d4d73
# ╟─cfdb51b7-71df-4da8-a11d-41b8ee1964fc
# ╟─0fb6d8c3-c1c2-4f6d-96c3-1a47a00458a5
# ╟─8960d456-a118-462d-be93-db4b884d85d6
# ╟─611b7d1b-991d-45eb-987e-c390b6498a04
# ╠═2bf2e225-251c-44c2-a310-fa151e5b4971
# ╠═05ecf689-b662-41ba-b3fd-2b39e3060317
# ╠═f909bf22-0953-4b41-9822-036f7042ccca
# ╠═5a45700e-ddd3-40d4-8504-4712118904e9
# ╠═daa1c63b-b988-4d38-ad53-c22630d11af2
# ╠═d5e4c316-fadc-4dbe-8f51-471929992f45
# ╟─f75080e9-6836-4f3e-839d-cbc53c637a93
# ╟─b44fc616-883b-47b7-a8d1-4d60879cc96d
# ╟─116e3f65-4c0f-4edb-a35f-a6abcf8697a5
# ╠═27a52b36-0fa9-4f86-bcd0-4cf5824daecb
# ╠═55a43d63-8f38-4962-b428-9c0f28e814c7
# ╠═ab72236e-a88a-4eb6-accf-8d67de2ac41f
# ╠═a40e8ed4-a6f0-4de7-8a9f-c0002a81dab2
# ╟─32162101-95fb-4558-a8d7-57590ad05179
# ╟─7ba30cb9-e856-4e5f-a50b-4c14a0c517cd
# ╟─56f6f8b4-669c-4bb8-9e3e-b693c342d67c
# ╟─d5da5d67-513f-49d9-a70b-2bfead33c5fd
# ╠═6817b779-c89b-47ce-b07e-ad48831b78cb
# ╠═45a3f241-cc85-4202-ad0e-00de0342f5fd
# ╠═f6016974-f880-4c69-90fc-37a14b4ff4ce
# ╠═e5b80139-a62d-4497-843b-77d4e1d2a6cf
# ╠═5e4a06c3-bed3-4a8e-9550-dd57bb4a33de
# ╟─9f8d5635-baf7-43f2-ba01-6d006b3d6a77
# ╟─b37414bb-2941-4a47-85c3-0c10c246cc0c
# ╠═df85dc89-82d9-4d4d-9665-8282155840fb
# ╠═4e169ac1-210c-49d8-8f5c-d7c62084fb18
# ╠═31928ce9-9fa0-402a-ba97-cfc37d5898e5
# ╠═6b9d2a4f-3e5c-4f7b-8a1d-2c4e6f8b0d13
# ╠═4b166f97-4217-4a97-95b2-edab3ff8aa46
# ╠═08f3aa33-2562-424a-8526-349c90fdd559
# ╠═b83e74d9-3d5e-4059-b8e6-0e6307b119bf
# ╟─ac7c95a2-969d-41f6-bf40-b299016a01ea
# ╠═e54d8f12-2cb3-4c9b-a1ef-a68e4ad5ec3a
# ╠═aeeda182-58d5-45f5-b269-bb9bb53a54d3
# ╠═75243d5b-b7b8-4b71-95c4-8648072543b4
# ╠═7d1e3a52-4b6c-4f0e-9a21-3c5e8f7b9d14
# ╠═9a5e0d2b-6717-4c81-9684-d2a0cb94718e
# ╠═01e67f8a-e865-40a9-979b-51da0a3e3f5d
# ╠═8811c92e-4cb3-4c0c-add3-6e7da17660a0
# ╠═8c0e3b5a-4f6d-4a8c-9b2e-3d5f7a9c1e24
# ╠═15d45948-dcea-4ac6-a552-7c23db9a5988
# ╠═e75df31a-d1ec-4690-8b58-f8177b34886c
# ╟─9a8d5350-20ea-4e6a-820a-0e4d8a3045d1
# ╠═69f6db9d-70cf-4c42-8ad4-684b798cdb3e
# ╠═0c29d1f3-3e97-49c3-9ccf-2ed17a205c4a
# ╠═70f5505d-1073-4b93-82b7-3c7674b657fd
# ╠═5fe0c8fa-fe29-4114-b0d1-3076eb730ede
# ╠═c238aa74-eace-4853-8b12-80a69ae1cd3a
# ╠═d99468dc-2991-44f5-becc-b8edf15f5778
# ╟─74ad599e-a3f3-423a-850b-ed84438f8dfe
# ╠═aa982a9c-9ff8-4e46-97d8-65c13a957f3c
# ╠═2abf3ae2-4e6a-4e90-9f2e-df9a68f29a56
# ╠═8d571f7f-8a17-4191-914e-be250d44f00e
# ╠═74e87670-dc56-477e-af66-4667a0502bc6
# ╠═8c221ce6-b5dd-4184-9271-da2f65e50f84
# ╠═a8018905-b9ec-4432-b672-40407e3d54d6
# ╠═1a1c84fa-19f7-4091-92f1-846bfa063dd9
# ╟─6510bc2c-7878-427b-9c90-b3fe8fe0884b
# ╠═2d09b1bf-4ba6-4e97-98d4-575d2bb8c69f
# ╠═9d1f4c6b-5a7e-4b9d-8c3f-4e6a8b0d2f35
# ╠═72e11792-5070-47b2-87d5-08c51895c710
# ╠═73bb4a14-b592-4da3-a025-a7acb725aa58
# ╟─41071a00-8b18-4c1c-9e9e-aba7a2cbd257
# ╠═71283248-3d0e-490e-97ff-6fed714ea769
# ╠═4ac5be16-921c-42c6-8d79-8d6445946ab2
# ╟─b44c76e6-18a2-465c-b293-722318aacc60
# ╠═5a62b1f2-d47f-400d-a669-a79c6cbcfd57
# ╠═ce3898c4-2b0e-48b8-a765-ee5fe8b69b5f
# ╠═d5172879-75ed-4750-87be-4920faa3718f
# ╠═2d80a57b-7c8e-4d38-a196-29df385aeb09
# ╟─ba7b7bc0-6d9f-4292-9c59-47ca6ba3d66c
# ╟─30746d3d-3c69-42aa-a58d-743161d21444
