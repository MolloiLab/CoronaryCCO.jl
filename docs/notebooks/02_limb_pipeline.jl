### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 329820de-86c9-4483-8daf-4cdb31dad119
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    # a writable font cache for CairoMakie (a notebook worker may have no HOME)
    ENV["XDG_CACHE_HOME"] = get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))
end;

# ╔═╡ b3d8dbc1-fce8-44ba-b4f4-062be2f1769e
using PlutoUI: TableOfContents

# ╔═╡ 7cd96038-c409-4a12-8367-f8f1c28c7ae8
using Printf: @sprintf

# ╔═╡ ed211fb3-e82c-4a06-ab74-8b8cf919e773
using Serialization: serialize, deserialize

# ╔═╡ 63e77379-a0ca-4b7f-9bda-a244fb1896be
using Statistics: mean, median, quantile

# ╔═╡ 2d88882b-057d-4f1d-a810-03f14c411fae
md"""
# Limb Pipeline

**From the raw XCAT leg to a vascular tree that reaches 10 µm terminal arterioles, grown in two passes.**

The leg is built the way the heart (`00_cardiac_pipeline`) and the lungs (`01_pulmonary_pipeline`) are, with one difference of scale: a leg needs about 60 million terminals, far more than one growth can hold. So it grows in two passes. The first grows the whole leg to 10 mL patches; the second grows each patch to 10 µm on its own, and the patches run in parallel.

This page runs the stages that take seconds or minutes: the geometry, the first pass, and the second pass on three sample patches, which time the full second pass. The full second pass, the flow solve, the contrast transport and the phantom follow the same steps as the heart and lungs, and are laid out at the end.

| Stage | What happens | Run here |
|---|---|---|
| **1. Load** | XCAT's leg: its skin, bones and imaged femoral artery and vein, as NURBS surfaces | yes |
| **2. Prepare** | the soft tissue to grow into, the bones as blockers, the imaged trees with penetrating stubs | yes |
| **3. Terminal count** | Murray's law fixes the count from the femoral artery and a 10 µm terminal arteriole | yes |
| **4. Pass 1** | the whole leg grown to 10 mL patches, arteries and veins | yes |
| **5. Patches** | each pass-1 terminal's own tissue | yes |
| **6. Pass 2** | three patches grown to 10 µm terminals; the full pass timed | three patches |
| **7. Flow, transport, phantom** | the leg's own physiology, the femoral AIF, the phantom | laid out |
| **8. Design choices** | what carries over from the heart and lungs | — |
| **9. Completing the build** | the remaining steps, in order | — |

Every heavy step is cached in `CACHE_DIR`, so a rerun takes seconds. The XCAT phantom is licensed: its files and anything made from them stay on lab storage.
"""

# ╔═╡ 4a662e6f-8a3c-4219-a3a6-03a7c848b45e
md"""
## Setup
"""

# ╔═╡ 939015ba-d0c0-4411-9223-131ae2d336e1
# ╠═╡ show_logs = false
import Vasculature as V

# ╔═╡ c91bff6a-765c-4bbd-93bd-b81f892a174e
# ╠═╡ show_logs = false
import CairoMakie as Mke

# ╔═╡ b24bcc05-e7c8-4d75-884e-876a5bf0b4b5
import GPUSelect

# ╔═╡ dcdca1e1-7737-4f0c-b019-17179dcd9530
TableOfContents()

# ╔═╡ 4e477957-3699-4072-99ee-464609f3d3f1
begin
    # licensed: this file and anything made from it stay on lab storage
    const LEG_NRB = get(
        ENV,
        "XCAT_LEG_NRB",
        joinpath(homedir(), "Documents", "dev", "MolloiLab", "xcat-data", "legs", "both_legs_1.nrb"),
    )
    # every cached stage (pass 1, the patches, each pass-2 patch); on a cluster, a shared scratch folder
    const CACHE_DIR = get(ENV, "LIMB_CACHE_DIR", joinpath(homedir(), "xcat-perfused", "limb_pipeline_cache"))
    mkpath(CACHE_DIR)
end;

# ╔═╡ c45cbbd6-6d31-4d1d-9d54-fc5dfb7870fb
# ╠═╡ show_logs = false
"Arrays to this machine's GPU (CUDA, Metal, ROCm, oneAPI), or the CPU if there is none."
const DEVICE = GPUSelect.Storage()

# ╔═╡ 77a6aaf7-fa17-470e-9e01-34bc52fbfd4f
"Run `compute` once and keep its result in `CACHE_DIR/name.jls`; afterwards load it."
function cached(compute, name)
    path = joinpath(CACHE_DIR, name * ".jls")
    isfile(path) && return deserialize(path)

    seconds = @elapsed result = compute()
    serialize(path, result)
    @info "$name: $(round(seconds; digits = 1)) s"
    return result
end

# ╔═╡ 8393831b-051f-4982-a98f-d28277891340
"The six face neighbours of a voxel."
const SIX_STEPS = (
    CartesianIndex(1, 0, 0),
    CartesianIndex(-1, 0, 0),
    CartesianIndex(0, 1, 0),
    CartesianIndex(0, -1, 0),
    CartesianIndex(0, 0, 1),
    CartesianIndex(0, 0, -1),
)

# ╔═╡ fc85e6c1-ab4a-43a8-aa91-91d6f74dba97
md"""
## 1. Load the Raw XCAT Leg

XCAT's legs come as one grouped NURBS file: every surface of both legs, each a set of patches under one name. The right leg is at negative x. This notebook reads its skin (`leg_right` and the foot), its bones, and its imaged arteries (`arteries_rleg`) and veins (`veins_rleg`). XCAT's femoral artery enters at the top of the leg with a radius of about 3.9 mm.

`Vasculature.XCAT.read_nrb` reads XCAT's single-object format, so this grouped one has its own reader below.
"""

# ╔═╡ f7984dcd-8660-4bb2-b99e-433d60b08960
begin
    const ENVELOPE = ["leg_right", "foot_right", "foot_right1", "foot_right2", "foot_right3", "foot_right4", "foot_right5"]
    const BONES = ["r_femur", "r_tibia", "r_fibula", "r_patella", "r_talus", "r_calcaneus", "rpelvis1", "rpelvis2", "rpelvis3"]
    const LONG_BONES = ["r_femur", "r_tibia", "r_fibula"] # the fine blockers
    const ARTERIES = "arteries_rleg"
    const VEINS = "veins_rleg"

    "A surface's object name, without the `__occNNN` a repeated name gets."
    group(name) = replace(name, r"__occ\d+$" => "")
    "Foot bones XCAT does not name by side: the right ones are told apart by x."
    side_bone(name) = startswith(name, "metatarsal") || startswith(name, "phalanx")
    "The surfaces this notebook reads."
    wanted(name) = name in ENVELOPE || name in BONES || name in (ARTERIES, VEINS) || (side_bone(name) && !endswith(name, "_inner"))
end;

# ╔═╡ 0928ebc1-1d0b-4b61-a361-637c87839d35
"""
    read_grouped_nrb(path, keep) -> Vector{NurbsSurface}

XCAT's grouped NRB: per patch a name, an object code line, a patch-scale line, `M :M`, `N :N`, then the knots and control
points of `Vasculature.XCAT.read_nrb`'s format. A name repeated over several patches gets `__occNNN`. Only names `keep`
accepts are parsed.
"""
function read_grouped_nrb(path::AbstractString, keep)
    L = readlines(path)
    n = length(L)
    out = V.XCAT.NurbsSurface[]
    seen = Dict{String, Int}()
    i = 1
    while i <= n - 4
        if isempty(strip(L[i])) || !endswith(rstrip(L[i + 3]), ":M") || !endswith(rstrip(L[i + 4]), ":N")
            i += 1
            continue
        end
        name = String(strip(L[i]))
        M = parse(Int, split(L[i + 3])[1])
        N = parse(Int, split(L[i + 4])[1])
        i += 5
        startswith(strip(L[i]), "U Knot") || error("$name: expected U knots at line $i")
        i += 1
        u = Float64[]
        while !startswith(strip(L[i]), "V Knot")
            push!(u, parse(Float64, strip(L[i])))
            i += 1
        end
        i += 1
        v = Float64[]
        while !startswith(strip(L[i]), "Control")
            push!(v, parse(Float64, strip(L[i])))
            i += 1
        end
        i += 1
        k = get(seen, name, 0) + 1
        seen[name] = k
        if keep(name)
            P = Array{Float64, 3}(undef, 3, N, M)
            for q in 1:(M * N)
                f = split(L[i + q - 1], (',', ' ', '\t'); keepempty = false)
                jj = mod1(q, N)
                ii = cld(q, N)
                for d in 1:3
                    P[d, jj, ii] = parse(Float64, f[d])
                end
            end
            push!(out, V.XCAT.NurbsSurface(k == 1 ? name : @sprintf("%s__occ%03d", name, k), u, v, P))
        end
        i += M * N
    end
    return out
end

# ╔═╡ 5776c7ba-8194-47d8-a28e-046ef1205abb
leg_surfaces = cached("leg_surfaces") do
    S = read_grouped_nrb(LEG_NRB, wanted)
    # the unsided foot bones of the right leg (negative x)
    filter!(s -> !side_bone(group(s.name)) || mean(@view s.P[1, :, :]) < 0, S)
end;

# ╔═╡ 2c4e27bb-61d1-4d72-9dfa-8e496e27c231
let count_of(f) = count(s -> f(group(s.name)), leg_surfaces)
    Markdown.parse(
        """
        | object | patches |
        |---|---|
        | skin (leg and foot) | $(count_of(in(ENVELOPE))) |
        | bones | $(count_of(n -> n in BONES || side_bone(n))) |
        | imaged arteries (`$ARTERIES`) | $(count_of(==(ARTERIES))) |
        | imaged veins (`$VEINS`) | $(count_of(==(VEINS))) |
        """,
    )
end

# ╔═╡ ff0ae1c4-2c23-4fee-b563-57ddb6d3e5ec
md"""
## 2. Prepare

**The tissue.** The leg's soft tissue: XCAT's skin surfaces minus its bones, on a 1 mm grid. It is muscle, fat and skin together, one tissue for growth and for flow (`LIMB_FLOW` in section 7 is a composite of the three). Separate muscle and fat masks (XCAT has the muscles) would let each take its own physiology later.

**The blockers.** No new vessel may cross bone: the femur, tibia and fibula at 0.4 mm, every bone at 1 mm out of the domain.

**The inlet trees.** XCAT's imaged artery and vein patches, as centrelines joined into one tree each, rooted where the femoral vessel enters the leg. Their nodes are spaced 2.8 mm × the organ's size scale, as for the heart and lungs. Every imaged vessel then gets a **penetrating stub** about every spacing where none is imaged (`V.XCAT.penetrating_stubs`), so growth starts along the whole vessel and not only from its few imaged ends. The imaged lumens join the domain as voxels vessels may cross but where no terminal ends (`V.NOTERM`).
"""

# ╔═╡ d57038e3-65cf-49fc-b365-1d3888b8b148
"A NURBS surface's triangle mesh, finer for surfaces with more control points."
mesh(s) = V.XCAT.tessellate(s; nu = clamp(4 * size(s.P, 2), 64, 480), nv = clamp(4 * size(s.P, 3), 64, 480))

# ╔═╡ 82d39e1f-f872-4449-838e-81678f93d077
"""
The leg on a 1 mm grid: its grid, the skin envelope, the bones, the soft tissue (envelope minus bone), and the fine
blockers (femur, tibia, fibula at 0.4 mm) as `(mask, origin, h)`.
"""
LEG = cached("leg") do
    skin = [mesh(s) for s in leg_surfaces if group(s.name) in ENVELOPE]
    bones = [s for s in leg_surfaces if group(s.name) in BONES || side_bone(group(s.name))]
    around = V.grid_around(reduce(hcat, first.(skin)), 1.0; margin = 3)
    grid = V.VoxelGrid(Float64.(round.(Int, around.origin)), 1.0, around.dims)

    envelope = reduce(.|, (V.voxelize(X, T, grid) for (X, T) in skin))
    bone = reduce(.|, (V.voxelize(mesh(s)..., grid) for s in bones))

    long = [mesh(s) for s in bones if group(s.name) in LONG_BONES]
    fine_grid = V.grid_around(reduce(hcat, first.(long)), 0.4; margin = 2)
    fine = reduce(.|, (V.voxelize(X, T, fine_grid) for (X, T) in long))

    (grid = grid, envelope = envelope, bone = bone, soft = envelope .& .!bone, blockers = (fine, fine_grid.origin, 0.4))
end;

# ╔═╡ d3bfa261-2daf-42e5-9bfe-4bc1d4a72f3c
"The soft tissue's volume (mL): the tissue every growth pass and the flow are sized by."
const SOFT_ML = count(LEG.soft) / 1000

# ╔═╡ ca00e0f2-712a-4f9e-8b82-2ba6543e3f18
"Whether the point `p` (mm) lies in the voxels of `mask` on the leg grid."
function inside(mask, p; grid = LEG.grid)
    k = round.(Int, (p .- grid.origin) ./ grid.h) .+ 1
    return all(1 .<= k .<= size(mask)) && mask[k...]
end

# ╔═╡ 2c77ebbd-5cc6-4081-aacf-1f6e10e9a78a
"""
Node spacing of the imaged centrelines (mm): the heart's 2.8 mm × the leg's size scale (its tissue over the 148.1 mL
reference, cube root). A branch leaves its host on the host's centreline, and the collision test exempts only the host,
so nodes closer than the host's final radius would make every new branch hit the next imaged piece.
"""
const STEP = 2.8 * cbrt(SOFT_ML / 148.1)

# ╔═╡ fa99ca18-b2a3-4a33-8cc1-54c84219aced
"A centreline resampled every `step` mm by arc length."
function resample(c::V.XCAT.Centerline, step::Real)
    arc = [0.0; cumsum([sqrt(sum(abs2, c.X[:, k + 1] .- c.X[:, k])) for k in 1:(size(c.X, 2) - 1)])]
    n = max(2, round(Int, arc[end] / step) + 1)
    X = zeros(3, n)
    r = zeros(n)
    for (q, t) in enumerate(range(0, arc[end], length = n))
        k = clamp(searchsortedlast(arc, t), 1, length(arc) - 1)
        w = (t - arc[k]) / max(arc[k + 1] - arc[k], eps())
        X[:, q] .= (1 - w) .* c.X[:, k] .+ w .* c.X[:, k + 1]
        r[q] = (1 - w) * c.r[k] + w * c.r[k + 1]
    end
    return V.XCAT.Centerline(c.name, X, r)
end

# ╔═╡ 6536e61a-f90e-4a09-973b-0f3f15bf2542
"""
    imaged_tree(surfaces) -> (tree, r_root, attach_tol)

The imaged patches of one vessel system joined into one tree: each patch's centreline (u and v swapped where the patch
runs the other way), the highest patch as the root, cut where it enters the leg's skin, resampled every `STEP` mm and
joined with the smallest attachment tolerance that joins them all. `r_root` is the root's mean radius over its first
10 mm.
"""
function imaged_tree(surfaces)
    lines = map(surfaces) do s
        c = V.XCAT.centerline(s; step = 1.0)
        # parametrized the other way round (u along the vessel): swap u and v
        median(c.r) > 15 ? V.XCAT.centerline(V.XCAT.NurbsSurface(s.name, s.v, s.u, permutedims(s.P, (1, 3, 2))); step = 1.0) : c
    end
    top = argmax([maximum(@view l.X[3, :]) for l in lines])
    c = lines[top]
    c.X[3, 1] < c.X[3, end] && (c = V.XCAT.Centerline(c.name, c.X[:, end:-1:1], reverse(c.r)))
    k0 = findfirst(k -> inside(LEG.envelope, c.X[:, k]), axes(c.X, 2))
    lines[top] = V.XCAT.Centerline(c.name, c.X[:, k0:end], c.r[k0:end])
    r_root = mean(lines[top].r[1:min(10, end)])
    pieces = [resample(l, STEP) for l in lines]
    for tol in (3.0, 5.0, 8.0, 12.0, 20.0)
        try
            return (tree = only(V.XCAT.base_forest(pieces; roots = [c.name], attach_tol = tol)), r_root, attach_tol = tol)
        catch e
            (e isa ErrorException && occursin("base_forest", e.msg)) || rethrow()
        end
    end
    return error("the imaged $(group(c.name)) patches do not join into one tree within 20 mm")
end

# ╔═╡ 4fb8cda2-0411-4a9e-a2f5-f4156b25a698
"""
The growth inputs of one vessel system: its imaged tree with penetrating stubs, the stub report, the root radius, and the
lumen voxels outside the tissue that vessels may cross but where no terminal ends.
"""
inlets = cached("inlets") do
    map((arteries = ARTERIES, veins = VEINS)) do name
        base = imaged_tree([s for s in leg_surfaces if group(s.name) == name])
        tree, report = V.XCAT.penetrating_stubs(
            base.tree,
            LEG.grid,
            LEG.soft;
            scale = STEP / 2.8,
            forbid = LEG.bone,
            blockers = (LEG.bone, LEG.grid.origin, LEG.grid.h),
        )
        pass = V.XCAT.vessel_tube(tree, LEG.grid; margin = 1.5) .& .!LEG.soft .& .!LEG.bone
        (tree = tree, report = report, r_root = base.r_root, attach_tol = base.attach_tol, pass = pass)
    end
end;

# ╔═╡ 6e3760ce-8bde-4094-b3ba-d19d5c6cb0cb
let a = inlets.arteries, v = inlets.veins
    Markdown.parse(
        @sprintf(
            """
            | | arteries | veins |
            |---|---|---|
            | root radius (mm) | %.2f | %.2f |
            | imaged nodes | %d | %d |
            | penetrating stubs | %d of %d sites | %d of %d sites |
            | lumen outside the tissue (mL) | %.0f | %.0f |

            Soft tissue %.0f mL (skin envelope %.0f mL, bone %.0f mL). Centrelines resampled every %.1f mm. XCAT's femoral artery is %.1f mm across; its femoral vein (%.1f mm) is narrower than the measured supine 13.6 ± 3.0 mm (Keiler et al. 2019).
            """,
            a.r_root, v.r_root,
            length(a.tree.parent), length(v.tree.parent),
            a.report.stubs, a.report.sites, v.report.stubs, v.report.sites,
            count(a.pass) / 1000, count(v.pass) / 1000,
            SOFT_ML, count(LEG.envelope) / 1000, count(LEG.bone .& LEG.envelope) / 1000, STEP,
            2a.r_root, 2v.r_root,
        ),
    )
end

# ╔═╡ 6193958a-729a-4632-b491-4e2c0ca5edf1
"The leg's front (x–z) and side (y–z) silhouettes, with segments drawn over them: one `linesegments!` per view."
function leg_figure(trees; title = "", size = (900, 900))
    fig = Mke.Figure(; size)
    g = LEG.grid
    xs = [g.origin[d] .+ (0:(g.dims[d] - 1)) .* g.h for d in 1:3]
    for (column, (view, label)) in enumerate(((1, "front (x–z)"), (2, "side (y–z)")))
        ax = Mke.Axis(fig[1, column]; title = label, aspect = Mke.DataAspect(), xlabel = view == 1 ? "x (mm)" : "y (mm)", ylabel = "z (mm)")
        silhouette = dropdims(any(LEG.soft; dims = 3 - view); dims = 3 - view)
        Mke.heatmap!(ax, xs[view], xs[3], silhouette; colormap = [:white, :gray88])
        for (segments, colour, scale) in trees
            a, b, r = segments
            points = Mke.Point2f[]
            for j in eachindex(r)
                push!(points, Mke.Point2f(a[view, j], a[3, j]), Mke.Point2f(b[view, j], b[3, j]))
            end
            widths = repeat(max.(0.4, scale .* r); inner = 2)
            Mke.linesegments!(ax, points; color = colour, linewidth = widths)
        end
        Mke.hidedecorations!(ax; label = false, ticklabels = false, ticks = false)
    end
    Mke.Label(fig[0, :], title; fontsize = 18, font = :bold)
    return fig
end

# ╔═╡ 2cc35552-467a-418b-9ebd-cd4a927e5b61
"An inlet tree's segments as `(a, b, r)`: each node with its parent."
function base_segments(tree)
    k = 2:length(tree.parent)
    return (tree.X[:, tree.parent[k]], tree.X[:, k], tree.r[k])
end

# ╔═╡ 826f47f1-d661-4ae1-83e7-a6a0fe1c1e62
leg_figure(
    [(base_segments(inlets.veins.tree), :royalblue, 0.6), (base_segments(inlets.arteries.tree), :firebrick, 0.6)];
    title = "Imaged arteries (red) and veins (blue), with their penetrating stubs",
)

# ╔═╡ 2a5c666d-73d6-49aa-b3c2-285b61e90e26
md"""
## 3. The Terminal Count

**Murray's law fixes the count.** With the bifurcation exponent γ = 3 (Murray) and every terminal carrying an equal share of the flow, a tree's terminal radius is

```math
r_\text{terminal} = r_\text{root} \, N^{-1/3}.
```

The femoral artery's radius is imaged. A terminal arteriole is about 10 µm in radius. So the leg needs N = (r\_root / 10 µm)³ terminals: about 59 million. The heart's terminal density (957,045 in 148.1 mL) over this leg would give 82 million, the same order: the leg's tissue is less densely supplied than the myocardium.

**Why γ = 3 and not the coronary 2.33.** A leg grown with the coronary exponent has conduits at 0.54–0.76 of XCAT's imaged radii. They take almost all of the arteriovenous pressure drop, and 39% of the capillary beds reverse. With γ = 3 the imaged vessels get 0.87–0.95 of their XCAT radii and no bed reverses. Lowering γ does give small terminals (9 µm at γ = 2.33), but only by starving them of pressure: a 10 mL block of thigh muscle grown that way loses about 800 mmHg at maximal flow, against 8 mmHg at γ = 3 with the full count.

**Why two passes.** 59 million terminals means 118 million segments, too many for one growth on a workstation. Growth is local, though: a terminal's patch of tissue does not depend on how the rest of the leg is drawn. So the leg grows in two passes:
- **Pass 1** grows the whole leg to about 10 mL patches, one terminal each.
- **Pass 2** grows each patch from its pass-1 terminal down to 10 µm, with the patch's share of the 59 million.

Every radius in both passes then follows Murray's law from the terminals it feeds, so the result is one Murray tree from the femoral artery to the terminals. The patches share no vessels, so they run in batches on a GPU or in parallel on a cluster.
"""

# ╔═╡ 413b0149-ca96-4678-956d-5e67b98f412a
begin
    "Radius (mm) of a terminal arteriole, the size growth must reach."
    const TERMINAL_RADIUS_MM = 0.01
    "Murray's exponent: the leg's conduits need it (γ = 2.33 starves the beds)."
    const GAMMA = 3.0
    "The terminal count Murray's law gives for the imaged femoral artery and a 10 µm terminal."
    const N_TOTAL = round(Int, (inlets.arteries.r_root / TERMINAL_RADIUS_MM)^GAMMA)
    "Tissue (mL) one pass-1 terminal feeds: one pass-2 patch."
    const PATCH_ML = 10.0
    "Pass-1 terminals: one per patch."
    const N_PASS1 = round(Int, SOFT_ML / PATCH_ML)
end;

# ╔═╡ 47c099f1-608a-42d3-8ea7-a4cc0b76c4bf
let heart_density = 957_045 / 148.1, per_patch = N_TOTAL / N_PASS1
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | femoral artery radius | %.2f mm |
            | terminal radius | %.0f µm |
            | terminals, Murray (γ = %.0f) | %.1f million |
            | terminals at the heart's density, for comparison | %.1f million |
            | pass 1: patches (%.0f mL each) | %d |
            | pass 2: terminals per patch | about %.0f thousand |
            | pass-1 terminal radius, Murray | %.0f µm |
            """,
            inlets.arteries.r_root,
            1000TERMINAL_RADIUS_MM,
            GAMMA,
            N_TOTAL / 1.0e6,
            heart_density * SOFT_ML / 1.0e6,
            PATCH_ML,
            N_PASS1,
            per_patch / 1000,
            1000 * inlets.arteries.r_root * N_PASS1^(-1 / GAMMA),
        ),
    )
end

# ╔═╡ 13d24d4f-32e5-4c8f-8ff2-9542863323b8
md"""
## 4. Pass 1: The Whole Leg to 10 mL Patches

HCT-CCO's staged growth (the only schedule the package grows with), from the imaged femoral artery and vein with their stubs, into the soft tissue, around the bones. The root radius is fixed at the imaged one; γ = 3. Veins grow the same way from the femoral vein, with the same count.
"""

# ╔═╡ 1357edb4-59f6-4f23-a522-3d67f4eb503e
"""
The leg's growth preset: `HEART_GROWTH`'s generic values with the leg's own tissue, imaged root and Murray's exponent.
A draft: once the leg is built it belongs in `src/organs.jl` (as `LIMB_GROWTH`), with its sources in
`test/reference/`.
"""
limb_growth(nterm, tissue_ml, root_radius) = V.GrowthSpec(
    V.HEART_GROWTH;
    nterm,
    flow_ml_min = nothing,
    tissue_ml,
    root_caps = [Float64(root_radius)],
    root_fixed = true,
    symmetry_floor = Tuple{Float64, Float64}[],
    gamma = GAMMA,
    device = DEVICE,
)

# ╔═╡ 63c202d8-e781-45bf-8efc-340e18870f5b
"""
Grow one tree into `tissue` (a mask on `grid`): `pass` voxels may be crossed but take no terminal. Returns the grown
tree (`V.vtree`) and the seconds it took.
"""
function grow_tree(spec, inlet, tissue, pass, grid, blockers, name)
    dir = mktempdir()
    allow = zeros(UInt16, size(tissue))
    allow[pass] .= V.NOTERM
    domain = V.write_vol(joinpath(dir, "$name.vol"), tissue .| pass, Tuple(round.(Int, grid.origin ./ grid.h)))
    gxl = V.XCAT.write_base_gxl(joinpath(dir, "$name.xml"), inlet)
    seconds = @elapsed engine = V.grow(spec, domain, [gxl]; blockers, territories = allow)
    return (tree = V.vtree(engine), seconds = seconds)
end

# ╔═╡ 1e7050c2-417a-4674-9f47-a887894d24bd
pass1 = cached("pass1") do
    map(inlets) do inlet
        spec = limb_growth(N_PASS1, SOFT_ML, inlet.r_root)
        grow_tree(spec, inlet.tree, LEG.soft, inlet.pass, LEG.grid, LEG.blockers, "pass1")
    end
end;

# ╔═╡ 2d36963b-6887-4e12-a57e-99d254c5c9cb
"The terminal segments of a grown tree: the segments without children."
function terminals(tree)
    has_child = falses(length(tree.r))
    foreach(p -> p > 0 && (has_child[p] = true), tree.parent)
    return findall(!, has_child)
end

# ╔═╡ 849eebe7-c721-44ca-b04b-8b75316b95be
let row(name, grown) = (
        t = grown.tree;
        ends = terminals(t);
        r = 1000 .* t.r[ends];
        @sprintf(
            "| %s | %d | %d | %.0f µm (%.0f–%.0f) | %.0f µm | %.1f s |",
            name, length(t.r), length(ends), median(r), quantile(r, 0.05), quantile(r, 0.95),
            1000 * inlets[name].r_root * length(ends)^(-1 / GAMMA), grown.seconds,
        )
    ),
        rows = [row(name, grown) for (name, grown) in pairs(pass1)]
    Markdown.parse(
        """
        | tree | segments | terminals | terminal radius, median (5–95%) | Murray: r\\_root · N^(-1/3) | growth |
        |---|---|---|---|---|---|
        $(join(rows, "\n"))

        Each pass-1 terminal is the root of one pass-2 patch. Pass 1 sets the tree's geometry; section 6 replaces its radii with Murray's, from the terminals each vessel feeds.
        """,
    )
end

# ╔═╡ dd3a7c54-c5ed-4139-9211-411c81f568f4
let segments(t) = (t.a, t.b, t.r)
    leg_figure(
        [(segments(pass1.veins.tree), (:royalblue, 0.6), 0.6), (segments(pass1.arteries.tree), (:firebrick, 0.8), 0.6)];
        title = "Pass 1: arteries (red) and veins (blue), one terminal per 10 mL patch",
    )
end

# ╔═╡ ff93d210-7eaa-41ca-979d-233217bc534c
md"""
## 5. The Patches

Each pass-1 terminal's patch is the tissue nearer to it, through the tissue, than to any other terminal: a flood from every terminal at once, one voxel step at a time, which goes around bone instead of through it. Staged growth gives every terminal an equal share of the tissue, so the patches come out near 10 mL.
"""

# ╔═╡ a1dd86e8-9921-431a-81c1-cc68ddb7311b
"""
    flood(mask, grid, points) -> (label, voxels, box)

Each voxel of `mask` labelled with the nearest of `points` (mm) through the mask (a breadth-first flood from all of them
at once, six neighbours per step), the voxel count per label, and each label's bounding box (`CartesianIndices`).
"""
function flood(mask, grid, points)
    label = zeros(Int32, size(mask))
    queue = CartesianIndex{3}[]
    for (k, p) in enumerate(eachcol(points))
        c = CartesianIndex(Tuple(clamp.(round.(Int, (p .- grid.origin) ./ grid.h) .+ 1, 1, size(mask))))
        mask[c] && label[c] == 0 || continue
        label[c] = k
        push!(queue, c)
    end
    head = 1
    while head <= length(queue)
        c = queue[head]
        head += 1
        for s in SIX_STEPS
            n = c + s
            checkbounds(Bool, mask, n) && mask[n] && label[n] == 0 || continue
            label[n] = label[c]
            push!(queue, n)
        end
    end
    n = size(points, 2)
    voxels = zeros(Int, n)
    lo = fill(CartesianIndex(size(mask)), n)
    hi = fill(CartesianIndex(1, 1, 1), n)
    for c in CartesianIndices(label)
        k = label[c]
        k == 0 && continue
        voxels[k] += 1
        lo[k] = min(lo[k], c)
        hi[k] = max(hi[k], c)
    end
    return (label = label, voxels = voxels, box = [lo[k]:hi[k] for k in 1:n])
end

# ╔═╡ 0a4968ed-9eec-4358-9e0e-9f3ffb03c846
"The arterial patches: each pass-1 arterial terminal's segment and its tissue."
PATCHES = cached("patches") do
    t = pass1.arteries.tree
    ends = terminals(t)
    f = flood(LEG.soft, LEG.grid, t.b[:, ends])
    (segment = ends, f...)
end;

# ╔═╡ 15144905-f653-42a0-97f5-67f3963ee428
let ml = PATCHES.voxels ./ 1000
    seeded = count(>(0), PATCHES.voxels)
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | patches | %d (%d with tissue) |
            | patch volume, median (5–95%%) | %.1f mL (%.1f–%.1f) |
            | soft tissue in a patch | %.1f%% |
            """,
            length(ml), seeded, median(ml[ml .> 0]), quantile(ml[ml .> 0], 0.05), quantile(ml[ml .> 0], 0.95),
            100 * sum(PATCHES.voxels) / count(LEG.soft),
        ),
    )
end

# ╔═╡ dd46fe50-0a2e-41cd-a507-bceab8d0e381
let g = LEG.grid, j = round(Int, (mean(pass1.arteries.tree.b[2, :]) - g.origin[2]) / g.h) + 1
    fig = Mke.Figure(size = (700, 900))
    ax = Mke.Axis(fig[1, 1]; title = "Patches, coronal slice at y = $(round(Int, g.origin[2] + (j - 1) * g.h)) mm", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    slice = Float32.(PATCHES.label[:, j, :])
    colour = Float32.(mod.(slice .* 0.6180339f0, 1.0f0))
    colour[slice .== 0] .= NaN32
    Mke.heatmap!(ax, g.origin[1] .+ (0:(g.dims[1] - 1)) .* g.h, g.origin[3] .+ (0:(g.dims[3] - 1)) .* g.h, colour; colormap = :tab20, nan_color = :white)
    fig
end

# ╔═╡ bf480e2e-8ebb-4030-81fd-6a5e062bdfd2
md"""
## 6. Pass 2: Patches to 10 µm, Stitched into One Tree

**Every radius follows Murray's law from the terminals it feeds.** Pass 1 gives the tree's geometry and its patches, not its radii. Growth balances each terminal's path resistance, so pass-1 terminals that carry the same flow come out up to a third thinner or thicker than each other. A patch grown from a thin one would end in terminals too thin for their flow. So every vessel takes the radius Murray's law gives for the 10 µm terminals below it, r = 10 µm · n^(1/3): each patch's root from the patch's own count, each pass-1 vessel from the counts of all the patches below it. The femoral artery comes out at its imaged radius, because the total count was chosen from it.

**Each patch grows on its own, then joins the tree.** Its tissue is the patch. Its inlet is the pass-1 terminal segment that feeds it, at its Murray radius. Its count is its share of the 59 million by tissue volume, so every terminal carries the same flow. The segment's lumen, and the neighbouring patches' tissue within 3 mm, may be crossed but take no terminal, so terminals reach the patch's border and neighbouring patches meet without a seam. The grown tree then replaces that pass-1 segment, so the two passes become one tree.

Here every patch under one pass-1 branch in the calf is grown and stitched in: a contiguous piece of the leg, complete to 10 µm. The checks:
- **terminal radius** near 10 µm;
- **pressure lost** from the femoral artery to each terminal at maximal flow (every vessel dilated, each terminal an equal share of the leg's flow, Poiseuille with the Pries–Secomb in vivo viscosity). The capillaries must keep a pressure above the femoral vein's 11 mmHg.
"""

# ╔═╡ ac61ef6f-c424-4d8f-90fa-41906373d29d
"The segments of a tree (its `parent` vector) in breadth-first order from the root, and each segment's children."
function walk(parent)
    children = [Int[] for _ in eachindex(parent)]
    order = Int[]
    for j in eachindex(parent)
        parent[j] == 0 ? push!(order, j) : push!(children[parent[j]], j)
    end
    for j in order # `order` grows as it is read
        append!(order, children[j])
    end
    return order, children
end

# ╔═╡ ca33ac06-932e-4f4f-a34d-d39d651aefce
"""
    downstream(parent, weight) -> (order, children, n)

The terminals below each segment: a terminal `j` counts `weight(j)` (the final terminals it stands for), every other
segment the sum over its children.
"""
function downstream(parent, weight)
    order, children = walk(parent)
    n = zeros(length(parent))
    for j in reverse(order)
        n[j] = isempty(children[j]) ? weight(j) : sum(n[c] for c in children[j])
    end
    return order, children, n
end

# ╔═╡ 5b216de6-c728-457d-a182-01c990235112
"Each patch's final terminal count: its share of `N_TOTAL` by tissue volume, so every terminal carries the same flow."
const PATCH_TERMINALS = round.(Int, N_TOTAL .* PATCHES.voxels ./ sum(PATCHES.voxels))

# ╔═╡ 6e378f0f-405d-48fe-ac44-36b1621c638f
"""
The pass-1 arteries at their Murray radii: each segment's final terminal count below it (`below`) and its radius
`TERMINAL_RADIUS_MM · below^(1/γ)` (mm). A pass-1 terminal counts its patch's terminals (at least one).
"""
PASS1_MURRAY = let t = pass1.arteries.tree, patch_of = Dict(j => k for (k, j) in enumerate(PATCHES.segment))
    _, _, below = downstream(t.parent, j -> max(PATCH_TERMINALS[patch_of[j]], 1))
    (below = below, radius = TERMINAL_RADIUS_MM .* below .^ (1 / GAMMA))
end;

# ╔═╡ 59b8c1c5-e00d-4e74-901c-1ac10ddf5f87
let t = pass1.arteries.tree, ends = PATCHES.segment, root = findfirst(==(0), t.parent)
    ratio = PASS1_MURRAY.radius ./ t.r
    µm(x) = 1000 .* x
    Markdown.parse(
        @sprintf(
            """
            | | value |
            |---|---|
            | femoral artery: imaged / Murray from the counts | %.2f / %.2f mm |
            | pass-1 terminals as grown, median (5–95%%) | %.0f µm (%.0f–%.0f) |
            | pass-1 terminals at their Murray radius, median (5–95%%) | %.0f µm (%.0f–%.0f) |
            | every pass-1 vessel: Murray ÷ grown, median (5–95%%) | %.2f (%.2f–%.2f) |
            """,
            inlets.arteries.r_root, PASS1_MURRAY.radius[root],
            median(µm(t.r[ends])), quantile(µm(t.r[ends]), 0.05), quantile(µm(t.r[ends]), 0.95),
            median(µm(PASS1_MURRAY.radius[ends])), quantile(µm(PASS1_MURRAY.radius[ends]), 0.05), quantile(µm(PASS1_MURRAY.radius[ends]), 0.95),
            median(ratio), quantile(ratio, 0.05), quantile(ratio, 0.95),
        ),
    )
end

# ╔═╡ 641271f3-d904-4d16-893c-848557c02b94
begin
    "Neighbouring tissue (mm) a patch may grow through but not end in, so patches meet without a seam."
    const PATCH_MARGIN_MM = 3.0

    "`mask` grown `steps` voxels through its six face neighbours."
    function dilate(mask, steps)
        out = copy(mask)
        for _ in 1:steps
            grown = copy(out)
            for c in CartesianIndices(out)
                out[c] || continue
                for s in SIX_STEPS
                    n = c + s
                    checkbounds(Bool, out, n) && (grown[n] = true)
                end
            end
            out = grown
        end
        return out
    end
end

# ╔═╡ eedb0097-f6eb-4984-bcd0-0a994f928cd2
"""
    grow_patch(k) -> (tree, seconds, nterm, ml, root_radius, segment)

Pass 2 for patch `k`: the patch cropped to its box (and the pass-1 segment that feeds it), grown from that segment at
its Murray radius with the patch's share of `N_TOTAL` terminals. The neighbouring patches' tissue within
`PATCH_MARGIN_MM` and the segment's lumen may be crossed but take no terminal, so the patch's terminals reach its border
and neighbouring patches meet without a seam.
"""
function grow_patch(k; tree = pass1.arteries.tree)
    j = PATCHES.segment[k]
    a = Float64.(tree.a[:, j])
    b = Float64.(tree.b[:, j])
    r = Float64(PASS1_MURRAY.radius[j])
    g = LEG.grid
    # the patch's box, grown to hold its feeding segment and the margin
    pad = ceil(Int, r / g.h) + ceil(Int, PATCH_MARGIN_MM / g.h) + 2
    ka = round.(Int, (a .- g.origin) ./ g.h) .+ 1
    lo = max.(Tuple(min(first(PATCHES.box[k]), CartesianIndex(Tuple(ka)))) .- pad, 1)
    hi = min.(Tuple(max(last(PATCHES.box[k]), CartesianIndex(Tuple(ka)))) .+ pad, g.dims)
    box = CartesianIndices(Tuple(UnitRange.(lo, hi)))
    grid = V.VoxelGrid(g.origin .+ (Tuple(first(box)) .- 1) .* g.h, g.h, size(box))

    tissue = PATCHES.label[box] .== k
    soft = LEG.soft[box]
    margin = dilate(tissue, round(Int, PATCH_MARGIN_MM / g.h)) .& soft .& .!tissue
    inlet = V.XCAT.BaseTree("patch_$k", [0, 1], hcat(a, b), [r, r], ["pass1", "pass1"])
    pass = (V.XCAT.vessel_tube(inlet, grid; margin = 1.5) .& .!LEG.bone[box] .| margin) .& .!tissue
    ml = count(tissue) / 1000
    nterm = PATCH_TERMINALS[k]
    grown = grow_tree(limb_growth(nterm, ml, r), inlet, tissue, pass, grid, LEG.blockers, "patch_$k")
    return (tree = grown.tree, seconds = grown.seconds, nterm = nterm, ml = ml, root_radius = r, segment = j)
end

# ╔═╡ b57241e2-449d-4e78-b18e-44d5b6bc793d
"""
The region grown here: the pass-1 branch with 8 to 16 patches below it whose end lies nearest 30% of the way up the leg
(the calf), and the patches below it.
"""
REGION = let t = pass1.arteries.tree, z = t.b[3, PATCHES.segment]
    _, _, patches = downstream(t.parent, j -> 1)
    height = minimum(z) + 0.3 * (maximum(z) - minimum(z))
    candidates = findall(n -> 8 <= n <= 16, patches)
    branch = candidates[argmin(abs.(t.b[3, candidates] .- height))]
    function under(j)
        while j != 0 && j != branch
            j = t.parent[j]
        end
        return j == branch
    end
    (branch = branch, patches = [k for (k, j) in enumerate(PATCHES.segment) if under(j) && PATCHES.voxels[k] > 0])
end

# ╔═╡ 78ceb802-8b99-480e-9767-bd24ca40f15b
"Pass 2 for every patch of the region, each cached on its own."
pass2 = [
    cached("patch_$k") do
        grow_patch(k)
    end for k in REGION.patches
];

# ╔═╡ 2bdfc5d1-1351-4b5a-b678-90995faf4dd1
"""
    stitch(tree, radius, grown) -> (parent, a, b, r, pass, from)

One tree from the pass-1 `tree` (at its Murray `radius`) and the pass-2 trees `grown`: each grown tree replaces the
pass-1 terminal segment it grew from and hangs from that segment's parent. `pass` is each segment's pass (1 or 2) and
`from` a pass-1 segment's index in `tree` (0 for a grown one).
"""
function stitch(tree, radius, grown)
    replaced = Set(p.segment for p in grown)
    keep = [j for j in eachindex(tree.r) if !(j in replaced)]
    index = zeros(Int, length(tree.r))
    index[keep] = eachindex(keep)
    parent = Int[p == 0 ? 0 : index[p] for p in tree.parent[keep]]
    a = [tree.a[:, keep]]
    b = [tree.b[:, keep]]
    r = [Float32.(radius[keep])]
    pass = [fill(1, length(keep))]
    from = [keep]
    n = length(keep)
    for p in grown
        t = p.tree
        up = index[tree.parent[p.segment]]
        append!(parent, [q == 0 ? up : q + n for q in t.parent])
        push!(a, t.a)
        push!(b, t.b)
        push!(r, t.r)
        push!(pass, fill(2, length(t.r)))
        push!(from, zeros(Int, length(t.r)))
        n += length(t.r)
    end
    return (
        parent = parent, a = reduce(hcat, a), b = reduce(hcat, b), r = reduce(vcat, r),
        pass = reduce(vcat, pass), from = reduce(vcat, from),
    )
end

# ╔═╡ 6ae0879c-c30e-450d-be22-f831d2ac317d
"The leg's arterial tree with the region complete to 10 µm: pass 1 at its Murray radii, the region's patches grafted in."
stitched = stitch(pass1.arteries.tree, PASS1_MURRAY.radius, pass2);

# ╔═╡ e5c89aca-c566-48d4-9ccf-59bbb4ca9e8b
"The final terminals a terminal of `stitched` stands for: one if it was grown in pass 2, its patch's count if not."
terminal_weight = let patch_of = Dict(j => k for (k, j) in enumerate(PATCHES.segment))
    j -> stitched.pass[j] == 2 ? 1.0 : Float64(max(PATCH_TERMINALS[patch_of[stitched.from[j]]], 1))
end

# ╔═╡ 5b657f30-24f2-4481-88c0-1c23e8732280
let fig = Mke.Figure(size = (1500, 860)), t = stitched, g = LEG.grid
    region = findall(==(2), t.pass)
    ends = let has_child = falses(length(t.r))
        foreach(p -> p > 0 && (has_child[p] = true), t.parent)
        filter(j -> !has_child[j], region)
    end
    # the region's box (its terminals) with 10 mm of surrounding leg, and an 8 mm detail at the region's middle
    box_lo = vec(minimum(t.b[:, ends]; dims = 2)) .- 10
    box_hi = vec(maximum(t.b[:, ends]; dims = 2)) .+ 10
    centre = vec(median(t.b[:, ends]; dims = 2))
    detail_lo = centre .- 4
    detail_hi = centre .+ 4
    # the pass-1 path from the femoral artery to the region's branch
    path = Int[findfirst(==(REGION.branch), t.from)]
    while t.parent[path[end]] != 0
        push!(path, t.parent[path[end]])
    end
    within(j, lo, hi) = lo[1] <= t.b[1, j] <= hi[1] && lo[3] <= t.b[3, j] <= hi[3]
    function draw!(ax, segments; colour = nothing, scale = 1.0)
        isempty(segments) && return
        segments = segments[sortperm(t.r[segments])] # thin first, so thick vessels lie on top
        points = Mke.Point2f[]
        for j in segments
            push!(points, Mke.Point2f(t.a[1, j], t.a[3, j]), Mke.Point2f(t.b[1, j], t.b[3, j]))
        end
        return Mke.linesegments!(
            ax, points; color = colour === nothing ? repeat(log10.(1000 .* t.r[segments]); inner = 2) : colour,
            colormap = :viridis, colorrange = (1, 3.6), linewidth = repeat(max.(0.4, scale .* t.r[segments]); inner = 2),
        )
    end
    xs = [g.origin[d] .+ (0:(g.dims[d] - 1)) .* g.h for d in 1:3]
    silhouette = dropdims(any(LEG.soft; dims = 2); dims = 2)
    outline!(ax, lo, hi) = Mke.lines!(ax, [lo[1], hi[1], hi[1], lo[1], lo[1]], [lo[3], lo[3], hi[3], hi[3], lo[3]]; color = :black, linewidth = 1.5)

    # 1. the leg: pass 1, the path to the region, the region's box
    leg = Mke.Axis(fig[1, 1]; title = "The leg", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    Mke.heatmap!(leg, xs[1], xs[3], silhouette; colormap = [:white, :gray92])
    draw!(leg, findall(==(1), t.pass); colour = :gray55, scale = 2.0)
    draw!(leg, path; colour = :firebrick, scale = 2.0)
    outline!(leg, box_lo, box_hi)

    # 2. the region and 10 mm around it: pass 1 nearby in grey, the region's vessels of 20 µm or more in colour
    zoom = Mke.Axis(fig[1, 2]; title = "The region: $(length(pass2)) patches stitched to pass 1 (vessels ≥ 20 µm)", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    Mke.heatmap!(zoom, xs[1], xs[3], silhouette; colormap = [:white, :gray92])
    draw!(zoom, filter(j -> t.pass[j] == 1 && within(j, box_lo .- 20, box_hi .+ 20), eachindex(t.r)); colour = :gray60, scale = 6.0)
    draw!(zoom, [path; filter(j -> t.r[j] >= 0.02f0, region)]; scale = 6.0)
    outline!(zoom, detail_lo, detail_hi)
    Mke.limits!(zoom, box_lo[1], box_hi[1], box_lo[3], box_hi[3])

    # 3. 8 mm at the region's middle: every vessel, down to the 10 µm terminals
    detail = Mke.Axis(fig[1, 3]; title = "8 mm detail: every vessel, down to 10 µm", aspect = Mke.DataAspect(), xlabel = "x (mm)", ylabel = "z (mm)")
    draw!(detail, filter(j -> within(j, detail_lo .- 2, detail_hi .+ 2), region); scale = 40.0)
    Mke.limits!(detail, detail_lo[1], detail_hi[1], detail_lo[3], detail_hi[3])
    Mke.Colorbar(fig[1, 4]; colormap = :viridis, limits = (1, 3.6), label = "log₁₀ radius (µm)")
    Mke.colsize!(fig.layout, 1, Mke.Relative(0.16))
    Mke.colgap!(fig.layout, 1, 60)
    Mke.colgap!(fig.layout, 2, 60)

    # lines from each box's right corners to the next panel's left corners, in figure pixels (they follow the layout),
    # on a scene made after the axes so they are drawn over them
    overlay = Mke.Scene(fig.scene; clear = false, camera = Mke.campixel!)
    function connect!(from, lo, hi, to)
        points = Mke.lift(from.scene.viewport, from.finallimits, to.scene.viewport) do vf, lf, vt
            pixel(x, z) = Mke.Point2f(vf.origin .+ (Mke.Point2f(x, z) .- lf.origin) ./ lf.widths .* vf.widths)
            [pixel(hi[1], hi[3]), Mke.Point2f(vt.origin[1], vt.origin[2] + vt.widths[2]), pixel(hi[1], lo[3]), Mke.Point2f(vt.origin...)]
        end
        return Mke.linesegments!(overlay, points; color = :gray30, linewidth = 1)
    end
    connect!(leg, box_lo, box_hi, zoom)
    connect!(zoom, detail_lo, detail_hi, detail)
    fig
end

# ╔═╡ 5f97677c-390b-41e6-85ab-ceff75a4c140
# Murray's law across the stitched tree: each segment's radius against the final terminals below it, pass 1 and pass 2
# together; the line is r = 10 µm · n^(1/3)
let fig = Mke.Figure(size = (760, 560)), t = stitched
    _, _, n = downstream(t.parent, terminal_weight)
    ax = Mke.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "terminals downstream", ylabel = "radius (µm)", title = "One Murray tree from the femoral artery to 10 µm")
    grown = findall(==(2), t.pass)
    keep = grown[1:max(1, length(grown) ÷ 40000):end]
    Mke.scatter!(ax, n[keep], 1000 .* t.r[keep]; markersize = 2, color = (:seagreen, 0.3), label = "pass 2 (the region)")
    first_pass = findall(==(1), t.pass)
    Mke.scatter!(ax, n[first_pass], 1000 .* t.r[first_pass]; markersize = 4, color = (:firebrick, 0.6), label = "pass 1")
    line = 10.0 .^ range(0, log10(N_TOTAL), length = 50)
    Mke.lines!(ax, line, 1000TERMINAL_RADIUS_MM .* line .^ (1 / GAMMA); color = :black, label = "r = 10 µm · n^(1/3)")
    Mke.axislegend(ax; position = :lt)
    fig
end

# ╔═╡ cc1cfb89-1c01-4fe9-ba4a-0abe7041e4ec
let seconds_per_terminal = sum(p.seconds for p in pass2) / sum(p.nterm for p in pass2), segments = 2N_TOTAL
    Markdown.parse(
        @sprintf(
            """
            **The full pass 2, from this region.** %.0f µs of growth per terminal, so the %d patches (%.1f million terminals) take about **%.1f hours** one after another on this machine, or about %.0f minutes spread over 40 GPUs. The arterial tree holds about %.0f million segments, %.1f GB at 40 bytes a segment; the veins as many again.
            """,
            1.0e6seconds_per_terminal, N_PASS1, N_TOTAL / 1.0e6, seconds_per_terminal * N_TOTAL / 3600,
            seconds_per_terminal * N_TOTAL / 60 / 40, segments / 1.0e6, 40 * segments / 1.0e9,
        ),
    )
end

# ╔═╡ ba9579ba-1796-47e6-b072-c4466c0e68cb
md"""
## 7. Flow, Transport and the Phantom

### The leg's own physiology

`LIMB_FLOW` takes every value from measurements in human legs where one exists, and derives the rest from cited measurements. The leg is supine and at rest. Its soft tissue is one composite tissue: muscle 0.62, fat 0.34 and skin 0.04 by volume. Fourteen fields differ from the heart's. These set the oxygen and the flow:

- **O₂ demand** `mvo2` 0.00146 mL O₂/min/g, from the measured resting leg flow (2.6 mL/min/100 mL, Koch Esteves et al. 2023) and the femoral venous saturation;
- **O₂ extraction** `oef` 0.315: femoral venous SO₂ 67.1% against arterial 97.9% (Pearson et al. 2011);
- **venous pressure** `p_out` 11 mmHg (common femoral vein, supine);
- **maximal dilation** `r_total_dilated` 300 mmHg·min·g/mL: 90 mmHg over the PET flow under adenosine (muscle 0.40, fat 0.105 mL/min/g; Heinonen et al. 2010, 2014);
- **tone** `d_active` 300 µm: skeletal-muscle feed arteries regulate;
- **blood and exchange** `fb` 0.037, `ve` 0.11, `ps` 0.014 mL/min/g (EDTA-derived), `rho` 1.02 g/mL.

Flow is pressure-driven: growth knows nothing of these values, and the leg's flow comes out of the arterial and venous pressures, the trees and the O₂-driven tone. The resting flow below is what the O₂ balance implies; the solve must find it, not be tuned to it.
"""

# ╔═╡ 51296634-0886-4732-a48f-f1313041b7da
"""
The leg's flow physiology: the fields that differ from the heart's, each with its source. It moves into `src/organs.jl`
as a preset once its sources are in `test/reference/`, as every preset's are. The heart-only fields take the leg's
no-change values: no measured blood-volume rise under adenosine (`fb_dilated = fb`) and a constant PS.
"""
const LIMB_FLOW = V.FlowSpec(
    V.HEART_FLOW;
    p_out = 11.0f0, # mmHg. Common femoral vein, supine rest, catheter (Pearson 2011)
    fb = 0.037f0, # mL/mL. Muscle 3.3 mL/100 g by PET (Raitakari 1995), fat and skin (Leggett 1991), + conduits
    fb_dilated = 0.037f0, # no measured rise under adenosine
    share_art = 0.15f0, # systemic arteries : veins = 13 : 64 (Guyton & Hall), applied to 1 - share_cap
    share_cap = 0.11f0, # muscle capillary volume 0.48 of muscle blood (Gidlöf 1988; Raitakari 1995)
    share_ven = 0.74f0, # the remainder
    cap_length_mm = 0.44f0, # functional capillary path 438 µm (Sarelius 1986)
    cap_diameter_um = 5.3f0, # human quadriceps capillary lumen (Gidlöf 1988)
    ve = 0.11f0, # muscle 0.07, fat 0.14, skin 0.46 by composition (Goodall 2020; Thompson 2005; Levitt 2003)
    ps = 0.014f0, # forearm EDTA PS (Keiding 1988) scaled to iodinated contrast (Levitt 2002)
    ps_slope = 0.0f0,
    rho = 1.02f0, # ICRU-44 muscle, adipose, skin by composition (Hubbell & Seltzer 1995)
    r_total_dilated = 300.0f0, # 90 mmHg / PET flow under adenosine (Heinonen 2010, 2014)
    d_active = 300.0f0, # muscle feed arteries regulate (Williams & Segal 1993; Ives 2012)
    oef = 0.315f0, # supine femoral venous SO2 67.1% vs arterial 97.9% (Pearson 2011)
    mvo2 = 0.00146f0, # 2.6 mL/min/100 mL (Koch Esteves 2023) × O2 content × oef / rho
)

# ╔═╡ 016fdbaa-8dd3-4459-b9a3-f141d21e2b05
"""
    pressure_drop(tree, q, weight) -> (terminals, drop)

Pressure lost (mmHg) from the root to each terminal when the tree carries `q` mL/min, each terminal taking a share in
proportion to `weight(j)`, the final terminals it stands for. Every vessel at its radius (fully dilated), Poiseuille
with the Pries–Secomb in vivo viscosity at the discharge hematocrit.
"""
function pressure_drop(tree, q, weight; phys = LIMB_FLOW)
    order, children, n = downstream(tree.parent, weight)
    unit = q / 60 / n[order[1]] # mL/s per final terminal
    P = zeros(length(tree.r))
    for j in order
        L = Float32(sqrt(sum(abs2, tree.b[:, j] .- tree.a[:, j])))
        r = tree.r[j]
        drop = V.Flow.poiseuille(L, r, phys.plasma_viscosity_mPas * V.Flow.mu_rel(2000 * r, phys.hd)) * unit * n[j]
        P[j] = (tree.parent[j] == 0 ? 0.0 : P[tree.parent[j]]) + drop
    end
    ends = findall(isempty, children)
    return (terminals = ends, drop = P[ends])
end

# ╔═╡ 45930316-b6ba-4992-b30b-ead047135d61
let f = LIMB_FLOW, q_max = LIMB_FLOW.p_in / LIMB_FLOW.r_total_dilated * LIMB_FLOW.rho * SOFT_ML # mL/min, the whole leg dilated
    rows = map(zip(REGION.patches, pass2)) do (k, p)
        r = 1000 .* p.tree.r[terminals(p.tree)]
        @sprintf(
            "| %d | %.1f mL | %.0f µm | %d | %.1f µm (%.1f–%.1f) | %.0f s |",
            k, p.ml, 1000p.root_radius, p.nterm, median(r), quantile(r, 0.05), quantile(r, 0.95), p.seconds,
        )
    end
    loss = pressure_drop(stitched, q_max, terminal_weight)
    grown = loss.drop[stitched.pass[loss.terminals] .== 2]
    Markdown.parse(
        @sprintf(
            """
            | patch | tissue | root radius (Murray) | terminals | terminal radius, median (5–95%%) | growth |
            |---|---|---|---|---|---|
            %s

            The region: %d patches, %.1f mL, %.2f million terminals, %.1f million segments, grown in %.0f s. At maximal flow (%.1f L/min for the whole leg), the pressure lost from the femoral artery to the region's 10 µm terminals is %.0f mmHg (median; 95%%: %.0f), which leaves the capillaries %.0f mmHg of the %.0f mmHg arterial pressure, above the femoral vein's %.0f mmHg.
            """,
            join(rows, "\n"),
            length(pass2), sum(p.ml for p in pass2), sum(p.nterm for p in pass2) / 1.0e6,
            count(==(2), stitched.pass) / 1.0e6, sum(p.seconds for p in pass2),
            q_max / 1000, median(grown), quantile(grown, 0.95), f.p_in - median(grown), f.p_in, f.p_out,
        ),
    )
end

# ╔═╡ 04cb85af-10a2-46b7-b415-f5879a214023
let f = LIMB_FLOW, content = 1.34 * f.hb / 100 * f.sao2 # mL O2 per mL blood
    rest = f.mvo2 / (content * f.oef) # mL/min/g
    maximal = f.p_in / f.r_total_dilated
    mass = f.rho * SOFT_ML
    fp(q) = q * (1 - f.hd)
    Markdown.parse(
        @sprintf(
            """
            | | the O₂ balance implies | measured |
            |---|---|---|
            | resting flow | %.4f mL/min/g, %.0f mL/min for this %.0f g leg | 300 ± 30 mL/min (Chiesa 2016), 353 mL/min (Klein 2003) |
            | leg VO₂ | %.1f mL/min | 18–28 mL/min |
            | maximal flow | %.2f mL/min/g, %.1f L/min | nutritive, under adenosine |
            | flow reserve | %.1f | — |
            | contrast extraction, rest / maximal | %.2f / %.2f | — |
            | mean blood transit at rest, fb / flow | %.0f s | — |

            The resting transit is long: at rest a leg holds iodine for over a minute, so very little leaves between two volume scans a few seconds apart. Under adenosine it is about %.0f s, and venous outflow matters as it did for the heart at stress.
            """,
            rest, rest * mass, mass, f.mvo2 * mass, maximal, maximal * mass / 1000, maximal / rest,
            1 - exp(-f.ps / fp(rest)), 1 - exp(-f.ps / fp(maximal)),
            60 * f.fb / (rest * f.rho), 60 * f.fb / (maximal * f.rho),
        ),
    )
end

# ╔═╡ ccb1198f-120e-49a1-aa97-0de84552442e
md"""
### Transport, compartments and the phantom

These follow the heart and lungs step for step, on the full pass-2 trees.

1. **Veins, pass 2.** The same patches from the femoral vein's pass-1 terminals (`grow_patch` with `pass1.veins.tree` and a venous flood).
2. **Flow.** `V.flow(LIMB_FLOW, arteries, veins, tissue_centres)`, at rest and at maximal dilation. With 118 million segments per tree, each patch is solved on its own and enters the pass-1 tree as one resistance. The patches share no vessels, so this is exact. The checks: resting inflow about 330 mL/min, femoral venous SO₂ 0.66, VO₂ 18–28 mL/min, no reversed beds, tone not saturated.
3. **Transport.** The injection (`INJECTION`, below) through the whole-body model. The leg is the body's right-leg bed, fed from the abdominal aorta and draining to the femoral vein. That bed's share of the cardiac output must be the grown leg's own flow, as `coronary_share` sets it for the heart and `organ_flow` for the lungs, so the body and the leg carry the same blood. `transport` does this for the heart and lungs today; a peripheral bed needs the same.
4. **The AIF.** A tracking point in the imaged **femoral artery**, at the top of the scanned slab. The package gives a tracking point in a vessel its own label at the segment's blood concentration, so its curve is the true arterial input. Put the AIF ROI inside that one segment's lumen: the other imaged segments are painted at their mass over XCAT's lumen, which differs from the grown radius.
5. **Compartments.** The muscles (or thigh, calf and foot), each also grouped with the conduit vessels that feed and drain it, as the heart's "myocardium + epicardial vessels" and the lungs' "lungs + hilar vessels": the ideal first-pass compartment counts the iodine in transit.
6. **The phantom.** XCAT has no 0.2 mm label volume of the leg, and the whole leg at 0.2 mm is 8 billion voxels. Build the phantom for one scanned slab (the 160 mm beam: the thigh, or the calf), voxelizing XCAT's leg surfaces there.
"""

# ╔═╡ f6355eaf-f55a-4d71-99f2-f68b2dd5e88d
"""
The injection, the same as the heart's and lungs' (Hadjiabdolhamid et al. 2024, Shojazadeh et al. 2026): 0.5 mL/kg of
contrast at 5 mL/s, then 0.5 mL/kg of a 30:70 contrast:saline chaser at the same rate, for the body model's 70 kg.
"""
const INJECTION = let kg = 70, conc = 350
    V.Body.Injection(; phases = [(0.5 * kg, 5, conc), (0.5 * kg, 5, 0.3 * conc)])
end

# ╔═╡ c6d2e2bb-6ace-4a49-b0de-ebfe038aab20
md"""
## 8. Design Choices Carried Over from the Heart and Lungs

**One source of truth**
- One code path for every organ: grow → flow → transport → phantom. Organs differ only in their spec values (`GrowthSpec`, `FlowSpec`, `TransportSpec`, `PhantomSpec`). The leg's presets live in this notebook until they join the others in `src/organs.jl`.
- Every organ has its own physiology. The leg takes none of the heart's values where a leg measurement exists: γ = 3, its own O₂ demand, extraction, pressures and blood volumes.
- Flow is pressure-driven. Ground truth comes out of the trees, the pressures and the tone; measured flow is only a comparator. Growth and flow are never tuned to make a perfusion measurement agree.

**Growth**
- HCT-CCO's staged schedule only.
- Imaged vessels get penetrating stubs, so growth starts along their whole length (the heart's and lungs' imaged trees had too few ends without them).
- The imaged lumens are crossable but take no terminal; blockers keep vessels out of bone (the lungs' fissures, the heart's chambers).
- Terminal size follows the count, r\_root · N^(-1/γ). Reach real terminal arterioles by growing more terminals (two passes), never by lowering γ.

**Anatomy and timing**
- Check XCAT against measurements. Its aorta is about three standard deviations narrow, so the aorta's timing was taken from the body model (`V.XCAT.scale_tree`). Its main pulmonary artery is large. Its femoral vein is narrower than measured.
- The AIF must be in time with the compartment's inlet. The ideal compartment counts the iodine still in transit in its conduit vessels (epicardial, hilar; here the femoral and popliteal branches).
- V2 about 1 s before the arterial peak loses the least iodine out of the veins. At rest the leg barely loses any; under adenosine it will.
- One injection for every organ: the lab's 0.5 mL/kg protocol.

**Painting**
- An imaged vessel's radius is set by its flow, not by XCAT's lumen. Its iodine is spread over XCAT's lumen, so its painted concentration can differ from its blood's.
- A tracking point's segment is painted at its true blood concentration. Keep the CT AIF ROI inside that segment.

**Notebooks**
- Loops over notebook globals are slow; they run inside functions.
- Every heavy step is cached.
- Many segments are drawn with one `linesegments!`.
"""

# ╔═╡ 333ce447-64ae-4ec5-898b-8a8d86427091
md"""
## 9. Completing the Build

1. **Pass 2 for every arterial patch**, in batches on a GPU or as parallel jobs. Each patch is independent and cached on its own; the patches are then grafted onto pass 1.
2. **Veins**: the same two passes from the femoral vein.
3. **Terminal size from leg data.** The 10 µm terminal is the size the heart's terminals reach. Skeletal-muscle (and fat) terminal arteriole sizes, or microvascular unit volumes, set `TERMINAL_RADIUS_MM` or the count for the leg.
4. **The leg's presets in the package**: `LIMB_GROWTH` and `LIMB_FLOW` in `src/organs.jl`, their sources in `test/reference/`, with tests.
5. **Flow**: each patch, then the leg, at rest and at maximal dilation, against the checks in section 7.
6. **Transport**: the peripheral bed's share of the cardiac output from the grown leg's flow, then the injection with the femoral AIF.
7. **Compartments and a phantom slab**, then idealized first-pass perfusion against ground truth, as for the heart and lungs.
"""


# ╔═╡ Cell order:
# ╟─2d88882b-057d-4f1d-a810-03f14c411fae
# ╟─4a662e6f-8a3c-4219-a3a6-03a7c848b45e
# ╠═329820de-86c9-4483-8daf-4cdb31dad119
# ╠═b3d8dbc1-fce8-44ba-b4f4-062be2f1769e
# ╠═7cd96038-c409-4a12-8367-f8f1c28c7ae8
# ╠═ed211fb3-e82c-4a06-ab74-8b8cf919e773
# ╠═63e77379-a0ca-4b7f-9bda-a244fb1896be
# ╠═939015ba-d0c0-4411-9223-131ae2d336e1
# ╠═c91bff6a-765c-4bbd-93bd-b81f892a174e
# ╠═b24bcc05-e7c8-4d75-884e-876a5bf0b4b5
# ╠═dcdca1e1-7737-4f0c-b019-17179dcd9530
# ╠═4e477957-3699-4072-99ee-464609f3d3f1
# ╠═c45cbbd6-6d31-4d1d-9d54-fc5dfb7870fb
# ╠═77a6aaf7-fa17-470e-9e01-34bc52fbfd4f
# ╠═8393831b-051f-4982-a98f-d28277891340
# ╟─fc85e6c1-ab4a-43a8-aa91-91d6f74dba97
# ╠═f7984dcd-8660-4bb2-b99e-433d60b08960
# ╠═0928ebc1-1d0b-4b61-a361-637c87839d35
# ╠═5776c7ba-8194-47d8-a28e-046ef1205abb
# ╟─2c4e27bb-61d1-4d72-9dfa-8e496e27c231
# ╟─ff0ae1c4-2c23-4fee-b563-57ddb6d3e5ec
# ╠═d57038e3-65cf-49fc-b365-1d3888b8b148
# ╠═82d39e1f-f872-4449-838e-81678f93d077
# ╠═d3bfa261-2daf-42e5-9bfe-4bc1d4a72f3c
# ╠═ca00e0f2-712a-4f9e-8b82-2ba6543e3f18
# ╠═2c77ebbd-5cc6-4081-aacf-1f6e10e9a78a
# ╠═fa99ca18-b2a3-4a33-8cc1-54c84219aced
# ╠═6536e61a-f90e-4a09-973b-0f3f15bf2542
# ╠═4fb8cda2-0411-4a9e-a2f5-f4156b25a698
# ╟─6e3760ce-8bde-4094-b3ba-d19d5c6cb0cb
# ╠═6193958a-729a-4632-b491-4e2c0ca5edf1
# ╠═2cc35552-467a-418b-9ebd-cd4a927e5b61
# ╟─826f47f1-d661-4ae1-83e7-a6a0fe1c1e62
# ╟─2a5c666d-73d6-49aa-b3c2-285b61e90e26
# ╠═413b0149-ca96-4678-956d-5e67b98f412a
# ╟─47c099f1-608a-42d3-8ea7-a4cc0b76c4bf
# ╟─13d24d4f-32e5-4c8f-8ff2-9542863323b8
# ╠═1357edb4-59f6-4f23-a522-3d67f4eb503e
# ╠═63c202d8-e781-45bf-8efc-340e18870f5b
# ╠═1e7050c2-417a-4674-9f47-a887894d24bd
# ╠═2d36963b-6887-4e12-a57e-99d254c5c9cb
# ╟─849eebe7-c721-44ca-b04b-8b75316b95be
# ╟─dd3a7c54-c5ed-4139-9211-411c81f568f4
# ╟─ff93d210-7eaa-41ca-979d-233217bc534c
# ╠═a1dd86e8-9921-431a-81c1-cc68ddb7311b
# ╠═0a4968ed-9eec-4358-9e0e-9f3ffb03c846
# ╟─15144905-f653-42a0-97f5-67f3963ee428
# ╟─dd46fe50-0a2e-41cd-a507-bceab8d0e381
# ╟─bf480e2e-8ebb-4030-81fd-6a5e062bdfd2
# ╠═ac61ef6f-c424-4d8f-90fa-41906373d29d
# ╠═ca33ac06-932e-4f4f-a34d-d39d651aefce
# ╠═5b216de6-c728-457d-a182-01c990235112
# ╠═6e378f0f-405d-48fe-ac44-36b1621c638f
# ╟─59b8c1c5-e00d-4e74-901c-1ac10ddf5f87
# ╠═eedb0097-f6eb-4984-bcd0-0a994f928cd2
# ╠═641271f3-d904-4d16-893c-848557c02b94
# ╠═b57241e2-449d-4e78-b18e-44d5b6bc793d
# ╠═78ceb802-8b99-480e-9767-bd24ca40f15b
# ╠═2bdfc5d1-1351-4b5a-b678-90995faf4dd1
# ╠═6ae0879c-c30e-450d-be22-f831d2ac317d
# ╠═e5c89aca-c566-48d4-9ccf-59bbb4ca9e8b
# ╠═016fdbaa-8dd3-4459-b9a3-f141d21e2b05
# ╟─45930316-b6ba-4992-b30b-ead047135d61
# ╟─5b657f30-24f2-4481-88c0-1c23e8732280
# ╟─5f97677c-390b-41e6-85ab-ceff75a4c140
# ╟─cc1cfb89-1c01-4fe9-ba4a-0abe7041e4ec
# ╟─ba9579ba-1796-47e6-b072-c4466c0e68cb
# ╠═51296634-0886-4732-a48f-f1313041b7da
# ╟─04cb85af-10a2-46b7-b415-f5879a214023
# ╟─ccb1198f-120e-49a1-aa97-0de84552442e
# ╠═f6355eaf-f55a-4d71-99f2-f68b2dd5e88d
# ╟─c6d2e2bb-6ace-4a49-b0de-ebfe038aab20
# ╟─333ce447-64ae-4ec5-898b-8a8d86427091
