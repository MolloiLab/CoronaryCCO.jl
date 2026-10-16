# SPDX-License-Identifier: GPL-3.0-or-later
# XCAT heart → growth inputs: myocardial domain (LV + RV), coronary base forest (LCA, RCA),
# territory split, and writers for the engine's `.vol` / GXL formats.

"""
    Centerline

Polyline `X` (3×K, mm) with radius `r` (K, mm), sampled from a tube surface, running from the
vessel's proximal end to its distal end.
"""
struct Centerline
    name::String
    X::Matrix{Float64}
    r::Vector{Float64}
end

"""
    centerline(s::NurbsSurface; step = 1.0) -> Centerline

Centerline of an XCAT tube surface: ring centroids and mean ring radii along v (u runs around
the tube), resampled every `step` mm. The collapsed cap rings at both ends are dropped.
"""
function centerline(s::NurbsSurface; step::Real = 1.0, nring::Int = 400, naround::Int = 32)
    tv = range(s.v[4], s.v[end - 3], length = nring); tu = range(s.u[4], s.u[end - 3], length = naround + 1)[1:naround]
    C = zeros(3, nring); R = zeros(nring)
    for (b, y) in enumerate(tv)
        ring = [evaluate(s, x, y) for x in tu]
        c = ntuple(d -> sum(p[d] for p in ring) / naround, 3)
        C[:, b] .= c
        R[b] = sum(sqrt(sum((p[d] - c[d])^2 for d in 1:3)) for p in ring) / naround
    end
    # drop the closed caps: rings narrower than half the vessel's median radius at either end
    med = sort(R)[(end + 1) ÷ 2]
    a = findfirst(>(0.5 * med), R); b = findlast(>(0.5 * med), R)
    C = C[:, a:b]; R = R[a:b]
    # arc-length resampling
    sarc = [0.0; cumsum([sqrt(sum(abs2, C[:, k + 1] .- C[:, k])) for k in 1:(size(C, 2) - 1)])]
    n = max(2, round(Int, sarc[end] / step) + 1)
    ts = range(0, sarc[end], length = n)
    X = zeros(3, n); r = zeros(n)
    for (q, t) in enumerate(ts)
        k = clamp(searchsortedlast(sarc, t), 1, length(sarc) - 1)
        w = (t - sarc[k]) / max(sarc[k + 1] - sarc[k], eps())
        X[:, q] .= (1 - w) .* C[:, k] .+ w .* C[:, k + 1]; r[q] = (1 - w) * R[k] + w * R[k + 1]
    end
    return Centerline(s.name, X, r)
end

"""
    BaseTree

Binary base tree for one coronary root: node 1 is the ostium; `parent[k]` (0 for the ostium),
positions `X` (3×n, mm), radii `r` and the vessel each node came from.
"""
struct BaseTree
    name::String
    parent::Vector{Int}
    X::Matrix{Float64}
    r::Vector{Float64}
    vessel::Vector{String}
end

nnodes(t::BaseTree) = length(t.parent)

"""
    base_forest(lines; roots, attach_tol = 3.0) -> Vector{BaseTree}

Joins centerlines into one binary tree per root. `roots` names the vessels that start at an
ostium (e.g. `["dias_lad1", "dias_rca1"]`). Every other vessel attaches at the node nearest to
its proximal end on the closest already-placed vessel (within `attach_tol` mm of that vessel's
surface); its proximal end is the end nearer to that tree. A node may carry at most two
children, so a third branch moves to the next node downstream.
"""
function base_forest(lines::Vector{Centerline}; roots::Vector{String}, attach_tol::Real = 3.0)
    byname = Dict(l.name => l for l in lines)
    forest = BaseTree[]
    pending = [l for l in lines if !(l.name in roots)]
    trees = Dict{String, Any}()
    for rn in roots
        l = byname[rn]
        t = (parent = collect(0:(size(l.X, 2) - 1)), X = copy(l.X), r = copy(l.r), vessel = fill(l.name, size(l.X, 2)), nchild = [fill(1, size(l.X, 2) - 1); 0])
        trees[rn] = t
    end
    # attach pending vessels repeatedly (a branch may hang off another branch)
    while !isempty(pending)
        best = (Inf, 0, "", 0, false)
        for (pi, l) in enumerate(pending), (rn, t) in trees
            for flip in (false, true)
                p0 = flip ? l.X[:, end] : l.X[:, 1]
                d2 = [sum(abs2, t.X[:, k] .- p0) for k in 1:size(t.X, 2)]
                k = argmin(d2); d = sqrt(d2[k]) - t.r[k]
                d < best[1] && (best = (d, pi, rn, k, flip))
            end
        end
        d, pi, rn, k, flip = best
        d <= attach_tol || error("base_forest: $(pending[pi].name) is $(round(d, digits = 2)) mm from any tree")
        l = pending[pi]; t = trees[rn]
        # a branch starting near its parent's distal end continues from that end node
        e = k
        while (c = findfirst(==(e), t.parent)) !== nothing && t.vessel[c] == t.vessel[k]
            e = c
        end
        if e != k && t.nchild[e] == 0
            p0 = flip ? l.X[:, end] : l.X[:, 1]
            sqrt(sum(abs2, t.X[:, e] .- p0)) <= sqrt(sum(abs2, t.X[:, k] .- p0)) + 2 * t.r[k] && (k = e)
        end
        while t.nchild[k] >= 2
            k = findfirst(==(k), t.parent)      # move downstream along the main child
        end
        X = flip ? l.X[:, end:-1:1] : l.X; r = flip ? reverse(l.r) : l.r
        n0 = size(t.X, 2); m = size(X, 2)
        # skip the branch's first point if it sits on the parent node
        first = sum(abs2, X[:, 1] .- t.X[:, k]) < 0.25 ? 2 : 1
        newpar = [k; (n0 + 1):(n0 + m - first)]
        t2 = (
            parent = [t.parent; newpar], X = hcat(t.X, X[:, first:end]), r = [t.r; r[first:end]],
            vessel = [t.vessel; fill(l.name, m - first + 1)], nchild = [t.nchild; fill(1, m - first); 0],
        )
        t2.nchild[k] += 1
        trees[rn] = t2
        deleteat!(pending, pi)
    end
    for rn in roots
        t = trees[rn]
        push!(forest, BaseTree(rn, t.parent, t.X, t.r, t.vessel))
    end
    return forest
end

"""
    write_base_gxl(path, t::BaseTree)

Writes the base tree in the GXL layout `import_gxl` reads: node 0 is the ostium ("root
node"), leaves are "terminal node", and each edge carries the child node's radius.
"""
function write_base_gxl(path::AbstractString, t::BaseTree)
    n = nnodes(t)
    nch = zeros(Int, n); foreach(p -> p > 0 && (nch[p] += 1), t.parent)
    f(x) = @sprintf("%.6f", x)
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" ?>\n<gxl>\n  <graph id=\"n0_root_first_edge_n0_n1\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">")
        println(io, "    <info_graph>\n      <attr name=\" pPerf\">\n        <float>13300</float>\n      </attr>\n      <attr name=\" pTerm\">\n        <float>8400</float>\n      </attr>\n    </info_graph>")
        for k in 1:n
            kind = t.parent[k] == 0 ? "root node" : nch[k] == 0 ? "terminal node" : "bifurication"
            println(io, "    <node id=\"n$(k - 1)\">\n      <attr name=\" nodeType\">\n        <string> $kind </string>\n      </attr>")
            println(io, "      <attr name=\" position\">\n        <tup>\n          <float>$(f(t.X[1, k]))</float>\n          <float>$(f(t.X[2, k]))</float>\n          <float>$(f(t.X[3, k]))</float>\n        </tup>\n      </attr>\n    </node>")
        end
        for k in 2:n
            println(io, "    <edge id=\"e$(k - 2)\" to=\"n$(k - 1)\" from=\"n$(t.parent[k] - 1)\">\n      <attr name=\" radius\">\n        <float>$(f(t.r[k]))</float>\n      </attr>\n    </edge>")
        end
        println(io, "  </graph>\n</gxl>")
    end
    return path
end

"""
    write_vol(path, mask::AbstractArray{Bool,3}, lower)

DGtal `.vol` (Version 2, uncompressed) with inside voxels 255, so voxel `(i, j, k)` is the
integer point `lower .+ (i, j, k) .- 1` (mm when the voxel size is 1 mm).
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
    heart_masks(S, g; phase = "dias", res = (160, 160)) -> NamedTuple of Array{Bool,3}

LV and RV myocardium and cavities of an XCAT heart on grid `g`:
LV myocardium = inside `lv_4` minus inside `lv_0`, plus the papillary muscles;
RV myocardium = inside `rvm` minus `rv` minus the LV;
LV cavity = inside `lv_0` minus papillary muscles; RV cavity = inside `rv` minus the LV.
"""
function heart_masks(S::Dict{String, NurbsSurface}, g::VoxelGrid; phase::AbstractString = "dias", res = (160, 160))
    vox(n) = voxelize(tessellate(S["$(phase)_$n"]; nu = res[1], nv = res[2])..., g)
    lv4 = vox("lv_4"); lv0 = vox("lv_0"); rv = vox("rv"); rvm = vox("rvm")
    pap = vox("lpap1") .| vox("lpap2") .| vox("upap")
    lvmyo = (lv4 .& .!lv0) .| (pap .& lv0)
    rvmyo = rvm .& .!rv .& .!lv4
    return (; lvmyo, rvmyo, lvcav = lv0 .& .!pap, rvcav = rv .& .!lv4)
end

"""
    label_masks(L, origin, h, g; roi = trues(g.dims), lv = 15, rv = 16, lvpool = 19, rvpool = 20,
                pools = (19, 20, 21, 22, 24, 25)) -> NamedTuple of Array{Bool,3}

Heart masks on the coarse grid `g` from a fine XCAT label volume `L` (voxel size `h` mm, voxel
(1,1,1) centred at `origin`). Each coarse voxel looks at the fine voxels whose centres fall
inside it:
- `lvmyo` / `rvmyo`: at least half of them are LV / RV myocardium (the larger count wins);
- `lvcav` / `rvcav`: at least half are LV / RV blood pool;
- `pool`: any of them is blood (the four chamber pools, great arteries and veins by default);
- `domain`: myocardium with no blood at all, so a vessel grown inside it cannot enter a chamber.
XCAT reuses some label values for distant organs (e.g. 16 for the kidneys), so `roi` limits the
muscle and cavity masks to the heart region; `pool` covers the whole grid. The default label values are those of the XCAT heart phantom.
"""
function label_masks(
        L::AbstractArray{UInt8, 3}, origin, h::Real, g::VoxelGrid; roi::AbstractArray{Bool, 3} = trues(g.dims),
        lv = 15, rv = 16, lvpool = 19, rvpool = 20, pools = (19, 20, 21, 22, 24, 25)
    )
    lvmyo = falses(g.dims); rvmyo = falses(g.dims); lvcav = falses(g.dims); rvcav = falses(g.dims); pool = falses(g.dims)
    isp = falses(256); for p in pools
        isp[p + 1] = true
    end
    Threads.@threads for k in 1:g.dims[3]
        for j in 1:g.dims[2], i in 1:g.dims[1]
            c = center(g, i, j, k)
            r = ntuple(d -> (ceil(Int, (c[d] - g.h / 2 - origin[d]) / h - 1.0e-9) + 1):(ceil(Int, (c[d] + g.h / 2 - origin[d]) / h - 1.0e-9)), 3)
            n = 0; nl = 0; nr = 0; npl = 0; npr = 0; anyp = false
            for kk in r[3], jj in r[2], ii in r[1]
                checkbounds(Bool, L, ii, jj, kk) || continue
                l = L[ii, jj, kk]; n += 1
                nl += l == lv; nr += l == rv; npl += l == lvpool; npr += l == rvpool; anyp |= isp[l + 1]
            end
            n == 0 && continue
            pool[i, j, k] = anyp            # blood anywhere on the grid (the aorta leaves the pericardium)
            roi[i, j, k] || continue
            if 2(nl + nr) >= n
                nl >= nr ? (lvmyo[i, j, k] = true) : (rvmyo[i, j, k] = true)
            end
            2npl >= n && (lvcav[i, j, k] = true)
            2npr >= n && (rvcav[i, j, k] = true)
        end
    end
    return (; lvmyo, rvmyo, lvcav, rvcav, pool, domain = (lvmyo .| rvmyo) .& .!pool)
end

"""
    vessel_tube(t::BaseTree, g; margin = 1.5) -> Array{Bool,3}

Voxels of `g` within `r + margin` mm of the base tree's centerline. Epicardial arteries run in
the fat outside the myocardium; adding this tube to a tree's growth domain lets new segments
connect to them without leaving the domain.
"""
function vessel_tube(t::BaseTree, g::VoxelGrid; margin::Real = 1.5)
    M = zeros(Bool, g.dims)
    for k in 2:nnodes(t)
        a = t.X[:, t.parent[k]]; b = t.X[:, k]; R = max(t.r[k], t.r[t.parent[k]]) + margin
        lo = [max(1, floor(Int, (min(a[d], b[d]) - R - g.origin[d]) / g.h) + 1) for d in 1:3]
        hi = [min(g.dims[d], ceil(Int, (max(a[d], b[d]) + R - g.origin[d]) / g.h) + 1) for d in 1:3]
        ab = b .- a; L2 = max(sum(abs2, ab), eps())
        for kk in lo[3]:hi[3], j in lo[2]:hi[2], i in lo[1]:hi[1]
            p = center(g, i, j, kk)
            w = clamp(sum((p[d] - a[d]) * ab[d] for d in 1:3) / L2, 0.0, 1.0)
            sum((p[d] - a[d] - w * ab[d])^2 for d in 1:3) <= R^2 && (M[i, j, kk] = true)
        end
    end
    return M
end

# Marks the capsule of radius R around segment a→b on grid g.
function _capsule!(M::AbstractArray{Bool, 3}, g::VoxelGrid, a, b, R)
    lo = [max(1, floor(Int, (min(a[d], b[d]) - R - g.origin[d]) / g.h) + 1) for d in 1:3]
    hi = [min(g.dims[d], ceil(Int, (max(a[d], b[d]) + R - g.origin[d]) / g.h) + 1) for d in 1:3]
    ab = b .- a; L2 = max(sum(abs2, ab), eps())
    for k in lo[3]:hi[3], j in lo[2]:hi[2], i in lo[1]:hi[1]
        p = center(g, i, j, k)
        w = clamp(sum((p[d] - a[d]) * ab[d] for d in 1:3) / L2, 0.0, 1.0)
        sum((p[d] - a[d] - w * ab[d])^2 for d in 1:3) <= R^2 && (M[i, j, k] = true)
    end
    return M
end

"""
    vessel_bridges(t::BaseTree, myo, g; every = 2.0, reach = 15.0, margin = 1.5) -> Array{Bool,3}

Epicardial arteries run in the fat a few mm outside the myocardium. Every `every` mm along
the centerline, a capsule (radius r + `margin`) joins the vessel to the nearest myocardial
voxel within `reach` mm, so new branches can pass from the vessel into the wall, as real
penetrating branches do, without the growth domain having to include the whole fat layer.
"""
function vessel_bridges(t::BaseTree, myo::AbstractArray{Bool, 3}, g::VoxelGrid; every::Real = 2.0, reach::Real = 15.0, margin::Real = 1.5)
    M = zeros(Bool, g.dims)
    nreach = ceil(Int, reach / g.h)
    sarc = 0.0; last = -Inf
    for k in 1:nnodes(t)
        k > 1 && (sarc += sqrt(sum(abs2, t.X[:, k] .- t.X[:, t.parent[k]])))
        sarc - last < every && k > 1 && continue
        last = sarc
        p = t.X[:, k]
        c = [round(Int, (p[d] - g.origin[d]) / g.h) + 1 for d in 1:3]
        best = (Inf, c)
        for kk in max(1, c[3] - nreach):min(g.dims[3], c[3] + nreach), j in max(1, c[2] - nreach):min(g.dims[2], c[2] + nreach),
                i in max(1, c[1] - nreach):min(g.dims[1], c[1] + nreach)
            myo[i, j, kk] || continue
            q = center(g, i, j, kk)
            d2 = sum((q[d] - p[d])^2 for d in 1:3)
            d2 < best[1] && (best = (d2, [i, j, kk]))
        end
        best[1] <= reach^2 || continue
        q = collect(center(g, best[2]...))
        # push one more mm into the wall so the bridge overlaps the myocardium
        dir = (q .- p) ./ max(sqrt(best[1]), eps())
        _capsule!(M, g, p, q .+ dir, t.r[k] + margin)
    end
    return M
end
