# SPDX-License-Identifier: GPL-3.0-or-later
# What a phantom classifier needs to name an XCAT voxel: the 1 mm tissue voxel it lies in, the imaged segment whose lumen
# it is, and, for a great-vessel voxel, its vessel and its place along it (the blood between it and the heart), so the
# phantom times it by volume over flow (`pool_delay`). Every organ's phantom uses these; the organ differs only in which
# labels are its tissue and its imaged vessels. XCAT label values: 6 lung, 10 liver, 19-22 chamber blood (LV, RV, LA, RA),
# 24 veins, 25 arteries, 28 aorta, 124 and 125 the pulmonary veins and arteries in the lungs (the full chest's).

# the voxel kinds (`VOX_POOL`, ...) are the phantom's, defined after this module: looked up in the parent at call time
const _PHANTOM = parentmodule(@__MODULE__)

"""
    tissue_numbers(mask) -> Array{Int32, 3}

Each voxel's number among the `true` voxels of `mask` (in `findall` order), 0 elsewhere: the tissue voxel numbers
`flow` and `phantom` use.
"""
function tissue_numbers(mask::AbstractArray{Bool, 3})
    numbers = zeros(Int32, size(mask))
    numbers[findall(mask)] .= Int32.(1:count(mask))
    return numbers
end

"""
    nearest_tissue_voxel(numbers, k) -> Int32

The number of the tissue voxel at grid index `k` (a tuple), or of the nearest one within 2 voxels (0: none).
"""
function nearest_tissue_voxel(numbers::AbstractArray{Int32, 3}, k)
    dims = size(numbers)
    all(1 .<= k .<= dims) && numbers[k...] > 0 && return numbers[k...]
    best = typemax(Int); number = Int32(0)
    for di in -2:2, dj in -2:2, dk in -2:2
        n = k .+ (di, dj, dk)
        all(1 .<= n .<= dims) || continue
        numbers[n...] > 0 || continue
        d = di^2 + dj^2 + dk^2
        d < best && (best = d; number = numbers[n...])
    end
    return number
end

"""
    imaged_segments(tree) -> (indices, starts, ends, radii)

The imaged (base) segments of a grown tree (`vtree`): their indices in the tree, end points (3 × n, mm) and radii (mm).
"""
function imaged_segments(tree)
    k = findall(tree.isbase)
    return (k, tree.a[:, k], tree.b[:, k], tree.r[k])
end

"""
    nearest_segment(segments, p) -> (index, distance)

Of `segments` ([`imaged_segments`](@ref)), the one whose surface is nearest to `p` (mm): its index in the tree and the
distance from `p` to its surface (mm; negative inside the lumen).
"""
function nearest_segment((indices, starts, ends, radii), p)
    best = Inf32; nearest = 1
    for k in eachindex(radii)
        axis = ends[:, k] .- starts[:, k]
        offset = Float32.(p) .- starts[:, k]
        length2 = sum(abs2, axis)
        s = length2 > 0 ? clamp(sum(offset .* axis) / length2, 0.0f0, 1.0f0) : 0.0f0
        distance = sqrt(sum(abs2, offset .- s .* axis)) - radii[k]
        distance < best && (best = distance; nearest = k)
    end
    return indices[nearest], best
end

"""
    SegmentIndex

A bucket grid over a tree's [`imaged_segments`](@ref) for fast nearest-lumen queries near the vessels: each cell lists the
segments whose lumen, widened by the index's `reach`, touches it. Build it with [`segment_index`](@ref).
"""
struct SegmentIndex{S}
    segments::S
    origin::NTuple{3, Float64}
    cell::Float64
    dims::NTuple{3, Int}
    first::Vector{Int}   # cell c's segments are members[first[c]:(first[c + 1] - 1)]
    members::Vector{Int32}
end

"""
    segment_index(segments; reach, cell = 4.0) -> SegmentIndex

An index over `segments` ([`imaged_segments`](@ref)) for [`nearest_segment`](@ref) queries at points within `reach` mm of
a lumen; `cell` is the bucket size (mm).
"""
function segment_index(segments; reach::Real, cell::Real = 4.0)
    _, starts, ends, radii = segments
    pad = maximum(radii; init = 0.0f0) + reach
    lo = Tuple(min.(vec(minimum(starts; dims = 2)), vec(minimum(ends; dims = 2))) .- pad .- cell)
    hi = Tuple(max.(vec(maximum(starts; dims = 2)), vec(maximum(ends; dims = 2))) .+ pad .+ cell)
    dims = Tuple(ceil.(Int, (hi .- lo) ./ cell))
    lists = [Int32[] for _ in 1:prod(dims)]
    L = LinearIndices(dims)
    for k in eachindex(radii)
        r = radii[k] + reach
        a = min.(starts[:, k], ends[:, k]) .- r; b = max.(starts[:, k], ends[:, k]) .+ r
        ia = clamp.(floor.(Int, (a .- lo) ./ cell) .+ 1, 1, dims); ib = clamp.(floor.(Int, (b .- lo) ./ cell) .+ 1, 1, dims)
        for z in ia[3]:ib[3], y in ia[2]:ib[2], x in ia[1]:ib[1]
            push!(lists[L[x, y, z]], Int32(k))
        end
    end
    first = cumsum([1; length.(lists)])
    return SegmentIndex(segments, Float64.(lo), Float64(cell), dims, first, reduce(vcat, lists; init = Int32[]))
end

"""
    nearest_segment(index::SegmentIndex, p) -> (index, distance)

As [`nearest_segment`](@ref) over the indexed segments, for a point within the index's `reach` of a lumen; elsewhere
`(0, Inf32)`.
"""
function nearest_segment(ix::SegmentIndex, p)
    indices, starts, ends, radii = ix.segments
    c = floor.(Int, (p .- ix.origin) ./ ix.cell) .+ 1
    all(1 .<= c .<= ix.dims) || return (0, Inf32)
    l = LinearIndices(ix.dims)[c...]
    best = Inf32; nearest = 0
    for m in ix.first[l]:(ix.first[l + 1] - 1)
        k = ix.members[m]
        axis = ends[:, k] .- starts[:, k]
        offset = Float32.(p) .- starts[:, k]
        length2 = sum(abs2, axis)
        s = length2 > 0 ? clamp(sum(offset .* axis) / length2, 0.0f0, 1.0f0) : 0.0f0
        distance = sqrt(sum(abs2, offset .- s .* axis)) - radii[k]
        distance < best && (best = distance; nearest = Int(k))
    end
    return nearest == 0 ? (0, Inf32) : (indices[nearest], best)
end

"""
    closing(mask, r) -> Array{Bool, 3}

Morphological closing of a 3-D mask by a ball of `r` voxels: dilate, then erode, each from an exact squared Euclidean
distance transform (Felzenszwalb & Huttenlocher, one line at a time).
"""
function closing(mask::AbstractArray{Bool, 3}, r::Real)
    dilated = _squared_distance(mask) .<= r^2
    return .!(_squared_distance(.!dilated) .<= r^2)
end

function _squared_distance(inside::AbstractArray{Bool, 3})
    dims = size(inside)
    distance = Array{Float32}(undef, dims)
    for axis in 1:3
        others = [dims[d] for d in 1:3 if d != axis]
        m = dims[axis]
        AK.foreachindex(1:prod(others); min_elems = 1) do line
            f = Vector{Float64}(undef, m); v = Vector{Int}(undef, m); z = Vector{Float64}(undef, m + 1)
            a = (line - 1) % others[1] + 1; b = (line - 1) ÷ others[1] + 1
            index(u) = axis == 1 ? (u, a, b) : axis == 2 ? (a, u, b) : (a, b, u)
            for u in 1:m
                f[u] = axis == 1 ? (inside[index(u)...] ? 0.0 : 1.0e20) : distance[index(u)...]
            end
            # the lower envelope of the parabolas rooted at every sample
            k = 1; v[1] = 1; z[1] = -Inf; z[2] = Inf
            for u in 2:m
                s = ((f[u] + u^2) - (f[v[k]] + v[k]^2)) / (2 * (u - v[k]))
                while s <= z[k]
                    k -= 1
                    s = ((f[u] + u^2) - (f[v[k]] + v[k]^2)) / (2 * (u - v[k]))
                end
                k += 1; v[k] = u; z[k] = s; z[k + 1] = Inf
            end
            k = 1
            for u in 1:m
                while z[k + 1] < u
                    k += 1
                end
                distance[index(u)...] = (u - v[k])^2 + f[v[k]]
            end
        end
    end
    return distance
end

"""
    great_vessel_transit(labels, origin, h; cell_mm = 1.0, svc_mm = 70) -> NamedTuple

Which great vessel each vessel voxel belongs to, and where it lies along it, on a `cell_mm` grid: its tree (0 off the
vessels) and the blood (mL) between it and the heart. `labels` is XCAT's label volume with voxel (1, 1, 1) centred at
`origin` (mm), spacing `h`. The trees grow breadth first through XCAT's vessel labels from where they meet the heart, by
connection rather than by mesh:

1. the aorta (label 28) and the systemic arteries joined to it, from the LV blood;
2. the pulmonary arteries (label 25 and those in the lungs, 125), from the RV blood, never stepping beside the aorta, so
   they cannot run into it where the two touch; grown first, so the aortic tree takes the rest;
3. the superior vena cava, from the upper half of the RA's veins: its first `svc_mm` (the SVC itself, about 70 mm long
   [mahlon2007_svc_ct_length]) is tree 3, and the veins beyond it are the right (6) or left (7) brachiocephalic tree, by
   side (the liver's is the right);
4. the inferior vena cava, from the lower half of the RA's veins;
5. the pulmonary veins (label 24 and those in the lungs, 124), from the LA blood.

The venous trees grow at once, each voxel joining the first to reach it. A voxel's volume is its tree's volume nearer the
heart: a cross-section of the pulmonary or venous trees cuts every parallel branch, which together carry the whole flow, so
all of them count; in the aortic tree only the aorta counts, since its branches run beside it, not in its path. Past the
arch's top the descending aorta carries its own, smaller flow, so it is timed from the arch's volume (`apex_ml`) on; a
brachiocephalic voxel is timed from the SVC's top (`confluence_ml`). Returns `(; origin_mm, h, tree, volume, descending,
apex_ml, confluence_ml, tree_ml, ostia_ml)`; `ostia_ml` is the aorta's volume where the coronary arteries (label 26) leave it.
"""
function great_vessel_transit(labels::AbstractArray{UInt8, 3}, origin, h::Real; cell_mm::Real = 1.0, svc_mm::Real = 70)
    step = round(Int, cell_mm / h)
    coarse = labels[1:step:end, 1:step:end, 1:step:end]
    hc = step * h
    sides = (CartesianIndex(1, 0, 0), CartesianIndex(-1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, -1, 0), CartesianIndex(0, 0, 1), CartesianIndex(0, 0, -1))
    touches(c, label) = any(s -> checkbounds(Bool, coarse, c + s) && coarse[c + s] == label, sides)
    tree = zeros(UInt8, size(coarse)); distance = fill(Int32(-1), size(coarse))
    order = Dict(t => CartesianIndex{3}[] for t in 0x01:0x07) # each tree's cells, nearest the heart first
    # grows trees from their seeds at once, breadth first: a cell joins the first tree to reach it
    function grow!(seeds, allowed)
        queue = CartesianIndex{3}[]
        for (c, t) in seeds
            tree[c] == 0 || continue
            tree[c] = t; distance[c] = 0; push!(queue, c)
        end
        head = 1
        while head <= length(queue)
            c = queue[head]; head += 1
            for s in sides
                n = c + s
                checkbounds(Bool, coarse, n) && tree[n] == 0 && allowed(n, tree[c]) || continue
                tree[n] = tree[c]; distance[n] = distance[c] + 1; push!(queue, n)
            end
        end
        foreach(c -> push!(order[tree[c]], c), queue)
        return nothing
    end
    seeds(label, blood, t) = [(c, t) for c in CartesianIndices(coarse) if coarse[c] == label && touches(c, blood)]
    grow!(seeds(0x19, 0x14, 0x02), (n, t) -> coarse[n] in (0x19, 0x7d) && !touches(n, 0x1c))
    grow!(seeds(0x1c, 0x13, 0x01), (n, t) -> coarse[n] in (0x1c, 0x19))
    caval = seeds(0x18, 0x16, 0x00)
    z_middle = sort([c[3] for (c, _) in caval])[(length(caval) + 1) ÷ 2]
    grow!(
        [[(c, c[3] >= z_middle ? 0x03 : 0x04) for (c, _) in caval]; seeds(0x18, 0x15, 0x05)],
        (n, t) -> coarse[n] == 0x18 || (t == 0x05 && coarse[n] == 0x7c),
    )
    # volumes nearer the heart, tree by tree
    cell_ml = hc^3 / 1000
    volume = fill(NaN32, size(coarse))
    for (t, cells) in order
        counted(c) = t != 0x01 || coarse[c] == 0x1c
        nearer = 0; before = 0
        for (k, c) in enumerate(cells)
            k > 1 && distance[c] > distance[cells[k - 1]] && (before = nearer)
            volume[c] = before * cell_ml
            counted(c) && (nearer += 1)
        end
    end
    # the descending aorta: past the arch's top
    aorta = [c for c in order[0x01] if coarse[c] == 0x1c]
    apex = aorta[argmax([c[3] for c in aorta])]
    descending = falses(size(coarse))
    for c in aorta
        descending[c] = distance[c] > distance[apex]
    end
    # the SVC's first svc_mm, then the brachiocephalic trees by side
    trunk = [c for c in order[0x03] if distance[c] * hc <= svc_mm]
    confluence_ml = maximum(volume[c] for c in trunk)
    svc_x = sum(c[1] for c in trunk) / length(trunk)
    liver = [c[1] for c in CartesianIndices(coarse) if coarse[c] == 0x0a]
    right_is_low_x = sum(liver) / length(liver) < svc_x
    for c in order[0x03]
        distance[c] * hc <= svc_mm && continue
        tree[c] = (c[1] < svc_x) == right_is_low_x ? 0x06 : 0x07
    end
    return (
        origin_mm = Tuple(Float64.(origin)), h = hc, tree = tree, volume = volume, descending = descending,
        apex_ml = Float32(volume[apex]), confluence_ml = Float32(confluence_ml),
        tree_ml = Dict(t => count(==(t), tree) * cell_ml for t in 0x01:0x07),
        ostia_ml = sort([volume[c] for c in aorta if touches(c, 0x1a)]),
    )
end

"""
    scale_tree(transit, tree, k) -> NamedTuple

`transit` ([`great_vessel_transit`](@ref)) with the blood along tree `tree` scaled by `k`: each of its voxels keeps its place
along the vessel, and the blood upstream of it is `k` times XCAT's. For a phantom whose drawn vessel is narrower than the
blood the body model that times it holds (XCAT vmale_50's aorta). The marks on that tree (`apex_ml`, `ostia_ml` on the aorta,
`confluence_ml` on the superior vena cava) scale with it; `tree_ml` stays XCAT's.
"""
function scale_tree(transit, tree::UInt8, k::Real)
    volume = copy(transit.volume)
    on = transit.tree .== tree
    volume[on] .*= Float32(k)
    aorta = tree == 0x01; svc = tree == 0x03
    return merge(
        transit, (;
            volume,
            apex_ml = aorta ? Float32(k) * transit.apex_ml : transit.apex_ml,
            ostia_ml = aorta ? Float32(k) .* transit.ostia_ml : transit.ostia_ml,
            confluence_ml = svc ? Float32(k) * transit.confluence_ml : transit.confluence_ml,
        )
    )
end

"""
    vessel_volume_at(transit, p, t) -> Float32

Where the great-vessel voxel at `p` (mm) lies along tree `t` of `transit` ([`great_vessel_transit`](@ref)): the blood (mL)
between it and the heart, interpolated trilinearly over the tree's cells around it (smooth along the vessel), else its
nearest cell on the tree's value.
"""
function vessel_volume_at(transit, p, t)
    volume = transit.volume; tree = transit.tree
    f = (p .- transit.origin_mm) ./ transit.h .+ 1
    k = floor.(Int, f); w = f .- k
    sum_v = 0.0f0; sum_w = 0.0f0
    for corner in CartesianIndices((0:1, 0:1, 0:1))
        c = CartesianIndex(k .+ Tuple(corner))
        checkbounds(Bool, volume, c) && tree[c] == t && !isnan(volume[c]) || continue
        wc = Float32(prod(ifelse.(Tuple(corner) .== 1, w, 1 .- w)))
        sum_v += wc * volume[c]; sum_w += wc
    end
    sum_w > 0 && return sum_v / sum_w
    best = Inf32
    for step in CartesianIndices((-1:1, -1:1, -1:1))
        c = CartesianIndex(round.(Int, f)) + step
        checkbounds(Bool, volume, c) && tree[c] == t && !isnan(volume[c]) && (best = min(best, volume[c]))
    end
    return isfinite(best) ? best : 0.0f0
end

"""
    vessel_tree_at(transit, p, candidates) -> UInt8

The great-vessel tree the voxel at `p` (mm) is on: of `candidates`, the one most of the cells around it are on
(trilinear weights), else its nearest cell's, else the first candidate.
"""
function vessel_tree_at(transit, p, candidates)
    tree = transit.tree
    f = (p .- transit.origin_mm) ./ transit.h .+ 1
    k = floor.(Int, f); w = f .- k
    best = 0x00; best_w = 0.0f0
    for t in candidates
        wt = 0.0f0
        for corner in CartesianIndices((0:1, 0:1, 0:1))
            c = CartesianIndex(k .+ Tuple(corner))
            checkbounds(Bool, tree, c) && tree[c] == t || continue
            wt += Float32(prod(ifelse.(Tuple(corner) .== 1, w, 1 .- w)))
        end
        wt > best_w && (best = t; best_w = wt)
    end
    best != 0x00 && return best
    for step in CartesianIndices((-1:1, -1:1, -1:1))
        c = CartesianIndex(round.(Int, f)) + step
        checkbounds(Bool, tree, c) && tree[c] in candidates && return tree[c]
    end
    return first(candidates)
end

"""
    great_vessel_pool(transit, p, t) -> (VOX_POOL, pool, ml)

The body pool ([`XCAT_BODY_POOLS`](@ref)) of a voxel at `p` (mm) on great-vessel tree `t`, and the blood (mL) along its
vessel that times it: the aorta (5) or, past the arch's top, the descending aorta (27, timed from the arch's end at its
own flow); the pulmonary arteries (6); the superior (8) and inferior (11) venae cavae; the pulmonary veins (7); the right
(9) and left (10) brachiocephalic veins, timed from the SVC's top.
"""
function great_vessel_pool(transit, p, t)
    ml = vessel_volume_at(transit, p, t)
    if t == 0x01
        c = CartesianIndex(round.(Int, (p .- transit.origin_mm) ./ transit.h) .+ 1)
        descending = checkbounds(Bool, transit.descending, c) && transit.descending[c]
        return descending ? (_PHANTOM.VOX_POOL, Int32(27), max(ml - transit.apex_ml, 0.0f0)) : (_PHANTOM.VOX_POOL, Int32(5), ml)
    end
    t == 0x02 && return (_PHANTOM.VOX_POOL, Int32(6), ml)
    t == 0x03 && return (_PHANTOM.VOX_POOL, Int32(8), ml)
    t == 0x04 && return (_PHANTOM.VOX_POOL, Int32(11), ml)
    t == 0x05 && return (_PHANTOM.VOX_POOL, Int32(7), ml)
    return (_PHANTOM.VOX_POOL, Int32(t == 0x06 ? 9 : 10), max(ml - transit.confluence_ml, 0.0f0))
end
