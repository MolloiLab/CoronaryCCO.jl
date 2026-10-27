# SPDX-License-Identifier: GPL-3.0-or-later
# XCAT's surface meshes and the vessel trees made from their tubes.
#
# Besides its NURBS heart, XCAT draws organs and vessels as triangle meshes (`*_organs.raw`, `*_vessels.raw`): one block of
# triangles per surface, each headed by its name, a name repeated for every tube of a vessel system (the lungs' arteries
# and veins are hundreds of short tubes, one calibre each, that meet only near each other). These functions read the
# blocks, take each tube's centreline, and join the tubes into inlet trees: joined by nearness into trees, the trees
# grafted onto a trunk (the main pulmonary artery) or carried to the chamber they drain into (the left atrium). Joining
# segments get their radii from growth, which sizes every vessel by its flow; only their geometry is imposed here.
# Host-only Float64 (XCAT precompute).

"""
    read_mesh_blocks(path; names = nothing, offset = (0, 0, 0)) -> Dict{String, Vector{Tuple{Matrix{Float64}, Matrix{Int}}}}

The triangle blocks of an XCAT `.raw` mesh file, by name: each block as `(points, triangles)`, points 3 × n (mm, plus
`offset`) and triangles 3 × m (indices into the points), the form [`voxelize`](@ref) takes. A name heads every block of
its surface, so a vessel system drawn as many tubes gives one entry per tube. `names` keeps only those surfaces.
`offset` moves the mesh into another frame (XCAT writes its meshes `offset` mm from its NURBS frame; its log gives it).
"""
function read_mesh_blocks(path::AbstractString; names = nothing, offset = (0.0, 0.0, 0.0))
    out = Dict{String, Vector{Tuple{Matrix{Float64}, Matrix{Int}}}}()
    name = ""; P = Float64[]
    function flush!()
        if !isempty(P) && (names === nothing || name in names)
            X = reshape(copy(P), 3, :) .+ collect(Float64, offset)
            push!(get!(out, name, Tuple{Matrix{Float64}, Matrix{Int}}[]), (X, reshape(collect(1:size(X, 2)), 3, :)))
        end
        return empty!(P)
    end
    for line in eachline(path)
        t = strip(line); isempty(t) && continue
        if isletter(t[1])
            flush!(); name = String(t)
        else
            append!(P, parse.(Float64, split(t)))
        end
    end
    flush!()
    return out
end

"""
    tube_centerline(points; name = "tube", step = 1.0, max_nodes = 400) -> Centerline

The centreline of one tube mesh, given as its vertices (3 × n, mm): slices across its principal axis about `step` mm
apart, each slice's centroid a node. XCAT's mesh tubes have one calibre, so every node takes the tube's median slice
radius (the end slices are partial rings). A tube too short to slice is its two ends at its axis' extremes.
"""
function tube_centerline(P::AbstractMatrix{<:Real}; name::AbstractString = "tube", step::Real = 1.0, max_nodes::Integer = 400)
    c = vec(sum(P; dims = 2)) ./ size(P, 2); Q = P .- c
    axis = _principal_axis(Q)
    s = vec(axis' * Q); lo, hi = extrema(s); L = hi - lo
    k = clamp(round(Int, L / step), 1, max_nodes - 1)
    edges = range(lo, hi; length = k + 2)
    X = Vector{Vector{Float64}}(); r = Float64[]
    for i in 1:(k + 1)
        sel = findall(x -> edges[i] <= x <= edges[i + 1], s); length(sel) < 3 && continue
        q = vec(sum(P[:, sel]; dims = 2)) ./ length(sel); push!(X, q)
        push!(r, _median([sqrt(sum(abs2, (P[:, j] .- q) .- ((P[:, j] .- q)' * axis) .* axis)) for j in sel]))
    end
    if length(X) < 2
        across = [sqrt(sum(abs2, Q[:, j] .- (Q[:, j]' * axis) .* axis)) for j in axes(Q, 2)]
        X = [c .+ lo .* axis, c .+ hi .* axis]; r = fill(_median(across), 2)
    end
    return Centerline(String(name), reduce(hcat, X), fill(_median(r), length(r)))
end

_median(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

# The unit direction of largest spread of the centred points Q (3 × n), by power iteration on Q Qᵀ.
function _principal_axis(Q)
    M = Q * Q'; v = [1.0, 0.7, 0.3]
    for _ in 1:200
        w = M * v; n = sqrt(sum(abs2, w)); n == 0 && return [1.0, 0.0, 0.0]
        w ./= n; sum(abs2, w .- v) < 1.0e-24 && return w
        v = w
    end
    return v
end

"""
    tube_forest(lines; toward, join_mm = 6.0, min_nodes = 12) -> Vector{BaseTree}

Joins tube centrelines into trees, largest tube first: a tube attaches by its nearer end to the nearest node of an already
placed tube at least half its radius, when that end lies within the two radii plus `join_mm` of it (a node that already
has a child is split at the midpoint of its parent segment, so the tree stays binary); otherwise it starts a new tree,
whose root is its end nearer `toward` (mm; the hilum or the heart). Trees of fewer than `min_nodes` nodes (stray single
tubes) are dropped. Each tree's nodes run from its root, parents first.
"""
function tube_forest(lines::Vector{Centerline}; toward, join_mm::Real = 6.0, min_nodes::Integer = 12)
    hub = collect(Float64, toward)
    X = Vector{Vector{Float64}}(); R = Float64[]; par = Int[]; tree = Int[]; vessel = String[]; roots = Int[]
    for t in sortperm([maximum(l.r) for l in lines]; rev = true)
        l = lines[t]; tx = l.X; tr = l.r; rt = maximum(tr)
        best = (Inf, 0, 1)
        for (e, p) in ((1, tx[:, 1]), (2, tx[:, end])), n in eachindex(X)
            R[n] >= 0.5rt || continue
            d = sqrt(sum(abs2, p .- X[n])); d < best[1] && (best = (d, n, e))
        end
        d, n, e = best
        if n == 0 || d > R[n] + rt + join_mm
            sum(abs2, tx[:, 1] .- hub) > sum(abs2, tx[:, end] .- hub) && (tx = tx[:, end:-1:1]; tr = reverse(tr))
            push!(roots, length(X) + 1); a = 0; id = length(roots)
        else
            e == 2 && (tx = tx[:, end:-1:1]; tr = reverse(tr))
            a = n; id = tree[n]
            if count(==(n), par) >= 1 # the node already has a child: split its parent segment
                p = par[n]
                push!(X, p == 0 ? X[n] : (X[n] .+ X[p]) ./ 2); push!(R, R[n]); push!(par, p); push!(tree, id); push!(vessel, vessel[n])
                j = length(X); par[n] = j; a = j
                p == 0 && (roots[id] = j) # the junction is now this tree's root
            end
        end
        for i in axes(tx, 2)
            push!(X, tx[:, i]); push!(R, tr[i]); push!(par, i == 1 ? a : length(X) - 1); push!(tree, id); push!(vessel, l.name)
        end
    end
    out = BaseTree[]
    for (id, r0) in enumerate(roots)
        nodes = findall(==(id), tree); length(nodes) < min_nodes && continue
        push!(out, _subtree("$(lines[1].name)_$(length(out) + 1)", r0, par, X, R, vessel))
    end
    return out
end

# The tree hanging from node r0, renumbered from it (parents first).
function _subtree(name, r0, par, X, R, vessel)
    kids = Dict{Int, Vector{Int}}()
    for i in eachindex(par)
        par[i] > 0 && push!(get!(kids, par[i], Int[]), i)
    end
    order = [r0]; h = 1
    while h <= length(order)
        append!(order, get(kids, order[h], Int[])); h += 1
    end
    number = Dict(i => q for (q, i) in enumerate(order))
    parent = [i == r0 ? 0 : number[par[i]] for i in order]
    return BaseTree(name, parent, reduce(hcat, X[order]), R[order], vessel[order])
end

"""
    graft(trunk::Centerline, trees; name = trunk.name) -> BaseTree

One tree: the `trunk` (e.g. the main pulmonary artery from its valve) with `trees` joined to it, in order of their roots'
distance from the trunk's distal end, each at the nearest node already in the tree (the trunk's or an earlier tree's). The
trunk runs from the end farther from every tree root; its rings narrower than half its widest (a valve sinus or dome
collapsing at an end) are dropped. A node that already has two children is split at the midpoint of its parent segment.
The root edge takes the trunk's median ring radius, robust to the sinus and the bifurcation dome at its two ends.
"""
function graft(trunk::Centerline, trees::Vector{BaseTree}; name::AbstractString = trunk.name)
    isempty(trees) && throw(ArgumentError("graft: no trees to join"))
    tips = reduce(hcat, [t.X[:, 1] for t in trees])
    nearest_root(p) = minimum(sqrt(sum(abs2, p .- tips[:, k])) for k in axes(tips, 2))
    keep = findall(>=(0.5 * maximum(trunk.r)), trunk.r)
    C = trunk.X[:, keep]; r = trunk.r[keep]
    length(r) >= 2 || throw(ArgumentError("graft: the trunk has fewer than two rings"))
    nearest_root(C[:, 1]) < nearest_root(C[:, end]) && (C = C[:, end:-1:1]; r = reverse(r))
    X = [C[:, i] for i in axes(C, 2)]; R = copy(r); par = collect(0:(length(R) - 1)); vessel = fill(trunk.name, length(R))
    R[2] = _median(r) # the GXL keeps each edge's radius on its child: the root edge carries the trunk's
    distal = C[:, end]
    for k in sortperm([sqrt(sum(abs2, t.X[:, 1] .- distal)) for t in trees])
        t = trees[k]
        n = argmin([sum(abs2, x .- t.X[:, 1]) for x in X])
        n = _free_node!(X, R, par, vessel, n)
        off = length(X)
        for i in 1:nnodes(t)
            push!(X, t.X[:, i]); push!(R, t.r[i]); push!(par, t.parent[i] == 0 ? n : off + t.parent[i]); push!(vessel, t.vessel[i])
        end
    end
    return BaseTree(String(name), par, reduce(hcat, X), R, vessel)
end

# A node of the tree that can take one more child: `n` itself, or a new node at the midpoint of its parent segment when
# `n` already has two children. Returns that node.
function _free_node!(X, R, par, vessel, n)
    count(==(n), par) < 2 && return n
    p = par[n]
    push!(X, p == 0 ? X[n] : (X[n] .+ X[p]) ./ 2); push!(R, R[n]); push!(par, p); push!(vessel, vessel[n])
    par[n] = length(X)
    return length(X)
end

"""
    drain(trees, surface; root_radius, name = "veins") -> BaseTree

One draining tree: the largest of `trees` takes the others, each joined by its root at its nearest node, and a straight
segment carries the joined tree from its root to the nearest point of `surface` (3 × n points, mm: the chamber it drains
into), which becomes the root, of radius `root_radius` (mm).
"""
function drain(trees::Vector{BaseTree}, surface::AbstractMatrix{<:Real}; root_radius::Real, name::AbstractString = "veins")
    isempty(trees) && throw(ArgumentError("drain: no trees to join"))
    main = argmax(nnodes.(trees)); m = trees[main]
    X = [m.X[:, i] for i in 1:nnodes(m)]; R = copy(m.r); par = copy(m.parent); vessel = copy(m.vessel)
    rest = [k for k in eachindex(trees) if k != main]
    for k in rest[sortperm([minimum(sum(abs2, trees[k].X[:, 1] .- x) for x in X) for k in rest])]
        t = trees[k]
        n = _free_node!(X, R, par, vessel, argmin([sum(abs2, x .- t.X[:, 1]) for x in X]))
        off = length(X)
        for i in 1:nnodes(t)
            push!(X, t.X[:, i]); push!(R, t.r[i]); push!(par, t.parent[i] == 0 ? n : off + t.parent[i]); push!(vessel, t.vessel[i])
        end
    end
    ostium = surface[:, argmin(vec(sum(abs2, surface .- X[1]; dims = 1)))]
    pushfirst!(X, ostium); pushfirst!(R, root_radius); pushfirst!(vessel, "ostium")
    par = [0; [p == 0 ? 1 : p + 1 for p in par]]
    R[2] = root_radius # the GXL keeps each edge's radius on its child: the root edge carries the ostium's
    return BaseTree(String(name), par, reduce(hcat, X), R, vessel)
end
