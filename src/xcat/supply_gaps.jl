# SPDX-License-Identifier: GPL-3.0-or-later
# Supply gaps of an imaged vessel tree.
#
# An imaged (base) tree is only as complete as the imaging. Where it resolves no vessel near some tissue, the growth has
# to supply that tissue from far away, one greedy connection at a time, and builds a long, tortuous, oversized trunk
# whose iodine arrives late (XCAT's left lung: 15.5% of its volume more than 20 mm from any imaged artery). Real organs
# have a resolved branch there. `fill_supply_gaps` adds one: tissue farther than `max_gap_mm` from the tree is a gap;
# each connected gap gets a straight branch from the nearest node that can reach it without leaving the organ, sized
# by coronary/pulmonary structure–function scaling from the tissue it supplies, Q ∝ D^qexp with qexp = 7/3 (Huo & Kassab
# 2012). Repeated until no gap is left. The grower then runs unchanged on the completed tree (and resizes every vessel
# by its flow), so this is a property of the input tree, organ-generic, and a no-op on a tree without gaps.
# Host-only Float64 (XCAT precompute, charter decision E2).

"""
    fill_supply_gaps(parent, X, r, domain, lower; r_resolved, h = 1.0, max_gap_mm = nothing, min_gap_ml = nothing, qexp = 7 / 3,
                     step_mm = 2.0, forbidden = nothing, region = nothing, free_mm = 3.0, max_iter = 200)
        -> (; parent, X, r, added, gaps)

`parent` (0 = root), `X` (3 × n node positions, mm) and `r` (radii, mm) describe the tree as nodes; segment k runs from
`parent[k]` to `k`. `domain` is the organ mask on a grid whose voxel (1, 1, 1) is centred at `lower` (mm) with spacing
`h`. A gap is a connected set (6-neighbour) of domain voxels farther than `max_gap_mm` from every segment, of at least
`min_gap_ml`. Both follow from the imaging's resolution by default: the most tissue a vessel thinner than `r_resolved`
(the smallest radius the imaging resolves, mm; required, since a grown or flow-resized tree's smallest radius says
nothing about the imaging) can supply is V_max = V_organ · (r_resolved / r_root)^qexp, so a region larger
than V_max that is farther than G = (3 V_max / 4π)^(1/3) from the imaged tree must be fed by a vessel the imaging would
have shown, and one is missing (`min_gap_ml` = V_max, `max_gap_mm` = G). Its branch ends at the gap voxel nearest the gap's centroid and starts at the nearest node whose radius is
at least the new branch's radius and whose straight path to it obeys the organ's blockers: no point of it lies in
`forbidden` (a Bool mask on the domain grid: chambers, airways, ...), and every point farther than `free_mm` from the start
lies in the domain and, with a `region` label map (left / right lung, lobes, ...), in the target's own region. Nodes are
placed every `step_mm`. A start node that already has two children is moved to the midpoint of its parent segment (no
node gets three children). Radius: r_root × (V_gap / V_domain)^(1/qexp). Returns the new arrays, one entry per added
branch in `added`, the gap volumes found in each round in `gaps`, and in `unfilled` the gaps no node can reach without
crossing a blocker or a region boundary (an organ region with no imaged vessel of its own).
"""
function fill_supply_gaps(
        parent::AbstractVector{<:Integer}, X::AbstractMatrix{<:Real}, r::AbstractVector{<:Real}, domain::AbstractArray{Bool, 3}, lower;
        r_resolved::Real, h::Real = 1.0, max_gap_mm = nothing, min_gap_ml = nothing, qexp::Real = 7 / 3, step_mm::Real = 2.0,
        forbidden::Union{Nothing, AbstractArray{Bool, 3}} = nothing, region::Union{Nothing, AbstractArray{<:Integer, 3}} = nothing,
        free_mm::Real = 3.0, max_iter::Integer = 200
    )
    par = Int.(collect(parent)); P = [Float64.(X[:, k]) for k in axes(X, 2)]; R = Float64.(collect(r))
    lo = Float64.(collect(lower)); dims = size(domain)
    vox = findall(domain); vml = h^3 / 1000; Vdom = length(vox) * vml
    root = findfirst(==(0), par); rroot = R[root]
    Vmax = Vdom * (r_resolved / rroot)^qexp                    # mL: the most an unresolved vessel can supply
    min_gap_ml = min_gap_ml === nothing ? Vmax : Float64(min_gap_ml)
    max_gap_mm = max_gap_mm === nothing ? cbrt(3 * Vmax * 1000 / (4π)) : Float64(max_gap_mm)
    centre(I) = lo .+ (Tuple(I) .- 1) .* h
    vidx(p) = (I = round.(Int, (p .- lo) ./ h) .+ 1; all(1 .<= I .<= dims) ? CartesianIndex(I...) : nothing)
    blocked(I) = I !== nothing && forbidden !== nothing && forbidden[I]
    # the straight path a → t: never through a blocker; past the first free_mm, inside the organ and the target's region
    function path_ok(a, t)
        L = sqrt(sum(abs2, t .- a)); ns = max(2, ceil(Int, 4L / h)); It = vidx(t); rt = region === nothing || It === nothing ? 0 : region[It]
        for s in 0:ns
            p = a .+ (s / ns) .* (t .- a); I = vidx(p); blocked(I) && return false
            s * L / ns <= free_mm && continue
            (I === nothing || !domain[I]) && return false
            region === nothing || region[I] == rt || return false
        end
        return true
    end
    added = NamedTuple[]; gaps = Vector{Float64}[]; unfilled = Float64[]
    for _ in 1:max_iter
        d = _distance_to_tree(par, P, vox, lo, h)
        far = falses(dims); for (i, I) in enumerate(vox)
            d[i] > max_gap_mm && (far[I] = true)
        end
        # gaps per region: a gap is filled from its own side (left / right lung, lobe, ...)
        comps = region === nothing ? _components(far) : reduce(vcat, [_components(far .& (region .== v)) for v in unique(region[vox])])
        sizes = [length(c) * vml for c in comps]; push!(gaps, sort(sizes; rev = true))
        big = [c for c in comps if length(c) * vml >= min_gap_ml]
        isempty(big) && break
        nadd = length(added)
        for c in sort(big; by = length, rev = true)
            Vg = length(c) * vml; rn = rroot * (Vg / Vdom)^(1 / qexp)
            cen = sum(centre(I) for I in c) ./ length(c)
            tgt = centre(c[argmin([sum(abs2, centre(I) .- cen) for I in c])])
            # nearest node that is thick enough and reaches the target through the organ
            cand = sort([k for k in eachindex(P) if R[k] >= rn]; by = k -> sum(abs2, P[k] .- tgt))
            # the start point: the node itself, or the midpoint of its parent segment when it already branches twice
            # (no node gets three children); the path is checked from the point the branch really starts at
            start(k) = count(==(k), par) >= 2 && par[k] != 0 ? (P[k] .+ P[par[k]]) ./ 2 : P[k]
            a = 0
            for k in cand
                count(==(k), par) >= 2 && par[k] == 0 && continue
                path_ok(start(k), tgt) && (a = k; break)
            end
            a == 0 && continue
            if count(==(a), par) >= 2
                p = par[a]; push!(P, (P[a] .+ P[p]) ./ 2); push!(R, R[a]); push!(par, p); par[a] = length(P); a = length(P)
            end
            L = sqrt(sum(abs2, tgt .- P[a])); m = max(1, round(Int, L / step_mm)); prev = a
            for i in 1:m
                push!(P, P[a] .+ (i / m) .* (tgt .- P[a])); push!(R, rn); push!(par, prev); prev = length(P)
            end
            push!(added, (; from = a, radius_mm = rn, length_mm = L, gap_ml = Vg, target = tgt))
        end
        length(added) == nadd && (unfilled = [length(c) * vml for c in big]; break)     # no gap reachable any more
    end
    return (; parent = par, X = reduce(hcat, P), r = R, added, gaps, unfilled, max_gap_mm, min_gap_ml)
end

# distance (mm) from each listed domain voxel to the nearest point of any segment, via a uniform bucket grid of
# segment samples (spacing ≤ h/2 along each segment)
function _distance_to_tree(par, P, vox, lo, h)
    pts = Vector{Float64}[]
    for k in eachindex(par)
        p = par[k]; p == 0 && (push!(pts, P[k]); continue)
        L = sqrt(sum(abs2, P[k] .- P[p])); ns = max(1, ceil(Int, 2L / h))
        for s in 0:ns
            push!(pts, P[p] .+ (s / ns) .* (P[k] .- P[p]))
        end
    end
    cell = 8.0; key(p) = (floor(Int, p[1] / cell), floor(Int, p[2] / cell), floor(Int, p[3] / cell))
    grid = Dict{NTuple{3, Int}, Vector{Int}}(); for (i, p) in enumerate(pts)
        push!(get!(grid, key(p), Int[]), i)
    end
    d = fill(Inf, length(vox))
    AK.foreachindex(d; min_elems = 1024) do i
        q = lo .+ (Tuple(vox[i]) .- 1) .* h; c = key(q); best = Inf; ring = 0
        while true
            for a in -ring:ring, b in -ring:ring, e in -ring:ring
                max(abs(a), abs(b), abs(e)) == ring || continue
                for j in get(grid, (c[1] + a, c[2] + b, c[3] + e), Int[])
                    best = min(best, sum(abs2, pts[j] .- q))
                end
            end
            sqrt(best) <= ring * cell && break      # every unvisited point is at least ring·cell away
            ring += 1; ring > 64 && break
        end
        d[i] = sqrt(best)
    end
    return d
end

# connected components (6-neighbour) of a Bool volume
function _components(M::AbstractArray{Bool, 3})
    seen = falses(size(M)); out = Vector{Vector{CartesianIndex{3}}}()
    nb = (CartesianIndex(1, 0, 0), CartesianIndex(-1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, -1, 0), CartesianIndex(0, 0, 1), CartesianIndex(0, 0, -1))
    for I in CartesianIndices(M)
        (M[I] && !seen[I]) || continue
        comp = CartesianIndex{3}[]; st = [I]; seen[I] = true
        while !isempty(st)
            J = pop!(st); push!(comp, J)
            for o in nb
                K = J + o; checkbounds(Bool, M, K) && M[K] && !seen[K] && (seen[K] = true; push!(st, K))
            end
        end
        push!(out, comp)
    end
    return out
end
