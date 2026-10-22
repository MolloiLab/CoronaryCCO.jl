# SPDX-License-Identifier: GPL-3.0-or-later
# Penetrating stubs of an imaged vessel tree.
#
# Imaged arteries rarely show their small side branches, but real ones give off branches into the tissue along their
# whole length: the pig LAD trunk carries ~18 daughters of orders 7-10 over its 48 mm, one every 2.8 mm (Estes et al.
# 1966; Farrer-Brown 1968; Kassab et al. 1993). HCT-CCO's staged growth attaches new vessels anywhere only in its first
# stage; later stages attach them near the listed ends of the tree. A bare imaged trunk has few ends, so the tissue is
# then reached from a few points by long chains of the smallest vessels. `penetrating_stubs` adds a short branch from
# the vessel into the tissue at every site where one is missing; each stub tip is an end of the inlet tree, so growth
# starts a subtree there, and the grower sizes every vessel by its flow. Host-only Float64 (XCAT precompute).

"""
    penetrating_stubs(t::BaseTree, g::VoxelGrid, tissues; scale = 1, spacing = 2.8scale, reach = 6scale, depth = 1.5scale,
                      domain = nothing, forbid = nothing, blockers = nothing, radius = 0.1) -> (tree, report)

Adds stubs, short branches from the vessels of `t` into the tissue they run beside, about every `spacing` mm along each
vessel. A stub runs from a node of `t` to the nearest voxel of `tissues` and `depth` mm on into it; its tip is a new end
of the tree (radius `radius` mm until growth sizes it by flow).

The step decides where stubs are needed, so it can be run on any inlet tree:
- Sites lie about `spacing` mm apart along each vessel: the node nearest one `spacing` past the previous site, a side
  branch of `t` (an imaged branch counts as a stub that is already there) or the root. The first node of every branch is
  a site. Sample the centrelines at `spacing` or finer (`centerline(s; step = spacing)`).
- A site gets no stub when it already has two children (the tree is binary); when no voxel of `tissues` lies within
  `reach` mm; when the tip lies outside the grid, outside `domain` or in `forbid`; or when the stub crosses `blockers`.
So a densely imaged tree gets few stubs, a bare trunk one per `spacing`, and stretches far from the tissue (e.g. a left
main running free in fat) none. To grow from the imaged tree alone, do not call it.

Where the vessel runs through the tissue (its nearest tissue voxel lies within its own radius, as for a femoral artery
inside a soft-tissue domain), the stub leaves sideways, perpendicular to the vessel and `radius_k + depth` mm long, turning
by the golden angle from one such stub to the next; of eight directions the first whose tip is inside and clear is taken.

`tissues` is one Bool mask on `g` or several (the heart's LV and RV walls); where several are within reach, successive
stubs alternate between them. `domain` (default: the union of `tissues`) and `forbid` are Bool masks on `g`;
`blockers` is `(mask, origin_mm, cell_mm)` as for [`grow`](@ref), a fine mask no stub may cross (chambers, great
vessels). `report` counts the sites and what became of them: `(; sites, stubs, sideways, no_tissue, outside, blocked)`.

Lengths: the heart's are the defaults at `scale = 1`: `spacing` 2.8 mm (the pig LAD's mean segment length; RCA 3.2,
Kassab et al. 1993), `reach` 6 mm, `depth` 1.5 mm. Like every other length set for the heart, they scale with the organ
(docs/deviations.md §16): pass `scale = size_scale(spec)` and sample the centrelines at `2.8 * scale` mm. An organ whose
feeding arteries' branch spacing is measured passes that as `spacing`.
"""
function penetrating_stubs(
        t::BaseTree, g::VoxelGrid, tissues; scale::Real = 1, spacing::Real = 2.8scale, reach::Real = 6scale,
        depth::Real = 1.5scale, domain = nothing,
        forbid::Union{Nothing, AbstractArray{Bool, 3}} = nothing, blockers = nothing, radius::Real = 0.1
    )
    T = tissues isa AbstractArray{Bool, 3} ? [tissues] : collect(tissues)
    isempty(T) && throw(ArgumentError("penetrating_stubs: give at least one tissue mask"))
    all(size(M) == g.dims for M in T) || throw(DimensionMismatch("penetrating_stubs: every tissue mask must have the grid's size $(g.dims)"))
    dom = domain === nothing ? reduce(.|, T) : domain
    size(dom) == g.dims || throw(DimensionMismatch("penetrating_stubs: domain must have the grid's size $(g.dims)"))
    nch = zeros(Int, nnodes(t))
    for q in t.parent
        q > 0 && (nch[q] += 1)
    end
    sites = _stub_sites(t, nch, spacing)
    par = copy(t.parent); X = copy(t.X); r = copy(t.r); ves = copy(t.vessel)
    inside(tip) = (
        ci = [round(Int, (tip[i] - g.origin[i]) / g.h) + 1 for i in 1:3];
        all(1 .<= ci .<= collect(g.dims)) && dom[ci...] && (forbid === nothing || !forbid[ci...])
    )
    clear(a, b) = blockers === nothing || _stub_clear(blockers, a, b)
    side = 0; nside = 0; nst = 0; nnot = 0; nout = 0; nblk = 0
    for k in sites
        p = t.X[:, k]
        cand = filter(!isnothing, [_nearest_voxel(p, M, g, reach) for M in T])
        isempty(cand) && (nnot += 1; continue)
        d, q = cand[1 + side % length(cand)]; side += 1
        if d <= t.r[k]
            # the vessel runs through the tissue: the stub leaves sideways, turning by the golden angle from one such
            # stub to the next; the first of eight directions whose tip is inside and clear is taken
            u, v = _perpendiculars(t, k); tip = nothing; why = :outside
            for j in 0:7
                θ = 2π * (nside * 0.3819660112501051 + j / 8)
                c = p .+ (t.r[k] + depth) .* (cos(θ) .* u .+ sin(θ) .* v)
                inside(c) || continue
                clear(p, c) || (why = :blocked; continue)
                tip = c; break
            end
            tip === nothing && (why === :blocked ? (nblk += 1) : (nout += 1); continue)
            nside += 1
        else
            dir = (q .- p) ./ max(d, 1.0e-6); tip = q .+ depth .* dir
            inside(tip) || (nout += 1; continue)
            clear(p, tip) || (nblk += 1; continue)
        end
        push!(par, k); X = hcat(X, tip); push!(r, radius); push!(ves, t.vessel[k]); nst += 1
    end
    report = (; sites = length(sites), stubs = nst, sideways = nside, no_tissue = nnot, outside = nout, blocked = nblk)
    return BaseTree(t.name, par, X, r, ves), report
end

# Sites about `spacing` apart along each vessel, in node order. `since` is the arc length back to the last site, imaged
# side branch or the root, so it depends only on a node's path from the root (walked from the root, whatever the node
# order); a node is a site when it is the first of its vessel or the node nearest one spacing past that point (the next
# node along the vessel would overshoot by more).
function _stub_sites(t::BaseTree, nch::Vector{Int}, spacing::Real)
    n = nnodes(t)
    count(==(0), t.parent) == 1 || throw(ArgumentError("penetrating_stubs: the tree must have exactly one root"))
    len(k) = sqrt(sum(abs2, t.X[:, k] .- t.X[:, t.parent[k]]))
    kids = [Int[] for _ in 1:n]; nxt = zeros(Int, n)
    for k in 1:n
        p = t.parent[k]; p == 0 && continue
        push!(kids[p], k); t.vessel[p] == t.vessel[k] && (nxt[p] = k)
    end
    since = zeros(n); site = falses(n); stack = [findfirst(==(0), t.parent)]
    while !isempty(stack)
        p = pop!(stack)
        for k in kids[p]
            e = len(k); s = since[p] + e
            if nch[k] >= 2
                since[k] = 0.0
            elseif t.vessel[p] != t.vessel[k] || s >= spacing - (nxt[k] > 0 ? len(nxt[k]) : e) / 2
                site[k] = true; since[k] = 0.0
            else
                since[k] = s
            end
            push!(stack, k)
        end
    end
    return findall(site)
end

# Two unit vectors perpendicular to the vessel at node k (its direction from the parent node) and to each other.
function _perpendiculars(t::BaseTree, k::Int)
    a = t.X[:, k] .- t.X[:, t.parent[k]]; n = sqrt(sum(abs2, a))
    a = n > 0 ? a ./ n : [0.0, 0.0, 1.0]
    e = zeros(3); e[argmin(abs.(a))] = 1.0                       # the axis least aligned with the vessel
    u = e .- (e' * a) .* a; u ./= sqrt(sum(abs2, u))
    v = [a[2] * u[3] - a[3] * u[2], a[3] * u[1] - a[1] * u[3], a[1] * u[2] - a[2] * u[1]]
    return u, v
end

# The voxel of M nearest to p within `reach` mm: (distance, centre) or nothing.
function _nearest_voxel(p, M::AbstractArray{Bool, 3}, g::VoxelGrid, reach::Real)
    c = [round(Int, (p[d] - g.origin[d]) / g.h) + 1 for d in 1:3]; nr = ceil(Int, reach / g.h); best = (Inf, c)
    for kk in max(1, c[3] - nr):min(g.dims[3], c[3] + nr), j in max(1, c[2] - nr):min(g.dims[2], c[2] + nr),
            i in max(1, c[1] - nr):min(g.dims[1], c[1] + nr)
        M[i, j, kk] || continue
        q = center(g, i, j, kk); d2 = sum((q[d] - p[d])^2 for d in 1:3)
        d2 < best[1] && (best = (d2, [i, j, kk]))
    end
    return best[1] <= reach^2 ? (sqrt(best[1]), collect(center(g, best[2]...))) : nothing
end

# Whether the straight stub a → b misses every cell of the fine blocker mask (sampled 8 times per cell).
function _stub_clear((B, O, h), a, b)
    n = max(1, ceil(Int, 8 * sqrt(sum(abs2, b .- a)) / h))
    for s in 0:n
        q = a .+ (s / n) .* (b .- a); I = [round(Int, (q[d] - O[d]) / h) + 1 for d in 1:3]
        all(1 .<= I .<= collect(size(B))) && B[I...] != 0 && return false
    end
    return true
end
