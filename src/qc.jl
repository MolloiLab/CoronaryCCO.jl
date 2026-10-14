# SPDX-License-Identifier: GPL-3.0-or-later
# Geometric quality checks for a grown tree: overlapping vessels and sharp turns.

"""
    tree_qc(E; rmin = 0.25, hops = 2) -> NamedTuple

Checks the vessels with radius ≥ `rmin` (mm) of a grown tree (or forest):

- `overlaps`: pairs of segments whose tubes intersect (centerline distance < sum of radii),
  with the deepest penetration in mm. Pairs that must touch at a shared junction are skipped:
  segments within `hops` steps of each other (a step joins segments sharing an end node:
  parent and child, or siblings), and segments joined through the tree by a centerline path
  shorter than their two radii plus the largest radius along the path (short, wide imaged
  segments around a bifurcation, such as a left main and its first branches).
- turns: joints between a segment and its parent whose direction changes by more than 90°
  or 135°, split into imaged base-tree vessels and grown vessels.
"""
function tree_qc(E::Engine{T}; rmin::Real = 0.25, hops::Integer = 2) where {T}
    m = E.m
    big = [j for j in 3:nsegments(m) if m.radius[j] >= rmin]
    rbig = isempty(big) ? zero(T) : maximum(m.radius[big])
    # segments sharing an end node with j (index 1 is the root sentinel, not a segment)
    function neighbours(j)
        out = Int32[]
        p = m.parent[j]
        if p > 1
            push!(out, p)
            for c in (m.left[p], m.right[p])
                (c != 0 && c != j) && push!(out, c)
            end
        end
        for c in (m.left[j], m.right[j])
            c != 0 && push!(out, c)
        end
        return out
    end
    # checked segments reachable from j through checked segments: (steps, length of the
    # segments strictly between, largest radius among them), explored while within `hops`
    # steps or while the path is shorter than r_j + twice the largest checked radius
    function near(j)
        seen = Dict{Int32, Tuple{Int, T, T}}(Int32(j) => (0, zero(T), zero(T)))
        front = [Int32(j)]; h = 0
        while !isempty(front)
            h += 1
            nxt = Int32[]
            for i in front
                ci = i == j ? zero(T) : seen[i][2] + m.glen[i]
                ri = i == j ? zero(T) : max(seen[i][3], m.radius[i])
                (h > hops && ci > m.radius[j] + 2 * rbig) && continue
                for k in neighbours(i)
                    (m.radius[k] >= rmin && !haskey(seen, k)) || continue
                    seen[k] = (h, ci, ri); push!(nxt, k)
                end
            end
            front = nxt
        end
        return seen
    end
    nb = Dict(j => near(j) for j in big)
    function adjacent(a, b)
        x = get(nb[a], b, nothing)
        x === nothing && return false
        # tubes meeting at a junction fill a ball as wide as the widest vessel there
        return x[1] <= hops || x[2] < m.radius[a] + m.radius[b] + x[3]
    end
    noverlap = 0; worst_overlap = 0.0; worst_pair = (Int32(0), Int32(0))
    for x in 1:length(big), y in (x + 1):length(big)
        a = big[x]; b = big[y]
        adjacent(a, b) && continue
        d = seg2seg_distance(segA(m, a), segB(m, a), segA(m, b), segB(m, b))
        pen = Float64(m.radius[a] + m.radius[b] - d)
        if pen > 0
            noverlap += 1
            pen > worst_overlap && (worst_overlap = pen; worst_pair = (Int32(a), Int32(b)))
        end
    end
    turns = Dict(:base90 => 0, :base135 => 0, :grown90 => 0, :grown135 => 0)
    worst_turn = 0.0
    for j in big
        p = m.parent[j]
        (p <= 2 || m.radius[p] < rmin) && continue
        u = sub3(m.pos[p], m.pos[m.parent[p]]); v = sub3(m.pos[j], m.pos[p])
        nu = norm3(u); nv = norm3(v)
        (nu == 0 || nv == 0) && continue
        ang = acosd(clamp(Float64(dot3(u, v) / (nu * nv)), -1.0, 1.0))
        b = E.isbase[j] ? :base : :grown
        ang > 90 && (turns[Symbol(b, 90)] += 1)
        ang > 135 && (turns[Symbol(b, 135)] += 1)
        worst_turn = max(worst_turn, ang)
    end
    return (
        vessels = length(big), overlaps = noverlap, worst_overlap_mm = worst_overlap, worst_pair,
        base_turns_90 = turns[:base90], base_turns_135 = turns[:base135],
        grown_turns_90 = turns[:grown90], grown_turns_135 = turns[:grown135],
        worst_turn_deg = worst_turn,
    )
end
