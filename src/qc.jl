# SPDX-License-Identifier: GPL-3.0-or-later
# Geometric quality checks for a grown tree: overlapping vessels and sharp turns.

"""
    tree_qc(E; rmin = 0.25) -> NamedTuple

Checks the vessels with radius ≥ `rmin` (mm) of a grown tree:

- `overlaps`: pairs of non-adjacent segments (not parent/child, not siblings) whose tubes
  intersect (centerline distance < sum of radii), with the deepest penetration in mm.
- turns: joints between a segment and its parent whose direction changes by more than 90°
  or 135°, split into imaged base-tree vessels and grown vessels.
"""
function tree_qc(E::Engine{T}; rmin::Real = 0.25) where {T}
    m = E.m
    big = [j for j in 3:nsegments(m) if m.radius[j] >= rmin]
    adjacent(a, b) = m.parent[a] == b || m.parent[b] == a || m.parent[a] == m.parent[b]
    noverlap = 0; worst_overlap = 0.0
    for x in 1:length(big), y in (x + 1):length(big)
        a = big[x]; b = big[y]
        adjacent(a, b) && continue
        d = seg2seg_distance(segA(m, a), segB(m, a), segA(m, b), segB(m, b))
        pen = Float64(m.radius[a] + m.radius[b] - d)
        if pen > 0
            noverlap += 1
            worst_overlap = max(worst_overlap, pen)
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
        vessels = length(big), overlaps = noverlap, worst_overlap_mm = worst_overlap,
        base_turns_90 = turns[:base90], base_turns_135 = turns[:base135],
        grown_turns_90 = turns[:grown90], grown_turns_135 = turns[:grown135],
        worst_turn_deg = worst_turn,
    )
end
