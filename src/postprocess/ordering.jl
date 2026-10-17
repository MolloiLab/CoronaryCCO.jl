# SPDX-License-Identifier: GPL-3.0-or-later
# Diameter-defined Strahler ordering against the literature (Kassab) thresholds, as run()
# in coronary_postprocess.py (HCT-CCO commit 9bc141f) performs it.

# LITERATURE_MU / LITERATURE_SD (Kassab 1993, upstream's MU/SD arrays) live in morphometry.jl.

# Boundaries between consecutive literature orders: midpoint of the lower order's mean + SD
# and the upper order's mean - SD, evaluated left to right as numpy does.
# Any organ's reference morphometry (per-order mean and SD diameter, same units as the tree's
# diameters) gives its own boundaries; Kassab's coronary arrays are upstream's.
_literature_boundaries(mu, sd) = [(((mu[k + 1] + sd[k + 1]) + mu[k]) - sd[k]) / 2.0 for k in 1:(length(mu) - 1)]
const _LITERATURE_DB = _literature_boundaries(LITERATURE_MU, LITERATURE_SD)

# Ports literature_order_of: 11 minus the number of boundaries above d, clipped to 1..11
# (K orders for K - 1 boundaries `db`).
function literature_order_of(d::Real, db = _LITERATURE_DB)
    K = length(db) + 1
    o = K - count(b -> d < b, db)
    return clamp(o, 1, K)
end

# Shared by the stack-based walks below: upstream's iterative post-order with a `seen` set.
function _postorder!(post::Vector{Int}, seen::BitVector, root::Int, children)
    st = [(root, false)]
    while !isempty(st)
        nd, done = pop!(st)
        if done
            push!(post, nd)
            continue
        end
        seen[nd] && continue
        seen[nd] = true
        push!(st, (nd, true))
        for c in children[nd]
            push!(st, (c, false))
        end
    end
    return post
end

# Ports strahler_orders_per_segment: leaves are 1; a parent takes its largest child order,
# plus one when at least two children share it. Segments unreachable from a root stay 0.
function strahler_orders_per_segment(conn, n::Int)
    strahler = zeros(Int, n)
    post = Int[]
    seen = falses(n)
    for root in 1:n
        isempty(conn.parent[root]) && _postorder!(post, seen, root, conn.children)
    end
    for s in post
        ch = conn.children[s]
        if isempty(ch)
            strahler[s] = 1
        else
            mx = maximum(strahler[c] for c in ch)
            strahler[s] = count(c -> strahler[c] == mx, ch) >= 2 ? mx + 1 : mx
        end
    end
    return strahler
end

# Ports group_strahler_elements: maximal connected runs of equal Strahler order, grown by a
# depth-first walk over parent (siblings included, as upstream) and children lists. The
# discovery order of each run is kept because numpy's pairwise mean depends on it.
function group_strahler_elements(strahler::Vector{Int}, conn, n::Int)
    seg2elem = fill(0, n)
    off = [1]; val = Int32[]                   # elements as one flat list (see Ragged), in discovery order
    stk = Int[]
    for s in 1:n
        seg2elem[s] == 0 || continue
        o = strahler[s]
        e = length(off)
        seg2elem[s] = e
        push!(stk, s); push!(val, s)
        while !isempty(stk)
            cur = pop!(stk)
            for nb in Iterators.flatten((conn.parent[cur], conn.children[cur]))
                if seg2elem[nb] == 0 && strahler[nb] == o
                    seg2elem[nb] = e
                    push!(stk, nb)
                    push!(val, nb)
                end
            end
        end
        push!(off, length(val) + 1)
    end
    return seg2elem, Ragged(off, val)
end

# Ports refine_root_element: if the root's element spans more than one literature order, split
# it into maximal connected runs of equal literature order (the first run keeps the element id).
function refine_root_element(seg2elem, elem_segs, elem_order, diameters, conn, db = _LITERATURE_DB)
    n = length(seg2elem)
    root = findfirst(s -> isempty(conn.parent[s]), 1:n)
    root === nothing && return seg2elem, elem_segs, elem_order
    e0 = seg2elem[root]
    segs0 = elem_segs[e0]
    lit = Dict(s => literature_order_of(diameters[s], db) for s in segs0)
    length(unique(values(lit))) <= 1 && return seg2elem, elem_segs, elem_order
    seg2elem = copy(seg2elem)
    elem_order = copy(elem_order)
    comps = [Int[]]
    visited = Set{Int}()
    first = true
    for s in segs0
        s in visited && continue
        o = lit[s]
        stk = [s]
        comp = Int[]
        while !isempty(stk)
            cur = pop!(stk)
            cur in visited && continue
            push!(visited, cur)
            push!(comp, cur)
            for nb in Iterators.flatten((conn.parent[cur], conn.children[cur]))
                if haskey(lit, nb) && !(nb in visited) && lit[nb] == o
                    push!(stk, nb)
                end
            end
        end
        if first
            comps[1] = comp
            elem_order[e0] = o
            seg2elem[comp] .= e0
            first = false
        else
            push!(comps, comp)
            push!(elem_order, o)
            seg2elem[comp] .= length(elem_segs) + length(comps) - 1
        end
    end
    # element e0 becomes the first run, the others are appended (same lists as upstream's)
    lists = Iterators.flatten(((e == e0 ? comps[1] : elem_segs[e] for e in 1:length(elem_segs)), comps[2:end]))
    off = [1]; val = Int32[]
    for l in lists
        append!(val, l); push!(off, length(val) + 1)
    end
    return seg2elem, Ragged(off, val), elem_order
end

# Ports element_connectivity_counts: for every (daughter order m, parent order n), one count per
# parent element of order n of its distinct daughter elements of order m.
function element_connectivity_counts(seg2elem, elem_order, conn)
    ne = length(elem_order)
    # segments of each element (counting sort), then each element's distinct daughter elements
    cnt = zeros(Int, ne); foreach(e -> cnt[e] += 1, seg2elem)
    off = cumsum(vcat(1, cnt)); segs = Vector{Int32}(undef, length(seg2elem)); nxt = off[1:(end - 1)]
    for s in eachindex(seg2elem)
        e = seg2elem[s]; segs[nxt[e]] = s; nxt[e] += 1
    end
    present = sort(unique(elem_order))
    slot = zeros(Int, maximum(present; init = 0)); for (i, m) in enumerate(present)
        slot[m] = i
    end
    counts = Dict((m, p) => Int32[] for m in present for p in present)
    buf = Int32[]; tally = zeros(Int32, length(present))
    for e in 1:ne
        empty!(buf)
        for k in off[e]:(off[e + 1] - 1), c in conn.children[segs[k]]
            b = seg2elem[c]; b == e || push!(buf, b)
        end
        unique!(sort!(buf)); fill!(tally, 0)
        for b in buf
            tally[slot[elem_order[b]]] += 1
        end
        p = elem_order[e]
        for (i, m) in enumerate(present)
            push!(counts[(m, p)], tally[i])
        end
    end
    return present, counts
end

"""
    diameter_defined_orders(t::PPTree, diameters = 2 .* t.radius; elements = :series, reference) -> NamedTuple

The ordering of coronary_postprocess.run(): topological Strahler order per segment, grouping
into maximal same-order elements, relabelling each element by its mean diameter against the
literature thresholds, and the root-element guard. `orders[s]` is the final order of segment `s`.

Elements are runs of same-order segments joined end to end (`elements = :series`).
`elements = :upstream` reproduces upstream, whose neighbour lists also join sibling segments
(two daughters of one bifurcation with equal Strahler order become a single element).

`reference = (mu, sd)` is the per-order mean and SD diameter of the morphometry the orders are
labelled against (same units as `diameters`; default Kassab 1993's coronary arrays, as upstream).
Another organ passes its own, e.g. Huang et al. 1996's pulmonary arteries.
"""
function diameter_defined_orders(
        t::PPTree, diameters::Vector{Float64} = 2 .* t.radius; elements::Symbol = :series,
        reference = (LITERATURE_MU, LITERATURE_SD)
    )
    elements in (:series, :upstream) || throw(ArgumentError("elements must be :series or :upstream, got :$elements"))
    n = nsegments(t)
    conn = segment_connections(t; siblings = elements === :upstream)
    strahler = strahler_orders_per_segment(conn, n)
    seg2elem, elem_segs = group_strahler_elements(strahler, conn, n)
    db = _literature_boundaries(reference...)
    elem_order = [literature_order_of(_np_mean(diameters[c]), db) for c in elem_segs]
    seg2elem, elem_segs, elem_order = refine_root_element(seg2elem, elem_segs, elem_order, diameters, conn, db)
    orders = [elem_order[seg2elem[s]] for s in 1:n]
    return (; orders, strahler, seg2elem, elem_segs, elem_order, conn)
end
