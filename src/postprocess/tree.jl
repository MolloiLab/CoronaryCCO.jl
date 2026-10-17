# SPDX-License-Identifier: GPL-3.0-or-later
# Reading the tree XML and the segment adjacency used by every ordering routine
# (coronary_postprocess.py at HCT-CCO commit 9bc141f). Segments are the XML edges in file
# order; Python's 0-based segment and point ids are 1-based here.

const _PP_NODE = r"<node id=\"n(\d+)\">(.*?)</node>"s
const _PP_POS = r"name=\"\s*position\">\s*<tup>\s*<float>([-\d.eE+]+)</float>\s*<float>([-\d.eE+]+)</float>\s*<float>([-\d.eE+]+)</float>"
const _PP_EDGE = r"<edge id=\"e\d+\" to=\"n(\d+)\" from=\"n(\d+)\">(.*?)</edge>"s
const _PP_RADIUS = r"name=\"\s*radius\">\s*<float>([-\d.eE+]+)</float>"

"""
    PPTree

The tree as coronary_postprocess.py sees it: `points[i]` is the position of node `n(i-1)`
(origin for ids without a position), and segment `s` runs from `points[from[s]]` to
`points[to[s]]` with radius `radius[s]`.
"""
struct PPTree
    points::Vector{NTuple{3, Float64}}
    from::Vector{Int}
    to::Vector{Int}
    radius::Vector{Float64}
end

nsegments(t::PPTree) = length(t.radius)

# Ports build_polydata_from_tree_xml (and the identical parsing in compute_hemodynamics):
# the same regular expressions, the same handling of missing fields.
function read_tree_xml(path::AbstractString)
    text = read(path, String)
    pos = Dict{Int, NTuple{3, Float64}}()
    for m in eachmatch(_PP_NODE, text)
        pm = match(_PP_POS, m.captures[2])
        pm === nothing && continue
        pos[parse(Int, m.captures[1])] = (
            parse(Float64, pm.captures[1]),
            parse(Float64, pm.captures[2]),
            parse(Float64, pm.captures[3]),
        )
    end
    np = isempty(pos) ? 0 : maximum(keys(pos)) + 1
    points = [get(pos, i, (0.0, 0.0, 0.0)) for i in 0:(np - 1)]
    from, to, radius = Int[], Int[], Float64[]
    for m in eachmatch(_PP_EDGE, text)
        rm = match(_PP_RADIUS, m.captures[3])
        a, b = parse(Int, m.captures[2]) + 1, parse(Int, m.captures[1]) + 1
        max(a, b) <= np || error("read_tree_xml: edge references node n$(max(a, b) - 1) without a position")
        push!(from, a)
        push!(to, b)
        push!(radius, rm === nothing ? 0.0 : parse(Float64, rm.captures[1]))
    end
    return PPTree(points, from, to, radius)
end

# Lists of integers stored flat (`val`) with offsets (`off`): `r[i]` is a view of list i. One Vector{Int} per segment
# costs ~60 bytes of overhead each, which a 68 M-segment lung tree cannot afford; the lists, their order and their
# contents are the same.
struct Ragged <: AbstractVector{SubArray{Int32, 1, Vector{Int32}, Tuple{UnitRange{Int}}, true}}
    off::Vector{Int}
    val::Vector{Int32}
end
Base.size(r::Ragged) = (length(r.off) - 1,)
Base.@propagate_inbounds Base.getindex(r::Ragged, i::Int) = view(r.val, r.off[i]:(r.off[i + 1] - 1))

# Lists from a count pass and a fill pass: `each(i, f)` calls f(x) for every member x of list i, in order.
function _ragged(n::Int, each)
    off = ones(Int, n + 1)
    for i in 1:n
        k = 0; each(i, _ -> (k += 1)); off[i + 1] = off[i] + k
    end
    val = Vector{Int32}(undef, off[end] - 1)
    for i in 1:n
        k = off[i]; each(i, x -> (val[k] = x; k += 1))
    end
    return Ragged(off, val)
end

# Ports build_segment_connections: neighbours through segment ENDPOINTS. `children[s]` are the
# other segments touching s's end point; `parent[s]` the segment ending at its start point.
# Upstream lists every other segment touching the start point, so s's siblings land in
# `parent[s]` as well; `siblings = true` reproduces that.
function segment_connections(t::PPTree; siblings::Bool = false)
    n = nsegments(t)
    # segments at each point, in segment order (from and to of segment s, as upstream appends them)
    cnt = zeros(Int, length(t.points))
    for s in 1:n
        cnt[t.from[s]] += 1; cnt[t.to[s]] += 1
    end
    off = cumsum(vcat(1, cnt)); val = Vector{Int32}(undef, off[end] - 1); nxt = off[1:(end - 1)]
    for s in 1:n
        val[nxt[t.from[s]]] = s; nxt[t.from[s]] += 1
        val[nxt[t.to[s]]] = s; nxt[t.to[s]] += 1
    end
    p2s = Ragged(off, val)
    parent = _ragged(
        n, (s, f) -> for x in p2s[t.from[s]]
            (x != s && (siblings || t.to[x] == t.from[s])) && f(x)
        end
    )
    children = _ragged(
        n, (s, f) -> for x in p2s[t.to[s]]
            x != s && f(x)
        end
    )
    return (parent = parent, children = children)
end
