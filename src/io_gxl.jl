# SPDX-License-Identifier: GPL-3.0-or-later
# Ports CoronaryTreeIO::readTreeFromXml / writeTreeToXml (src/helpers/CoronaryTreeIO.h,
# HCT-CCO commit 9bc141f), including the state upstream carries over from the constructor.

using Printf

const _RE_NODE = r"^\s*<node id=\"n(\d+)\">"
const _RE_ATTR = r"^\s*<attr name=\" ([^\"]{1,63})\">"
const _RE_FLOAT = r"^\s*<float>([^<]*)</float>"
const _RE_EDGE = r"^\s*<edge id=\"e\d+\" to=\"n(\d+)\" from=\"n(\d+)\">"

# sscanf("%f") into a float, then widened (or kept) to T.
function _scanf_float(::Type{T}, line) where {T}
    m = match(_RE_FLOAT, line)
    m === nothing && return nothing
    v = tryparse(Float32, strip(m.captures[1]))
    return v === nothing ? nothing : T(v)
end

"""
    import_gxl(T, path, params) -> CCOModel{T}

Builds the model exactly as upstream does for `-i base_tree.xml`: the constructor state
(a two-segment tree whose terminal list holds upstream index 1) followed by
readTreeFromXml, which replaces the segments but keeps that terminal entry and the first two
parent/children slots, then reconciles flows and runs a full hydraulic update.
"""
function import_gxl(::Type{T}, path::AbstractString, params::CCOParams) where {T}
    m = CCOModel{T}(params)
    # Constructor leftovers: parent/children slots for upstream 0 and 1, terminal list [1].
    parent = Int32[1, 1]; left = Int32[0, 0]; right = Int32[0, 0]
    terminals = Int32[2]
    segs = Vector{NamedTuple{(:id, :pos, :radius, :flow, :resistance), Tuple{Int, P3{T}, T, T, T}}}()
    lines = readlines(path)
    cur = (id = 0, pos = (zero(T), zero(T), zero(T)), radius = one(T), flow = zero(T), resistance = zero(T))
    to = 0
    i = 1
    while i <= length(lines)
        line = lines[i]
        if occursin("<node", line)
            mm = match(_RE_NODE, line)
            if mm !== nothing
                cur = (
                    id = parse(Int, mm.captures[1]), pos = (zero(T), zero(T), zero(T)),
                    radius = one(T), flow = zero(T), resistance = zero(T),
                )
            end
        elseif occursin("<attr name=\" nodeType\">", line)
            i += 1
            occursin("<string> terminal node </string>", lines[i]) && push!(terminals, Int32(cur.id + 1))
        elseif occursin("<attr name=\"", line)
            mm = match(_RE_ATTR, line)
            if mm !== nothing
                name = mm.captures[1]
                if name == "pPerf"
                    i += 1; v = _scanf_float(T, lines[i]); v === nothing || (m.pperf = v)
                elseif name == "pTerm"
                    i += 1; v = _scanf_float(T, lines[i]); v === nothing || (m.pterm = v)
                elseif name == "position"
                    i += 2; x = _scanf_float(T, lines[i])
                    i += 1; y = _scanf_float(T, lines[i])
                    i += 1; z = _scanf_float(T, lines[i])
                    p = cur.pos
                    p = (x === nothing ? p[1] : x, y === nothing ? p[2] : y, z === nothing ? p[3] : z)
                    cur = merge(cur, (pos = p,))
                elseif name in ("flow", "resistance", "radius")
                    i += 1; v = _scanf_float(T, lines[i])
                    if v !== nothing
                        k = findfirst(s -> s.id == to, segs)
                        k === nothing || (segs[k] = merge(segs[k], NamedTuple{(Symbol(name),)}((v,))))
                    end
                end
            end
        elseif occursin("</node>", line)
            if cur.id > 1
                push!(parent, Int32(1)); push!(left, Int32(0)); push!(right, Int32(0))
            end
            push!(segs, cur)
        elseif occursin("<edge", line)
            mm = match(_RE_EDGE, line)
            if mm !== nothing
                to = parse(Int, mm.captures[1]); from = parse(Int, mm.captures[2])
                parent[to + 1] = from + 1
                if left[from + 1] == 0
                    left[from + 1] = to + 1
                else
                    right[from + 1] = to + 1
                end
            end
        end
        i += 1
    end
    # myVectSegments is rebuilt from the nodes in file order (ids are assumed 0..n-1).
    n = length(segs)
    length(parent) == n || error("import_gxl: node ids must be 0..n-1 in order")
    parent[1] = 2   # tree.myVectParent[0] = 1
    m.pos = [s.pos for s in segs]
    m.radius = [s.radius for s in segs]
    m.resistance = [s.resistance for s in segs]
    m.flow = zeros(T, n)
    m.kterm = zeros(Int32, n)
    m.beta = ones(T, n)
    m.parent = parent; m.left = left; m.right = right
    m.terminals = terminals
    m.listed = falses(n); m.listed[terminals] .= true
    m.listed = Vector{Bool}(m.listed)
    m.glen = zeros(T, n)
    for j in 2:n
        refresh_glen!(m, j)
    end
    # Reconcile: unit terminal flow up each terminal's path, then a full hydraulic update.
    for t in terminals
        update_flow!(m, t, m.qterm)
    end
    m.kterm_total = length(terminals)
    update_length_factor!(m)
    order = postorder(m)
    update_resistance_from_root!(m, order)
    update_root_radius!(m, order)
    return m
end

"""
    import_gxl(T, paths::AbstractVector, params) -> CCOModel{T}

A forest: each file's tree is imported as above and appended, so the model has one root
segment per tree (parent = the shared sentinel). Terminal flow, pressures and the length
factor are shared, as for a single tree.
"""
function import_gxl(::Type{T}, paths::AbstractVector{<:AbstractString}, params::CCOParams) where {T}
    m = import_gxl(T, first(paths), params)
    for path in paths[2:end]
        t = import_gxl(T, path, params)
        o = Int32(nsegments(m) - 1)                   # t's index i >= 2 becomes i + o
        sh(i) = i == 0 ? Int32(0) : i == 1 ? Int32(1) : Int32(i + o)
        for f in (:pos, :radius, :kterm, :resistance, :flow, :beta, :listed, :glen)
            append!(getfield(m, f), getfield(t, f)[2:end])
        end
        append!(m.parent, sh.(t.parent[2:end]))
        append!(m.left, sh.(t.left[2:end])); append!(m.right, sh.(t.right[2:end]))
        append!(m.terminals, sh.(t.terminals))
        m.kterm_total += t.kterm_total
    end
    update_length_factor!(m)
    return m
end

# C++ `ostream << double` with default flags: %g with precision 6.
_g(x) = @sprintf("%g", Float64(x))

"""
    export_gxl(m, path)

CoronaryTreeIO::writeTreeToXml: nodes (type from the terminal list; a node is the "root" if
its parent is upstream index 0), then one edge per segment with flow, resistance and radius.
"""
function export_gxl(m::CCOModel, path::AbstractString)
    open(path, "w") do io
        println(io, "<gxl><graph id=\"", path, "\" edgeids=\" true\" edgemode=\" directed\" hypergraph=\" false\">")
        println(io, "<info_graph>")
        println(io, "    <attr name=\" pPerf\">\n      <float>", _g(m.pperf), "</float>\n    </attr>")
        println(io, "    <attr name=\" pTerm\">\n      <float>", _g(m.pterm), "</float>\n    </attr>")
        println(io, "</info_graph>")
        for j in eachindex(m.pos)
            kind = m.parent[j] == 1 ? "root node" : m.listed[j] ? "terminal node" : "bifurication"
            p = m.pos[j]
            println(
                io, "  <node id=\"n", j - 1, "\">\n    <attr name=\" nodeType\">\n      <string> ", kind,
                " </string>\n    </attr>\n    <attr name=\" position\">\n      <tup>"
            )
            for c in p
                println(io, "        <float>", _g(c), "</float>")
            end
            println(io, "      </tup>\n    </attr>\n  </node>")
        end
        for j in 2:length(m.pos)
            sv(v) = isfinite(v) ? v : zero(v)
            println(io, "  <edge id=\"e", j - 1, "\" to=\"n", j - 1, "\" from=\"n", m.parent[j] - 1, "\">")
            println(io, "    <attr name=\" flow\">\n      <float>", _g(sv(m.flow[j])), "</float>\n    </attr>")
            println(io, "    <attr name=\" resistance\">\n      <float>", _g(sv(m.resistance[j])), "</float>\n    </attr>")
            println(io, "    <attr name=\" radius\">\n      <float>", _g(sv(m.radius[j])), "</float>\n    </attr>")
            println(io, "  </edge>")
        end
        println(io, "</graph></gxl>")
    end
    return path
end
