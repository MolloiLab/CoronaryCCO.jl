# SPDX-License-Identifier: GPL-3.0-or-later
# The tree's VTK XML PolyData file and its OBJ tube mesh, written without VTK. They hold what
# coronary_postprocess.run() (HCT-CCO commit 9bc141f) writes with vtkXMLPolyDataWriter, and
# with vtkCellDataToPointData + vtkTubeFilter + vtkOBJWriter.

using Base64: base64encode

# One inline binary DataArray: base64 of a UInt32 byte count followed by the little-endian data.
function _vtk_array(io::IO, type, name, data::Vector{T}; ncomp::Int = 1, indent = "        ") where {T}
    buf = IOBuffer()
    write(buf, htol(UInt32(sizeof(data))))
    for v in data
        write(buf, htol(v))
    end
    print(io, indent, "<DataArray type=\"", type, "\" Name=\"", name, "\"")
    ncomp > 1 && print(io, " NumberOfComponents=\"", ncomp, "\"")
    println(io, " format=\"binary\">")
    println(io, indent, "  ", base64encode(take!(buf)))
    println(io, indent, "</DataArray>")
    return nothing
end

"""
    write_tree_vtp(path, t::PPTree, orders)

VTK XML PolyData with one point per node id (`Float32`, as vtkPoints stores them), one
two-point line per segment, and cell arrays `radius` (`Float64`) and `final_order` (`Int32`,
the active scalars), the content run() writes to `<prefix>_tree.vtp`.
"""
function write_tree_vtp(path::AbstractString, t::PPTree, orders::Vector{Int})
    n = nsegments(t)
    pts = Float32[c for p in t.points for c in p]
    conn = Int64[v - 1 for s in 1:n for v in (t.from[s], t.to[s])]
    offs = Int64[2s for s in 1:n]
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">")
        println(io, "  <PolyData>")
        println(
            io, "    <Piece NumberOfPoints=\"", length(t.points), "\" NumberOfVerts=\"0\" NumberOfLines=\"", n,
            "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">"
        )
        println(io, "      <PointData>\n      </PointData>")
        println(io, "      <CellData Scalars=\"final_order\">")
        _vtk_array(io, "Float64", "radius", t.radius)
        _vtk_array(io, "Int32", "final_order", Int32.(orders))
        println(io, "      </CellData>")
        println(io, "      <Points>")
        _vtk_array(io, "Float32", "Points", pts; ncomp = 3)
        println(io, "      </Points>")
        println(io, "      <Lines>")
        _vtk_array(io, "Int64", "connectivity", conn)
        _vtk_array(io, "Int64", "offsets", offs)
        println(io, "      </Lines>")
        println(io, "    </Piece>")
        println(io, "  </PolyData>")
        println(io, "</VTKFile>")
    end
    return path
end

# vtkMath::Normalize: divide by the Euclidean norm unless it is zero; returns (v, norm).
function _vnormalize(v::NTuple{3, Float64})
    d = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3])
    return d == 0.0 ? (v, d) : ((v[1] / d, v[2] / d, v[3] / d), d)
end
_vcross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
# VTK calls the platform libm's cos/sin (cf. `cpow`); Julia's own differ in the last bit for
# a few of the tube angles.
_ccos(x::Float64) = ccall(:cos, Float64, (Float64,), x)
_csin(x::Float64) = ccall(:sin, Float64, (Float64,), x)
_f32(v::NTuple{3, Float64}) = (Float64(Float32(v[1])), Float64(Float32(v[2])), Float64(Float32(v[3])))

# vtkPolyLine::GenerateSlidingNormals for a two-point line: a vector orthogonal to the
# segment built from its first non-zero component, normalized and stored as Float32.
function _sliding_normal(snext::NTuple{3, Float64})
    i = findfirst(!=(0.0), snext)
    nrm = zeros(3)
    nrm[mod(i + 2, 1:3)] = 0.0
    nrm[mod(i + 1, 1:3)] = 1.0
    nrm[i] = -snext[mod(i + 1, 1:3)] / snext[i]
    return _f32(_vnormalize((nrm[1], nrm[2], nrm[3]))[1])
end

"""
    write_tree_obj(path, t::PPTree; sides = 24)

Capped tube mesh of the tree, built as run() builds it: the cell `radius` averaged onto the
points (vtkCellDataToPointData), then for every segment a vtkTubeFilter tube (absolute
radius, `sides` sides sharing vertices, triangle-strip caps), written like vtkOBJWriter
(positions, normals, strips split into triangles, `f v//n`).
"""
function write_tree_obj(path::AbstractString, t::PPTree; sides::Int = 24)
    n = nsegments(t)
    np = length(t.points)
    # vtkCellDataToPointData: each point gets the mean radius of the cells that use it.
    rsum = zeros(np)
    rcnt = zeros(Int, np)
    for s in 1:n, v in (t.from[s], t.to[s])
        rsum[v] += t.radius[s]
        rcnt[v] += 1
    end
    prad = [c == 0 ? 0.0 : r / c for (r, c) in zip(rsum, rcnt)]
    theta = 2.0 * pi / sides
    cs = [_ccos(k * theta) for k in 0:(sides - 1)]
    sn = [_csin(k * theta) for k in 0:(sides - 1)]
    pts = NTuple{3, Float32}[]
    nrms = NTuple{3, Float32}[]
    strips = Vector{Int}[]
    sizehint!(pts, 4 * sides * n)
    sizehint!(nrms, 4 * sides * n)
    for s in 1:n
        ids = (t.from[s], t.to[s])
        p0, p1 = _f32(t.points[ids[1]]), _f32(t.points[ids[2]])
        d = (p1[1] - p0[1], p1[2] - p0[2], p1[3] - p0[3])
        unit, len = _vnormalize(d)
        len == 0.0 && continue          # vtkTubeFilter skips lines with coincident points
        nrm = _sliding_normal(unit)
        offset = length(pts)
        startcap = _vnormalize((-d[1], -d[2], -d[3]))[1]
        endcap = _vnormalize(unit)[1]
        for j in 1:2
            # vtkTubeFilter::GeneratePoints bevels with the mean of the previous and next
            # directions; on a two-point line those are the raw and the normalized direction at
            # the first point, and the normalized and the renormalized direction at the second.
            sprev, snext = j == 1 ? (d, unit) : (unit, endcap)
            p = j == 1 ? p0 : p1
            sv = _vnormalize(((sprev[1] + snext[1]) / 2.0, (sprev[2] + snext[2]) / 2.0, (sprev[3] + snext[3]) / 2.0))[1]
            w = _vnormalize(_vcross(sv, nrm))[1]
            np_ = _vnormalize(_vcross(w, sv))[1]
            r = prad[ids[j]]
            for k in 1:sides
                nk = (w[1] * cs[k] + np_[1] * sn[k], w[2] * cs[k] + np_[2] * sn[k], w[3] * cs[k] + np_[3] * sn[k])
                push!(pts, Float32.((p[1] + 1.0 * r * nk[1], p[2] + 1.0 * r * nk[2], p[3] + 1.0 * r * nk[3])))
                push!(nrms, Float32.(nk))
            end
        end
        for (ring, capn) in ((0, startcap), (sides, endcap)), k in 1:sides
            push!(pts, pts[offset + ring + k])
            push!(nrms, Float32.(capn))
        end
        # vtkTubeFilter::GenerateStrips (0-based point ids, offset by this tube's first point)
        for k in 0:(sides - 1)
            i1, i2 = k, (k + 1) % sides
            push!(strips, [offset + i2, offset + i1, offset + sides + i2, offset + sides + i1])
        end
        start = offset + 2 * sides
        cap = [start, start + 1]
        i1, i2 = sides - 1, 2
        for k in 0:(sides - 3)
            if isodd(k)
                push!(cap, start + i2); i2 += 1
            else
                push!(cap, start + i1); i1 -= 1
            end
        end
        push!(strips, cap)
        start += sides
        cap = [start, start + sides - 1]
        i1, i2 = sides - 2, 1
        for k in 0:(sides - 3)
            if isodd(k)
                push!(cap, start + i1); i1 -= 1
            else
                push!(cap, start + i2); i2 += 1
            end
        end
        push!(strips, cap)
    end
    open(path, "w") do io
        println(io, "# Generated by Vasculature.jl (vtkTubeFilter-equivalent tube mesh)")
        buf = UInt8[]
        dg = UInt8[]
        flush!(buf) = length(buf) > 1 << 20 ? (write(io, buf); empty!(buf)) : buf
        for (tag, vs) in (("v ", pts), ("vn ", nrms)), v in vs
            append!(buf, codeunits(tag))
            for c in 1:3
                c > 1 && push!(buf, UInt8(' '))
                _put_to_chars!(buf, dg, Float64(v[c]))
            end
            push!(buf, UInt8('\n'))
            flush!(buf)
        end
        # vtkOBJWriter splits each strip into triangles, flipping every other one.
        for st in strips, i in 3:length(st)
            append!(buf, codeunits("f"))
            for v in (iseven(i) ? (st[i - 1], st[i - 2], st[i]) : (st[i - 2], st[i - 1], st[i]))
                push!(buf, UInt8(' '))
                _put_int!(buf, v + 1)
                push!(buf, UInt8('/'), UInt8('/'))
                _put_int!(buf, v + 1)
            end
            push!(buf, UInt8('\n'))
            flush!(buf)
        end
        write(io, buf)
    end
    return path
end
