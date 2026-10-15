# SPDX-License-Identifier: GPL-3.0-or-later
# XCAT NURBS surfaces: read and write the `.nrb` text format, evaluate the B-spline surfaces
# exactly, and tessellate them into triangle meshes.
#
# Format (one block per surface): name, "<M> :M", "<N> :N", "U Knot Vector", N+4 knots,
# "V Knot Vector", M+4 knots, "Control Points", M·N lines "x,y,z" (mm). The u direction (N
# points) runs fastest; for XCAT vessels it goes around the tube and v runs along it. XCAT surfaces are non-rational cubics (all weights 1), so this is a B-spline evaluator.

"""
    NurbsSurface

Cubic B-spline surface with control points `P[:, j, i]`: `j` along u (fastest in the file,
`nu` points) and `i` along v (`nv` points); knot vectors `u` (length nu + 4) and `v`
(length nv + 4). Coordinates in mm.
"""
struct NurbsSurface
    name::String
    u::Vector{Float64}
    v::Vector{Float64}
    P::Array{Float64, 3}    # (3, nu, nv)
end

nu(s::NurbsSurface) = size(s.P, 2)
nv(s::NurbsSurface) = size(s.P, 3)

"""
    read_nrb(path) -> Dict{String, NurbsSurface}

Reads every surface in an XCAT `.nrb` file (for example `vmale_50_heart.nrb`), keyed by name
(`"dias_lad1"`, `"sys_lv_0"`, ...).
"""
function read_nrb(path::AbstractString)
    L = readlines(path)
    out = Dict{String, NurbsSurface}()
    i = 1
    num(s) = parse(Float64, strip(s))
    while i <= length(L)
        name = strip(L[i]); i += 1
        isempty(name) && continue
        M = parse(Int, split(L[i])[1]); i += 1
        N = parse(Int, split(L[i])[1]); i += 1
        startswith(L[i], "U Knot") || error("read_nrb: expected U knots for $name at line $i")
        i += 1
        u = Float64[]
        while !startswith(L[i], "V Knot")
            push!(u, num(L[i])); i += 1
        end
        i += 1
        v = Float64[]
        while !startswith(L[i], "Control")
            push!(v, num(L[i])); i += 1
        end
        i += 1
        length(u) == N + 4 && length(v) == M + 4 || error("read_nrb: $name has $(length(u))/$(length(v)) knots for $M×$N points")
        P = Array{Float64, 3}(undef, 3, N, M)
        for k in 1:(M * N)
            f = split(L[i], (',', ' ', '\t'); keepempty = false)
            jj = mod1(k, N); ii = cld(k, N)
            for d in 1:3
                P[d, jj, ii] = parse(Float64, f[d])
            end
            i += 1
        end
        out[String(name)] = NurbsSurface(String(name), u, v, P)
    end
    return out
end

"""
    write_nrb(path, surfaces)

Writes surfaces in XCAT `.nrb` format, in the given order (a vector of `NurbsSurface`, or a
Dict, written sorted by name).
"""
function write_nrb(path::AbstractString, surfaces)
    S = surfaces isa AbstractDict ? [surfaces[k] for k in sort(collect(keys(surfaces)))] : collect(surfaces)
    open(path, "w") do io
        for s in S
            println(io, s.name); println(io, nv(s), " :M"); println(io, nu(s), " :N")
            println(io, "U Knot Vector"); foreach(x -> println(io, _g(x)), s.u)
            println(io, "V Knot Vector"); foreach(x -> println(io, _g(x)), s.v)
            println(io, "Control Points")
            for i in 1:nv(s), j in 1:nu(s)
                println(io, _g(s.P[1, j, i]), ",", _g(s.P[2, j, i]), ",", _g(s.P[3, j, i]))
            end
        end
    end
    return path
end

_g(x::Real) = (y = round(Float64(x); digits = 4); isinteger(y) ? string(Int(y)) : string(y))

# Index of the knot span containing t (clamped cubic), and the four nonzero cubic basis values.
function _span(knots::Vector{Float64}, n::Int, t::Float64)
    t >= knots[n + 1] && return n
    lo, hi = 4, n + 1
    while hi - lo > 1
        mid = (lo + hi) >> 1
        t < knots[mid] ? (hi = mid) : (lo = mid)
    end
    return lo
end

function _basis(knots::Vector{Float64}, s::Int, t::Float64)
    # Cox–de Boor, degree 3 (Piegl & Tiller, The NURBS Book, A2.2); nonzero functions s-3..s
    N = zeros(4); left = zeros(4); right = zeros(4)
    N[1] = 1.0
    for j in 1:3
        left[j + 1] = t - knots[s + 1 - j]
        right[j + 1] = knots[s + j] - t
        saved = 0.0
        for r in 0:(j - 1)
            den = right[r + 2] + left[j - r + 1]
            tmp = den == 0 ? 0.0 : N[r + 1] / den
            N[r + 1] = saved + right[r + 2] * tmp
            saved = left[j - r + 1] * tmp
        end
        N[j + 1] = saved
    end
    return N
end

"""
    evaluate(s, t_u, t_v) -> NTuple{3, Float64}

Point on the surface at parameters (`t_u`, `t_v`), each within its knot range.
"""
function evaluate(s::NurbsSurface, tu::Real, tv::Real)
    iu = _span(s.u, nu(s), Float64(tu)); iv = _span(s.v, nv(s), Float64(tv))
    Bu = _basis(s.u, iu, Float64(tu)); Bv = _basis(s.v, iv, Float64(tv))
    x = 0.0; y = 0.0; z = 0.0
    for a in 1:4, b in 1:4
        w = Bu[a] * Bv[b]
        w == 0 && continue
        j = iu - 4 + a; i = iv - 4 + b
        x += w * s.P[1, j, i]; y += w * s.P[2, j, i]; z += w * s.P[3, j, i]
    end
    return (x, y, z)
end

"""
    tessellate(s; nu = 64, nv = 64) -> (points::Matrix{Float64} (3×K), triangles::Matrix{Int} (3×T))

Samples the surface on a regular parameter grid and triangulates it. Degenerate triangles
(collapsed caps) are dropped.
"""
function tessellate(s::NurbsSurface; nu::Int = 64, nv::Int = 64)
    tu = range(s.u[4], s.u[end - 3], length = nu); tv = range(s.v[4], s.v[end - 3], length = nv)
    X = Matrix{Float64}(undef, 3, nu * nv)
    for (a, x) in enumerate(tu), (b, y) in enumerate(tv)
        p = evaluate(s, x, y)
        k = b + nv * (a - 1)
        X[1, k] = p[1]; X[2, k] = p[2]; X[3, k] = p[3]
    end
    T = Int[]
    id(a, b) = b + nv * (a - 1)
    area2(i, j, k) = (
        e1 = X[:, j] .- X[:, i]; e2 = X[:, k] .- X[:, i];
        sum(abs2, (e1[2] * e2[3] - e1[3] * e2[2], e1[3] * e2[1] - e1[1] * e2[3], e1[1] * e2[2] - e1[2] * e2[1]))
    )
    for a in 1:(nu - 1), b in 1:(nv - 1)
        p, q, r, t = id(a, b), id(a + 1, b), id(a + 1, b + 1), id(a, b + 1)
        area2(p, q, r) > 1.0e-12 && append!(T, (p, q, r))
        area2(p, r, t) > 1.0e-12 && append!(T, (p, r, t))
    end
    return X, reshape(T, 3, :)
end
