# SPDX-License-Identifier: GPL-3.0-or-later
# Ports CoronaryGeometry (src/helpers/CoronaryGeometry.h) and HctccoModel::getProjDistance
# from HCT-CCO commit 9bc141f, operation for operation.

# projectOnStraightLine: projection of C on line AB, and upstream's "inside" test (true only
# if some axis has A < B with the projection not beyond B on that axis).
@inline function project_on_line(A::P3{T}, B::P3{T}, C::P3{T}) where {T}
    A == C && return A, true
    B == C && return B, true
    vAB = sub3(B, A)
    vABn = div3(vAB, norm3(vAB))
    d = dot3(sub3(C, A), vABn)
    proj = add3(A, mul3(vABn, d))
    res = (A[1] < B[1] && proj[1] <= B[1]) || (A[2] < B[2] && proj[2] <= B[2]) ||
        (A[3] < B[3] && proj[3] <= B[3])
    return proj, d >= zero(T) && res
end

# HctccoModel::getProjDistance(p0, p1, p)
@inline function proj_distance(p0::P3{T}, p1::P3{T}, p::P3{T}) where {T}
    proj, inside = project_on_line(p0, p1, p)
    inside && return norm3(sub3(proj, p))
    return cmin(norm3(sub3(p0, p)), norm3(sub3(p1, p)))
end

# CoronaryGeometry::segment2segmentDistance
@inline function seg2seg_distance(sA::P3{T}, sB::P3{T}, sC::P3{T}, sD::P3{T}) where {T}
    A, B, C, D = sA, sB, sC, sD
    if B[1] - A[1] == 0 && B[2] - A[2] == 0 && B[3] - A[3] == 0   # sameAB: swap AB <-> CD
        A, B, C, D = sC, sD, sA, sB
    end
    d1 = sub3(B, A); d2 = sub3(D, C); d12 = sub3(C, A)
    z = zero(T)
    D1 = ((z + d1[1] * d1[1]) + d1[2] * d1[2]) + d1[3] * d1[3]
    D2 = ((z + d2[1] * d2[1]) + d2[2] * d2[2]) + d2[3] * d2[3]
    R = ((z + d1[1] * d2[1]) + d1[2] * d2[2]) + d1[3] * d2[3]
    S1 = ((z + d1[1] * d12[1]) + d1[2] * d12[2]) + d1[3] * d12[3]
    S2 = ((z + d2[1] * d12[1]) + d2[2] * d12[2]) + d2[3] * d12[3]
    sAB = d1[1] == 0 && d1[2] == 0 && d1[3] == 0
    sCD = d2[1] == 0 && d2[2] == 0 && d2[3] == 0
    den = D1 * D2 - R * R
    u = z; t = z
    step = if !sAB && sCD
        4
    elseif sAB && sCD
        5
    elseif !sAB && !sCD && den == 0
        3
    else
        2
    end
    if step <= 2
        t = (S1 * D2 - S2 * R) / den
        t = t < 0 ? z : (t > 1 ? one(T) : t)
    end
    if step <= 3
        u = (t * R - S2) / D2
        if u < 0
            u = z
        elseif u > 1
            u = one(T)
        else
            step = 5
        end
    end
    if step <= 4
        t = (u * R + S1) / D1
        t = t < 0 ? z : (t > 1 ? one(T) : t)
    end
    DD = z
    for i in 1:3
        dd = (d1[i] * t - d2[i] * u) - d12[i]
        DD += dd * dd
    end
    return sqrt(DD)
end

# CoronaryGeometry::isIntersecting: thick segments AB (radius rAB) and CD (radius rCD).
@inline function is_intersecting(A::P3{T}, B::P3{T}, rAB::T, C::P3{T}, D::P3{T}, rCD::T) where {T}
    two = T(2)
    cAB = div3(add3(A, B), two)
    cCD = div3(add3(C, D), two)
    lAB = norm3(sub3(A, B))
    lCD = norm3(sub3(C, D))
    dC = norm3(sub3(cAB, cCD))
    dC > (((lAB + lCD) / two + rAB) + rCD) && return false
    return !(seg2seg_distance(A, B, C, D) > rAB + rCD)
end
