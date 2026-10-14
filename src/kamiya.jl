# SPDX-License-Identifier: GPL-3.0-or-later
# Ports HctccoModel::kamyiaOptimization (src/HctccoModel.ih, HCT-CCO commit 9bc141f): the
# Kamiya-Togawa bifurcation position and radius update with one inner solve
# (`kamiya_ceres`, which reproduces upstream's Ceres call).

"""
    kamiya_optimization(gamma, lf, pCurrent, pParent, rCurrent, pL, rL, fL, pR, rR, fR)
        -> (ok, pOpt, r0, r1, r2)

One call of upstream `kamyiaOptimization(..., nbIter = 1, ...)`. `(pL, rL, fL)` is the new
left segment (the old segment's distal part) and `(pR, rR, fR)` the new terminal.
"""
@inline function kamiya_optimization(
        gamma::T, lf::T, pCurrent::P3{T}, pParent::P3{T}, rCurrent::T,
        pL::P3{T}, rL::T, fL::T, pR::P3{T}, rR::T, fR::T
    ) where {T}
    o = one(T)
    R0 = rCurrent * rCurrent
    f1 = fL; f2 = fR; f0 = f1 + f2
    ratioSeed = o
    (f2 > 0 && f1 > 0) && (ratioSeed = upow(f1 / f2, T(0.25)))
    bLs = upow(o / (o + upow(o / ratioSeed, -gamma)), o / gamma)
    bRs = upow(o / (o + upow(ratioSeed, -gamma)), o / gamma)
    rLs = cmax(T(1.0e-8), rCurrent * bLs)
    rRs = cmax(T(1.0e-8), rCurrent * bRs)
    R1 = rLs * rLs
    R2 = rRs * rRs
    den = cmax(T(1.0e-8), T(2) * f0)
    w1 = ((f0 * pParent[1] + f1 * pL[1]) + f2 * pR[1]) / den
    w2 = ((f0 * pParent[2] + f1 * pL[2]) + f2 * pR[2]) / den
    w3 = ((f0 * pParent[3] + f1 * pL[3]) + f2 * pR[3]) / den
    pb = (isfinite(w1) ? w1 : pCurrent[1], isfinite(w2) ? w2 : pCurrent[2], isfinite(w3) ? w3 : pCurrent[3])
    l0 = norm3(sub3(pParent, pb)) * lf
    l1 = norm3(sub3(pL, pb)) * lf
    l2 = norm3(sub3(pR, pb)) * lf
    dP1 = (f0 * l0) / (R0 * R0) + (f1 * l1) / (R1 * R1)
    dP2 = (f0 * l0) / (R0 * R0) + (f2 * l2) / (R2 * R2)
    ok, rr1, rr2, _ = kamiya_ceres(gamma, dP1, dP2, f0, f1, f2, l0, l1, l2, R1, R2)
    # Equation 27, evaluated whether or not the solve succeeded (as upstream does).
    R0 = upow(f0 * (upow(rr1, gamma) / f1 + upow(rr2, gamma) / f2), o / gamma)
    R1 = rr1; R2 = rr2
    q0 = R0 / l0; q1 = R1 / l1; q2 = R2 / l2
    dsum = (q0 + q1) + q2
    pb = (
        (((pParent[1] * R0) / l0 + (pL[1] * R1) / l1) + (pR[1] * R2) / l2) / dsum,
        (((pParent[2] * R0) / l0 + (pL[2] * R1) / l1) + (pR[2] * R2) / l2) / dsum,
        (((pParent[3] * R0) / l0 + (pL[3] * R1) / l1) + (pR[3] * R2) / l2) / dsum,
    )
    return ok, pb, sqrt(R0), sqrt(R1), sqrt(R2)
end
