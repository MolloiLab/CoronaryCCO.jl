# SPDX-License-Identifier: GPL-3.0-or-later
# Bit-exact replica of the single Ceres Solver call in upstream `CoronaryGeometry::kamyiaOpt`.
#
# Upstream solves Kamiya's two radius equations (CostOptPos, "Equation 28") with
# `AutoDiffCostFunction<CostOptPos, 2, 1, 1>` and default `Solver::Options`. This file
# reproduces, operation for operation, what Ceres Solver 2.2.0 does on that problem when it
# is built with Eigen 3.4.0, -ffp-contract=off and only Eigen's sparse Cholesky
# (test/oracle/build.sh): trust-region Levenberg-Marquardt with Jacobi scaling, solving each
# step with SPARSE_NORMAL_CHOLESKY + EIGEN_SPARSE (SimplicialLDLT, natural ordering). With
# T = Float64 on the CPU it returns the same bits as that build for every recorded call
# (test/test_kamiya.jl).
#
# It mirrors the numerical procedure of Ceres Solver (Copyright 2023 Google Inc., BSD
# 3-Clause; http://ceres-solver.org) and of Eigen's SimplicialLDLT (MPL-2.0) so that our
# results match upstream's. No source is copied; the comments name the Ceres/Eigen file and
# function each block follows. Ceres's license notice: "Redistribution and use in source and
# binary forms, with or without modification, are permitted provided that ... Neither the
# name of Google Inc. nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission."
#
# Everything is scalar and allocation free, generic in T, so the same code can run in a GPU
# kernel or in Float32. Rules that keep the Float64 path bitwise exact:
#   * every pow goes through `upow` (the C library `pow` that Ceres calls);
#   * no fma/muladd/@fastmath, and each expression keeps the C++ evaluation order;
#   * the `zero(T) + ...` terms are the zero-initialised accumulators of Ceres's loops.
#     They only change the sign of a zero, but cost nothing and keep the replica literal.

# ---------------------------------------------------------------------------------------
# Jets (include/ceres/jet.h): value `a` and the two partial derivatives (d/dx1, d/dx2).

struct _Jet2{T}
    a::T
    v1::T
    v2::T
end

# Binary + and - (operator+, operator-).
@inline _jadd(f::_Jet2, g::_Jet2) = _Jet2(f.a + g.a, f.v1 + g.v1, f.v2 + g.v2)
@inline _jsub(f::_Jet2, g::_Jet2) = _Jet2(f.a - g.a, f.v1 - g.v1, f.v2 - g.v2)
# Jet * Jet: (f.a * g.a, f.a * g.v + f.v * g.a).
@inline _jmul(f::_Jet2, g::_Jet2) =
    _Jet2(f.a * g.a, f.a * g.v1 + f.v1 * g.a, f.a * g.v2 + f.v2 * g.a)
# s * Jet and Jet * s: (f.a * s, f.v * s).
@inline _jscale(f::_Jet2, s) = _Jet2(f.a * s, f.v1 * s, f.v2 * s)
# Jet / s multiplies by the reciprocal: s_inverse = 1 / s; (f.a * s_inverse, f.v * s_inverse).
@inline function _jdiv(f::_Jet2{T}, s::T) where {T}
    sinv = one(T) / s
    return _Jet2(f.a * sinv, f.v1 * sinv, f.v2 * sinv)
end
# pow(Jet, double): tmp = g * pow(f.a, g - 1); (pow(f.a, g), tmp * f.v).
@inline function _jpow(f::_Jet2{T}, g::T) where {T}
    tmp = g * upow(f.a, g - one(T))
    return _Jet2(upow(f.a, g), tmp * f.v1, tmp * f.v2)
end

# ---------------------------------------------------------------------------------------
# Residuals: upstream CostOptPos::operator(). The functor recomputes the same pow
# expression four times; it is deterministic, so it is evaluated once here.

# T = double instantiation, used when Ceres needs only the cost (AutoDiffCostFunction
# calls the functor directly). Note `pow(x, e) / f1` is a true division here.
@inline function _lm_residuals(
        gamma::T, dp1::T, dp2::T, f0::T, f1::T, f2::T, l0::T, l1::T,
        l2::T, x1::T, x2::T
    ) where {T}
    e = (T(3) + gamma) / T(2)
    P = upow(f0 * ((upow(x1, e) / f1) + (upow(x2, e) / f2)), T(2) / T(3))
    r1 = ((dp1 * x1) * x1) * P - ((f0 * l0) * x1) * x1 - (f1 * l1) * P
    r2 = ((dp2 * x2) * x2) * P - ((f0 * l0) * x2) * x2 - (f2 * l2) * P
    return r1, r2
end

# T = Jet<double, 2> instantiation, used when Ceres needs the Jacobian (autodiff.h seeds
# x1 = (x1, [1, 0]) and x2 = (x2, [0, 1])). Its residual values can differ from the double
# path in the last bit because Jet / double multiplies by the reciprocal.
@inline function _lm_jets(
        gamma::T, dp1::T, dp2::T, f0::T, f1::T, f2::T, l0::T, l1::T,
        l2::T, x1::T, x2::T
    ) where {T}
    o, z = one(T), zero(T)
    X1 = _Jet2(x1, o, z)
    X2 = _Jet2(x2, z, o)
    e = (T(3) + gamma) / T(2)
    S = _jadd(_jdiv(_jpow(X1, e), f1), _jdiv(_jpow(X2, e), f2))
    W = _jpow(_jscale(S, f0), T(2) / T(3))
    R1 = _jsub(
        _jsub(_jmul(_jmul(_jscale(X1, dp1), X1), W), _jmul(_jscale(X1, f0 * l0), X1)),
        _jscale(W, f1 * l1)
    )
    R2 = _jsub(
        _jsub(_jmul(_jmul(_jscale(X2, dp2), X2), W), _jmul(_jscale(X2, f0 * l0), X2)),
        _jscale(W, f2 * l2)
    )
    return R1, R2
end

# array_utils.cc IsArrayValid: finite and not the "uninitialised" marker kImpossibleValue.
@inline _lm_valid(x::T) where {T} = isfinite(x) & (x != T(1.0e302))

# residual_block_utils.cc IsEvaluationValid: residuals and Jacobian entries all valid.
@inline _lm_jets_valid(R1::_Jet2, R2::_Jet2) =
    _lm_valid(R1.a) & _lm_valid(R2.a) & _lm_valid(R1.v1) & _lm_valid(R2.v1) &
    _lm_valid(R1.v2) & _lm_valid(R2.v2)

# Eigen's NaN-propagating packet max (NEON fmax), used by lpNorm<Infinity>.
@inline _lm_pmax(a, b) = (isnan(a) | isnan(b)) ? a + b : ifelse(a < b, b, a)

# ---------------------------------------------------------------------------------------
# Evaluator (program_evaluator.h ProgramEvaluator::Evaluate + residual_block.cc).

# Cost only: ResidualBlock::Evaluate with jacobians == nullptr. Returns (ok, cost).
@inline function _lm_cost(
        gamma::T, dp1::T, dp2::T, f0::T, f1::T, f2::T, l0::T, l1::T,
        l2::T, x1::T, x2::T
    ) where {T}
    r1, r2 = _lm_residuals(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, x1, x2)
    (_lm_valid(r1) & _lm_valid(r2)) || return false, zero(T)
    # cost = 0.5 * squaredNorm(residuals), summed into zero-initialised accumulators.
    cost = zero(T) + (zero(T) + T(0.5) * (r1 * r1 + r2 * r2))
    return isfinite(cost), cost
end

# ---------------------------------------------------------------------------------------
# The solver: Solver::Solve -> TrustRegionMinimizer::Minimize with LevenbergMarquardtStrategy
# and SparseNormalCholeskySolver. Returns (termination_type, x1, x2, iterations) where
# termination_type uses Ceres's codes: 0 CONVERGENCE, 1 NO_CONVERGENCE, 2 FAILURE.
function _kamiya_ceres_full(
        gamma::T, dp1::T, dp2::T, f0::T, f1::T, f2::T, l0::T, l1::T,
        l2::T, x1in::T, x2in::T
    ) where {T <: AbstractFloat}
    CONVERGENCE, NO_CONVERGENCE, FAILURE = 0, 1, 2
    z = zero(T)
    # Solver::Options defaults (include/ceres/solver.h).
    max_num_iterations = 50
    max_num_consecutive_invalid_steps = 5
    function_tolerance = T(1.0e-6)
    gradient_tolerance = T(1.0e-10)
    parameter_tolerance = T(1.0e-8)
    min_relative_decrease = T(1.0e-3)
    min_trust_region_radius = T(1.0e-32)
    max_trust_region_radius = T(1.0e16)
    min_lm_diagonal = T(1.0e-6)
    max_lm_diagonal = T(1.0e32)
    dmax = floatmax(T)  # std::numeric_limits<double>::max()

    # trust_region_preprocessor.cc IsProgramValid: parameters must be finite. On failure the
    # summary keeps its default termination (FAILURE) and no iterations.
    (_lm_valid(x1in) & _lm_valid(x2in)) || return FAILURE, x1in, x2in, 0
    # reorder_program.cc: AMD on the 2x2 block Hessian is the identity, so the parameter
    # order stays (x1, x2) and SimplicialLDLT runs with natural ordering.

    # TrustRegionMinimizer::IterationZero -> EvaluateGradientAndJacobian.
    x1, x2 = x1in, x2in
    R1, R2 = _lm_jets(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, x1, x2)
    _lm_jets_valid(R1, R2) || return FAILURE, x1in, x2in, 0
    r1, r2 = R1.a, R2.a
    x_cost = z + (z + T(0.5) * (r1 * r1 + r2 * r2))
    isfinite(x_cost) || return FAILURE, x1in, x2in, 0
    # Jacobian, row-major per column block: J[row, col].
    J11, J21, J12, J22 = R1.v1, R2.v1, R1.v2, R2.v2
    # Gradient J' r (small_blas.h MatrixTransposeVectorMultiply).
    g1 = z + (z + ((z + J11 * r1) + J21 * r2))
    g2 = z + (z + ((z + J12 * r1) + J22 * r2))
    # Jacobi scaling, fixed at iteration 0: scale = 1 / (1 + sqrt(diag(J'J)))
    # (BlockSparseMatrix::SquaredColumnNorm), then ScaleColumns.
    s1 = one(T) / (one(T) + sqrt(z + (J11 * J11 + J21 * J21)))
    s2 = one(T) / (one(T) + sqrt(z + (J12 * J12 + J22 * J22)))
    J11, J21, J12, J22 = J11 * s1, J21 * s1, J12 * s2, J22 * s2
    # gradient_max_norm = |x - Plus(x, -gradient)|_inf.
    gmax = _lm_pmax(abs(x1 - (x1 + (-g1))), abs(x2 - (x2 + (-g2))))

    # Minimizer state (trust_region_minimizer.cc Init).
    p1, p2 = x1, x2                      # parameters_: best point, returned to the user
    minimum_cost = dmax
    num_consecutive_invalid_steps = 0
    atleast_one_successful_step = false
    # LevenbergMarquardtStrategy state.
    radius = T(1.0e4)
    decrease_factor = T(2)
    reuse_diagonal = false
    d1 = d2 = z
    # TrustRegionStepEvaluator state (use_nonmonotonic_steps = false -> max steps 0).
    ev_minimum_cost = ev_current_cost = ev_reference_cost = ev_candidate_cost = x_cost
    acc_reference_mcc = acc_candidate_mcc = z
    num_consecutive_nonmonotonic_steps = 0

    iteration = 0
    step_is_successful = true            # iteration 0 counts as a successful step
    niter = 0                            # summary.iterations.size()
    while true
        # FinalizeIterationAndCheckIfMinimizerCanContinue.
        if step_is_successful && x_cost < minimum_cost
            minimum_cost = x_cost
            p1, p2 = x1, x2
        end
        niter += 1
        iteration >= max_num_iterations && return NO_CONVERGENCE, p1, p2, niter
        (step_is_successful && !(gmax > gradient_tolerance)) && return CONVERGENCE, p1, p2, niter
        !(radius > min_trust_region_radius) && return CONVERGENCE, p1, p2, niter

        iteration += 1
        step_is_successful = false

        # ---- ComputeTrustRegionStep -> LevenbergMarquardtStrategy::ComputeStep ----------
        if !reuse_diagonal
            # diagonal = SquaredColumnNorm(J), clamped to [min_lm_diagonal, max_lm_diagonal].
            d1 = z + (J11 * J11 + J21 * J21)
            d2 = z + (J12 * J12 + J22 * J22)
            d1 = ifelse(d1 < min_lm_diagonal, min_lm_diagonal, d1)
            d2 = ifelse(d2 < min_lm_diagonal, min_lm_diagonal, d2)
            d1 = ifelse(max_lm_diagonal < d1, max_lm_diagonal, d1)
            d2 = ifelse(max_lm_diagonal < d2, max_lm_diagonal, d2)
        end
        D1 = sqrt(d1 / radius)
        D2 = sqrt(d2 / radius)
        reuse_diagonal = true

        # SparseNormalCholeskySolver::SolveImpl: rhs = J' r; lhs = J'J + D'D, built by
        # InnerProductComputer::Compute over the Jacobian row block then the appended
        # diagonal row blocks (lower triangle).
        b1 = z + ((z + J11 * r1) + J21 * r2)
        b2 = z + ((z + J12 * r1) + J22 * r2)
        A11 = (z + ((z + J11 * J11) + J21 * J21)) + (z + D1 * D1)
        A21 = z + ((z + J12 * J11) + J22 * J21)
        A22 = (z + ((z + J12 * J12) + J22 * J22)) + (z + D2 * D2)

        # Eigen SimplicialLDLT<Upper, NaturalOrdering>: factorize_preordered (LDLT, shift
        # scale 1 and offset 0) then _solve_impl (L, D^-1, L').
        solved = false
        q1 = q2 = z
        y1 = z + A11
        dd1 = y1 * one(T) + z
        if dd1 != z
            y1 = z + A21
            y2 = z + A22
            dd2 = y2 * one(T) + z
            l21 = y1 / dd1
            dd2 = dd2 - l21 * y1
            if dd2 != z
                q1, q2 = b1, b2
                if q1 != z
                    q2 = q2 - q1 * l21
                end
                q1 = (one(T) / dd1) * q1
                q2 = (one(T) / dd2) * q2
                q1 = q1 - l21 * q2
                solved = _lm_valid(q1) & _lm_valid(q2)
            end
        end

        step_is_valid = false
        st1 = st2 = model_cost_change = z
        if solved
            # step = -solution; model_residuals = J * step (BlockSparseMatrix
            # RightMultiplyAndAccumulate); model_cost_change = -(r + m/2) . m.
            st1, st2 = -q1, -q2
            m1 = (z + (z + J11 * st1)) + (z + J12 * st2)
            m2 = (z + (z + J21 * st1)) + (z + J22 * st2)
            model_cost_change = -(z + (z + (m1 * (r1 + m1 / T(2)) + m2 * (r2 + m2 / T(2)))))
            step_is_valid = model_cost_change > z
        end

        if !step_is_valid
            # HandleInvalidStep; LevenbergMarquardtStrategy::StepIsInvalid == StepRejected.
            num_consecutive_invalid_steps += 1
            num_consecutive_invalid_steps >= max_num_consecutive_invalid_steps &&
                return FAILURE, x1in, x2in, niter
            radius = radius / decrease_factor
            decrease_factor *= T(2)
            reuse_diagonal = true
            continue
        end
        num_consecutive_invalid_steps = 0
        delta1, delta2 = st1 * s1, st2 * s2

        # ComputeCandidatePointAndEvaluateCost: Plus is x + delta; a failed evaluation
        # counts as an infinitely bad step.
        c1, c2 = x1 + delta1, x2 + delta2
        ok, candidate_cost = _lm_cost(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, c1, c2)
        ok || (candidate_cost = dmax)

        # ParameterToleranceReached: |x - candidate| <= tol * (|x| + tol).
        if atleast_one_successful_step
            x_norm = sqrt(x1 * x1 + x2 * x2)
            e1, e2 = x1 - c1, x2 - c2
            step_norm = sqrt(e1 * e1 + e2 * e2)
            step_norm > parameter_tolerance * (x_norm + parameter_tolerance) ||
                return CONVERGENCE, p1, p2, niter
        end
        # FunctionToleranceReached: |cost change| <= function_tolerance * cost.
        abs(x_cost - candidate_cost) > function_tolerance * x_cost ||
            return CONVERGENCE, p1, p2, niter

        # IsStepSuccessful -> TrustRegionStepEvaluator::StepQuality.
        if candidate_cost >= dmax
            relative_decrease = -dmax
        else
            rd = (ev_current_cost - candidate_cost) / model_cost_change
            hrd = (ev_reference_cost - candidate_cost) / (acc_reference_mcc + model_cost_change)
            relative_decrease = ifelse(rd < hrd, hrd, rd)
        end

        if relative_decrease > min_relative_decrease
            # HandleSuccessfulStep: move, then re-evaluate residuals, Jacobian, gradient.
            atleast_one_successful_step = true
            x1, x2 = c1, c2
            R1, R2 = _lm_jets(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, x1, x2)
            _lm_jets_valid(R1, R2) || return FAILURE, x1in, x2in, niter
            r1, r2 = R1.a, R2.a
            x_cost = z + (z + T(0.5) * (r1 * r1 + r2 * r2))
            isfinite(x_cost) || return FAILURE, x1in, x2in, niter
            J11, J21, J12, J22 = R1.v1, R2.v1, R1.v2, R2.v2
            g1 = z + (z + ((z + J11 * r1) + J21 * r2))
            g2 = z + (z + ((z + J12 * r1) + J22 * r2))
            J11, J21, J12, J22 = J11 * s1, J21 * s1, J12 * s2, J22 * s2
            gmax = _lm_pmax(abs(x1 - (x1 + (-g1))), abs(x2 - (x2 + (-g2))))
            step_is_successful = true

            # LevenbergMarquardtStrategy::StepAccepted.
            t = one(T) - upow(T(2) * relative_decrease - one(T), T(3))
            third = one(T) / T(3)
            radius = radius / ifelse(third < t, t, third)
            radius = ifelse(radius < max_trust_region_radius, radius, max_trust_region_radius)
            decrease_factor = T(2)
            reuse_diagonal = false

            # TrustRegionStepEvaluator::StepAccepted.
            ev_current_cost = candidate_cost
            acc_candidate_mcc += model_cost_change
            acc_reference_mcc += model_cost_change
            if ev_current_cost < ev_minimum_cost
                ev_minimum_cost = ev_current_cost
                num_consecutive_nonmonotonic_steps = 0
                ev_candidate_cost = ev_current_cost
                acc_candidate_mcc = z
            else
                num_consecutive_nonmonotonic_steps += 1
                if ev_current_cost > ev_candidate_cost
                    ev_candidate_cost = ev_current_cost
                    acc_candidate_mcc = z
                end
            end
            if num_consecutive_nonmonotonic_steps == 0
                ev_reference_cost = ev_candidate_cost
                acc_reference_mcc = acc_candidate_mcc
            end
        else
            # LevenbergMarquardtStrategy::StepRejected.
            radius = radius / decrease_factor
            decrease_factor *= T(2)
            reuse_diagonal = true
        end
    end
    return
end

"""
    kamiya_ceres(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, x1, x2) -> (usable, x1, x2, iterations)

Solve upstream's Kamiya radius equations (CostOptPos) from the start point `(x1, x2)`
exactly as upstream `kamyiaOpt` does with Ceres Solver 2.2.0 and default options.
`usable` is `Solver::Summary::IsSolutionUsable()`; when it is false the start point is
returned unchanged, as Ceres restores the user's parameters. `iterations` is
`summary.iterations.size()`. Allocation free and generic in `T`.
"""
@inline function kamiya_ceres(
        gamma::T, dp1::T, dp2::T, f0::T, f1::T, f2::T, l0::T, l1::T,
        l2::T, x1::T, x2::T
    ) where {T <: AbstractFloat}
    term, y1, y2, n = _kamiya_ceres_full(gamma, dp1, dp2, f0, f1, f2, l0, l1, l2, x1, x2)
    return term != 2, y1, y2, n
end
