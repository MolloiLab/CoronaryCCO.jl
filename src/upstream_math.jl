# SPDX-License-Identifier: GPL-3.0-or-later
# Scalar and 3-vector math that must round exactly like the upstream C++ build.

# Julia's `^(::Float64, ::Float64)` is a pure-Julia implementation, not the C library's
# `pow`. Upstream calls `pow`/`std::pow` from the platform libm, so the Float64 CPU path
# calls the same symbol. Float32 runs use exp(y·log x) in Float32 for positive x: Julia's
# Float32 `^` widens to a Float64 pow, which was the largest single cost of candidate
# evaluation; the fastmath forms call libm expf/logf on CPUs and native instructions on GPUs. The relative error is a few
# ulp times |y·log x|, far below every tolerance the growth uses. Other bases (negative
# bases with integer exponents occur in the trust-region update) fall back to `^`. A negative
# base with a non-integer exponent gives NaN, as C's `pow` does (Julia's `^` throws): the
# Levenberg–Marquardt evaluation then sees an invalid step and rejects it, as Ceres does.
@inline cpow(x::Float64, y::Float64) = ccall(:pow, Float64, (Float64, Float64), x, y)
@inline upow(x::Float64, y::Float64) = cpow(x, y)
@inline upow(x::Float32, y::Float32) = (x > 0.0f0 && isfinite(x)) ? @fastmath(exp(y * log(x))) :
    (x < 0.0f0 && !isinteger(y)) ? NaN32 : x^y

# Points are plain 3-tuples. DGtal's PointVector does every operation componentwise; its
# dot product and L2 norm accumulate left to right starting from +0.0, which fixes the sign
# of a zero result, so the sums below start from zero(T) as well. Julia never fuses a*b+c
# into an FMA unless asked, matching an upstream build compiled with -ffp-contract=off.
const P3{T} = NTuple{3, T}

@inline add3(a::P3, b::P3) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline sub3(a::P3, b::P3) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline mul3(a::P3, s) = (a[1] * s, a[2] * s, a[3] * s)
@inline div3(a::P3, s) = (a[1] / s, a[2] / s, a[3] / s)
@inline dot3(a::P3{T}, b::P3{T}) where {T} = ((zero(T) + a[1] * b[1]) + a[2] * b[2]) + a[3] * b[3]
@inline norm3(a::P3{T}) where {T} = sqrt(((zero(T) + a[1] * a[1]) + a[2] * a[2]) + a[3] * a[3])
@inline norm3(x::T, y::T, z::T) where {T} = norm3((x, y, z))

# static_cast<int>(double): truncation toward zero (aarch64 fcvtzs semantics for NaN/overflow).
@inline itrunc(x::AbstractFloat) = unsafe_trunc(Int32, x)
@inline itrunc3(p::P3) = (itrunc(p[1]), itrunc(p[2]), itrunc(p[3]))
