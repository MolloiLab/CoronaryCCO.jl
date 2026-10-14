# SPDX-License-Identifier: GPL-3.0-or-later
# Scalar math that must round exactly like the upstream C++ build.

# Julia's `^(::Float64, ::Float64)` is a pure-Julia implementation, not the C library's
# `pow`. Upstream calls `pow`/`std::pow` from the platform libm, so on the CPU we call the
# same symbol. A portable, GPU-compilable pow that rounds like the reference libm is a
# separate work item (see docs/porting-notes.md).
@inline cpow(x::Float64, y::Float64) = ccall(:pow, Float64, (Float64, Float64), x, y)

# DGtal `PointVector::norm()` for a 3D double point: sqrt(x*x + y*y + z*z), summed left to
# right. Julia never contracts a*b+c into an FMA unless asked, matching an upstream build
# compiled without contraction (e.g. x86-64 baseline, or clang with -ffp-contract=off).
@inline norm3(x::Float64, y::Float64, z::Float64) = sqrt(x * x + y * y + z * z)
