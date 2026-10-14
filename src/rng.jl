# SPDX-License-Identifier: GPL-3.0-or-later
# Portable stand-in for C rand(), identical to cco_rand() in oracle/cco_oracle.h.
# Upstream seeds rand() from the wall clock; both sides use this generator instead so a
# seed fixes the candidate stream on every platform.

mutable struct CRand
    state::UInt64
end
CRand(seed::Integer = 0) = CRand(UInt64(seed))

"SplitMix64; returns the top 31 bits, i.e. an Int32 in [0, 2^31 - 1] like glibc/macOS rand()."
@inline function crand!(r::CRand)
    r.state += 0x9E3779B97F4A7C15
    z = r.state
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    z = z ⊻ (z >> 31)
    return Int32(z >> 33)
end

reseed!(r::CRand, seed::Integer) = (r.state = UInt64(seed); r)
