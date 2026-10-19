# SPDX-License-Identifier: GPL-3.0-or-later
# Bitwise parity of `kamiya_ceres` with upstream's Ceres call (src/kamiya_lm.jl).
#
# Reference records come from the pinned Ceres 2.2.0 build (test/oracle/build.sh):
#   test/oracle/out/ceres/ceres_n2000.txt, ceres_n4000.txt  every kamyiaOpt call made by
#       HctCco (CCO_SEED=1, case0, -n 2000/4000, CCO_CERES_DUMP=...);
#   ceres_fuzz.txt, ceres_wild.txt  synthetic inputs solved by test/oracle/build/bin/kamiya_ref
#       (`kamiya_ref fuzz 200000 7`, `kamiya_ref wild 200000 11`), covering rejected steps,
#       failed evaluations, diagonal clamping, iteration and trust-region limits.
# Each line: gamma dP1 dP2 f0 f1 f2 l0 l1 l2 x1_in x2_in x1_out x2_out usable iterations
# termination_type, floats as C99 hex. The files are regenerable and not committed; the
# test is skipped when they are absent.
using Vasculature
using Test

const CERES_DIR = joinpath(@__DIR__, "oracle", "out", "ceres")

function check_ceres_file(path; show = 5)
    total = exact = 0
    for line in eachline(path)
        f = split(line)
        length(f) == 16 || continue
        v = ntuple(i -> parse(Float64, f[i]), 13)
        usable = f[14] == "1"
        iters = parse(Int, f[15])
        term = parse(Int, f[16])
        ok, y1, y2, it = Vasculature.kamiya_ceres(v[1:11]...)
        t, = Vasculature._kamiya_ceres_full(v[1:11]...)
        total += 1
        if y1 === v[12] && y2 === v[13] && ok == usable && it == iters && t == term
            exact += 1
        elseif total - exact <= show
            @info "Ceres mismatch" path line = total expected = (v[12], v[13], usable, iters, term) got = (y1, y2, ok, it, t)
        end
    end
    return exact, total
end

@testset "kamiya_ceres matches Ceres 2.2.0 bitwise" begin
    files = isdir(CERES_DIR) ? sort(
            filter(
                f -> startswith(basename(f), "ceres_") && endswith(f, ".txt"),
                readdir(CERES_DIR; join = true)
            )
        ) : String[]
    if isempty(files)
        @info "No Ceres reference records in $CERES_DIR; run test/oracle/build.sh and regenerate them."
        @test_skip false
    end
    for path in files
        exact, total = check_ceres_file(path)
        println("  ", rpad(basename(path), 22), lpad(exact, 7), " / ", total, " exact")
        @test total > 0
        @test exact == total
    end
end

# Measured through top-level functions: a closure defined inside @testset may not be inferred,
# which boxes the returned tuple at the call site and reports a false allocation.
kamiya_run(a) = Vasculature.kamiya_ceres(a...)
kamiya_alloc(a) = @allocated kamiya_run(a)

@testset "kamiya_ceres is allocation free and generic" begin
    args = (3.0, 290196.27, 193949.4, 12.495, 8.33, 4.165, 1.7, 3.73, 3.19, 0.008, 0.0113)
    kamiya_alloc(args)
    @test kamiya_alloc(args) == 0
    args32 = Float32.(args)
    kamiya_alloc(args32)
    @test kamiya_alloc(args32) == 0
    ok64, a64, b64, _ = kamiya_run(args)
    ok32, a32, b32, _ = kamiya_run(args32)
    @test ok64 && ok32
    @test a32 isa Float32 && b32 isa Float32
    @test isapprox(a32, a64; rtol = 1.0e-4) && isapprox(b32, b64; rtol = 1.0e-4)
    # A non-finite start fails in preprocessing and leaves the start point untouched.
    @test Vasculature.kamiya_ceres(args[1:9]..., NaN, 0.01) === (false, NaN, 0.01, 0)
end
