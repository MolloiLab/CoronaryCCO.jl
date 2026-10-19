# SPDX-License-Identifier: GPL-3.0-or-later
# Bitwise parity against the instrumented upstream build (test/oracle/build.sh). Each dump line is
# one segment: index parent left right kterm x y z radius resistance flow beta, floats in
# C99 hex. Generate with, e.g.,
#   cd test/oracle/out/parity && CCO_SEED=1 CCO_DUMP=tree_n500.txt ../../build/bin/HctCco \
#     -d ../../build/src/samples/case0.vol -i ../../build/src/samples/case0.xml \
#     -n 500 -a 10000 -m 1 -p -26 41 13 -x final_tree.xml --no-postprocess

const PARITY_DIR = joinpath(@__DIR__, "oracle", "out", "parity")
# Committed reference dump. Bits depend on the platform's libm `pow`, so it is only checked on
# the platform it was generated on (Apple Silicon macOS, with test/oracle/build.sh).
const PARITY_FIXTURES = joinpath(@__DIR__, "fixtures", "parity", "darwin-aarch64")

"Compare a grown Float64 model with an oracle dump; returns (mismatching segments, first mismatch)."
function compare_dump(m, path)
    rows = [split(l) for l in eachline(path)]
    length(rows) == nsegments(m) || return (abs(length(rows) - nsegments(m)), "segment count $(length(rows)) vs $(nsegments(m))")
    hx(s) = parse(Float64, s)
    up(j) = j == 0 ? 0 : j - 1
    bad = 0; first = ""
    for (k, r) in enumerate(rows)
        j = parse(Int, r[1]) + 1
        got = (
            up(m.parent[j]), up(m.left[j]), up(m.right[j]), Int(m.kterm[j]), m.pos[j]..., m.radius[j],
            m.resistance[j], m.flow[j], m.beta[j],
        )
        want = (parse(Int, r[2]), parse(Int, r[3]), parse(Int, r[4]), parse(Int, r[5]), hx.(r[6:12])...)
        if !all(got .=== want)
            bad += 1
            bad == 1 && (first = "segment $(j - 1): got $got want $want")
        end
    end
    return bad, first
end

@testset "bitwise parity with upstream HCT-CCO (Float64)" begin
    dumps = String[]
    for dir in (PARITY_DIR, (Sys.isapple() && Sys.ARCH === :aarch64) ? PARITY_FIXTURES : "")
        isdir(dir) || continue
        for f in readdir(dir)
            occursin(r"^tree_n\d+\.txt$", f) && !any(endswith(d, "/" * f) for d in dumps) && push!(dumps, joinpath(dir, f))
        end
    end
    if isempty(dumps)
        @info "no oracle dumps in $PARITY_DIR; skipping parity (see test/oracle/build.sh)"
        @test_skip false
    end
    for f in dumps
        n = parse(Int, match(r"tree_n(\d+)", f).captures[1])
        E = grow_case0(n; T = Float64, mode = :exact, upstream_attach = true)   # upstream's isAddable rule (deviations §15)
        bad, first = compare_dump(E.m, f)
        bad == 0 || @info "parity n=$n: $bad segments differ; first: $first"
        @test bad == 0
    end
end
