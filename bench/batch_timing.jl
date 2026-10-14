# SPDX-License-Identifier: GPL-3.0-or-later
# Wall time of batched growth on case0 (corrected options) on CPU threads and, optionally, a GPU.
#
#   julia -t 8 --project=bench bench/batch_timing.jl 100000,800000 [n_seq] [metal | metalonly]
#
# `gpu` is `metal` (Apple, unified memory) here; any GPU array constructor works the same way.
using Vasculature, Printf, Metal   # Metal is only a bench dependency

const D = joinpath(pkgdir(Vasculature), "test", "data")
grow(n; kw...) = generate(
    joinpath(D, "case0.vol"), joinpath(D, "case0.xml"); nterm = n, aperf = 10000, min_dist = 1,
    root = (-26, 41, 13), seed = 1, T = Float32, gamma = 2.33, volume = :full, metric = :euclidean,
    sampling = :continuous, fix_base = true, collision = :final, max_branch_angle = 90,
    growth = :batched, kw...
)

function main()
    ns = parse.(Int, split(get(ARGS, 1, "100000"), ","))
    nseq = parse(Int, get(ARGS, 2, "2000"))
    gpu = get(ARGS, 3, "")
    devs = Pair{String, Any}[]
    gpu == "metalonly" || push!(devs, "cpu" => identity)
    if gpu in ("metal", "metalonly")
        push!(devs, "metal" => (x -> Metal.MtlArray{eltype(x), ndims(x), Metal.SharedStorage}(x)))
    end
    for (lab, dev) in devs
        grow(3000; n_seq = 1000, device = dev)          # compile
        for n in ns
            t = @elapsed E = grow(n; n_seq = nseq, device = dev)
            s = E.stats
            @printf(
                "%-6s n=%7d: %7.1f s | locations %.1f, candidates %.1f (snapshot %.1f), commit %.1f (refresh %.1f), %d rounds\n",
                lab, n, t, s.t_location, s.t_candidates, s.t_view, s.t_commit, s.t_refresh, s.rounds
            )
            flush(stdout)
        end
    end
    return
end

main()
