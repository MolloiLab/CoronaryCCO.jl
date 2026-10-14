# SPDX-License-Identifier: GPL-3.0-or-later
# Batched vs sequential growth on case0 with the corrected options.
#
#   julia -t 8 --project=bench bench/batch_validate.jl [nterm] [nseeds] [n_seq] [batch_frac] [cpu | metal]
#
# CCO trees grown from different seeds differ in structure but are equally valid (Schreiner
# et al. 1997, doi:10.1006/jtbi.1997.0413). Batched growth passes if every metric of its trees
# falls inside the range spanned by sequential trees over the same seeds. Results go to
# bench/out/batch_validate_<nterm>.md.
using Vasculature, Printf, Statistics, Serialization, Metal   # Metal is only a bench dependency

const D = joinpath(pkgdir(Vasculature), "test", "data")
grow(n, seed; kw...) = generate(
    joinpath(D, "case0.vol"), joinpath(D, "case0.xml"); nterm = n, aperf = 10000, min_dist = 1,
    root = (-26, 41, 13), seed = seed, T = Float32, gamma = 2.33, volume = :full,
    metric = :euclidean, sampling = :continuous, fix_base = true, collision = :final,
    max_branch_angle = 90, kw...
)

# Per-order element count, mean diameter (µm) and mean length (mm).
function orders(m; maxiter)
    R = dds_orders(m; maxiter = maxiter)
    K = maximum(R.order)
    return [
        (
                count = count(==(o), R.order), d = 1000 * mean(R.diameter[R.order .== o]),
                l = mean(R.length[R.order .== o]),
            ) for o in 1:K
    ]
end

function metrics(E, t)
    m = E.m; st = E.stats; n = nsegments(m)
    leaves = [j for j in 2:n if m.left[j] == 0 && m.right[j] == 0]
    dl = sort(2000 .* Float64.(m.radius[leaves]))
    depth = map(leaves) do j
        k = 0
        while j != 1
            k += 1; j = m.parent[j]
        end
        k
    end
    q25 = tree_qc(E; rmin = 0.25); q10 = tree_qc(E; rmin = 0.1)
    return (
        time = t, t_loc = st.t_location, t_cand = st.t_candidates, t_commit = st.t_commit,
        volume = pi * sum(Float64(m.radius[j])^2 * Float64(m.glen[j]) for j in 2:n),
        root_d = 2 * Float64(m.radius[2]), leaf_med = dl[(end + 1) ÷ 2], leaf_p5 = dl[max(1, length(dl) ÷ 20)],
        depth = mean(depth), overlaps25 = q25.overlaps, worst25 = q25.worst_overlap_mm,
        turns90 = q25.grown_turns_90, overlaps10 = q10.overlaps, worst10 = q10.worst_overlap_mm,
        turns90_10 = q10.grown_turns_90, worst_turn = q25.worst_turn_deg,
        strahler = orders(m; maxiter = 0), dds = orders(m; maxiter = 100),
    )
end

function main()
    n = parse(Int, get(ARGS, 1, "100000"))
    ns = parse(Int, get(ARGS, 2, "5"))
    nseq = parse(Int, get(ARGS, 3, "10000"))
    frac = parse(Float64, get(ARGS, 4, "0.01"))
    gpu = get(ARGS, 5, "cpu")
    device = identity
    gpu == "metal" && (device = x -> MtlArray{eltype(x), ndims(x), Metal.SharedStorage}(x))
    grow(300, 1); grow(3000, 1; growth = :batched, n_seq = 1000, device = device)   # compile
    S = []; B = []
    out = joinpath(@__DIR__, "out"); mkpath(out)
    for seed in 1:ns
        cache = joinpath(out, "seq_$(n)_$(seed)_t$(Threads.nthreads()).jls")
        if isfile(cache)
            push!(S, deserialize(cache))
        else
            t = @elapsed E = grow(n, seed); push!(S, metrics(E, t)); serialize(cache, S[end])
        end
        @printf("seed %d sequential %.1f s vol %.2f\n", seed, S[end].time, S[end].volume)
        t = @elapsed E = grow(n, seed; growth = :batched, n_seq = nseq, batch_frac = frac, device = device); push!(B, metrics(E, t))
        @printf("seed %d batched    %.1f s vol %.2f\n", seed, t, B[end].volume)
    end
    io = IOBuffer()
    println(io, "# Batched vs sequential growth, case0, n = $n, $(Threads.nthreads()) threads\n")
    println(io, "Batched ($gpu): n_seq = $nseq, batch_frac = $frac. Seeds 1–$ns. Range = min–max over seeds.\n")
    println(io, "| Metric | Sequential range | Batched range | Batched mean inside sequential range |")
    println(io, "|---|---|---|---|")
    fails = String[]
    row(name, f; fmt = "%.3f", tol = 0.0) = begin
        s = [f(x) for x in S]; b = [f(x) for x in B]
        lo, hi = minimum(s), maximum(s)
        ok = lo - tol <= mean(b) <= hi + tol
        ok || push!(fails, name)
        F = Printf.Format(fmt)
        println(
            io, "| $name | ", Printf.format(F, lo), "–", Printf.format(F, hi), " | ", Printf.format(F, minimum(b)),
            "–", Printf.format(F, maximum(b)), " | ", ok ? "yes" : "**no**", " |"
        )
    end
    row("wall time (s)", x -> x.time; fmt = "%.1f")
    row("total volume (mm³)", x -> x.volume; fmt = "%.2f")
    row("root diameter (mm)", x -> x.root_d)
    row("leaf diameter median (µm)", x -> x.leaf_med; fmt = "%.2f")
    row("leaf diameter p5 (µm)", x -> x.leaf_p5; fmt = "%.2f")
    row("mean terminal depth (segments)", x -> x.depth; fmt = "%.2f")
    row("overlaps r ≥ 0.25 mm", x -> x.overlaps25; fmt = "%.0f", tol = 2.0)
    row("worst overlap r ≥ 0.25 mm (mm)", x -> x.worst25; tol = 0.05)
    row("grown turns > 90°, r ≥ 0.25 mm", x -> x.turns90; fmt = "%.0f")
    row("overlaps r ≥ 0.1 mm", x -> x.overlaps10; fmt = "%.0f", tol = 2.0)
    row("grown turns > 90°, r ≥ 0.1 mm", x -> x.turns90_10; fmt = "%.0f", tol = 1.0)
    for (key, label) in ((:strahler, "Strahler"), (:dds, "DDS"))
        K = minimum(length(getfield(x, key)) for x in vcat(S, B))
        for o in 1:K
            row("$label order $o count", x -> getfield(x, key)[o].count; fmt = "%.0f")
            row("$label order $o diameter (µm)", x -> getfield(x, key)[o].d; fmt = "%.2f")
            row("$label order $o length (mm)", x -> getfield(x, key)[o].l)
        end
        println(
            io, "| $label orders (per tree) | ", join(sort(unique(length(getfield(x, key)) for x in S)), ", "), " | ",
            join(sort(unique(length(getfield(x, key)) for x in B)), ", "), " | |"
        )
    end
    # Kassab (1993) comparison, orders aligned at the top (left main = order 11)
    println(io, "\n## DDS mean diameter vs Kassab et al. 1993 (orders aligned at the root)\n")
    println(io, "| Kassab order | Kassab µm (±SD) | Sequential µm | Batched µm |")
    println(io, "|---|---|---|---|")
    for kk in 11:-1:1
        vals(X) = begin
            v = Float64[]
            for x in X
                K = length(x.dds); o = K - (11 - kk)
                o >= 1 && push!(v, x.dds[o].d)
            end
            isempty(v) ? "—" : @sprintf("%.1f", mean(v))
        end
        @printf(io, "| %d | %.1f (±%.1f) | %s | %s |\n", kk, 1000Vasculature.LITERATURE_MU[kk], 1000Vasculature.LITERATURE_SD[kk], vals(S), vals(B))
    end
    println(
        io, "\nTimings (mean, s): sequential location/candidates/commit ",
        join([@sprintf("%.1f", mean(getfield.(S, f))) for f in (:t_loc, :t_cand, :t_commit)], " / "),
        "; batched ", join([@sprintf("%.1f", mean(getfield.(B, f))) for f in (:t_loc, :t_cand, :t_commit)], " / "), "."
    )
    println(io, "\nMetrics outside the sequential range: ", isempty(fails) ? "none" : join(fails, "; "))
    path = joinpath(out, "batch_validate_$(n)_$(nseq)_$(frac)_$(gpu).md")
    write(path, take!(io))
    println(read(path, String))
    return
end

main()
