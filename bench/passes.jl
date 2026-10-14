# SPDX-License-Identifier: GPL-3.0-or-later
# Feasibility benchmark: cost of the irreducible per-candidate whole-tree passes that
# dominate HCT-CCO at full scale, on CPU threads and (optionally) one GPU backend.
#
# Per Kamiya iteration of every candidate bifurcation, upstream recomputes
#   (1) radii top-down   r[i] = r1 * B[i],  B[i] = B[parent] * beta[i]  (unary: B[parent])
#   (2) volume bottom-up v[i] = (r[i]*r[i])*l[i] + (v[L] + v[R])        (unary: no add)
#   (3) three capsule-vs-tree intersection scans.
# Sibling betas change along the whole path to the root, so (1)-(2) touch every segment;
# keeping upstream's operation order makes them O(N) per evaluation.
#
# Usage: julia -t auto --project=bench bench/passes.jl [nseg] [ncand] [gpu]
#        gpu ∈ {none, metal}; GPU work is capped to a few seconds.

import AcceleratedKernels as AK
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @synchronize, @groupsize
using Random
const GPU = get(ARGS, 3, "none")
GPU == "metal" && @eval using Metal

# Synthetic tree grown with upstream's split operation: segment s gains children
# (newLeft, newRight); newLeft inherits s's old children. The first 10% of terminals
# attach anywhere, then only to leaves or parents of leaves (upstream stage 3).
function grow_tree(nterm; rng = Xoshiro(1))
    parent = Int32[1, 1]; left = Int32[0, 0]; right = Int32[0, 0]
    isleaf = Bool[false, true]; leaves = Int32[2]
    for k in 2:nterm
        s = if k < nterm ÷ 10
            rand(rng, 2:length(parent))
        else
            t = rand(rng, leaves); rand(rng, Bool) && parent[t] > 1 ? parent[t] : t
        end
        nl, nr = Int32(length(parent) + 1), Int32(length(parent) + 2)
        push!(parent, s, s); push!(left, left[s], 0); push!(right, right[s], 0)
        left[nl] != 0 && (parent[left[nl]] = nl); right[nl] != 0 && (parent[right[nl]] = nl)
        left[s], right[s] = nl, nr
        push!(isleaf, isleaf[s], true); isleaf[s] = false
        isleaf[nl] && (leaves[findfirst(==(s), leaves)] = nl)
        push!(leaves, nr)
    end
    return parent, left, right
end

# Breadth-first order (parents before children) and level offsets.
function levels(parent, left, right)
    order = Int32[2]; offs = Int32[1]; lo = 1
    while lo <= length(order)
        hi = length(order)
        for j in lo:hi
            i = order[j]
            left[i] != 0 && push!(order, left[i]); right[i] != 0 && push!(order, right[i])
        end
        push!(offs, hi + 1); lo = hi + 1
    end
    return order, offs
end

# ---------- CPU: one candidate per task, serial passes in BFS order ----------
function eval_candidate!(r, cum, v, order, parent, left, right, beta, len, r1::T) where {T}
    cum[2] = one(T); r[2] = r1
    @inbounds for j in 2:length(order)
        i = order[j]; p = parent[i]
        unary = (left[p] == 0) ⊻ (right[p] == 0)
        cum[i] = unary ? cum[p] : cum[p] * beta[i]
        r[i] = r1 * cum[i]
    end
    @inbounds for j in length(order):-1:1
        i = order[j]; L = left[i]; R = right[i]
        vi = (r[i] * r[i]) * len[i]
        v[i] = (L != 0 && R != 0) ? vi + (v[L] + v[R]) : vi
    end
    return @inbounds v[2]
end

function cpu_batch!(out, R, CUM, V, order, parent, left, right, BETA, len, r1)
    AK.foreachindex(out; max_tasks = Threads.nthreads()) do c
        out[c] = eval_candidate!(view(R, :, c), view(CUM, :, c), view(V, :, c), order,
                                 parent, left, right, view(BETA, :, c), len, r1)
    end
    return out
end

# ---------- GPU: one workgroup per candidate, levels swept with barriers ----------
@kernel function gpu_eval!(out, R, CUM, V, @Const(order), @Const(offs), nlev,
                           @Const(parent), @Const(left), @Const(right), @Const(BETA),
                           @Const(len), r1)
    c = @index(Group, Linear); t = @index(Local, Linear); G = @groupsize()[1]
    T = eltype(R)
    if t == 1
        CUM[2, c] = one(T); R[2, c] = r1
    end
    @synchronize
    for lev in 2:nlev
        a = offs[lev]; b = offs[lev+1] - 1
        j = a + t - 1
        while j <= b
            i = order[j]; p = parent[i]
            unary = (left[p] == 0) ⊻ (right[p] == 0)
            CUM[i, c] = unary ? CUM[p, c] : CUM[p, c] * BETA[i, c]
            R[i, c] = r1 * CUM[i, c]
            j += G
        end
        @synchronize
    end
    for lev in nlev:-1:1
        a = offs[lev]; b = offs[lev+1] - 1
        j = a + t - 1
        while j <= b
            i = order[j]; L = left[i]; Rr = right[i]
            vi = (R[i, c] * R[i, c]) * len[i]
            V[i, c] = (L != 0 && Rr != 0) ? vi + (V[L, c] + V[Rr, c]) : vi
            j += G
        end
        @synchronize
    end
    if t == 1
        out[c] = V[2, c]
    end
end

# ---------- capsule-vs-tree scan (upstream CoronaryGeometry::isIntersecting) ----------
@inline function seg2seg_dist(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)
    T = typeof(ax)
    d1x, d1y, d1z = bx - ax, by - ay, bz - az
    d2x, d2y, d2z = dx - cx, dy - cy, dz - cz
    rx, ry, rz = cx - ax, cy - ay, cz - az
    D1 = d1x * d1x + d1y * d1y + d1z * d1z
    D2 = d2x * d2x + d2y * d2y + d2z * d2z
    Rr = d1x * d2x + d1y * d2y + d1z * d2z
    S1 = d1x * rx + d1y * ry + d1z * rz
    S2 = d2x * rx + d2y * ry + d2z * rz
    den = D1 * D2 - Rr * Rr
    t = den > 0 ? clamp((S1 * D2 - S2 * Rr) / den, zero(T), one(T)) : zero(T)
    u = D2 > 0 ? (t * Rr - S2) / D2 : zero(T)
    if u < 0 || u > 1
        u = clamp(u, zero(T), one(T))
        t = D1 > 0 ? clamp((u * Rr + S1) / D1, zero(T), one(T)) : zero(T)
    end
    ex = d1x * t - d2x * u - rx; ey = d1y * t - d2y * u - ry; ez = d1z * t - d2z * u - rz
    return sqrt(ex * ex + ey * ey + ez * ez)
end

function scan!(hit, Q, x, y, z, parent, rad)
    fill!(hit, false)
    AK.foreachindex(x) do i
        i == 1 && return
        p = parent[i]
        for q in 1:size(Q, 2)
            d = seg2seg_dist(Q[1, q], Q[2, q], Q[3, q], Q[4, q], Q[5, q], Q[6, q],
                             x[i], y[i], z[i], x[p], y[p], z[p])
            d <= Q[7, q] + rad[i] && (hit[q] = true)
        end
    end
    return hit
end

# ---------- CPU, preorder layout: parents precede children, subtrees contiguous ----------
function preorder(parent, left, right)
    n = length(parent); perm = Int32[]; stack = Int32[2]
    while !isempty(stack)
        i = pop!(stack); push!(perm, i)
        right[i] != 0 && push!(stack, right[i]); left[i] != 0 && push!(stack, left[i])
    end
    inv = zeros(Int32, n); inv[perm] .= 1:length(perm)
    m(j) = j == 0 ? Int32(0) : inv[j]
    p2 = Int32[i == 1 ? Int32(1) : inv[parent[perm[i]]] for i in eachindex(perm)]
    p2[1] = 1
    l2 = Int32[m(left[perm[i]]) for i in eachindex(perm)]
    r2 = Int32[m(right[perm[i]]) for i in eachindex(perm)]
    # 1 if this node's parent is unary (no beta multiply), precomputed per node
    un = Bool[i == 1 ? false : ((l2[p2[i]] == 0) ⊻ (r2[p2[i]] == 0)) for i in eachindex(perm)]
    return perm, p2, l2, r2, un
end

function eval_seq!(r, cum, v, parent, left, right, unary, beta, len, r1::T) where {T}
    n = length(parent)
    cum[1] = one(T); r[1] = r1
    @inbounds for i in 2:n
        p = parent[i]
        cum[i] = unary[i] ? cum[p] : cum[p] * beta[i]
        r[i] = r1 * cum[i]
    end
    @inbounds for i in n:-1:1
        L = left[i]; R = right[i]
        vi = (r[i] * r[i]) * len[i]
        v[i] = (L != 0 && R != 0) ? vi + (v[L] + v[R]) : vi
    end
    return @inbounds v[1]
end

function bench(f, reps)
    f()                                   # compile / warm up
    t = @elapsed for _ in 1:reps; f(); end
    return t / reps
end

function main(nseg_target = 80_000, ncand = 20, gpu = "none")
    T = Float32
    parent, left, right = grow_tree(nseg_target ÷ 2)
    n = length(parent)
    order, offs = levels(parent, left, right)
    nlev = length(offs) - 1
    rng = Xoshiro(2)
    BETA = rand(rng, T, n, ncand) .* T(0.3) .+ T(0.6)
    len = rand(rng, T, n) .+ T(0.1)
    R = zeros(T, n, ncand); CUM = similar(R); V = similar(R); out = zeros(T, ncand)
    widths = diff(offs)
    println("segments=$n  candidates=$ncand  depth(levels)=$nlev  max level width=$(maximum(widths))  threads=$(Threads.nthreads())")

    t1 = bench(() -> eval_candidate!(view(R, :, 1), view(CUM, :, 1), view(V, :, 1), order,
                                     parent, left, right, view(BETA, :, 1), len, T(1.5)), 20)
    println("CPU 1 candidate, radius+volume passes:        $(round(t1 * 1e6, digits = 1)) µs")
    tb = bench(() -> cpu_batch!(out, R, CUM, V, order, parent, left, right, BETA, len, T(1.5)), 20)
    println("CPU $ncand candidates in parallel (threads):     $(round(tb * 1e6, digits = 1)) µs")

    perm, p2, l2, r2, un = preorder(parent, left, right)
    B2 = BETA[perm, :]; len2 = len[perm]
    t2 = bench(() -> eval_seq!(view(R, 1:length(perm), 1), view(CUM, 1:length(perm), 1), view(V, 1:length(perm), 1),
                               p2, l2, r2, un, view(B2, :, 1), len2, T(1.5)), 50)
    same = eval_seq!(view(R, 1:length(perm), 1), view(CUM, 1:length(perm), 1), view(V, 1:length(perm), 1),
                     p2, l2, r2, un, view(B2, :, 1), len2, T(1.5)) ===
           eval_candidate!(view(R, :, 2), view(CUM, :, 2), view(V, :, 2), order, parent, left, right, view(BETA, :, 1), len, T(1.5))
    println("CPU 1 candidate, preorder layout:             $(round(t2 * 1e6, digits = 1)) µs  (bitwise same volume: $same)")
    out2 = zeros(T, ncand)
    t3 = bench(() -> AK.foreachindex(out2; max_tasks = Threads.nthreads()) do c
            m = length(perm)
            out2[c] = eval_seq!(view(R, 1:m, c), view(CUM, 1:m, c), view(V, 1:m, c), p2, l2, r2, un, view(B2, :, c), len2, T(1.5))
        end, 50)
    println("CPU $ncand candidates, preorder, threads:        $(round(t3 * 1e6, digits = 1)) µs")

    x, y, z = rand(rng, T, n) .* 100, rand(rng, T, n) .* 100, rand(rng, T, n) .* 100
    rad = rand(rng, T, n) .* T(0.2)
    Q = rand(rng, T, 7, 3ncand) .* 100; Q[7, :] .= T(0.3)
    hit = zeros(Bool, 3ncand)
    ts = bench(() -> scan!(hit, Q, x, y, z, parent, rad), 20)
    println("CPU capsule scan, $(3ncand) capsules × $n segs:  $(round(ts * 1e6, digits = 1)) µs")

    if gpu == "metal"
        backend = Metal.MetalBackend()
        dev(a) = (b = KA.allocate(backend, eltype(a), size(a)); copyto!(b, a); b)
        dR, dC, dV, dout = dev(R), dev(CUM), dev(V), dev(out)
        dorder, doffs, dpar, dl, dr = dev(order), dev(offs), dev(parent), dev(left), dev(right)
        dB, dlen = dev(BETA), dev(len)
        for G in (256, 1024)
            k = gpu_eval!(backend, G)
            run() = (k(dout, dR, dC, dV, dorder, doffs, Int32(nlev), dpar, dl, dr, dB, dlen,
                       T(1.5); ndrange = G * ncand); KA.synchronize(backend))
            tg = bench(run, 10)
            ok = Array(dout) ≈ cpu_batch!(out, R, CUM, V, order, parent, left, right, BETA, len, T(1.5))
            println("Metal $ncand candidates, 1 workgroup each (G=$G): $(round(tg * 1e6, digits = 1)) µs  (matches CPU: $ok)")
        end
        dx, dy, dz, drad, dQ, dhit = dev(x), dev(y), dev(z), dev(rad), dev(Q), dev(hit)
        tgs = bench(() -> (scan!(dhit, dQ, dx, dy, dz, dpar, drad); KA.synchronize(backend)), 10)
        println("Metal capsule scan:                            $(round(tgs * 1e6, digits = 1)) µs")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse(Int, get(ARGS, 1, "80000")), parse(Int, get(ARGS, 2, "20")), get(ARGS, 3, "none"))
end
