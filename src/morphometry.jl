# SPDX-License-Identifier: GPL-3.0-or-later
# First-principles morphometry: diameter-defined Strahler (DDS) ordering computed from the
# tree's own diameters (Jiang, Kassab and Fung, J. Appl. Physiol. 76:882-892, 1994; HCT-CCO
# paper Eq. 12), with no literature thresholds and no diameter rescaling.

"Kassab et al. 1993 porcine coronary arterial diameters by order 1..11 (mm): mean and SD
(the MU/SD arrays of upstream's coronary_postprocess.py)."
const LITERATURE_MU = [0.009, 0.013, 0.0187, 0.0346, 0.0716, 0.15, 0.303, 0.467, 0.715, 1.492, 3.176]
const LITERATURE_SD = [0.0009, 0.0017, 0.0026, 0.0074, 0.0172, 0.035, 0.054, 0.056, 0.13, 0.3, 0.6]

"""
    dds_orders(m; maxiter = 100, tol = 0.01) -> NamedTuple

Diameter-defined Strahler ordering (Jiang, Kassab and Fung 1994), from the tree's own
diameters:

1. Start from topological Strahler orders (leaves 1; a parent takes the largest child order,
   plus one when two children share it; a single child passes its order through).
2. Elements are maximal connected runs of segments with the same order; an element's diameter
   is the mean of its segments' diameters. Per order, take the mean D̄ₙ and SD σₙ of element
   diameters; the upper boundary of order n is ½[(D̄ₙ + σₙ) + (D̄ₙ₊₁ − σₙ₊₁)] (Eq. 12).
3. Recompute orders from the leaves up: where children meet, the parent starts from the
   largest child order k (k is reached by two or more children, or by one child), and is
   promoted to k + 1 only if its diameter exceeds the upper boundary of order k.
4. Repeat 2-3 until fewer than `tol` of the segments change order.

Returns per-element orders, diameters and lengths, plus per-order count, mean and SD of
element diameters.
"""
function dds_orders(m::CCOModel; maxiter::Int = 100, tol::Real = 0.01)
    n = nsegments(m)
    post = postorder(m)
    d = [j == 1 ? 0.0 : 2 * Float64(m.radius[j]) for j in 1:n]
    order = zeros(Int, n)
    for j in post                      # 1. topological Strahler
        a = m.left[j]; b = m.right[j]
        order[j] = a == 0 && b == 0 ? 1 : a == 0 || b == 0 ? order[a == 0 ? b : a] :
            order[a] == order[b] ? order[a] + 1 : max(order[a], order[b])
    end
    elems = nothing; iters = 0
    while true
        elems = _elements(m, order)
        iters >= maxiter && break
        iters += 1
        up = _upper_bounds(elems, d, order)
        new = zeros(Int, n)
        for j in post                  # 3. diameter-gated promotion at confluences
            a = m.left[j]; b = m.right[j]
            if a == 0 && b == 0
                new[j] = 1
            else
                k = a == 0 ? new[b] : b == 0 ? new[a] : max(new[a], new[b])
                new[j] = (k <= length(up) && d[j] > up[k]) ? k + 1 : k
            end
        end
        changed = count(j -> new[j] != order[j], 2:n) / (n - 1)
        order = new
        changed < tol && (elems = _elements(m, order); break)
    end
    eorder, D, len = elems
    K = maximum(eorder)
    stats = map(1:K) do o
        sel = D[eorder .== o]
        μ = isempty(sel) ? NaN : sum(sel) / length(sel)
        (
            order = o, count = length(sel), mean = μ,
            sd = isempty(sel) ? NaN : sqrt(sum(abs2, sel .- μ) / length(sel)),
        )
    end
    return (order = eorder, diameter = D, length = len, stats = stats, iterations = iters)
end

# Elements: maximal connected runs of segments with the same order; returns per-element
# (order, mean diameter, total geometric length).
function _elements(m::CCOModel, order)
    n = nsegments(m)
    elem = zeros(Int, n); ne = 0
    for j in 2:n
        elem[j] != 0 && continue
        ne += 1; elem[j] = ne; stack = [j]
        while !isempty(stack)
            i = pop!(stack)
            for k in (m.parent[i], m.left[i], m.right[i])
                (k <= 1 || elem[k] != 0 || order[k] != order[i]) && continue
                elem[k] = ne; push!(stack, k)
            end
        end
    end
    dsum = zeros(ne); cnt = zeros(Int, ne); len = zeros(ne); eo = zeros(Int, ne)
    for j in 2:n
        e = elem[j]
        dsum[e] += 2 * Float64(m.radius[j]); cnt[e] += 1; len[e] += Float64(m.glen[j]); eo[e] = order[j]
    end
    return eo, dsum ./ cnt, len
end

# Upper boundary of each order from element statistics (Eq. 12); the top order is unbounded.
function _upper_bounds(elems, d, order)
    eo, D, _ = elems
    K = maximum(eo)
    μ = zeros(K); σ = zeros(K)
    for o in 1:K
        sel = D[eo .== o]
        isempty(sel) && continue
        μ[o] = sum(sel) / length(sel); σ[o] = sqrt(sum(abs2, sel .- μ[o]) / length(sel))
    end
    return [((μ[o] + σ[o]) + (μ[o + 1] - σ[o + 1])) / 2 for o in 1:(K - 1)]
end
