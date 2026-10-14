# Feasibility: can 40k terminals go from hours to ≤ 10 minutes?

Short answer: **yes, on CPU alone, with margin.** The projected full 40k-terminal run is
about 1–3 minutes on an Apple M4 laptop CPU, versus roughly 3 hours for an optimized build
of upstream on the same machine. This is a projection from measured upstream work counts
and measured costs of the ported kernels at full tree size, not yet an end-to-end timing
of the finished port.

Everything below can be reproduced with `oracle/build.sh` (upstream reference build) and
`bench/passes.jl` (kernel benchmark).

## 1. Where upstream spends its time

Upstream (`HctccoModel.ih`, `CoronaryGrowth.h`) adds one terminal per step. In every step
it samples a point, takes its 20 nearest connectable segments, and runs `isAddable` on
each. Each Kamiya iteration inside `isAddable`:

- copies the whole tree twice (`cTree1 = *this`, `cTreeCurr = cTree1`);
- recomputes every segment's resistance from the root, with four `pow` calls per node;
- recomputes every radius top-down and the total volume bottom-up;
- scans every segment three times for capsule intersections.

The nearest-segment search sorts all segments and then calls `std::find` over the terminal
list for each one, which is O(N·T) per sampled point in the late stage.

Measured on the `case0` sample (oracle build: Release, `-ffp-contract=off`, seed 1):

| n terminals | wall (s) | Kamiya iters / step | isAddable / step | neighbors / step | max depth |
|---:|---:|---:|---:|---:|---:|
| 500 | 2.4 | 15.9 | 21.9 | 30.8 | 35 |
| 1000 | 6.1 | 17.0 | 18.1 | 26.1 | 46 |
| 2000 | 19.9 | 17.7 | 15.3 | 23.3 | 64 |
| 4000 | 80.9 | 19.6 | 14.3 | 21.9 | 68 |
| 8000 | 374.3 | 20.6 | 13.8 | 21.3 | 79 |

The work per step is roughly constant, about 20 full-tree evaluations, while the cost of
each evaluation grows with tree size. Upstream's per-step time is 4.2–5.7 µs per segment
and rising (11.5 ms per step at about 2.7k segments, 60.7 ms at about 10.7k, 93.7 ms at
about 14.7k). Integrated
to 40k terminals (about 80k segments), that is about 3 hours on an M4. Fitting the wall
times directly gives the same answer: 4000 → 8000 terminals costs 4.6× (exponent ≈ 2.2),
so 8000 → 40000 costs 35–48×, i.e. **≈ 3.6–5 hours** for the full run. This matches the
upstream README's "order of hours" and the ~5 h reported on other hardware.

A `sample` profile at n = 8000 puts most of the time in the recursive whole-tree updates.
`pow` alone takes 23% of the samples. Ceres does not appear as a hotspot.

## 2. What the port changes, without changing any decision

| Upstream cost | Port | Same bits? |
|---|---|---|
| Two whole-tree copies per Kamiya iteration | Candidate state as a small delta over the shared tree | Yes: no arithmetic involved |
| Full resistance recompute with 4 `pow` per node | Recompute only the path from the new bifurcation to the root | Yes: every other node has identical inputs, so recomputing it reproduces the same bits |
| Radius top-down and volume bottom-up via recursion | Same formulas and operation order, swept over a preorder array layout | Yes: verified, `===` on the volume |
| Capsule scans over all segments | Conservative uniform grid, then the exact upstream predicate | Yes: the any-hit boolean is unchanged |
| Sort all segments plus O(N·T) `std::find` | Grid k-nearest search with connectability flags | Yes, with one caveat below |
| Candidates evaluated one after another | Candidates evaluated in parallel; their results are reduced in upstream order | Yes: the minimum is taken in upstream's candidate order |

Radii and volume really do need a full pass per evaluation. When a new bifurcation is
tried, every sibling along the path to the root gets a new β. So every radius in the tree
changes, and the total volume is a sum over all of them. These passes are what remains
irreducible.

## 3. Measured kernel costs at full size (80k segments, 20 candidates, Float32)

From `bench/passes.jl` on an Apple M4 (8 threads, 10-core GPU):

| Kernel | Time |
|---|---:|
| Radius + volume, 1 candidate, BFS order | 1.18 ms |
| Radius + volume, 1 candidate, **preorder layout** | **0.68 ms** (bitwise same volume) |
| Radius + volume, 20 candidates, preorder, CPU threads | **1.97 ms** |
| Radius + volume, 20 candidates, Metal, one workgroup each | 3.8 ms |
| Capsule scan, 60 capsules × 80k segments, brute force, CPU | 16.7 ms |
| Capsule scan, same, brute force, Metal | 3.0 ms |

## 4. Projection for the full run

Integrated over the run (tree size grows from about 250 to 80k segments, so about 1.6×10⁹
segment-steps in total):

- **Radius + volume passes:** 20 evaluations per step × 8.5 ns per segment, with a
  measured 6.9× parallel speedup on 8 threads. **≈ 40 s.**
- **Path-only resistance:** about 70 nodes × 4 `pow` per evaluation. **≈ 1–5 s.**
- **Capsule, nearest-neighbor and distance queries on a grid:** each query touches a few
  hundred segments instead of 80k. **≈ 5–10 s.** Without a grid these would cost about
  130 s even on the GPU, so the grid is required, not optional.
- **Kamiya solve:** about 67 Levenberg–Marquardt iterations per step on a 2×2 system.
  **≈ 1–3 s.**
- **Tree mutation in preorder** (insert two segments): O(N) memmove per step. **≈ 5 s.**

**Total ≈ 1–3 minutes on the M4 CPU.** That leaves room for the work per step to grow
past the measured ~20 evaluations at 40k terminals, and for the overheads not modelled
here. A data-center GPU is not needed to hit the target. A GPU can still take the
brute-force scans and wide candidate batches if the grid or CPU becomes the bottleneck.

## 5. Risks and open items

1. **Float32 convergence margin.** Kamiya stops when the volume change between
   iterations falls below an absolute 0.01. The total volume is about 5×10³ at n = 4000.
   A Float32 ulp there is about 6×10⁻⁴, so the margin is 16×, but it shrinks as the tree
   grows. Plan: in Float32 runs, compute the volume change of each candidate directly
   rather than differencing two large totals. The Float64 parity path keeps upstream's
   exact formulation.
2. **Implementation-defined tie order.** Upstream sorts segments by distance with
   `std::sort`, which is unstable. Exact distance ties happen when a point's nearest
   location is a node shared by several segments. libc++ and libstdc++ break such ties
   differently, so upstream itself is not reproducible across platforms. The oracle
   should switch to `std::stable_sort` (ties broken by index), and the port will match
   that.
3. **Ceres parity.** The Float64 parity path must reproduce Ceres' Levenberg–Marquardt
   (default options, Jacobi scaling, dense QR, autodiff of `pow`) step for step on the 2×2
   Kamiya system. This bounds how exact the Float64 check can be; it does not affect
   speed. Plan: dump Ceres inputs and outputs from the oracle and test the Julia solver
   against millions of recorded cases.
4. **`pow` on GPU.** Vendor `pow` implementations differ in their last bits. For
   cross-backend bitwise identity in Float32, `pow`, `exp` and `log` will be implemented
   in plain Julia arithmetic, which every backend compiles the same way.
5. **Metal has no Float64.** On the local M4 GPU only Float32 runs are possible, which is
   all production needs. The Float64 parity tests run on the CPU.
