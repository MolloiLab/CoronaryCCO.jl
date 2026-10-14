# CoronaryCCO.jl

[![CI](https://github.com/MolloiLab/CoronaryCCO.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolloiLab/CoronaryCCO.jl/actions/workflows/CI.yml)
[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3+-blue.svg)](LICENSE)

A pure-Julia, vendor-neutral port of **HCT-CCO** (Vajire, Choy, Kassab and Lee), which
grows subject-specific coronary arterial trees from a segmented epicardial base tree
inside a CT-derived myocardial domain.

> **This is a port, not a new method.** The algorithm belongs to the HCT-CCO authors, and
> HCT-CCO itself extends **OpenCCO** (Kerautret, Ngo, Passat, Talbot and Jaquet). If you
> use this package, please cite both of their articles (see [How to cite](#how-to-cite)).
> Upstream: <https://github.com/Sujal-vajire/HCT-CCO>.

## Why a port

Upstream needs on the order of hours to grow a full 40,000-terminal tree. This port aims
for minutes on the same problem. It keeps the algorithm's decisions and changes how the
work is executed:

- The whole-tree hydraulic and volume updates run in a cache-friendly preorder layout,
  across candidate bifurcations in parallel.
- Resistance updates touch only the path from the new bifurcation to the root. The rest
  of the tree has unchanged inputs, so recomputing it would reproduce the same bits.
- Brute-force scans over every segment are replaced by conservative spatial queries
  that return the same answers.
- Kernels use [AcceleratedKernels.jl](https://github.com/JuliaGPU/AcceleratedKernels.jl)
  (and [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) where a
  custom kernel is needed), so the same code runs on multithreaded CPUs and on CUDA, ROCm,
  oneAPI and Metal GPUs.

## Numerical contract

- Every kernel is generic over the float type `T`. Production runs use `Float32` on every
  backend.
- `Float32` runs are bitwise identical across backends (CPU threads, Metal, CUDA, ...).
- The test suite also runs the same code at `Float64` on the CPU and compares it bit for
  bit against an instrumented build of upstream HCT-CCO (see [`oracle/`](oracle)). This
  proves the port reproduces upstream's logic exactly, decision by decision.
- Upstream draws candidate points from `rand()`. Both sides replace it with the same
  portable SplitMix64 generator, so a seed defines one candidate stream everywhere.

## Status

Early development. Package skeleton, parameter and tree layout, the upstream reference
build, and the feasibility benchmark are in place. The growth loop port is next. See
[`docs/feasibility.md`](docs/feasibility.md) for the measured speed-up budget.

## How to cite

Please cite the upstream work:

1. S. L. Vajire, J. S. Choy, G. S. Kassab and L.-C. Lee, "Hierarchy CT-Driven Constrained
   Constructive Optimization: A Framework for Coronary Network Generation." (Citation to
   be updated upon publication.) Code: <https://github.com/Sujal-vajire/HCT-CCO>.
2. B. Kerautret, P. Ngo, N. Passat, H. Talbot and C. Jaquet, "OpenCCO: An Implementation
   of Constrained Constructive Optimization for Generating 2D and 3D Vascular Trees,"
   *Image Processing On Line*, 13, 258–279, 2023.
   <https://doi.org/10.5201/ipol.2023.477>

and, as appropriate, the methods they build on: Schreiner and Buxbaum (1993), Karch et al.
(1999) and Kamiya and Togawa (1972). Full references are in [`NOTICE.md`](NOTICE.md) and
[`CITATION.cff`](CITATION.cff).

## License

GPL-3.0-or-later, the same terms as HCT-CCO. See [`LICENSE`](LICENSE) and
[`NOTICE.md`](NOTICE.md) for provenance and the licenses of upstream dependencies.
