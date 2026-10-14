# Vasculature.jl

[![CI](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml)
[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3+-blue.svg)](LICENSE)

Pure-Julia, vendor-neutral tools that take an organ mesh or mask to a grown vascular network
and a voxelized material phantom: **grow → flow → mesh + voxelize + materialize**. The phantom
is the hand-off point. CT simulation (for example with
[BasisSimulator.jl](https://github.com/MolloiLab/BasisSimulator.jl)) and perfusion analysis
stay in their own packages and appear here only as worked examples.

The growth engine is a port of **HCT-CCO** (Vajire, Choy, Kassab and Lee), which grows
subject-specific coronary arterial trees from a segmented epicardial base tree inside a
CT-derived myocardial domain. The first target is coronary perfusion; brain and other organs
are planned.

> **The growth engine is a port, not a new method.** The algorithm belongs to the HCT-CCO
> authors, and HCT-CCO itself extends **OpenCCO** (Kerautret, Ngo, Passat, Talbot and Jaquet).
> If you use tree growth from this package, please cite both of their articles (see
> [How to cite](#how-to-cite)). Upstream: <https://github.com/Sujal-vajire/HCT-CCO>.

## What it does

```julia
using Vasculature
E = generate("case0.vol", "case0.xml"; nterm = 40_000, aperf = 10_000, min_dist = 1,
             root = (-26, 41, 13), seed = 1, export_xml = "tree.xml")
postprocess("tree.xml", "tree"; domain_vol = "case0.vol")   # orders, morphometry, hemodynamics, VTP/OBJ
dds_orders(E.m)                                              # orders computed from the tree itself
```

On an Apple M4 laptop CPU (8 threads), case0 with HCT-CCO's sample settings:

| terminals | upstream C++ | Vasculature exact mode | Vasculature (default) |
|---:|---:|---:|---:|
| 4,000 | 81 s | 2.3 s | 1.6 s |
| 40,000 | ≈ 3.6–5 h (extrapolated) | 109 s | 14 s |
| 200,000 | – | – | ≈ 100 s |

- **Exact mode** (test-only) reproduces upstream bit for bit in Float64: on case0, every
  segment's topology, coordinates, radius, resistance, flow and radius ratio match the C++
  reference build at 500 and 2,000 terminals. Ceres's Levenberg–Marquardt solve is replicated
  bit for bit on 514,000 recorded and synthetic cases.
- **The default engine** updates resistances, radii and volume only along the path from each
  new bifurcation to the root (the incremental update the HCT-CCO paper describes), using a
  spatial grid for all geometric queries. It equals exact mode to rounding.
- Kernels are written with AcceleratedKernels.jl / KernelAbstractions.jl, so they run on
  multithreaded CPUs and on CUDA, ROCm, oneAPI and Metal GPUs. Production runs use `Float32`.
- Post-processing (`postprocess`) is a pure-Julia port of HCT-CCO's Python script and
  reproduces its CSVs byte for byte; figures come from a CairoMakie extension.

## Departures from upstream

Upstream behaviour is the default. Where upstream has a bug or differs from the HCT-CCO
paper, an option gives the corrected behaviour: full-tree volume objective (`volume = :full`),
Euclidean wall distance (`metric = :euclidean`), continuous terminal sampling
(`sampling = :continuous`), imaged vessels kept on their centerlines (`fix_base = true`), and
in post-processing the decoded myocardial mass and end-to-end element grouping (defaults).
Each one is documented, with its origin and the supporting literature, in
[`docs/deviations.md`](docs/deviations.md).

With the paper's bifurcation exponent γ = 2.33 and these corrections, growth reaches all 11
Strahler orders without rescaling diameters after growth; see the deviations document for
the terminal counts this requires.

## How to cite

Please cite the upstream work:

1. S. L. Vajire, J. S. Choy, G. S. Kassab and L.-C. Lee, "Hierarchy Aware CT-Driven
   Constrained Constructive Optimization: A Framework for Coronary Network Generation," *Annals of Biomedical Engineering*, 2026.
   <https://doi.org/10.1007/s10439-026-04383-4>. Code: <https://github.com/Sujal-vajire/HCT-CCO>.
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
