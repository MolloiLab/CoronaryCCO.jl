# Vasculature.jl

[![CI](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-Getting%20Started-c4302b.svg)](https://molloilab.github.io/Vasculature.jl/getting-started/)
[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3+-blue.svg)](LICENSE)

| Heart: contrast through the grown coronary tree | Lung: contrast through the grown pulmonary trees |
|:---:|:---:|
| ![heart](docs/assets/heart_flow.webp) | ![lung](docs/assets/lung_flow.webp) |

A CT scan shows an organ's large vessels and nothing smaller. Vasculature.jl grows the rest of the tree, computes how
blood flows through it from the pressures alone, sends a contrast injection through it, and writes a CT phantom in which
the true perfusion of every voxel is known. One code path serves every organ; an organ is only the numbers you give it.

**Start with the [Getting Started guide](https://molloilab.github.io/Vasculature.jl/getting-started/).** It walks from
install to a first phantom, explains where the flow comes from, and shows how to build a new organ or a new disease,
with every rough edge we have hit along the way.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/MolloiLab/Vasculature.jl")
```

Growth runs on the CPU or on any GPU: add `Metal`, `CUDA` or `AMDGPU` for yours.

## The four stages

```julia
using Vasculature

tree = grow(HEART_GROWTH, "heart.vol", ["lca.xml", "rca.xml"])        # 1. grow the tree from the imaged vessels
F = flow(HEART_FLOW, vtree(tree), vtree(veins), tissue)              # 2. blood flow from the pressures
T = transport(HEART_TRANSPORT, F.rest, injection; times)             # 3. a contrast injection through the body
phantom(HEART_PHANTOM, scans, classify, grid; out = "healthy", ...)  # 4. the CT phantom, with its ground truth
```

## Documentation

- [Getting Started](https://molloilab.github.io/Vasculature.jl/getting-started/): the guide, from install to building
  your own organ.
- [Examples](https://molloilab.github.io/Vasculature.jl/examples/): the heart (`00_cardiac_pipeline`), the lungs
  (`01_pulmonary_pipeline`) and the leg (`02_limb_pipeline`), each a Pluto notebook in
  [docs/notebooks](docs/notebooks).
- [API](https://molloilab.github.io/Vasculature.jl/api/): every exported function and spec.

## How to cite

The growth algorithm is HCT-CCO, extended from OpenCCO. If you use tree growth, please cite:

- S. L. Vajire, J. S. Choy, G. S. Kassab, L.-C. Lee. Hierarchy aware CT-driven constrained constructive optimization:
  a framework for coronary network generation. *Ann. Biomed. Eng.*, 2026. doi:10.1007/s10439-026-04383-4
- B. Kerautret, P. Ngo, N. Passat, H. Talbot, C. Jaquet. OpenCCO: an implementation of constrained constructive
  optimization for generating 2D and 3D vascular trees. *Image Process. On Line* 13:258–279, 2023.
  doi:10.5201/ipol.2023.477

## License

GPL-3.0-or-later. The growth engine is a port of HCT-CCO; its provenance, attribution and upstream licenses are
in [docs/NOTICE.md](docs/NOTICE.md).
