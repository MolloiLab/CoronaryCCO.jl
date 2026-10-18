# Vasculature.jl

[![CI](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolloiLab/Vasculature.jl/actions/workflows/CI.yml)
[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3+-blue.svg)](LICENSE)

| Heart: contrast through the grown coronary tree | Lung: contrast through the grown pulmonary trees |
|:---:|:---:|
| ![heart](docs/assets/heart_flow.webp) | ![lung](docs/assets/lung_flow.webp) |

Grow a complete vascular tree inside any organ from the few vessels a scan shows, compute its blood flow from the
pressures alone, send a contrast injection through it, and export a voxelized perfusion phantom for CT simulation.
One algorithm for every organ; an organ differs only in the numbers you give it.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/MolloiLab/Vasculature.jl")
```

## 1. Grow

![imaged input → grown tree](docs/assets/step1_grow.png)

You give four things; everything else is the same for every organ:

1. **The organ**: a closed surface mesh of the tissue to perfuse (`mesh_domain` turns it into a 1 mm voxel domain), or a
   ready voxel mask (`.vol`). Terminals are placed only inside it.
2. **The inlets**: the vessels the scan shows, one tree per inlet, as GXL centrelines with radii. Growth starts from
   them and never moves them.
3. **The regions to avoid** (optional): anything inside or next to the organ that a vessel must not cross, such as the
   heart's chambers and great vessels, the lung's airways or the brain's ventricles. Give them as a fine voxel mask
   (`true` = forbidden) with its origin and cell size; every new vessel segment is checked against it.
4. **A `GrowthSpec`**: how many terminals, the tissue volume, the imaged inlet radii.

```julia
using Vasculature

domain = mesh_domain("heart.vol", vertices, triangles) # organ surface → 1 mm domain
# an organ is an instance of the one GrowthSpec struct (HEART_GROWTH, LUNG_GROWTH, ...); a run varies it
spec = GrowthSpec(HEART_GROWTH; device = MtlArray)
avoid = (chambers, origin_mm, 0.4) # forbidden mask, centre of its first cell (mm), cell size (mm)
tree = grow(spec, domain, ["lca.xml", "rca.xml"]; blockers = avoid)
```

On a laptop GPU (Apple M4) the heart grows 10⁶ terminals in about 40 seconds.

`flow_ml_min` sets how many terminals grow: every organ is resolved down to terminals of the same size (one carries
0.14 µL/min, about 12 µm across, the resolution of the heart tree scored against Kassab's morphometry), so the count is
the organ's flow over that. For tissue the organ feeds, the flow is its O₂ demand by Fick; for the lungs, the cardiac
output. You can also give `nterm` directly.

| | `flow_ml_min` | → terminals | `tissue_ml` | inlets |
|---|---|---|---|---|
| heart (adult, rest) | 135.6 | 957,045 | 148 (myocardium) | LCA, RCA |
| lungs (adult) | 4,800 (cardiac output) | 33.9 million | 2,800 | main pulmonary artery |
| brain (older adult, rest) | 510 (van Es et al. 2010) | 3.6 million | 1,016 (parenchyma) | both internal carotids, both vertebrals |
| one leg (adult, rest) | 290 (Koch Esteves et al. 2023) | 2.0 million | 11,000 | common femoral artery |
| any organ | VO₂ ÷ (arterial O₂ content × O₂ extraction) | flow ÷ 0.14 µL/min | its segmented volume | each imaged inlet (femoral artery, carotids + basilar, ...) |

More inlets are just more files: `grow(spec, brain, ["lica.xml", "rica.xml", "basilar.xml"])`.
Veins grow the same way from their imaged outlet trees (`grow(spec, domain, ["left_pulmonary_veins.xml", "right_pulmonary_veins.xml"])`).

## 2. Flow

![flow from pressure](docs/assets/step2_flow.png)

```julia
arteries, veins = vtree(tree), vtree(veintree) # grown trees as vessel lists
F = flow(HEART_FLOW, arteries, veins, tissue) # tissue: 3 × n voxel centres (mm)
F.rest.resist.mbf, F.stress.resist.mbf # mL/min/g at rest and at maximal dilation
```

No flow is prescribed: each vessel's flow follows from the inlet and outlet pressures, the vessel sizes and the
arterioles' tone. `HEART_FLOW` and `LUNG_FLOW` are instances of the one `FlowSpec` struct, every value
measured and cited; a variant is `FlowSpec(HEART_FLOW; p_in = 70)`.

## 3. Contrast

![contrast transport](docs/assets/step3_transport.png)

A clinical injection (dose, rate, saline flush, injection site) travels through a whole-body circulation into the
organ's grown tree, then through every capillary bed and vein, with exchange into the tissue and washout by the kidneys.
(Being folded into one call, `transport(spec, F, injection)`, like `grow` and `flow`.)

## 4. Phantom

![phantom](docs/assets/step4_phantom.png)

The result is one compartment volume (one integer per voxel) and one CSV per scenario: each compartment's base material
and its iodine concentration every second, plus the ground-truth blood flow. Load it at any time point as material
maps for a CT simulator such as [BasisSimulator.jl](https://github.com/MolloiLab/BasisSimulator.jl).
(Being folded into one call, `phantom(spec, ...)`, like `grow` and `flow`.)

![lung CT](docs/assets/step5_ct.png)

## How to cite

The growth algorithm is HCT-CCO, extended from OpenCCO. If you use tree growth, please cite:

- S. L. Vajire, J. S. Choy, G. S. Kassab, L.-C. Lee. Hierarchy aware CT-driven constrained constructive optimization:
  a framework for coronary network generation. *Ann. Biomed. Eng.*, 2026. doi:10.1007/s10439-026-04383-4
- B. Kerautret, P. Ngo, N. Passat, H. Talbot, C. Jaquet. OpenCCO: an implementation of constrained constructive
  optimization for generating 2D and 3D vascular trees. *Image Process. On Line* 13:258–279, 2023.
  doi:10.5201/ipol.2023.477

## License

GPL-3.0-or-later.
