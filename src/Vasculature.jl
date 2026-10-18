# SPDX-License-Identifier: GPL-3.0-or-later
# Vasculature.jl. The tree growth engine is a Julia port of HCT-CCO
# (Vajire, Choy, Kassab, Lee), which extends OpenCCO (Kerautret, Ngo, Passat, Talbot, Jaquet).
# See NOTICE.md for provenance.

"""
    Vasculature

Vendor-neutral vascular network tools: grow → flow → mesh + voxelize + materialize. Tree
growth is a port of HCT-CCO (Hierarchy CT-Driven Constrained Constructive Optimization).

Upstream: <https://github.com/Sujal-vajire/HCT-CCO> (GPL-3.0-or-later), itself built on
OpenCCO <https://github.com/kerautret/OpenCCO>. Please cite both; see `NOTICE.md`.
"""
module Vasculature

import AcceleratedKernels as AK
import KernelAbstractions as KA
import Mmap

include("upstream_math.jl")
include("params.jl")
include("rng.jl")
include("io_vol.jl")
include("domain.jl")
include("geometry.jl")
include("model.jl")
include("io_gxl.jl")
include("kamiya_lm.jl")
include("kamiya.jl")
include("layout.jl")
include("grid.jl")
include("grow.jl")
include("fast.jl")
include("view.jl")
include("batched.jl")
include("morphometry.jl")
include("qc.jl")
include("generate.jl")
include("postprocess.jl")
include("xcat/voxelize.jl")
include("spec.jl")
# perfusion model: flow from pressures, territories, whole-body contrast, coronary flow capacity, CT frames, phantom
include("perfusion/flow.jl")
include("perfusion/territory.jl")
include("perfusion/pipeline.jl")
include("perfusion/body.jl")
include("perfusion/transport.jl")
include("perfusion/cfc.jl")
include("perfusion/organct.jl")
include("perfusion/dynamic_phantom.jl")

export CCOParams, CCOModel, MaskDomain, Engine, nsegments
export GrowthSpec, grow, size_scale, terminal_count, mesh_domain, write_vol
export FlowSpec, flow, couple_veins, vtree, Unified
export TransportSpec, transport, circulation
export read_vol, import_gxl, export_gxl, grow!, dds_orders, tree_qc, set_territories!, set_forbidden!

end
