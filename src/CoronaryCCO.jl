# SPDX-License-Identifier: GPL-3.0-or-later
# CoronaryCCO.jl: Julia port of HCT-CCO (Vajire, Choy, Kassab, Lee), which extends
# OpenCCO (Kerautret, Ngo, Passat, Talbot, Jaquet). See NOTICE.md for provenance.

"""
    CoronaryCCO

Vendor-neutral Julia port of HCT-CCO (Hierarchy CT-Driven Constrained Constructive
Optimization) for growing subject-specific coronary arterial trees.

Upstream: <https://github.com/Sujal-vajire/HCT-CCO> (GPL-3.0-or-later), itself built on
OpenCCO <https://github.com/kerautret/OpenCCO>. Please cite both; see `NOTICE.md`.
"""
module CoronaryCCO

import AcceleratedKernels as AK
import KernelAbstractions as KA

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
include("morphometry.jl")
include("generate.jl")
include("postprocess.jl")

export CCOParams, CCOModel, MaskDomain, Engine, nsegments
export read_vol, import_gxl, export_gxl, grow!, generate, dds_orders

end
