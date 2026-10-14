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
include("tree.jl")

export CCOParams, CCOTree, nsegments

end
