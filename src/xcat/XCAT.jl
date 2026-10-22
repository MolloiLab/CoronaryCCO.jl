# SPDX-License-Identifier: GPL-3.0-or-later
# XCAT inputs: read NURBS surfaces (.nrb), take centrelines of imaged vessels and join them into inlet trees, add their
# penetrating stubs, voxelize organ and label masks, fill supply gaps. `Vasculature.XCAT.read_nrb(path)`, `centerline`,
# `base_forest`, `penetrating_stubs`, ...
module XCAT

using ..Vasculature: VoxelGrid, voxelize, grid_around, center, write_vol
import AcceleratedKernels as AK
using Printf

include("nurbs.jl")
include("heart.jl")
include("stubs.jl")
include("supply_gaps.jl")

end
