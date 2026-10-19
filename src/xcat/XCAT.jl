# SPDX-License-Identifier: GPL-3.0-or-later
# XCAT inputs: read NURBS surfaces (.nrb), take centrelines of imaged vessels and join them into inlet trees, voxelize
# organ and label masks, fill supply gaps. `Vasculature.XCAT.read_nrb(path)`, `centerline`, `base_forest`, ...
module XCAT

using ..Vasculature: VoxelGrid, voxelize, grid_around, write_vol
import AcceleratedKernels as AK
using Printf

include("nurbs.jl")
include("heart.jl")
include("supply_gaps.jl")

end
