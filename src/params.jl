# SPDX-License-Identifier: GPL-3.0-or-later
# Ports the parameters of HctccoModel (src/HctccoModel.h) and the HctCco CLI defaults
# (bin/HctCco.cpp) from HCT-CCO commit 9bc141f.

"""
    CCOParams(; kwargs...)

Run parameters. Defaults match the upstream `HctCco` command line and `HctccoModel`.

| Field | Upstream | Default |
|---|---|---|
| `nterm` | `-n, --nbTerm` | 500 |
| `aperf` | `-a, --aPerf` | 20000.0 |
| `gamma` | `-g, --gamma` | 3.0 |
| `root_radius_max` | `--rootRadiusMax` | 1.588 |
| `term_flow_cv` | `--termFlowCV` | 0.0 |
| `side_branch_bias` | `--sideBranchBias` | 0.0 |
| `min_distance_to_border` | `-m, --minDistanceToBorder` | 5.0 |
| `extension_mode` | `--extensionModel` | false |
| `anywhere_mode` | `--anywhereModel` | false |
| `num_neighbors` | `myNumNeighbor` | 20 |
| `nb_max_search` | `expandTree(nbMaxSearch)` | 100 |
| `nb_try_candidate` | `expandTree(nbTryCandite)` | 100 |
| `kamiya_iters` | `isAddable(nbIter)` | 100 |
| `kamiya_tol` | `isAddable(tolerance)` | 0.01 |
"""
Base.@kwdef struct CCOParams
    nterm::Int = 500
    aperf::Float64 = 20000.0
    gamma::Float64 = 3.0
    root_radius_max::Float64 = 1.588
    term_flow_cv::Float64 = 0.0
    side_branch_bias::Float64 = 0.0
    min_distance_to_border::Float64 = 5.0
    extension_mode::Bool = false
    anywhere_mode::Bool = false
    num_neighbors::Int = 20
    nb_max_search::Int = 100
    nb_try_candidate::Int = 100
    kamiya_iters::Int = 100
    kamiya_tol::Float64 = 0.01
end

# Physiological constants hard-coded in HctccoModel.h.
const MU = 3.6e-3      # my_mu, blood viscosity
const QPERF = 8.33e3   # my_qPerf
const PPERF = 1.33e4   # my_pPerf
