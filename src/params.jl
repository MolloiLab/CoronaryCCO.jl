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
| `upstream_attach` | (parity only; see below) | false |
| `max_connection` | `getMaxBranchLength()` (fixed) | 25.0 |
| `num_neighbors` | `myNumNeighbor` | 20 |
| `nb_max_search` | `expandTree(nbMaxSearch)` | 100 |
| `nb_try_candidate` | `expandTree(nbTryCandite)` | 100 |
| `kamiya_iters` | `isAddable(nbIter)` | 100 |
| `kamiya_tol` | `isAddable(tolerance)` | 0.01 |
| `symmetry_floor` | (Vasculature only) | none |

`symmetry_floor` is a list of `(mother diameter µm, minimum d_small/d_large)` knots, at most 12,
interpolated in log diameter and held flat beyond the ends. A candidate junction whose daughters
are more asymmetric than the floor at the mother's final-scale diameter is rejected. Measured
coronary junctions are more symmetric in small vessels than in large ones (VanBavel & Spaan
1992; Kalsho & Kassab 2004; Schwarz et al. 2020), while minimum-volume growth attaches terminal
twigs to vessels of any size. Empty (the default) is upstream behaviour.

Growth is always staged (HCT-CCO's hierarchy-aware schedule; there is no unrestricted mode):
for the first 10% of terminals a new vessel may attach to any segment, up to 35% to a terminal
or one of its three nearest ancestors, then only to a terminal or a terminal's parent
(`extension_mode`: 25% / 100%). The schedule is applied where candidates are chosen, as
upstream's `getN_NearestSegments` does. Upstream's `isAddable` additionally rejects every
non-terminal host at all stages, which empties the first two stages and stalls growth from an
imaged tree with few listed terminals (docs/deviations.md §15). `upstream_attach = true`
restores that check; it exists only for the bitwise parity test.
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
    upstream_attach::Bool = false
    max_connection::Float32 = 25.0f0     # mm; upstream getMaxBranchLength (a coronary value; set per organ)
    num_neighbors::Int = 20
    nb_max_search::Int = 100
    nb_try_candidate::Int = 100
    kamiya_iters::Int = 100
    kamiya_tol::Float64 = 0.01
    symmetry_floor::Vector{Tuple{Float64, Float64}} = Tuple{Float64, Float64}[]
end

# Lowest daughter symmetry d_small/d_large allowed under a mother of diameter d (µm): knots
# (log d, S) interpolated linearly, flat beyond the ends; 0 when there are no knots.
@inline function _symmetry_min(ld::NTuple{12, T}, sv::NTuple{12, T}, n::Int, d::T) where {T}
    n == 0 && return zero(T)
    x = log(d)
    x <= ld[1] && return sv[1]
    @inbounds for k in 2:n
        if x <= ld[k]
            w = (x - ld[k - 1]) / (ld[k] - ld[k - 1])
            return sv[k - 1] + w * (sv[k] - sv[k - 1])
        end
    end
    return sv[n]
end
function _symmetry_knots(::Type{T}, sf) where {T}
    n = min(length(sf), 12)
    length(sf) > 12 && error("symmetry_floor: at most 12 knots")
    issorted(first.(sf)) || error("symmetry_floor: diameters must increase")
    return n, ntuple(i -> i <= n ? T(log(sf[i][1])) : zero(T), 12), ntuple(i -> i <= n ? T(sf[i][2]) : zero(T), 12)
end
function symmetry_min(p::CCOParams, d::T) where {T}
    sf = p.symmetry_floor; n = length(sf)
    n == 0 && return zero(T)
    x = log(Float64(d))
    x <= log(sf[1][1]) && return T(sf[1][2])
    @inbounds for k in 2:n
        x1 = log(sf[k][1])
        if x <= x1
            x0 = log(sf[k - 1][1]); w = (x - x0) / (x1 - x0)
            return T(sf[k - 1][2] + w * (sf[k][2] - sf[k - 1][2]))
        end
    end
    return T(sf[n][2])
end

# Physiological constants hard-coded in HctccoModel.h.
const MU = 3.6e-3      # my_mu, blood viscosity
const QPERF = 8.33e3   # my_qPerf
const PPERF = 1.33e4   # my_pPerf
