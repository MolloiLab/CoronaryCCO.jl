# SPDX-License-Identifier: GPL-3.0-or-later
# Coronary flow capacity (Johnson & Gould): every region classified by its stress flow and its
# flow reserve, so the model's normal, reduced, ischemic and steal tissue can be counted the way
# PET reads patients. Thresholds come from test/reference (Rb-82 dipyridamole PET).

module CFC


export CFCThresholds, cfc_thresholds, cfc_class, CFC_NAMES

const CFC_NAMES = ("normal", "minimally reduced", "mildly reduced", "moderately reduced", "severely reduced (ischemia)", "steal", "scar")

"""
    CFCThresholds

Upper bounds of CFR and stress flow (mL/min/g) for the reduced classes, from minimally to
severely reduced. A region belongs to the most severe class whose two bounds both contain it
(nested rectangles, Johnson & Gould 2012 Fig. 1; the severe class is stated this way explicitly
in Gould et al. 2019). CFR < 1 is steal; scar is rest flow ≤ `scar` (2012) or rest and stress
both ≤ `scar` (2018 onward).
"""
struct CFCThresholds
    cfr::NTuple{4, Float32}
    stress::NTuple{4, Float32}
    scar::Float32
    scar_rest_only::Bool         # 2012: rest flow alone below the scar line; 2018+: rest and stress
    name::String
end

"""
    cfc_thresholds(set = :current) -> CFCThresholds

`:current` is the 2018-onward set (Gould et al. 2018; test/reference/gould2018_cfc_revised_thresholds.toml),
`:original` the 2012 map (test/reference/johnson2012_coronary_flow_capacity.toml).
"""
function cfc_thresholds(set::Symbol = :current)
    # the category boundaries (CFR and stress flow, mL/min/g, normal -> severe), as the reference files give them
    # (test/reference; the test suite checks these against them)
    set === :current && return CFCThresholds((2.9f0, 2.38f0, 1.6f0, 1.27f0), (2.17f0, 1.82f0, 1.09f0, 0.83f0), 0.3f0, false, "Gould 2018")
    # scar: rest flow below the lowest iso-rest-flow line, 0.2 cc/min/g (2012 Fig. 1)
    set === :original && return CFCThresholds((3.37f0, 2.7f0, 2.03f0, 1.74f0), (2.39f0, 1.76f0, 1.12f0, 0.91f0), 0.2f0, true, "Johnson & Gould 2012")
    throw(ArgumentError("set must be :current or :original"))
end

"""
    cfc_class(stress, rest, th) -> Int

1 normal, 2 minimally, 3 mildly, 4 moderately, 5 severely reduced, 6 steal, 7 scar.
"""
function cfc_class(stress::Real, rest::Real, th::CFCThresholds)
    rest <= th.scar && (th.scar_rest_only || stress <= th.scar) && return 7
    cfr = stress / max(rest, 1.0e-9)
    cfr < 1 && return 6
    for k in 4:-1:1
        cfr <= th.cfr[k] && stress <= th.stress[k] && return k + 1
    end
    return 1
end

end # module
