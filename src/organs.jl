# SPDX-License-Identifier: GPL-3.0-or-later
# Organs as instances of the one GrowthSpec and TransportSpec structs: each organ is only different values (the FlowSpec
# instances, HEART_FLOW and LUNG_FLOW, are in perfusion/organs.jl). Subject values are XCAT vmale_50's. A run varies an
# instance, e.g. a seed and a device such as `MtlArray`: GrowthSpec(LUNG_GROWTH; seed = 2, device = dev). [ids] are data/reference files.

# heart arteries: the scored v16 forest (LCA + RCA from XCAT's imaged dias_lad1/dias_rca1 trees)
const HEART_GROWTH = GrowthSpec(;
    # the reference heart's count: its resting Fick flow (135.6 mL/min) resolved to the terminals Kassab's morphometry is scored at
    nterm = 957_045, flow_ml_min = nothing, terminal_flow_ml_min = 135.6 / 957_045,
    # mL, XCAT vmale_50 myocardium (LV + RV); also the size HCT-CCO's constants were set for
    tissue_ml = 148.1, reference_ml = 148.1,
    wall_clearance = 0.5,          # mm, HCT-CCO default
    gamma = 2.33,                  # 7/3 from coronary structure-function scaling (Huo & Kassab 2012)
    root_caps = [2.54, 1.59],      # mm, LCA and RCA ostial radii (decision A3, QUEUE)
    root_fixed = false,
    # human coronary minimum daughter symmetry by mother diameter: geometric mean of each bin, lower whisker
    # [schwarz2020_human_coronary_symmetry]
    symmetry_floor = [
        (36.742346141747674, 0.278), (51.96152422706632, 0.411), (67.08203932499369, 0.347), (82.15838362577492, 0.314),
        (97.21111047611791, 0.272), (112.24972160321825, 0.242), (134.16407864998737, 0.175), (173.20508075688772, 0.131),
        (282.842712474619, 0.068), (565.685424949238, 0.019),
    ],
    seed = 1, T = Float32, device = identity,
)

# heart veins: the scored v9 tree into the coronary sinus ostium (grown with a request of 100,000; it holds 99,989)
const HEART_VEIN_GROWTH = GrowthSpec(;
    nterm = 100_000, flow_ml_min = nothing, terminal_flow_ml_min = NaN, tissue_ml = 148.1, reference_ml = 148.1,
    wall_clearance = 1.0, gamma = 2.33,
    root_caps = [2.7397784811568675],   # mm, coronary sinus ostium (veins/ostium.txt)
    root_fixed = false, symmetry_floor = [], seed = 1, T = Float32, device = identity,
)

# lung arteries: one tree from the main pulmonary artery (XCAT's imaged radius). nterm 200,000 is the pipeline's test size
# (QUEUE: H-ACINUS replaces it by the lung's own unit count)
const LUNG_GROWTH = GrowthSpec(;
    nterm = 200_000, flow_ml_min = nothing, terminal_flow_ml_min = NaN,
    tissue_ml = 2813.438,               # mL, XCAT's two lungs (lung.vol)
    reference_ml = 148.1, wall_clearance = 0.5, gamma = 2.33,
    root_caps = [15.411707],            # mm, main pulmonary artery (pulmonary_arteries.xml)
    root_fixed = true, symmetry_floor = [], seed = 1, T = Float32, device = identity,
)

# lung veins: one tree per lung into the left atrium (156,000: the pipeline's test size)
const LUNG_VEIN_GROWTH = GrowthSpec(
    LUNG_GROWTH; nterm = 156_000,
    root_caps = [8.170588, 8.170588]
)   # mm, the two pulmonary vein ostia (pulmonary_veins_left/right.xml)

# a bolus through the whole body with the heart as its coronary bed (rest: XCAT's heart, 81.6 mL/s at 60/min)
const HEART_TRANSPORT = TransportSpec(;
    phys = HEART_FLOW, tissue_ml = 148.052f0, bed = :coronary, input = :lv,   # mL, XCAT myocardium voxels
    co = 81.6f0, beat = 1.0f0, f_coronary = 0.033f0, tmax = 70.0f0, dt = 0.002f0, every = 0.1f0, coronary = nothing,
)

# a bolus with the lungs in their own body slot, fed by the right ventricle; the body's coronary bed is the heart's (lumped)
const LUNG_TRANSPORT = TransportSpec(
    HEART_TRANSPORT; phys = LUNG_FLOW, tissue_ml = 2813.438f0, bed = :lung, input = :rv,
    coronary = (; ml = 148.052f0, fb = HEART_FLOW.fb, ve = HEART_FLOW.ve, ps = HEART_FLOW.ps, rho = HEART_FLOW.rho)
)
