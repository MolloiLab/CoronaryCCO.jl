# SPDX-License-Identifier: GPL-3.0-or-later
# Organs as instances of the one GrowthSpec and TransportSpec structs: each organ is only different values (the FlowSpec
# instances, HEART_FLOW and LUNG_FLOW, are in perfusion/organs.jl). Subject values are XCAT vmale_50's. A run varies an
# instance, e.g. a seed and a device such as `MtlArray`: GrowthSpec(LUNG_GROWTH; seed = 2, device = dev). [ids] are test/reference files.

# heart arteries: the scored v16 forest (LCA + RCA from XCAT's imaged dias_lad1/dias_rca1 trees)
const HEART_GROWTH = GrowthSpec(;
    # the reference heart's count: its resting Fick flow (135.6 mL/min) resolved to the terminals Kassab's morphometry is scored at
    nterm = 957_045, flow_ml_min = nothing, terminal_flow_ml_min = 135.6 / 957_045,
    # mL, XCAT vmale_50 myocardium (LV + RV); also the size HCT-CCO's constants were set for
    tissue_ml = 148.1, reference_ml = 148.1,
    wall_clearance = 0.5,          # mm, HCT-CCO default
    gamma = 2.33,                  # 7/3 from coronary structure-function scaling (Huo & Kassab 2012)
    # mm, LCA and RCA ostial radii: XCAT's left main (2.54) and its proximal RCA lumen (11.3 mm², 1.90), past the tubes'
    # closing end caps
    root_caps = [2.54, 1.9],
    root_fixed = false,
    # the two coronaries are daughters of the aortic root and split its flow by the same law as every other branch
    # (RCA share 1.90^2.33 / (1.90^2.33 + 2.54^2.33) = 34%)
    split = :bifurcation,
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
    root_fixed = false, split = :free, symmetry_floor = [], seed = 1, T = Float32, device = identity,
)

# lung arteries: one tree from the main pulmonary artery (XCAT's imaged radius), grown to the lung's own unit, one terminal
# per acinus (the vessels inside each acinus are LUNG_FLOW's sub-terminal tables). The lung carries the whole cardiac
# output, so its terminals are sized by its tissue's gas-exchange units, not by a flow per terminal.
const LUNG_GROWTH = GrowthSpec(;
    # acini (transitional bronchioles) in both lungs, 2 x 14,956 per lung by micro-CT stereology; one terminal artery and
    # vein per acinus as in Clark 2011 [vasilescu2020_human_terminal_transitional_bronchiole_count, clark2011_acinar_terminal_definition]
    nterm = 29_912, flow_ml_min = nothing, terminal_flow_ml_min = NaN,
    tissue_ml = 2813.438,               # mL, XCAT's two lungs (lung.vol)
    reference_ml = 148.1, wall_clearance = 0.5,
    # per-junction exponent of human pulmonary arteries: 2.3 +/- 0.1 over 1,937 cast bifurcations, 2.31 +/- 0.60 by
    # CT/MRA per bifurcation; 7/3 lies in that range [horsfield1989_human_pulmonary_symmetry, altieri2025_pa_murray_ct_mra]
    gamma = 2.33,
    root_caps = [15.411707],            # mm, main pulmonary artery (pulmonary_arteries.xml)
    root_fixed = true, split = :free, symmetry_floor = [], seed = 1, T = Float32, device = identity,
)

# lung veins: one tree per lung into the left atrium, one venule per acinus (no human per-junction venous exponent is
# measured; the arteries' value)
const LUNG_VEIN_GROWTH = GrowthSpec(
    LUNG_GROWTH;
    root_caps = [8.170588, 8.170588]
)   # mm, the two pulmonary vein ostia (pulmonary_veins_left/right.xml)

const XCAT_MYOCARDIUM_ML = 148.052f0 # mL, XCAT vmale_50's myocardium voxels (LV + RV)

# a bolus through the whole body with the heart as its coronary bed (rest: XCAT's heart, 81.6 mL/s at 60/min)
const HEART_TRANSPORT = TransportSpec(;
    phys = HEART_FLOW, tissue_ml = XCAT_MYOCARDIUM_ML, bed = :coronary, input = :lv,
    co = 81.6f0, beat = 1.0f0, f_coronary = 0.033f0, tmax = 70.0f0, dt = 0.002f0, every = 0.1f0, coronary = nothing,
)

# a bolus with the lungs in their own body slot, fed by the right ventricle; the body's coronary bed is the heart's (lumped)
const LUNG_TRANSPORT = TransportSpec(
    HEART_TRANSPORT; phys = LUNG_FLOW, tissue_ml = 2813.438f0, bed = :lung, input = :rv,
    coronary = (; ml = XCAT_MYOCARDIUM_ML, fb = HEART_FLOW.fb, ve = HEART_FLOW.ve, ps = HEART_FLOW.ps, rho = HEART_FLOW.rho)
)

# ---------------------------------------------------------------------------------------------------------- phantoms
# The XCAT body around an organ, as pools of the whole-body model (shared by every organ's phantom; a classifier names them
# by index). XCAT labels: 19-22 chamber blood, 25/28 arteries, 24 veins, 6 lung, 10 liver, 14 spleen, 4/13 muscle, 23 skin,
# 1 body (fat), 5 rib, 8 cortical bone, 9 spine, 12 stomach, 3/29 soft tissue. Organ beds: blood + extracellular iodine over
# the body model's tissue volumes (ICRP 89 masses; the thorax holds 30% of each peripheral tissue). The great vessels are
# positional (`:out` pools): the arteries downstream of their valve, the descending aorta past the arch at its own flow
# (its `lead_ml`, the blood from the aortic valve to the arch's end, is the subject's: the caller sets it), and the veins
# upstream of the chamber they drain into.
const XCAT_BODY_POOLS = [
    Pool("right atrium blood", "ncat_blood", :ra, :mean), Pool("right ventricle blood", "ncat_blood", :rv, :mean), # 1, 2
    Pool("left atrium blood", "ncat_blood", :la, :mean), Pool("left ventricle blood", "ncat_blood", :lv, :mean), # 3, 4
    Pool("aorta and systemic arteries", "ncat_blood", :lv, :out), Pool("pulmonary arteries", "ncat_blood", :rv, :out), # 5, 6
    Pool("pulmonary veins", "ncat_blood", :lung, :out; upstream = true), # 7
    Pool("superior vena cava", "ncat_blood", :svc, :out; upstream = true), # 8
    Pool("right brachiocephalic vein", "ncat_blood", :brachiocephalic_r, :out; through = :svc, upstream = true), # 9
    Pool("left brachiocephalic vein", "ncat_blood", :brachiocephalic_l, :out; through = :svc, upstream = true), # 10
    Pool("inferior vena cava", "ncat_blood", :ivc_hepatic, :out; upstream = true), # 11
    Pool("atrial myocardium", "ncat_heart", :coronary, :organ; ml = XCAT_MYOCARDIUM_ML), # 12
    Pool("lung", "ncat_lung", :lung, :mean; scale = 0.1), # 13
    Pool("liver", "ncat_liver", :liver, :organ; ml = 1800), Pool("spleen", "ncat_spleen", :gut, :organ; ml = 1500), # 14, 15
    Pool("thoracic muscle", "ncat_muscle", :muscle_upper, :organ; ml = 0.3 * 28000), # 16
    Pool("thoracic skin", "skin", :skin_upper, :organ; ml = 0.3 * 3300), Pool("thoracic fat", "ncat_body", :fat_upper, :organ; ml = 0.3 * 15000), # 17, 18
    Pool("rib", "ncat_rib", :bone_upper, :organ; ml = 0.3 * 10000), Pool("cortical bone", "corticalbone", :bone_upper, :organ; ml = 0.3 * 10000), # 19, 20
    Pool("spine", "ncat_spine", :bone_upper, :organ; ml = 0.3 * 10000), Pool("stomach", "stomach", :other_upper, :organ; ml = 0.3 * 8000), # 21, 22
    Pool("muscle (other)", "ncat_muscle", :other_upper, :organ; ml = 0.3 * 8000), Pool("soft tissue", "softtissue", :other_upper, :organ; ml = 0.3 * 8000), # 23, 24
    Pool("ventricular myocardium", "ncat_heart", :coronary, :organ; ml = XCAT_MYOCARDIUM_ML), # 25 (when the heart is not the organ)
    Pool("organ tissue mean", "", :coronary, :tissue), # 26: the organ's own mean (its material is the spec's tissue_material)
    Pool("descending aorta", "ncat_blood", :lv, :out; through = :desc_aorta), # 27 (set lead_ml to the subject's arch)
]
# 7 is XCAT's trach_bronch, the airway wall (cartilage rings and mucosa; 3 mm in the trachea); 11 is the air inside it
const XCAT_STATIC = Dict(0 => "air", 2 => "cartilage", 7 => "cartilage", 11 => "air", 30 => "softtissue", 31 => "softtissue", 32 => "softtissue", 70 => "air")

# heart: the myocardium in 64 kinetic classes, grown (radius >= 0.1 mm) and imaged coronary vessels in 64 classes,
# great-artery transit binned by 1 mL of upstream blood, frames every second for 60 s
const HEART_PHANTOM = PhantomSpec(;
    times = collect(0.0f0:1.0f0:60.0f0), tissue_material = "ncat_heart", blood_material = "ncat_blood", tissue_classes = 64, vessel_classes = 64,
    vessel_rmin_mm = 0.1, delay_bin_ml = 1, energy_keV = 70, materials = XCAT_STATIC, pools = XCAT_BODY_POOLS,
)

# lung: lung parenchyma in 64 kinetic classes; the pulmonary trees' grown lumens of radius >= 0.25 mm drawn as blood
const LUNG_PHANTOM = PhantomSpec(HEART_PHANTOM; tissue_material = "ncat_lung", tissue_classes = 64, vessel_classes = 64, vessel_rmin_mm = 0.25)
