# Reference data for validating the coronary model

These files hold published experimental data for checking the synthetic coronary model. The model's
parts are the CCO arterial trees, the lumped capillary beds, the mirrored venous tree, Poiseuille
hemodynamics and contrast transit. Comparisons are made per vessel order or per diameter, as in
Kassab's morphometry plots. Every number comes from a cited, published source. There are no
unpublished data in this folder.

Each `*.toml` file covers one source paper and parses with Julia's `TOML` stdlib:

```julia
using TOML
d = TOML.parsefile(joinpath(@__DIR__, "kassab1993_arterial_morphometry.toml"))
d["data"]["rca_elements"]["diameter_mean_um"]
```

## File layout

- `[meta]` describes the study:
  - `id`: the file stem.
  - `data_type`: `experimental`, `experimental (review summary)`, `model` or `morphometry-derived ...`.
  - `validates`, `species`, `preparation` and `conditions`.
  - `verification`: the overall flag.
  - `access`: how the text was actually read.
  - `notes`, and `use_in_model` (how to compare the model to the data).
- `[source]` gives authors, year, title, journal, volume, issue, pages and DOI. `crossref_verified`
  records the date of the Crossref lookup and whether the title matched. A `[secondary_source]`
  block appears only when numbers come from a later paper that reproduces them.
- `[data.<name>]` tables hold parallel arrays. The unit is in the key suffix (`_um`, `_mm`,
  `_mmHg`, `_ml_per_s`, ...). Each table also has `location` (table, figure or text in the source),
  `provenance`, `statistic` (mean ± SD or SE, and what `n` counts) and `missing_note`.
- A missing number is stored as `nan` and explained in `missing_note`. No number is invented,
  interpolated or back-computed.

## Provenance flags

| Flag | Meaning |
|---|---|
| `FT` | Read from the full text (a table or the text). Tables were checked against 200–400 dpi page renders, because the OCR of these 1980s–90s scans mangles ± signs and zeros. |
| `abs` | Read from the abstract only. The full text was not accessible. |
| `fig` | Digitized from a figure. `digitization_error` gives the approximate error. |
| `secondary` | Copied from a later paper that reproduces the original numbers. Both papers are cited. |
| `missing` | A placeholder: the values exist in the paper but could not be read. |

In the table below, a file-level status of `mixed` means the file combines FT tables with
digitized figures or abstract statements. Each `[data.*]` table carries its own flag.

## Files

### Morphometry (pig): the per-order targets for the trees

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `kassab1993_arterial_morphometry.toml` | RCA, LAD and LCx per diameter-defined Strahler order (1–11): segment and element diameter and length (mean ± SD, n), S/E, connectivity matrices, total element counts, Horton fits | Pig, 5 hearts; silicone cast at 80 mmHg, arrested and vasodilated | FT (Tables 1–9) |
| `kassab1994_venous_morphometry.toml` | Sinusal (orders −1 to −12) and thebesian (−1 to −10) veins: major axis, axis ratio, length, S/E, connectivity matrices including the capillary row, element counts, Horton fits, arcades, anastomoses, blood-volume split | Pig; cast at 80 mmHg arterial inlet and 0 mmHg venous outlet, KCl arrest with adenosine | FT (Tables 1–9, text) |
| `kassab1994_capillary_network.toml` | Capillary segment diameter and length by segment type (C0a/C00/C0v/Ccc) for RV and LV, functional length, path length, cross-connections per path, capillaries per arteriole and venule order | Pig, 4 hearts; cast, arrested, vasodilated | FT (Tables 1–5, text) |
| `kassab1995_arteriolar_bifurcations.toml` | Murray ratio (D1³+D2³)/D0³ at bifurcations of 9–50 µm arterioles | Pig, 4 normal and 3 RVH; vasodilated cast | FT |
| `huo2012_intraspecific_scaling.toml` | Stem-crown scaling exponents fitted to the pig trees and to other organs: volume–diameter, flow–diameter, flow–length; Horton ε, γ | Pig and others; anatomy only | FT |
| `seiler1992_lumen_area_vs_bed_size.toml` | Human epicardial lumen area vs distal branch length and distal mass; area ratio at bifurcations | Human angiography | FT |
| `kassab1999_capillary_hemodynamics.toml` | Capillary distensibility and pressure–radius data used as model input | Pig, epicardial capillaries | secondary (full text not available) |

### Microvascular pressure and resistance distribution

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `chilian1986_microvascular_pressure.toml` | Systolic, diastolic and mean pressure in subepicardial arterioles (~100–325 µm) and venules (~110–225 µm); % of aortic pressure; % of resistance upstream | Cat, beating LV; resting tone vs papaverine | FT (Tables 1, 3–6, text) |
| `chilian1989_dipyridamole_resistance.toml` | Pressure by class (arteries, arterioles, venules, veins); % and absolute resistance in artery (>170 µm), microvessel and venous compartments; individual pressure–diameter points | Cat, beating LV; rest vs dipyridamole | mixed (Tables FT; Fig. 3 points and Fig. 4 venous bars fig) |
| `chilian1991_epi_endo_pressure.toml` | Subepicardial vs subendocardial pressures in ~100 µm arterioles and venules; resistance split | Pig, arrested, maximally dilated | mixed |
| `fujii1992_epi_endo_pressure.toml` | Epi vs endo arteriolar and venular pressure vs perfusion pressure | Dog, arrested, adenosine | mixed |
| `tillmanns1981_microvascular_pressure.toml` | Servo-null pressures in epimyocardial arterioles and venules by size class | Rat and cat, beating heart | FT |
| `nellis1981_microvascular_pressure.toml` | Pressures in small RV arteries and veins above and below 140 µm | Rabbit, beating RV | mixed |
| `mittal2005_coronary_flow_model.toml` | **Model benchmark, not data.** Flow–diameter exponent, inlet flow, arterial conductance, arterial transit times and volume on the Kassab anatomy | Pig anatomy; Poiseuille with Pries viscosity, 100 → 26 mmHg | FT |

### Velocity and flow

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `stepp1999_microvessel_velocity.toml` | Blood velocity vs diameter in 40–430 µm subepicardial arteries and arterioles; wall shear | Dog; rest vs intracoronary adenosine | mixed |
| `kanatsuka1990_microvessel_velocity.toml` | Red-cell velocity vs driving pressure in an arteriole, a venule and capillaries; zero-flow pressure | Dog; maximal dilation, long diastoles | mixed |
| `ashikawa1986_microvessel_velocity.toml` | Phasic red-cell velocity pattern in arterioles, capillaries and venules | Dog | abs |
| `toyota2005_microvessel_velocity.toml` | Systolic retrograde vs forward velocity, endo vs epi arterioles | Dog | abs |
| `kiyooka2005_microvessel_velocity.toml` | Capillary red-cell velocity at control, reactive hyperemia and adenosine | Dog | abs |
| `wieneke2005_flow_velocity.toml` | Proximal APV per artery at rest and with adenosine, IVUS lumen area, CFVR, total coronary flow | Human, normal arteries | FT |
| `anderson2000_velocity_area_mass.toml` | Resting proximal LAD APV vs lumen-area-to-mass ratio; LAD flow | Human | FT |
| `wilson1985_doppler_catheter_cfr.toml` | Hyperemic/rest velocity ratio after contrast and dipyridamole | Human | FT |
| `kern1996_cfvr_normal.toml` | Normal coronary flow velocity reserve by artery and sex | Human, intracoronary adenosine | abs |
| `ofili1993_jacc_flow_velocity.toml` | Proximal peak diastolic velocity and reserve in normal arteries | Human | abs |
| `ofili1993_ajc_proximal_distal.toml` | Proximal/distal velocity ratio in normal arteries | Human | abs |
| `ofili1995_doppler_apv.toml` | Proximal vs distal APV, flow and resistance per artery | Human | missing (the abstract has no numbers; the arrays are placeholders) |
| `corban2021_lad_flow.toml` | Resting and hyperemic mid-LAD flow (ml/min), CFVR, by sex | Human | abs |

### Contrast transit (angiographic)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `gibson1996_timi_frame_count.toml` | TIMI frame count per artery in normal flow (LAD 36.2 ± 2.6, RCA 20.4 ± 3.0, LCx 22.2 ± 4.1 frames at 30 fps), CTFC and the 1.7 LAD correction, landmarks | Human, normal arteries | FT (full text via Web Archive) |
| `abaci1999_tfc_confounders.toml` | Sensitivity of frame count to injection rate, catheter size, nitrate, heart rate and phase | Human | FT |

### Contrast transit (whole-body bolus, CT and MRI)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `leiknes2023_test_bolus_ccta.toml` | Ascending-aorta test-bolus time to peak (18.3 ± 2.6 s) and peak (194 ± 62 HU), main-scan aortic HU, contrast and saline volumes, MRI cardiac output (5.05 ± 1.13 L/min); Table S1 injection protocol (6 mL/s; test bolus 20 mL + 40 mL saline, 12-s delay, 2-s sampling); one digitized test-bolus curve | Human, 40 CCTA patients, iohexol 350, 100 kV; injection arm not stated | FT (Table 1, Table S1; Fig. S1 fig) |
| `borges2022_bolus_arrival_times.toml` | Arrival times (+50 HU) arm → SVC → pulmonary trunk → ascending aorta; PT→AA vs cardiac index | Human, n = 43; right antecubital, 4 mL/s, 1.5 mL/kg of 350 mgI/mL + 30 mL saline | FT |
| `nelsson2021_pulmonary_transit_time.toml` | Pulmonary transit time and pulmonary blood volume (PA→LA and RV→LV methods) | Human, CMR, n = 58 (21 healthy), rest | FT |
| `pelgrim2017_ctp_timing.toml` | Delay from aortic peak to non-ischemic myocardial peak; ischemia contrast window | Human, adenosine dynamic CTP | abs |
| `lee2022_dynamic_ctp.toml` | Dynamic CTP injection protocol and timing; stress myocardial HU (non-ischemic vs ischemic); CT MBF | Human, adenosine, 50 mL iohexol 350 at 5 mL/s + 30 mL saline, 100 kVp | FT |

### Contrast-bolus pathway (injection site to measurement site)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `deng2026_port_vs_peripheral_ct.toml` | Aortic HU and subclavian streak artifact, arm cannula vs chest port (same patients, same dose and rate) | Human, 195 | FT |
| `wang2025_arm_position_cta.toml` | Arterial HU with arms raised vs at the side (thoracic outlet) | Human, 185, 80 kV | FT |
| `nishiyama2024_arm_position_contrast.toml` | Collateral venous diversion and subclavian vein area, conventional vs crossed raised arms | Human, 92 | abs |
| `you2007_injection_arm_side.toml` | Right vs left arm injection: arterial timing | Human, 60 randomized | abs |
| `tseng2007_left_arm_reflux.toml` | Left brachiocephalic reflux vs retrosternal space | Human, 364 | abs |
| `skrok2012_mri_first_pass_bolus_kinetics.toml` | RV→LV pulmonary transit (6.4 s), LV bolus FWHM and rise time after a very short bolus | Human, 18 healthy, MRI | FT |
| `lehmonen2022_cmr_chamber_bolus_passage.toml` | Bolus arrival, peak and width (in cardiac cycles) in SVC, LA, LAA, LV; SVC→LA/LV transit | Human, 30 young controls, CMR | FT |
| `dock1961_pulmonary_blood_volume.toml` | PA→LA mean transit time and pulmonary blood volume by dye dilution | Human, 4 normal + 45 valve patients | FT |
| `bae2010_ct_contrast_timing_review.toml` | Arrival times by segment, injection-site shifts, recirculation, saline chaser effect, HU per mgI/mL at 80/100/120 kVp | Human, review | FT (values secondary) |
| `hazirolan2009_test_bolus_pa_aorta.toml` | PA and aortic test-bolus peak ΔHU and width above 100 HU (20 mL at 5 mL/s, 120 kV), contrast warming | Human, 64 | FT |
| `sakai2010_cardiac_function_aortic_peak.toml` | Descending-aorta peak time and peak HU vs cardiac index | Human, 29 | abs |
| `moradi2016_ef_aortic_test_bolus.toml` | Ascending-aorta test-bolus time to peak and peak HU vs ejection fraction | Human, 120 | FT |
| `tomita2019_kvp_iodine_sensitivity.toml` | Iodine HU ratio vs 120 kVp at 70–140 kVp | Phantom | FT |
| `mcdermott2024_phantom_idr_coronary_iodine.toml` | Peak coronary iodine concentration vs iodine delivery rate (10-s injection, CO 4.8 L/min) | Circulation phantom | FT |
| `salinaro2024_basilic_vein_diameter.toml` | Basilic vein diameter (4.4 ± 1.4 mm) and velocity, mid upper arm | Human, 28, 30° head-up, ultrasound | FT |
| `berk2026_subclavian_vein_sonoanatomy.toml` | Right subclavian vein diameter at three sites (9.5–10.6 mm) | Human, 20, supine, ultrasound | FT |
| `kim2008_ijv_subclavian_csa.toml` | Subclavian and internal jugular CSA, supine | Human, 20 men, ultrasound | abs |
| `kawalecrutkowska2025_ijv_csa_posture.toml` | Right vs left internal jugular CSA vs posture | Human, 25, ultrasound | FT |
| `mahlon2007_svc_ct_length.toml` | SVC length (70.7 ± 14.1 mm) | Human, 112, CT angiography | abs |
| `lin2009_great_vessels_mdct.toml` | SVC CSA percentiles | Human, 103 healthy, gated CT | abs |
| `wang2024_svc_brachiocephalic_cadaver.toml` | Brachiocephalic and SVC circumference and length (cadaver; low weight) | Human cadaver | FT |
| `dymond1982_svc_bolus_fwhm.toml` | SVC bolus FWHM after basilic-vein injection (first-pass radionuclide), right vs left arm RV→LV transit | Human, cardiac patients | FT |
| `claussen1984_bolus_geometry.toml` | Peripheral vs vena-cava catheter bolus geometry (qualitative) | Human | abs |
| `vancauteren2018_injector_dead_space.toml` | Injector-to-port tubing dead space (12 mL); cited venous holding volume | Minipig set-up; cited human value | FT (venous value secondary) |
| `cowan2017_injector_tubing_volume.toml` | Example connecting-tube fluid volume (2.78 mL) and syringe dead space | Patent US9694131B2 (device example) | FT |
| `halliwill1999_limb_venous_compliance.toml` | Forearm venous pressure-volume fit (ΔV = β0 + β1·P + β2·P²), compliance vs pressure; antecubital venous pressure 5.7 ± 1.1 mmHg | Human, healthy, supine, plethysmography | FT |
| `oue2019_forearm_venous_capacitance.toml` | Forearm venous compliance at 20 mmHg by sex, capacitance at 60 mmHg, outflow on deflation | Human, healthy, supine | FT |
| `buffa2015_upper_limb_volume.toml` | Arm tissue volume (wrist to mid upper arm) to scale per-100-mL compliance | Human, 30, water displacement | FT |
| `coskun2025_pvp_rap_pah.toml` | Arm venous pressure vs right atrial pressure | Human, 138 PAH patients, supine | FT |
| `amar2001_pvp_cvp_surgical.toml` | CVP minus peripheral venous pressure, ventilated and spontaneous | Human, surgical patients | abs |
| `munis2001_pvp_cvp_neurosurgical.toml` | Peripheral minus central venous pressure (~3 mmHg) | Human, 15 anesthetized | abs |
| `zhang2025_limb_basal_arterial_flow.toml` | Resting brachial, radial and ulnar artery flow by Doppler; BSA/BMI regression | Human, 304, supine | FT |
| `harvey2014_forearm_blood_flow.toml` | Resting forearm blood flow (mL/100 mL/min) and forearm volume | Human, 25 women, plethysmography | FT |

### Flow regulation, autoregulation and oxygen use

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `smith1993_coronary_autoregulation.toml` | Circumflex pressure–flow autoregulation: lower break point (45 ± 3 mmHg; 61 ± 2 after L-NAME), closed-loop gain (0.84 ± 0.09; 0.78 ± 0.07), slope below the break point, reactive hyperemia | Dog, conscious, paced | FT |
| `gould1974_critical_stenosis.toml` | Resting flow held to 85–90% diameter stenosis; hyperemic flow falls from about 50% | Reproduced figure caption (species not verified) | secondary |
| `kanatsuka1989_adenosine_microvessels.toml` | Size range of arterial microvessel dilation: pacing 40–380 µm (inverse to size) vs adenosine/dipyridamole < 150 µm only | Dog, beating heart, matched perfusion | FT |
| `mccommis2010_myocardial_oxygenation.toml` | Global OEF, MBF, MBV and MVO2 at rest and adenosine hyperemia | Human, healthy, CMR | FT |
| `bassingthwaighte2001_metabolic_heterogeneity.toml` | Flow heterogeneity summary, temporal stability of regional flow, regional O2 use rising with regional flow | Review of baboon, sheep, dog, rabbit and human data | FT (review summary) |

### Coronary flow capacity (PET stress perfusion and CFR categories)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `johnson2012_coronary_flow_capacity.toml` | Original CFC map: nested-rectangle categories on CFR vs stress flow (normal > 3.37 / 2.39; minimal, mild, moderate < 2.03 / 1.12; definite ischemia < 1.74 / 0.91; steal CFR < 1; scar rest < 0.2), and % of LV per category by group, including 241 young-volunteer scans (median 99% normal) | Human, 1,500 Rb-82 dipyridamole PET scans | mixed (text and tables FT; geometry from Fig. 1) |
| `johnson2011_pet_ischemic_thresholds.toml` | Worst-pixel ischemic thresholds (0.91 cc/min/g, CFR 1.74; indeterminate 1.12, 2.03) and steal frequency (CFR < 1 in 40% of definite ischemia) | Human, 1,674 PET studies (mostly Rb-82), dipyridamole | FT |
| `gould2018_cfc_revised_thresholds.toml` | Current Houston CFC cutoffs (normal > 2.9 / 2.17; minimal > 2.38 / 1.82; mild > 1.6 / 1.09; moderate > 1.27 / 0.83; severe ≤ 1.27 / ≤ 0.83; steal CFR < 1; scar ≤ 0.3) | Human, Rb-82 dipyridamole | FT |
| `gould2019_artery_specific_cfc.toml` | Prevalence and size of severe CFC; minimum-quadrant and global stress flow and CFR | Human, 3,774 Rb-82 PET scans | FT |
| `gould2013_anatomic_vs_physiologic.toml` | Literature-pooled rest, stress flow and CFR by population (normal controls 0.82 / 2.86 / 3.55); ischemic cutoffs by tracer; steal with collaterals | Human, 14,962 subjects from 252 papers | FT (Table 2 secondary) |
| `dewinter2022_o15_water_cfc.toml` | O-15 water CFC cutoffs (normal hMBF > 3.01 or CFR > 2.95; severe ≤ 1.5 / ≤ 1.5; steal CFR < 1), normal reference (hMBF 3.66 ± 1.30, CFR 3.69 ± 1.48), category distributions | Human, [15O]H2O PET, adenosine | FT |
| `murai2020_invasive_cfc_stenosis.toml` | Invasive Doppler CFC grades vs % diameter stenosis and FFR | Human, 133 vessels, adenosine | FT |

### Vessel-wall and regulation model parameters (model input, not data)

| File | What it provides | Species / conditions | Status |
|---|---|---|---|
| `fry2013_vessel_wall_parameters.toml` | Passive, active, myogenic, shear and conducted-metabolic parameters as functions of passive diameter D0; reference-state boundary conditions | Hamster cremaster network model | FT (Table 1) |
| `carlson2008_autoregulation_model.toml` | Large- and small-arteriole wall and tone parameters; control state (activation 0.5, M0 8.28); predicted flow-increase factors | Generic arteriole model | FT |
| `zhao2021_regulation_parameters.toml` | Collateral and arteriole wall parameters; ATP release and conducted-response constants; O2 transport and remodeling constants; autoregulation ratios | Rat hindlimb model (author manuscript) | FT |

### Flow heterogeneity

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `bassingthwaighte1989_flow_heterogeneity.toml` | RD at 1 g and fractal D per animal (Table 2) and per species; RD(m) = RD(1 g)·m^(1−D) | Baboon (awake), sheep, rabbit; microspheres | mixed (tables FT, Figs. 1, 6, 7 fig) |
| `king1985_flow_heterogeneity.toml` | RD at ~0.17 g, split into temporal, methodological and spatial parts | Baboon, awake | FT |
| `yipintsoi2016_flow_heterogeneity.toml` | RD at ~0.2 g and Hurst H, at rest, under catecholamines and under adenosine | Dog | FT |
| `bauer2001_flow_heterogeneity.toml` | Fractal D at sub-millimetre resolution vs perfusion pressure and nitroglycerin | Rat, isolated heart | FT |
| `stapleton1995_flow_heterogeneity.toml` | Heterogeneity at near-cellular voxels | Hamster | FT |
| `chareonthaitawee2001_human_flow_heterogeneity.toml` | Human PET flow at rest and hyperemia; regional CV | Human, healthy | FT |

### Blood volume and whole-bed transit time

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `knopp1976_coronary_transit.toml` | Shape of the aortic-root-to-coronary-sinus transport function: σ/t̄, skewness, kurtosis, appearance time/t̄, peak time/t̄, t̄ range 3–7 s; normalized curves | Dog, 60 curves in 3 dogs | mixed (text FT, Fig. 3 fig) |
| `judd1991_myocardial_blood_volume.toml` | Intramyocardial blood volume, split into vessels above and below 100 µm, transmural, diastolic vs systolic arrest | Rat | FT |
| `crystal1981_myocardial_blood_volume.toml` | Change in total vs small-vessel blood volume with adenosine | Dog | abs |
| `wu1992_blood_volume_flow.toml` | Blood volume fraction ∝ F^½ (epi and endo) | Dog | abs |
| `liu1995_blood_volume_flow.toml` | Blood volume = A·F + B·√F | Pig | abs |
| `mccommis2007_myocardial_blood_volume.toml` | MRI blood volume at rest, under dipyridamole and under dobutamine | Dog | abs |
| `nickander2019_cmr_myocardial_blood_volume.toml` | Human MBV and MBF at rest and adenosine stress (MBV/MBF gives the mean transit time) | Human, healthy | FT |
| `yang2024_cmr_myocardial_transit_time.toml` | Human aortic-root-to-coronary-sinus arrival delay on CMR | Human, rest | FT |
| `chatzantonis2020_cmr_myocardial_transit_time.toml` | The same MyoTT measure in controls | Human, rest | abs |
| `canty1991_contrast_extraction.toml` | First-pass extraction of nonionic contrast (ioversol) 33 ± 2% at vasodilated flow vs 50 ± 3% at ~half flow; intravascular-emulsion blood volume and appearance times. Flow only as circumflex mL/min (rest 20, vasodilated 123, stenotic 67), no per-gram flow or PS | Dog, closed chest, carbocromen vasodilation, fast CT | mixed (text and Table 1 FT; stenotic appearance times fig) |

## Caveats that matter for comparisons

- **Kassab orders are diameter-defined Strahler orders** (Eq. 1 of Kassab 1993), not generations.
  Order the CCO tree with the same rule before comparing. Orders 1–3 of the LCX are printed
  identical to the LAD in the source.
- **Venous diameters are major axes** of elliptical cross sections. Divide by the major-to-minor
  ratio stored in the file if you need an equivalent circular diameter.
- **The morphometry is from vasodilated casts** at 80 mmHg perfusion. The beating-heart pressure
  data (Chilian, Tillmanns, Nellis) are mostly from cats, rabbits and rats, with resting tone
  unless stated. A passive Poiseuille tree is closest to the papaverine and dipyridamole columns.
- **Connectivity matrices** use C(m, n) = (number of daughter elements of order m arising from
  parents of order n) / (number of elements of order n). Rows are daughters, columns are parents,
  and row 0 is capillaries. Read the `layout` or `convention` key in each file.
- **Several printed inconsistencies are kept as printed** and listed in each file's `notes`. They
  include Huo & Kassab's two volume–diameter exponents, Mittal's resistance unit, the rabbit N in
  Bassingthwaighte, and the RV order-4 value of 0.40 vs 0.040 in Kassab & Fung.

## Not included, and why

- **Per-order pressure, flow and velocity on the Kassab anatomy.** These are tabulated in Kassab,
  Berkley & Fung 1997, Ann Biomed Eng 25:204 (doi 10.1007/BF02738551), which is paywalled.
  Mittal 2005 shows them only as density plots.
- **Ofili 1995 proximal vs distal APV.** The Am Heart J full text is paywalled and the abstract
  has no numbers.
- **Abstract-only entries** (`abs`) can be upgraded once the full texts are available.
