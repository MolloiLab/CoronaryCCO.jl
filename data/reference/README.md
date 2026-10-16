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

### Human back pressure and hyperemic resistance (full dilation)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `nanto2001_human_zero_flow_pressure.toml` | Zero-flow pressure: rest 47 ± 15, IC ATP normal beats 36 ± 9, IC ATP long diastole 21 ± 7 mmHg (n = 15) | Human, LAD Doppler wire, no significant stenosis | FT |
| `vanherck2007_hyperemic_zero_flow_pressure.toml` | Hyperemic late-diastolic Pzf 14.6 ± 8.0 (stented) and ~12.7 ± 7 mmHg (reference vessel) in stable angina; rise after MI tracks LVEDP | Human, IC adenosine, Doppler wire | mixed (reference bar fig) |
| `wustmann2003_coronary_wedge_pressure.toml` | Wedge pressure after 1-min occlusion of a normal artery 23 ± 8 mmHg vs CVP 5.5, LVEDP 11; CFI 0.19 | Human, 100 normal arteries, no adenosine | FT |
| `verhoeff2005_hmr_reference_vessel.toml` | HMR (h-Pd/h-v) 1.73 ± 0.38 mmHg·s/cm in normal reference vessels; before/after PCI | Human, IC adenosine, dual-sensor wire | FT |
| `fournier2021_absolute_flow_resistance.toml` | Absolute hyperemic flow (whole heart 668 ± 185 mL/min; 4.9–5.9 mL/min/g) and Rmicro × mass 15.5–21.7 mmHg·min·g/mL in normals | Human, continuous thermodilution + CT mass, saline hyperemia | mixed (supplement not read) |
| `mahendiran2025_mass_indexed_microvascular_resistance.toml` | LAD Qhyper 4.56 ± 2.20 mL/min/g and Rμ,hyper × mass 20.7 ± 8.0 mmHg·min·g/mL (controls); resting and CMD values | Human ANOCA, continuous thermodilution + CT mass | FT |
| `melikian2010_imr_normal_controls.toml` | IMR 19 ± 5 U in patients without atherosclerosis | Human, bolus thermodilution | abs |

### Clinical hyperemic stimuli: does adenosine reach full dilation?

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `wilson1990_adenosine_vs_papaverine.toml` | CBFV ratio: IC papaverine 4.5, IC adenosine 16 µg 4.6 ± 0.7, IV adenosine 140 µg/kg/min 4.4 ± 0.9 (papaverine level in 84% of patients) | Human, normal arteries, Doppler catheter, same patients | abs |
| `rossen1991_dipyridamole_adenosine_papaverine.toml` | Velocity ratio papaverine 3.9 ± 1.1, IV adenosine 3.4 ± 1.2, IV dipyridamole 3.1 ± 1.2; resistance index 0.25 / 0.26 / 0.31 | Human, n = 13–15, non-stenotic artery, same patients | abs |
| `nitenberg2007_postocclusion_vs_adenosine.toml` | After PTCA, CFR 2.53 ± 0.58 (IC adenosine 60 µg) vs 3.31 ± 0.67 (30-s occlusion); CRR 2.58 vs 3.46 | Human, n = 10, Doppler wire, low IC dose | abs |
| `debruyne2017_saline_induced_hyperemia.toml` | IC saline 15–20 mL/min gives 111–112% of IC adenosine hyperemia; saline on top of IV adenosine adds nothing. Notes: pig mechanism (Adjedj 2020) and haemolysis/ATP (Gallinoro 2022) | Human, n = 24, Doppler wire, side-hole infusion catheter | abs |
| `everaars2019_thermodilution_vs_o15_pet.toml` | Invasive vs PET absolute flow in the same territories: 179 vs 185 and 167 vs 166 mL/min (bias −0.9 ± 35); saline-only 184 ± 75 vs adenosine 174 ± 72 mL/min | Human, n = 25, CCTA Voronoi mass, IV adenosine | FT |
| `jones1995_adenosine_arteriole_vs_small_artery.toml` | IC adenosine dilates arterioles < 100 µm (+25 ± 3%) far more than small arteries > 100 µm (+7 ± 3%); NO dependence | Dog, beating heart, intravital microscopy, submaximal doses | abs |

### Systolic compression of the dilated bed (beating vs arrested, heart rate, transmural)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `downey1975_arrested_vasodilated_transmural.toml` | Normalized flow by wall quarter: arrested and dilated, inner/outer ~1.07/0.91 (endo ~20% above epi); paced 200/min, 0.81/1.09 (endo ~20% below). Notes: Downey & Kirk 1975 beating vs arrested lines parallel above peak LVP (abs); Spaan 1981 intramyocardial pump p_im 53 ± 7 mmHg | Dog, vagal arrest or RV pacing, reactive hyperemia or IC adenosine, 42K/86Rb | fig |
| `cobb1974_awake_adenosine_transmural.toml` | Epi/endo flow at rest (0.75/0.83) and IV adenosine (4.98/4.49 mL/min/g); endo/epi 1.06–1.16 rest, 0.81–0.99 adenosine (HR 111), 0.63–0.78 after thoracotomy (HR 176) | Dog, awake then anesthetized, 7–10 µm microspheres | FT |
| `bache1977_vasodilated_tachycardia_transmural.toml` | Four-layer flow at maximal dilation vs pacing 100–250/min: mean 5.49 → 4.52 mL/min/g (−18%), endo/epi 1.00 → 0.40, outer half unchanged; diastolic time fraction 65 → 52% | Dog, awake, IV adenosine, microspheres | FT |
| `fokkema2005_diastolic_time_fraction_conductance.toml` | Layer conductance vs diastolic time fraction at full dilation: C_endo = −0.103 + 0.197·DTF + 0.00074·Pc; C_mid = −0.048 + 0.126·DTF + 0.00049·Pc; C_epi DTF-independent. Notes: Duncker 1994 beating-heart Pzf 12.6–23 mmHg | Goat, perfused left main, adenosine, microspheres | abs |
| `mcginn1990_cfr_heart_rate_human.toml` | CFR 4.5 → 3.8 → 3.2 with pacing 76 → 100 → 120/min, with hyperemic velocity unchanged; handgrip and volume effects | Human (mostly transplant), IC papaverine, Doppler catheter | abs |
| `park2000_lad_dsvr_hyperemia.toml` | LAD diastolic/systolic peak-velocity ratio 2.1 at rest → 1.7–1.8 at hyperemia (n = 9). Notes: resting DSVR 1.8 ± 0.5 (Segal 1992), 2.08 ± 0.47 (de Waard 2019) | Human, normal LAD, IC adenosine, Doppler wire | abs |

### Human PET myocardial blood flow (healthy volunteers)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `kaufmann1999_o15_water_reproducibility.toml` | Hyperemic MBF 3.51 ± 0.45 / 3.83 ± 0.49 mL/min/g, CVR ~4, MAP/MBF 25 ± 4 mmHg·min·g/mL at hyperemia | Human, 21 healthy men 45 ± 8 y, [15O]H2O, adenosine | FT |
| `uren1995_o15_water_age_bands.toml` | Rest and hyperemic MBF and reserve by decade (20–29 y: 2.95 ± 1.04; 30–39 y: 3.66 ± 0.89 mL/min/g) | Human, 56 normals, [15O]H2O, mostly dipyridamole | FT |

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
| `dodge1998_injection_rate_tfc.toml` | A 1 mL/s rise in injection rate (10th to 90th percentile of hand injections) shortens the frame count by < 2 frames (< 7%) | Human, power injector | abs |
| `stoel2009_dye_injection_pressure.toml` | Manual injection raises ostial pressure 6.0 ± 4.2 mmHg (6.7%) and distal pressure 3.9 mmHg; the implied TFC shortening is 1–2 frames. Pressure rise at 0, 3, 6 and 9 cm | Human, n = 25, 6F, nitroglycerin | FT |
| `hodgson1985_contrast_flow_first_6s.toml` | Coronary flow in the first 6 s after IC contrast: little rise during injection, nadir about 45% of baseline at 1.9 s | Dog, n = 16, flow probe, 1–6 mL/s | abs (one value secondary) |
| `stoel2003_frame_count_reserve.toml` | Rest vs IV adenosine frame count (LAD 35.1 → 10.3, LCx 25.1 → 8.3 frames at 25 fps), guidewire landmark lengths (LAD 14.3, LCx 11.4 cm), front velocity = 0.63 × Doppler APV, flow about 66–72 mL/min per artery at rest | Human, n = 38, manual injection, nitroglycerin | FT (Web Archive) |
| `stoel2003_frame_count_velocity.toml` | Guidewire landmark path lengths LAD 14.5, LCx 12.8, RCA 11.3 cm and frame count velocity per artery | Human, 119 vessels after PCI | abs |

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
| `tune2004_coronary_venous_po2_regulation.toml` | Coronary venous PO₂ is the regulated quantity (not SO₂ or extraction); resting LV MVO₂ ~60 µL/min/g, extraction ~75%, venous PO₂ ~18 Torr; strong dilation recruited below 15–20 Torr; β-feedforward ~25% of exercise flow rise | Review; canine exercise data | FT (Web Archive) |
| `crapo1999_arterial_blood_gas_reference.toml` | Inlet O₂: sea-level SaO₂ 96.9 ± 0.4% (18–24 y) to 95.5 ± 1.4% (> 64 y), PaO₂ by age; arterial Hb men 14.0 ± 1.0, women 12.3 ± 0.9 g/dL; 1,400 m values | Human, 339 healthy nonsmokers, CO-oximetry | FT (Web Archive) |
| `beutler2006_hemoglobin_reference.toml` | Hb by sex, age decade and race (white men 20–29 y 15.43 ± 0.77 / 15.51 ± 0.96 g/dL); lower limit of normal 13.7 g/dL for white men 20–59 y | Human, Scripps-Kaiser and NHANES-III, iron deficiency and inflammation excluded | FT |

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
| `arciero2008_atp_conducted_response.toml` | RBC ATP release R(S) = R0(1 − R1·S): R0 1.4 × 10⁻⁹ mol/s/cm³, R1 0.891; kd 2 × 10⁻⁴ cm/s; conducted length constant L0 = 1 cm; C_meta 30 µM⁻¹cm⁻¹; S(0) 0.97, P50 26 mmHg, nH 2.7 | Generic seven-segment pathway model | FT (Table 1) |
| `vinnakota2004_myocardial_density.toml` | Myocardial tissue density 1.053 g/mL (rat synthesis), ~1.055 g/mL stated for human; dog 1.063 ± 0.002; plasma 1.023, RBC 1.097 g/mL | Rat-based composition model | FT |
| `kesmarky2008_plasma_viscosity.toml` | Normal plasma viscosity 1.10–1.30 mPa·s at 37 °C, independent of age and sex | Human, review | abs |

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
| `ritman1994_iohexol_capillary_permeability.toml` | Iohexol PS during coronary adenosine (fast CT vs catheter reference), no concentration effect, bradykinin dose response | Pig | abs |
| `daghini2007_ct64_myocardial_permeability.toml` | Iopamidol Patlak permeability-surface product (0.37 mL/min/g) and fractional vascular volume at rest, 64-slice CT vs EBCT | Pig, 12 | abs |
| `malyar2011_ebct_contrast_extraction_microembolization.toml` | EBCT flow, blood volume and an extraction index (arbitrary units) at rest and intracoronary adenosine | Pig, 3 per group | FT |
| `bindschadler2014_iodine_exchange_model_ps.toml` | **Model assumptions, not data.** PS-vs-flow relations used for iodinated contrast (PS = 0.3853·F + 0.3; Alessio 2019: 0.35·F + 0.22), model volumes and dispersions | Pig-fitted / human | FT |
| `lassen2023_rb82_adenosine_cardiac_output.toml` | Cardiac output, HR and SV at rest and during adenosine 140 µg/kg/min; pulmonary transit and volume | Human, 33 healthy, Rb-82 PET | FT |
| `thomas2017_adenosine_cmr_hemodynamics.toml` | HR, blood pressure, SV index and EF after adenosine 140 µg/kg/min; CO +28% (abstract of the same group) | Human, 25 healthy, CMR | mixed |

### Lung (pulmonary circulation)

| File | What it validates | Species / conditions | Status |
|---|---|---|---|
| `huang1996_human_pulmonary_morphometry.toml` | Per-order element and segment diameter, length, S/E, connectivity, corrected counts and volume, arteries and veins, diameter-defined Strahler orders 1–15 (order-1 arteriole 20 ± 3 µm; 5.1 × 10⁷ per left lung). Cited as 80:1773 in the request; correct is 81:2123 | Human, 1 left lung (arteries) + 1 right lung (veins), silicone casts | FT |
| `singhal1973_human_pulmonary_arterial_tree.toml` | Whole-lung arterial tree by classical Strahler order 1–17: counts, diameter, length, area, volume (150.9 mL), velocity at 80 mL/s; order-1 13 µm, 3 × 10⁸ (extrapolated) | Human, 1 resin cast | FT |
| `bentley2020_healthy_exercise_pulmonary_hemodynamics.toml` | mPAP, PAWP, RAP, CO, TPG, Rp, Cp, RC time at semiupright rest and submaximal exercise | Human, 36 healthy, RHC | FT |
| `naeije2012_pulmonary_circulation_exercise.toml` | Limits of normal at rest (Ppa 13, Ppao 9, Pcp 10 mmHg, PVR 55 dyn·s·cm⁻⁵); pressure–flow slopes young vs old; 60/40 arterial vs capillary+venous resistance (Gaar) | Human, review | FT (values secondary) |
| `kovacs2009_normal_pulmonary_pressure.toml` | Normal mPAP at rest (14.0 ± 3.3 mmHg) by posture and age, and at mild exercise by age | Human, 1,187 healthy, review | abs |
| `kovacs2012_normal_pvr_exercise.toml` | Supine resting PVR and TPR by age; −12% PVR and −25% TPR with an 85% CO rise | Human, review | abs |
| `bhattacharya1980_lung_microvascular_pressures.toml` | Micropuncture resistance split: 45% in alveolar capillaries, most of the rest in arterioles, none in venules > 20 µm | Dog, isolated perfused lung | abs |
| `hakim1982_pulmonary_resistance_partition.toml` | Occlusion resistance split: middle segment < 16%, arteries and veins equal shares of the rest | Dog, perfused lobe, zone 3 | abs |
| `bouwsema2017_capillary_blood_volume_exercise.toml` | Resting capillary blood volume Vc (89 / 76 mL), DLCO, Dm; CO rest to 90% VO₂peak; capillary transit Vc/Q stays > 0.3 s | Human, 32 healthy young, multiple-FIO₂ DLCO | FT |
| `ebrahimi2022_pulmonary_model_parameters.toml` | **Model input, not data.** Sheet-flow capillary parameters (H₀ 3.5 µm, compliance 1.3 × 10⁻⁹ m/Pa, from cat), viscosities, 30,676 acinar units, tuned Strahler ratios | Human lung model | FT |
| `glenny1990_pulmonary_flow_fractal.toml` | Regional pulmonary flow heterogeneity is fractal, Ds = 1.09 ± 0.02; gravity a minor factor | Dog, 6, MAA | abs |
| `molina2012_lung_weight.toml` | Lung weight at autopsy (right 445 g, left 395 g, with blood) | Human, 232 men 18–35 y | abs |
| `ochs2004_alveolar_number.toml` | 480 × 10⁶ alveoli (274–790 × 10⁶), ~170 per mm³ parenchyma | Human, 6 lungs, stereology | abs |
| `haefelibleuer1988_human_acinus.toml` | Acinar volume 187 ± 79 mm³, ~9 intra-acinar generations, path length 8.8 mm | Human, 2 lungs, casts | abs |
| `huang2011_human_sheet_flow_parameters.toml` | **Model input, not data.** Human sheet-flow constants (h₀ 3.5 µm, α 0.127 µm/cmH₂O, S 0.88, A 126 m², µ 1.92 cP, k 12, f 1.8, L̄ 1188 µm scaled from cat); per-order human arteries/veins; model vs measured mPAP in 3 subjects. DOI registered with ISTIC, not Crossref | Human lung model | FT |
| `clark2010_acinar_sheet_flow_parameters.toml` | **Model input, not data.** Acinar sheet flow: α_C 1.30 × 10⁻⁹ m/Pa, H′ 3.5 µm, f 21.6 (= k·f), l_C 1186 µm, S 0.86, µ_C 1.92 cP, P_CU 15 cmH₂O; parameter sensitivities | Human lung model | FT |
| `clark2011_capillary_compliance_lung_volume.toml` | Sheet compliance vs transpulmonary pressure (dog 0.122 → 0.079 µm/cmH₂O, Glazier 1969; human scaled fit); vessel α 1.49 × 10⁻⁴ /Pa (Krenz & Dawson); sheet area 73 m² at TLC | Dog data, human model | FT (Glazier secondary) |
| `reeves2005_human_pulmonary_distensibility.toml` | Pulmonary vascular distensibility α 0.02 ± 0.002 /mmHg (n = 59 normal), ≈ 2 %/mmHg in isolated arteries; changes with hypoxia | Human, review of RHC data | abs |
| `yen1988_human_arteriole_venule_elasticity.toml` | Compliance of human arterioles and venules orders 1–3 (0.13–0.25 µm/cmH₂O) | Human, postmortem lungs | abs |
| `goresky1975_extravascular_lung_water.toml` | Accessible extravascular lung water 2.16 g/kg at rest, 2.55 g/kg at low exercise; central blood volume +50% at 3× cardiac output | Human, 11 normal adults, triple indicator dilution | abs |
| `brigham1983_lung_urea_ps.toml` | Closest analogue for lung contrast exchange: [¹⁴C]urea PS 0.50 ± 0.16 (survivors) and 3.44 ± 0.36 mL/s per L TLC (deaths); states no normal-human lung PSu exists. No lung PS for iodinated contrast found | Human, 14 with acute respiratory failure, multiple-indicator dilution | FT |

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
| `vanbavel1992_coronary_branching_symmetry.toml` | Bifurcation symmetry d_s/d_L by mother diameter (quartiles), junction exponents 2.82 / 2.50 / 2.35, area expansion; Strahler ratios of trees built from these local statistics | Pig, 1,663 nodes, casts | FT |
| `kalsho2004_bifurcation_asymmetry.toml` | Diameter asymmetry D_large/D_small and length asymmetry per Strahler order, RCA / LAD / LCx | Pig casts | FT |
| `kaimovitz2008_diameter_asymmetry.toml` | Smaller daughter ÷ mother diameter per (daughter order, mother order), LAD / LCx / RCA | Pig | FT (one table missing) |
| `schwarz2020_human_coronary_symmetry.toml` | Bifurcation symmetry d_S/d_L per mother-diameter class (whiskers, quartiles, median), histograms at 75–90 µm and ≥ 400 µm. Lower whisker per class is the growth symmetry floor | Human left coronary cast, ≈ 88,000 bifurcations | fig |
