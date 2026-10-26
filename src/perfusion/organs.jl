# SPDX-License-Identifier: GPL-3.0-or-later
# Organs as instances of the one `FlowSpec` struct: each organ is only different values. test/reference/ holds the
# measurements they come from and that the results are tested against. Each value's comment gives its unit, how it was
# obtained and its sources [test/reference ids]. A variant: `FlowSpec(LUNG_FLOW; p_in = 20)`.

# Human (and dog/cat where human data are missing) coronary circulation and myocardium.
const HEART_FLOW = FlowSpec(;
    # mmHg. Mean aortic pressure. Normal human mean arterial pressure; the PET subjects of the whole-bed resistance
    # anchor (r_total_dilated) had 86-89 mmHg (Kaufmann 1999). Chilian 1989 Table 2 gives 90 mmHg in the cat
    # (dipyridamole ~91). [kaufmann1999_o15_water_reproducibility, chilian1989_dipyridamole_resistance]
    p_in = 90.0f0,
    # mmHg. Effective outflow pressure: zero-flow pressure of the maximally dilated beating heart. Lumps
    # coronary-sinus pressure and the intramyocardial waterfall into one steady back pressure. Human values agree
    # (Nanto 2001, van Herck 2007). [kanatsuka1990_microvessel_velocity, nanto2001_human_zero_flow_pressure,
    # vanherck2007_hyperemic_zero_flow_pressure]
    p_out = 13.6f0,
    # Coronary flow is matched to O2 demand; coronary venous O2 tension is the regulated quantity (Tune, Gorman &
    # Feigl 2004). [tune2004_coronary_venous_po2_regulation]
    regulation = :o2,
    # Capillary bed resistance from parallel functional capillary tubes (Kassab & Fung 1994 dimensions).
    # [kassab1994_capillary_network]
    capillary = :tube,
    # Grown coronary veins take Kassab's venous/arterial cast volume ratio (the venous tree has no imaged anatomical
    # scale of its own). [kassab1994_venous_morphometry]
    venous_scale = :volume_share,
    # 1. Discharge (bulk) hematocrit, as used by the iodine-exchange model of Bindschadler 2014 (Alessio 2019: 0.45).
    # Tube hematocrit from the Fahraeus law (Pries 1992). Cite: ['Pries, Neuhaus & Gaehtgens 1992 Am J Physiol
    # 263:H1770'] [bindschadler2014_iodine_exchange_model_ps]
    hd = 0.45f0,
    # mL/mL. Intramyocardial blood volume fraction, healthy human at rest by CMR: 9.0 +/- 0.7 % (Nickander 2019).
    # McCommis 2007 (CMR, 6.2 mL/100 g) and Kassab's dilated cast (12.2 mL/100 g LV) bracket it.
    # [nickander2019_cmr_myocardial_blood_volume, mccommis2007_myocardial_blood_volume, kassab1994_venous_morphometry]
    fb = 0.09f0,
    # The CMR blood volume is intramyocardial: the epicardial (imaged) arteries and veins hold their own lumen volume on top.
    # [nickander2019_cmr_myocardial_blood_volume]
    fb_scope = :intramural,
    # 1. Derived: arterial 3.0 of 10.4 mL coronary cast volume (Kassab 1994, Discussion).
    # [kassab1994_venous_morphometry]
    share_art = 0.288f0,
    # 1. Derived: capillary 3.3 of 10.4 mL (Kassab 1994). [kassab1994_venous_morphometry]
    share_cap = 0.317f0,
    # 1. Derived: venous 4.1 of 10.4 mL (Kassab 1994); 0.394 rounded so the three shares sum to 1.
    # [kassab1994_venous_morphometry]
    share_ven = 0.395f0,
    # mm. Functional capillary length. [kassab1994_capillary_network]
    cap_length_mm = 0.5f0,
    # um. Functional capillary diameter. [kassab1994_capillary_network]
    cap_diameter_um = 6.0f0,
    # 1. Derived: tissue tube / discharge hematocrit 0.3825 / 0.45 (Bindschadler 2014).
    # [bindschadler2014_iodine_exchange_model_ps]
    cap_ht_ratio = 0.85f0,
    # mL/mL. Derived: interstitial fraction = ECV - plasma volume = 0.273 - 0.09 x (1 - 0.45) (Nickander 2019 ECV 27.3
    # % and MBV 9.0 %, same healthy cohort). [nickander2019_cmr_myocardial_blood_volume]
    ve = 0.224f0,
    # mL/min/g. Capillary permeability-surface product for iodinated contrast: Patlak influx constant of iopamidol in
    # normal pig myocardium, by 64-slice CT and EBCT, similar during adenosine.
    # [daghini2007_ct64_myocardial_permeability]
    ps = 0.37f0,
    # 1. Imaged (epicardial) coronaries are 13% wider in diameter at maximal hyperemia than at rest: adenosine, proximal
    # LAD +13.1 +/- 1.4% and mid LAD +15.8 +/- 2.5% (n = 15); +9.2 +/- 1.5% proximal (n = 12); distal +15% (n = 24);
    # IVUS area +17% (n = 19). Diseased segments do not dilate (papaverine: atherosclerotic area -7.4%).
    # [kang1998_flow_mediated_epicardial_dilation, shiode1996_adenosine_flow_mediated_dilation,
    # lupi1997_adenosine_dipyridamole_epicardial_diameter, reis1999_adenosine_epicardial_area, nabel1990_papaverine_flow_dilation]
    epicardial_dilation = 0.13f0,
    # g/mL. Myocardial tissue density (Vinnakota & Bassingthwaighte 2004, Table 3: 1.053; text: about 1.055 for
    # mammals including human). [vinnakota2004_myocardial_density]
    rho = 1.053f0,
    # mPa s. Derived: midpoint of the normal human plasma viscosity range at 37 C, 1.10-1.30 mPa s (Kesmarky 2008).
    # [kesmarky2008_plasma_viscosity]
    plasma_viscosity_mPas = 1.2f0,
    # 1. Derived: median SD/mean of element diameter over the 29 RCA, LAD and LCx orders with an SD (Kassab 1993
    # Tables): 0.124 (IQR 0.11-0.19). [kassab1993_arterial_morphometry]
    radius_cv = 0.12f0,
    # um. Arterial microcirculation cut-off of Chilian's resistance classes. [chilian1989_dipyridamole_resistance]
    d_art_micro = 170.0f0,
    # um. Venous microcirculation cut-off of Chilian's resistance classes. [chilian1989_dipyridamole_resistance]
    d_ven_micro = 150.0f0,
    # mmHg min g/mL. Not used: the scale kappa is anchored to the human whole-bed resistance (r_total_dilated).
    # Chilian 1989's dilated cat microvascular resistance (4.0) is now an emergent check.
    # [chilian1989_dipyridamole_resistance]
    r_micro_dilated = 0.0f0,
    # mmHg min g/mL. Whole-bed coronary resistance at maximal dilation in healthy men, MAP / MBF with no back pressure
    # (Kaufmann 1999 PET: 25 +/- 4 and 23 +/- 4 on two adenosine scans; mean 24). Anchors the scale kappa of the grown
    # microvessels (accepted calibration; Dale's decision A2, 2026-10-16); the stress-MBF row is therefore calibrated.
    # [kaufmann1999_o15_water_reproducibility]
    r_total_dilated = 24.0f0,
    # um. Arteries below this passive diameter carry smooth-muscle tone; adenosine and dipyridamole dilate only these.
    # [kanatsuka1989_adenosine_microvessels]
    d_active = 150.0f0,
    # dyn/cm per um. Passive tension strength C_pass = slope x D0. Carlson & Secomb 2005 / Carlson 2008 fits as
    # tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary refits:
    # namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cpass_slope = 6.666f0,
    # 1. Passive tension sensitivity C'_pass = intercept + slope x D0: intercept. Carlson & Secomb 2005 / Carlson 2008
    # fits as tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary refits:
    # namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cpass1_0 = 12.52f0,
    # 1/um. C'_pass slope. Carlson & Secomb 2005 / Carlson 2008 fits as tabulated by Fry et al. 2013 Table 1 (hamster
    # cremaster origin; coronary refits: namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cpass1_slope = -0.027f0,
    # dyn/cm. Maximal active tension C_act = coef x D0^exp: coefficient (refit by Fry 2013 on small vessels). Carlson
    # & Secomb 2005 / Carlson 2008 fits as tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary
    # refits: namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact_coef = 1.3f0,
    # 1. C_act exponent. Carlson & Secomb 2005 / Carlson 2008 fits as tabulated by Fry et al. 2013 Table 1 (hamster
    # cremaster origin; coronary refits: namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact_exp = 1.48f0,
    # 1. Active length dependence C'_act = intercept + slope x D0: intercept. Carlson & Secomb 2005 / Carlson 2008
    # fits as tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary refits:
    # namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact1_0 = 0.773f0,
    # 1/um. C'_act slope. Carlson & Secomb 2005 / Carlson 2008 fits as tabulated by Fry et al. 2013 Table 1 (hamster
    # cremaster origin; coronary refits: namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact1_slope = -0.00059f0,
    # 1. Active tension range C''_act = intercept + slope x D0: intercept. Carlson & Secomb 2005 / Carlson 2008 fits
    # as tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary refits:
    # namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact2_0 = 0.415f0,
    # 1/um. C''_act slope. Carlson & Secomb 2005 / Carlson 2008 fits as tabulated by Fry et al. 2013 Table 1 (hamster
    # cremaster origin; coronary refits: namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cact2_slope = -0.0008f0,
    # cm/dyn x um. Myogenic sensitivity C_myo = value / D0 (cm/dyn). Carlson & Secomb 2005 / Carlson 2008 fits as
    # tabulated by Fry et al. 2013 Table 1 (hamster cremaster origin; coronary refits:
    # namani2018_coronary_wall_parameters). [fry2013_vessel_wall_parameters]
    wall_cmyo = 1.369f0,
    # 1. Resting coronary-sinus O2 extraction: sets the single tone constant at rest (accepted calibration).
    # [mccommis2010_myocardial_oxygenation]
    oef = 0.73f0,
    # 1. Metabolic sensitivity of the control state: C_meta 30 /uM/cm (Carlson 2008 Table) x S_CR about 0.5 uM cm
    # (Carlson 2008 text; not in the TOML). Cite: ['Carlson, Arciero & Secomb 2008 Am J Physiol 295:H1572 (S_CR)']
    # [carlson2008_autoregulation_model]
    g_meta = 15.0f0,
    # mL O2/min/g. Derived: myocardial O2 use 5.2 umol/g/min in healthy volunteers (McCommis 2010) x 22.4 uL/umol.
    # [mccommis2010_myocardial_oxygenation]
    mvo2 = 0.116f0,
    # 1. Arterial O2 saturation by CO-oximetry at sea level, adults 25-44 y: 96.7 % (Crapo 1999 Table 4).
    # [crapo1999_arterial_blood_gas_reference]
    sao2 = 0.967f0,
    # g/dL. Arterial hemoglobin, healthy men at sea level: 14.0 +/- 1.0 g/dL (Crapo 1999, the same arterial samples as
    # sao2). Beutler 2006 venous reference: 15.4 (white men 20-29 y). [crapo1999_arterial_blood_gas_reference,
    # beutler2006_hemoglobin_reference]
    hb = 14.0f0,
    # 1. Slope of red-cell ATP release with saturation, R = R0(1 - R1 S), as used by Zhao 2021. Arciero 2008 printed
    # R1 = 0.891 (fitted to two release rates); the two papers disagree and the regulation is anchored by oef either
    # way. [zhao2021_regulation_parameters, arciero2008_atp_conducted_response]
    atp_r1 = 1.026f0,
    # cm. Length constant of the conducted metabolic signal, about 1 cm (Arciero 2008; Zhao 2021).
    # [arciero2008_atp_conducted_response, zhao2021_regulation_parameters]
    atp_l0 = 1.0f0,
    # mmHg. Not used (capillary = tube): pressure outside the capillary sheet.
    p_ext = 0.0f0,
    # um. Not used (capillary = tube).
    sheet_h0_um = 0.0f0,
    # um/cmH2O. Not used (capillary = tube).
    sheet_alpha_um_cmH2O = 0.0f0,
    # 1. Not used (capillary = tube).
    sheet_s = 0.0f0,
    # m2/mL. Not used (capillary = tube).
    sheet_area_m2_per_ml = 0.0f0,
    # cP. Not used (capillary = tube).
    sheet_mu_cP = 0.0f0,
    # 1. Not used (capillary = tube).
    sheet_kf = 0.0f0,
    # um. Not used (capillary = tube).
    sheet_path_um = 0.0f0,
    # cmH2O. Not used (capillary = tube).
    sheet_pmax_cmH2O = 0.0f0,
    # 1/cmH2O. Not used (rigid vessels in the heart model).
    beta_art_small = 0.0f0,
    # 1/cmH2O. Not used (rigid vessels in the heart model).
    beta_art_large = 0.0f0,
    # 1/cmH2O. Not used (rigid vessels in the heart model).
    beta_ven_small = 0.0f0,
    # 1/cmH2O. Not used (rigid vessels in the heart model).
    beta_ven_large = 0.0f0,
    # um. Not used (rigid vessels in the heart model).
    d_alveolus_um = 0.0f0,
    # cmH2O. Not used (rigid vessels in the heart model).
    p_pleural_cmH2O = 0.0f0,
    # cmH2O. Not used (rigid vessels in the heart model).
    p_ref_cmH2O = 0.0f0,
    # Not used (the heart's vessels below its terminals are carried by the capillary-bed resistance and the
    # microvascular scale kappa).
    sub_art = nothing,
    sub_ven = nothing,
)

# Human pulmonary circulation at rest, normoxia, end-expiration (zone 3 throughout: gravity not modelled).
const LUNG_FLOW = FlowSpec(;
    # mmHg. Mean pulmonary arterial pressure at rest, semiupright (17 +/- 3). Kovacs 2009 pooled supine: 14.0 +/- 3.3.
    # [bentley2020_healthy_exercise_pulmonary_hemodynamics, kovacs2009_normal_pulmonary_pressure]
    p_in = 17.0f0,
    # mmHg. Pulmonary arterial wedge pressure (left-atrial proxy) at rest, same cohort as p_in (11 +/- 2).
    # [bentley2020_healthy_exercise_pulmonary_hemodynamics]
    p_out = 11.0f0,
    # Normoxic pulmonary vessels carry little tone (Naeije 2012; distensibility ~2 %/mmHg, Reeves 2005): no active tone;
    # vessels distend passively (beta_* below) and the capillary sheet with its own law. Hypoxic vasoconstriction is not
    # modelled. [naeije2012_pulmonary_circulation_exercise, reeves2005_human_pulmonary_distensibility]
    regulation = :passive,
    # Alveolar-wall capillaries as Fung-Sobin sheet flow. [huang2011_human_sheet_flow_parameters,
    # clark2010_acinar_sheet_flow_parameters]
    capillary = :sheet,
    # Grown pulmonary veins keep their own radii: their root is anchored at two Huang 1996 order-15 veins per lung
    # (lung/pulm_veins.jl), not at a volume share. [huang1996_human_pulmonary_morphometry]
    venous_scale = :grown,
    # 1. Same blood as the heart (Bindschadler 2014 bulk hematocrit). [bindschadler2014_iodine_exchange_model_ps]
    hd = 0.45f0,
    # mL/mL. Derived: arteries 319 mL (Huang 1996 left lung 150.1 mL / 0.47 left-lung share of weight 395/840 g,
    # Molina 2012) + veins 163 mL (86.6 mL right lung / 0.53) + capillaries 83 mL (Bouwsema 2017, mean of 89 and 76) =
    # 565 mL over 2810 mL. Dock 1961 PBV ~470 mL (246 mL/m2) is the same order.
    # [huang1996_human_pulmonary_morphometry, molina2012_lung_weight, bouwsema2017_capillary_blood_volume_exercise,
    # dock1961_pulmonary_blood_volume]
    fb = 0.201f0,
    # Huang's volumes count every vessel, the hilar (imaged) ones included. [huang1996_human_pulmonary_morphometry]
    fb_scope = :whole,
    # 1. Derived: 319 / 565 (see fb). [huang1996_human_pulmonary_morphometry, molina2012_lung_weight]
    share_art = 0.565f0,
    # 1. Derived: 83 / 565 (see fb). [bouwsema2017_capillary_blood_volume_exercise]
    share_cap = 0.146f0,
    # 1. Derived: 163 / 565 (see fb). [huang1996_human_pulmonary_morphometry, molina2012_lung_weight]
    share_ven = 0.289f0,
    # mm. Not used (capillary = sheet).
    cap_length_mm = 0.0f0,
    # um. Not used (capillary = sheet).
    cap_diameter_um = 0.0f0,
    # 1. Tube/discharge hematocrit in unresolved microvessels (Fahraeus effect; coronary value, no lung measurement in
    # the repo). [bindschadler2014_iodine_exchange_model_ps]
    cap_ht_ratio = 0.85f0,
    # mL/mL. Derived: accessible extravascular lung water 2.16 g/kg in normal humans at rest (Goresky 1975) x 70 kg =
    # 151 mL over XCAT's 2810 mL lung. An extracellular solute reaches at most this water, so it is an upper bound
    # (Levitt 2003: 86 mL; pig extravascular ECV 21% of non-air tissue ~240 mL, Seemann 2025). The body model (body.jl
    # e_lung) uses the same 151 mL. [goresky1975_extravascular_lung_water, levitt2003_interstitial_pharmacokinetics,
    # seemann2025_porcine_lung_extravascular_ecv]
    ve = 0.054f0,
    # mL/min/g. Derived (2026-10-16): human lung DCE-MRI with gadobutrol (17 healthy volunteers, extended Tofts)
    # Ktrans 0.234 /min per mL of non-air tissue (Cheriyan 2022); about 1146 mL of non-air tissue around 510 mL of
    # lung blood (vp 0.445), and Ktrans ~ PS at the lung's high flow, gives 4.47 mL/s for gadobutrol (MW 605); scaled
    # to iohexol (MW 821) by PS ~ MW^-0.51 (human forearm EDTA vs inulin, Keiding 1988): 3.83 mL/s = 230 mL/min over
    # the 840 g lung. No direct measurement for iodinated contrast exists (Brigham 1983 urea PS is in acute
    # respiratory failure). The body model (body.jl ps_lung) uses the same value.
    # [cheriyan2022_lung_dce_mri_ve_ktrans, keiding1988_human_forearm_inulin_edta_ps, brigham1983_lung_urea_ps]
    ps = 0.274f0,
    # 1. The imaged pulmonary arteries keep their size (no hyperemic state is modelled for the lung).
    epicardial_dilation = 0.0f0,
    # g/mL. Derived: lung weight 840 g (Molina 2012, right 445 + left 395, with blood) over XCAT's lung volume 2810 mL
    # (air included). [molina2012_lung_weight]
    rho = 0.299f0,
    # mPa s. Same blood as the heart: midpoint of 1.10-1.30 mPa s (Kesmarky 2008). [kesmarky2008_plasma_viscosity]
    plasma_viscosity_mPas = 1.2f0,
    # 1. Diameter variability within one order: order-1 arterioles 20 +/- 3 um.
    # [huang1996_human_pulmonary_morphometry]
    radius_cv = 0.15f0,
    # um. Not used: no microvascular resistance anchor (r_micro_dilated = 0), so pure Poiseuille plus sheet flow.
    d_art_micro = 0.0f0,
    # um. Not used (see d_art_micro).
    d_ven_micro = 0.0f0,
    # mmHg min g/mL. Not used: 0 = no kappa anchor, so the lung's resistance is a prediction, not a fit.
    r_micro_dilated = 0.0f0,
    # mmHg min g/mL. Not used: 0 = no kappa anchor, so the lung's resistance is a prediction, not a fit.
    r_total_dilated = 0.0f0,
    # um. Not used: no vessel carries active tone (regulation = passive).
    d_active = 0.0f0,
    # dyn/cm per um. Not used: no vessel carries active tone (regulation = passive).
    wall_cpass_slope = 0.0f0,
    # 1. Not used: no vessel carries active tone (regulation = passive).
    wall_cpass1_0 = 0.0f0,
    # 1/um. Not used: no vessel carries active tone (regulation = passive).
    wall_cpass1_slope = 0.0f0,
    # dyn/cm. Not used: no vessel carries active tone (regulation = passive).
    wall_cact_coef = 0.0f0,
    # 1. Not used: no vessel carries active tone (regulation = passive).
    wall_cact_exp = 0.0f0,
    # 1. Not used: no vessel carries active tone (regulation = passive).
    wall_cact1_0 = 0.0f0,
    # 1/um. Not used: no vessel carries active tone (regulation = passive).
    wall_cact1_slope = 0.0f0,
    # 1. Not used: no vessel carries active tone (regulation = passive).
    wall_cact2_0 = 0.0f0,
    # 1/um. Not used: no vessel carries active tone (regulation = passive).
    wall_cact2_slope = 0.0f0,
    # cm/dyn x um. Not used: no vessel carries active tone (regulation = passive).
    wall_cmyo = 0.0f0,
    # 1. Not used (regulation = passive).
    oef = 0.0f0,
    # 1. Not used (regulation = passive).
    g_meta = 0.0f0,
    # mL O2/min/g. Not used (regulation = passive).
    mvo2 = 0.0f0,
    # 1. Not used (regulation = passive).
    sao2 = 0.0f0,
    # g/dL. Same blood as the heart (Crapo 1999). [crapo1999_arterial_blood_gas_reference]
    hb = 14.0f0,
    # 1. Not used (regulation = passive).
    atp_r1 = 0.0f0,
    # cm. Not used (regulation = passive); kept nonzero because the solver's signal decay divides by it.
    atp_l0 = 1.0f0,
    # mmHg. Definition: alveolar pressure at end-expiration is atmospheric, and all pressures are gauge.
    p_ext = 0.0f0,
    # um. Sheet thickness at zero transmural pressure, human (Sobin et al. 1979 as used by Huang 2011 and Clark 2010).
    # [huang2011_human_sheet_flow_parameters, clark2010_acinar_sheet_flow_parameters]
    sheet_h0_um = 3.5f0,
    # um/cmH2O. Sheet compliance, human (Clark 2010: 1.30e-9 m/Pa, the same value).
    # [huang2011_human_sheet_flow_parameters, clark2010_acinar_sheet_flow_parameters]
    sheet_alpha_um_cmH2O = 0.127f0,
    # 1. Vascular-space-tissue ratio: the vascular fraction of the sheet's face-on area (Sobin 1970), human value via
    # Huang 2011, paired with his A = 126 m2 (below). [huang2011_human_sheet_flow_parameters, sobin1970_vstr_definition]
    sheet_s = 0.88f0,
    # m2/mL. Huang 2011's sheet area, 126 m2 (Gehr 1978's capillary surface) over XCAT's lung volume 2810 mL; the lung flow
    # row's pass rests on it (73 m2: 3.09 L/min, below the 2-SD bound; 126 m2: 3.51). The solver
    # uses only the sheet conductance S A / (4 mu k f Lbar^2 alpha), and Huang's whole set (S 0.88, A 126 m2, h0, alpha,
    # k f, Lbar, all in this instance) reproduces his validated model's flows in three normal subjects within 1% (Table 4,
    # Harris 1968 data; scratchpad replication, QUEUE K3 revisited 2026-10-18), while the one-face area (73 m2, Clark
    # 2010/2011) with the same S gives Clark's PVR of 2.65-2.9 WU against ~1 WU measured. By Fung & Sobin's definition A
    # is one face, so one factor of Huang's set must offset ~1.7x; the likeliest is Lbar, scaled from the cat, not
    # measured in humans. Caveat: S A h0 = 0.88 x 126 m2 x 3.5 um = 388 mL is twice Weibel's morphometric capillary volume
    # (194 mL); here the sheet sets only conductance (volumes come from fb x share_*). [huang2011_human_sheet_flow_parameters,
    # fung1972_alveolar_sheet_elasticity, clark2010_acinar_sheet_flow_parameters, clark2011_capillary_compliance_lung_volume,
    # weibel2009_human_capillary_volume_morphometry]
    sheet_area_m2_per_ml = 0.04484f0,
    # cP. Apparent blood viscosity in the sheet. [huang2011_human_sheet_flow_parameters,
    # clark2010_acinar_sheet_flow_parameters]
    sheet_mu_cP = 1.92f0,
    # 1. Product k*f of the sheet-flow friction factors (Huang 2011: k 12, f 1.8; Clark 2010: f 21.6).
    # [huang2011_human_sheet_flow_parameters, clark2010_acinar_sheet_flow_parameters]
    sheet_kf = 21.6f0,
    # um. Mean arteriole-to-venule path length in the sheet: cat 556 um scaled by the human/cat alveolar diameter
    # ratio (Huang 2011; Clark 2010 uses 1186). [huang2011_human_sheet_flow_parameters,
    # clark2010_acinar_sheet_flow_parameters]
    sheet_path_um = 1188.0f0,
    # cmH2O. Transmural pressure above which the sheet no longer thickens: 32 cmH2O for human lungs (Clark 2011, from
    # Sobin et al. 1979 human septa; the reported 7.7 um maximum sheet height with h0 3.5 um and 0.127 um/cmH2O gives
    # about 33). Clark 2010 used 15 with no measurement behind it; cat sheets stay linear to 27 cmH2O (Sobin 1972,
    # Fung & Sobin 1972). [clark2011_capillary_compliance_lung_volume, fung1972_alveolar_sheet_elasticity,
    # sobin1972_cat_sheet_elasticity]
    sheet_pmax_cmH2O = 32.0f0,
    # 1/cmH2O. Derived: median of Huang 2011 Table 2 arterial beta for orders 1-5 (below the alveolar diameter):
    # 0.698, 0.482, 0.361, 0.361, 0.361 %/cmH2O -> 0.361 %/cmH2O (beta in %/cmH2O, see yen1990 notes).
    # [huang2011_human_sheet_flow_parameters, yen1990_human_pulmonary_vessel_elasticity]
    beta_art_small = 0.00361f0,
    # 1/cmH2O. Derived: median of Huang 2011 Table 2 arterial beta for orders 6-16: 1.135, 0.932, 0.877, 1.238 and
    # 1.281 (x7) %/cmH2O -> 1.281 %/cmH2O. [huang2011_human_sheet_flow_parameters,
    # yen1990_human_pulmonary_vessel_elasticity]
    beta_art_large = 0.01281f0,
    # 1/cmH2O. Derived: median of Huang 2011 Table 1 venous beta for orders 1-5: 1.010, 0.625, 0.433, 1.798, 1.260
    # %/cmH2O -> 1.010 %/cmH2O. [huang2011_human_sheet_flow_parameters, yen1990_human_pulmonary_vessel_elasticity]
    beta_ven_small = 0.0101f0,
    # 1/cmH2O. Derived: median of Huang 2011 Table 1 venous beta for orders 6-15: 1.260, 0.887 and 0.708 (x8) %/cmH2O
    # -> 0.708 %/cmH2O. [huang2011_human_sheet_flow_parameters, yen1990_human_pulmonary_vessel_elasticity]
    beta_ven_large = 0.00708f0,
    # um. Human alveolar diameter used by Huang 2011 (250 um, the human/cat ratio 250/117 of the path length): vessels
    # below it see alveolar pressure, larger ones pleural pressure (Huang 2011 Eq. 1).
    # [huang2011_human_sheet_flow_parameters]
    d_alveolus_um = 250.0f0,
    # cmH2O. Derived: mean of end-expiratory pleural (esophageal) pressure upright/sitting, -3.7 (Washko 2006, from
    # transpulmonary 3.7 at relaxation volume) and -3.3 cmH2O (Owens 2012, lean adults sitting).
    # [washko2006_esophageal_pressure_posture, owens2012_esophageal_pressure_sitting_supine]
    p_pleural_cmH2O = -3.5f0,
    # cmH2O. Transmural pressure at which the grown radii hold; distension D = D0 (1 + beta (P_tm - p_ref)). 0, as in
    # Huang 2011's validated whole-lung model (its "PA - Ppl = 10 cmH2O" is the casts' inflation, not a vessel pressure;
    # with 0 this solver reproduces his Table 4 flows within 1%). Caveat: the grown radii scale from XCAT's in vivo main-PA
    # radius, so the true reference of the largest vessels lies between their resting pressure and zero; at 0 the main PA
    # distends ~34% at rest. [huang2011_human_sheet_flow_parameters]
    p_ref_cmH2O = 0.0f0,
    # Vessels inside each acinus (below a grown terminal / above its venule), in series with its sheet: Huang 2011's whole-
    # lung tables per order (smallest first; number, diameter and length in mm of his casts of inflated lungs at zero
    # transmural pressure; beta per cmH2O, printed in %/cmH2O); every order narrower than the bed's terminal / venule.
    # These orders are inputs, so the lung morphometry rows score only the orders the trees resolve. Chosen after
    # comparison with Horsfield 1978's independent arteries (his orders 1-5 are extrapolated precapillaries counted strictly
    # in series; 1.1-1.3 vs 3.3-4.0 L/min on whole anatomy, QUEUE A42-7); Huang's set is validated against measured
    # pressures [huang2011_human_sheet_flow_parameters, horsfield1978_human_small_pulmonary_arteries]
    sub_art = OrderTable(
        [102411624, 28114394, 10203806, 4513692, 1348338, 571544, 172040, 44008, 12450, 3448, 900, 254, 86, 14, 4, 1],
        [0.02, 0.036, 0.056, 0.097, 0.15, 0.22, 0.34, 0.51, 0.77, 1.16, 1.75, 2.71, 4.16, 7.34, 14.8, 30],
        [0.22, 0.26, 0.36, 0.45, 0.68, 1.08, 1.92, 2.81, 3.73, 6.58, 12.35, 18.07, 25.97, 35.69, 25.3, 90.5],
        [0.00698, 0.00482, 0.00361, 0.00361, 0.00361, 0.01135, 0.00932, 0.00877, 0.01238, 0.01281, 0.01281, 0.01281, 0.01281, 0.01281, 0.01281, 0.01281],
    ),
    sub_ven = OrderTable(
        [79647106, 16988490, 4332666, 906102, 146050, 31134, 7442, 2134, 870, 286, 126, 62, 22, 8, 4],
        [0.018, 0.031, 0.067, 0.13, 0.23, 0.38, 0.62, 0.9, 1.42, 1.99, 2.88, 4, 5.86, 8.65, 12.97],
        [0.13, 0.21, 0.38, 1.06, 1.5, 2.92, 4.79, 6.78, 11.24, 14.78, 17.9, 26.49, 19.49, 34.99, 35.68],
        [0.0101, 0.00625, 0.00433, 0.01798, 0.0126, 0.0126, 0.00887, 0.00708, 0.00708, 0.00708, 0.00708, 0.00708, 0.00708, 0.00708, 0.00708],
    ),
)
