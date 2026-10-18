# SPDX-License-Identifier: GPL-3.0-or-later
# Validates the bolus model (body.jl) against human first-pass data, each study simulated with its
# own injection protocol, injection site and cardiac output, and scored against the values stored
# in data/reference. Usage: julia --project=. bolus_validation.jl [out.json]
import Vasculature; (isdefined(Main, :Body) || @eval(import Vasculature: Body)); (isdefined(Main, :Flow) || @eval(import Vasculature: Flow)); using .Body, .Flow, Printf, TOML
# the coronary bed: XCAT vmale_50 myocardium (LV + RV labels, 1 mm grid, 148.1 mL) with the heart's physiology, and one
# blood for the whole body (the heart's hematocrit), through the one transport spec (Vasculature.TransportSpec)
const HEART = HEART_FLOW
human_circulation(; kw...) = Vasculature.circulation(Vasculature.TransportSpec(Vasculature.HEART_TRANSPORT; tissue_ml = 148.1f0, kw...))

const REF = joinpath(@__DIR__, "..", "..", "data", "reference")
ref(id) = TOML.parsefile(joinpath(REF, id * ".toml"))["data"]
const HU = Dict(80 => 41.12f0, 100 => 31.74f0, 120 => 26.18f0)   # HU per mgI/mL of blood (Bae 2010)
const BSA = 1.97f0                                                 # XCAT vmale_50 (Du Bois, 176 cm, 77 kg)
const CO0 = human_circulation().co

rows = Dict{String, Any}[]
function row!(id, quantity, model, mean, sd; unit = "", note = "", score = true, calibrated = false)
    # one rule page-wide (Dale's decision B2, 2026-10-16): one model subject passes inside the population's ±2 SD
    pass = calibrated ? "calibrated" : score && isfinite(sd) ? abs(model - mean) <= 2 * sd : missing
    # a pass within 1% of the 2-SD bound flips on noise far inside the reference's own spread: say so (audit A16-3)
    pass === true && 2 * sd - abs(model - mean) < 0.01 * 2 * sd &&
        (note = (isempty(note) ? "" : note * "; ") * @sprintf("edge pass: %.3g inside the 2-SD bound %.4g, not counted as a gain", 2 * sd - abs(model - mean), mean + sign(model - mean) * 2 * sd))
    push!(rows, Dict("ref" => id, "quantity" => quantity, "model" => model, "mean" => mean, "sd" => sd, "unit" => unit, "pass" => pass, "note" => note))
    return @printf("%-40s %-44s model %7.2f  ref %7.2f ± %-6.2f %s %s\n", id, quantity, model, mean, sd, unit, pass === missing ? "info" : pass == "calibrated" ? "CAL" : pass ? "PASS" : "FAIL")
end
circ(co_Lmin; kw...) = human_circulation(; co = Float32(co_Lmin * 1000 / 60), kw...)
beat_for(co_Lmin) = Float32(60 * 81.6 / (co_Lmin * 1000))           # stroke volume fixed, rate follows output
# time of the first sample at or above a threshold, sampled every `dt` s from `t0`
function sampled(r, c, dt; t0 = 0.0f0)
    ts = t0:dt:r.t[end]; return ts, [c[clamp(round(Int, t / (r.t[2] - r.t[1])) + 1, 1, length(c))] for t in ts]
end
function first_at(r, c, thr; dt = 1.0f0)
    ts, y = sampled(r, c, dt); k = findfirst(>=(thr), y)
    return k === nothing ? NaN : ts[k]
end
ttp(r, c) = r.t[argmax(c)]
function fwhm(r, c)
    m = maximum(c); k = argmax(c); a = findlast(<(m / 2), c[1:k]); b = findfirst(<(m / 2), c[k:end])
    return (a === nothing || b === nothing) ? NaN : r.t[k + b - 1] - r.t[a]
end
dur_above(r, c, thr) = count(>=(thr), c) * (r.t[2] - r.t[1])

# ---- Dymond et al. 1982 (first-pass radionuclide): basilic-vein hand push of a small bolus with a
# rapid 20 mL saline flush (rate not stated: 1 mL at 2 mL/s, flush at 5 mL/s); SVC width at half
# height in the good-quality studies (sets the large veins' dispersion, so marked calibrated)
let d = ref("dymond1982_svc_bolus_fwhm")["svc_fwhm_interobserver_sample"]
    C = human_circulation(); r = simulate(C, Injection(; rate = 2, volume = 1, saline = 20, saline_rate = 5); dt = 0.001, every = 0.02, tmax = 30)
    row!(
        "dymond1982", "SVC bolus width at half height (arm veins)", fwhm(r, curve(r, :svc; what = :mean)), d["group_mean_fwhm_s"][1], d["group_sd_fwhm_s"][1]; unit = "s",
        note = "this measurement sets the large veins' dispersion", calibrated = true
    )
end
# ---- Borges et al. 2022: first arrival (+50 HU, 1-s sampling, ~100 kVp), right antecubital,
# 4 mL/s, 1.5 mL/kg (XCAT 77 kg → 115 mL) + 30 mL saline; cardiac index 2.6 L/min/m²
let d = ref("borges2022_bolus_arrival_times")["arrival_times"], co = 2.6f0 * BSA
    C = circ(co; beat = beat_for(co)); r = simulate(C, Injection(; site = :antecubital_r, rate = 4, volume = 115, saline = 30))
    thr = 50 / HU[100]
    a = [first_at(r, curve(r, s; what = :mean), thr) for s in (:svc, :pulm_trunk)]; aa = first_at(r, curve(r, :lv), thr)
    m = [a[1], a[2] - a[1], aa - a[2], aa]
    for (k, s) in enumerate(d["segment"])
        row!("borges2022", "first arrival, " * replace(s, "_to_" => " → ", "_" => " ", "svc" => "SVC"), m[k], d["time_mean_s"][k], d["time_sd_s"][k]; unit = "s")
    end
end
# ---- Leiknes et al. 2023 test bolus: 20 mL at 6 mL/s + 40 mL saline, ascending aorta sampled
# every 2 s from 12 s, MRI cardiac output 5.05 L/min (injection arm not stated: right)
let d = ref("leiknes2023_test_bolus_ccta")["test_bolus_timing"], co = 5.05f0
    C = circ(co; beat = beat_for(co)); r = simulate(C, Injection(; site = :antecubital_r, rate = 6, volume = 20, saline = 40))
    ts, y = sampled(r, curve(r, :lv), 2.0f0; t0 = 12.0f0)
    row!("leiknes2023", "test-bolus TTP (2-s samples from 12 s)", ts[argmax(y)], d["mean"][1], d["sd"][1]; unit = "s")
    row!(
        "leiknes2023", "test-bolus peak ΔHU at 100 kV", HU[100] * maximum(curve(r, :lv)), d["mean"][2] - 40, d["sd"][2]; unit = "HU", score = false,
        note = "absolute HU minus ~40 HU blood; not scored: the same curves gave Stewart-Hamilton output 9.16 vs MRI 5.05 L/min, so their area (and peak) is ~55% of what the injected iodine requires"
    )
end
# ---- Hazirolan et al. 2009: 20 mL at 5 mL/s + 20 mL saline, right antecubital, 120 kV, baseline
# subtracted; warmed and room-temperature groups pooled as a range; cardiac output not given
let d = ref("hazirolan2009_test_bolus_pa_aorta")["test_bolus_curves"]
    C = circ(CO0 * 60 / 1000); r = simulate(C, Injection(; site = :antecubital_r, rate = 5, volume = 20, saline = 20))
    for (lbl, c, i) in (("ascending aorta", curve(r, :lv), 1), ("main pulmonary artery", curve(r, :rv), 5))
        pk = d["peak_enhancement_HU"][i:(i + 1)]; sd = d["peak_enhancement_sd_HU"][i:(i + 1)]
        row!("hazirolan2009", "peak ΔHU 120 kV, " * lbl, HU[120] * maximum(c), sum(pk) / 2, (pk[1] - pk[2]) / 2 + sum(sd) / 2; unit = "HU")
        du = d["duration_above_100HU_s"][i:(i + 1)]; dsd = d["duration_above_100HU_sd_s"][i:(i + 1)]
        row!("hazirolan2009", "time above 100 HU, " * lbl, dur_above(r, HU[120] .* c, 100), sum(du) / 2, (du[1] - du[2]) / 2 + sum(dsd) / 2; unit = "s")
    end
    row!("hazirolan2009", "PA / aorta peak ratio", maximum(curve(r, :rv)) / maximum(curve(r, :lv)), (318 / 210 + 283 / 178) / 2, 0.15; note = "ratio of the group means; SD set from the two groups' spread")
end
# ---- Moradi et al. 2016: 15 mL test bolus (rate not stated: 5 mL/s, + 40 mL saline), metoprolol
# to < 65/min, 120 kV, peak printed as attenuation (absolute: ~40 HU blood subtracted)
let d = ref("moradi2016_ef_aortic_test_bolus")["test_bolus_and_main_scan"]
    C = circ(CO0 * 60 / 1000); r = simulate(C, Injection(; site = :antecubital_r, conc = 370, rate = 5, volume = 15, saline = 40))
    row!("moradi2016", "test-bolus TTP", ttp(r, curve(r, :lv)), d["mean"][1], d["sd"][1]; unit = "s")
    row!("moradi2016", "test-bolus peak ΔHU 120 kV", HU[120] * maximum(curve(r, :lv)), d["mean"][2] - 40, d["sd"][2]; unit = "HU")
end
# ---- Skrok et al. 2012 (MRI): 0.8 s bolus at 5 mL/s (4 mL) + 20 mL saline, healthy controls,
# cardiac index 3.01 L/min/m²; RV → LV peak-to-peak and LV FWHM (median, IQR → SD ≈ IQR/1.35)
let d = ref("skrok2012_mri_first_pass_bolus_kinetics")["bolus_kinetics"], co = 3.01f0 * BSA
    C = circ(co; beat = beat_for(co)); r = simulate(C, Injection(; site = :antecubital_r, rate = 5, volume = 4, saline = 20))
    rv = curve(r, :rv; what = :mean); lv = curve(r, :lv; what = :mean)
    row!("skrok2012", "RV → LV peak-to-peak transit", ttp(r, lv) - ttp(r, rv), d["ptt_rv_to_lv_peak_to_peak_median_s"][1], (d["ptt_q3_s"][1] - d["ptt_q1_s"][1]) / 1.35; unit = "s")
    row!("skrok2012", "LV FWHM", fwhm(r, lv), d["lv_fwhm_median_s"][1], (d["lv_fwhm_q3_s"][1] - d["lv_fwhm_q1_s"][1]) / 1.35; unit = "s")
end
# ---- Lehmonen et al. 2022 (MRI): gadoterate 0.1 mmol/kg (first half of 0.2 mmol/kg; 15 mL of
# 0.5 M for 75 kg) at 4.5 mL/s, heart rate 70/min; times in cardiac cycles (IQR → SD ≈ IQR/1.35)
let d = ref("lehmonen2022_cmr_chamber_bolus_passage"), co = 81.6f0 * 70 / 1000
    C = circ(co; beat = 60.0f0 / 70); b = 60.0f0 / 70
    r = simulate(C, Injection(; site = :antecubital_r, rate = 4.5, volume = 15, saline = 20, saline_rate = 4.5))
    sv = curve(r, :svc; what = :mean); la = curve(r, :la; what = :mean); lv = curve(r, :lv; what = :mean)
    tt = d["transit_times"]; tv = d["time_variables_controls"]
    row!("lehmonen2022", "SVC → LA peak-to-peak (cycles)", (ttp(r, la) - ttp(r, sv)) / b, tt["controls_median_cc"][2], tt["controls_iqr_width_cc"][2] / 1.35)
    row!("lehmonen2022", "SVC → LV peak-to-peak (cycles)", (ttp(r, lv) - ttp(r, sv)) / b, tt["controls_median_cc"][3], tt["controls_iqr_width_cc"][3] / 1.35)
    for (k, c) in ((1, sv), (3, la), (4, lv))
        # not scored: at 0.1 mmol/kg the saturation-recovery blood-pool signal is far past its linear range (linear only to
        # ≈ 0.01 mmol/kg, Utz 2007; peaks clipped, Li 2023, Kellman 2017), so the signal curve's width overstates the
        # concentration curve's; the CT width rows (iodine is linear in HU) test the same thing
        row!(
            "lehmonen2022", "FWHM (cycles), " * tv["site"][k], fwhm(r, c) / b, tv["fwhm_median_cc"][k], tv["fwhm_iqr_width_cc"][k] / 1.35;
            score = false, note = "not scored: MRI signal saturates at this dose (blood pool linear only to ≈ 0.01 mmol/kg), so the measured width is a signal width"
        )
    end
end
# ---- Dock et al. 1961 (dye dilution): pulmonary artery → left atrium mean transit, normal dynamics
let d = ref("dock1961_pulmonary_blood_volume")["normal_dynamics_subjects"], x = d["mtt_pa_to_la_s"], ci = d["cardiac_index_avg_L_per_min_per_m2"]
    co = sum(ci) / length(ci) * BSA; C = circ(co)
    tt = Dict(t.name => t.mtt for t in transit_table(C))
    m = sum(x) / length(x); s = sqrt(sum(abs2, x .- m) / (length(x) - 1))
    row!("dock1961", "PA → LA mean transit (at their CI)", tt[:pulm_trunk] / 2 + tt[:lung], m, s; unit = "s")
end
# ---- Bae 2010: arrival windows and the injection-site shift (statements of a review)
let d = ref("bae2010_ct_contrast_timing_review")
    C = human_circulation(); inj(site) = Injection(; site, rate = 5, volume = 100, saline = 40)
    r = simulate(C, inj(:antecubital_r)); thr = 50 / HU[120]
    w = d["normal_arrival_times"]
    for (k, c) in ((1, curve(r, :rv)), (2, curve(r, :lv)), (3, curve(r, :desc_aorta)))
        lo, hi = w["arrival_time_min_s"][k], w["arrival_time_max_s"][k]
        row!(
            "bae2010", "arrival (+50 HU), " * w["target"][k], first_at(r, c, thr; dt = 0.1f0), (lo + hi) / 2, (hi - lo) / 2; unit = "s", score = false,
            note = "review statement, not a measurement"
        )
    end
    rc = simulate(C, inj(:svc_top)); sh = d["injection_site_shift"]
    row!(
        "bae2010", "central venous vs antecubital shift", first_at(rc, curve(rc, :lv), thr; dt = 0.1f0) - first_at(r, curve(r, :lv), thr; dt = 0.1f0),
        (sh["timing_shift_min_s"][2] + sh["timing_shift_max_s"][2]) / 2, (sh["timing_shift_max_s"][2] - sh["timing_shift_min_s"][2]) / 2; unit = "s", score = false,
        note = "review statement, not a measurement"
    )
    rp = simulate(C, Injection(; site = :svc_top, rate = 3.3, volume = 75, saline = 30)); ra = simulate(C, Injection(; rate = 3.3, volume = 75, saline = 30))
    d2 = ref("deng2026_port_vs_peripheral_ct")["aortic_attenuation_by_route"]
    row!(
        "deng2026", "aortic peak, central port ÷ arm cannula", maximum(curve(rp, :lv)) / maximum(curve(ra, :lv)), d2["aorta_all_median_HU"][2] / d2["aorta_all_median_HU"][1], 0.05;
        note = "ratio of medians, 195 paired patients; SD set to 0.05 (paired, no SD of the ratio given); bolus-tracked scan, model compares peaks"
    )
end
# ---- conservation
let C = human_circulation(), r = simulate(C, Injection(; rate = 6, volume = 80, saline = 100); tmax = 120)
    row!("mass balance", "iodine unaccounted, % of injected", 100 * (r.injected - r.held - r.excreted) / r.injected, 0.0, 0.5; unit = "%")
end
# ---- injection sites: aortic and pulmonary-trunk curves for one protocol (50 mL at 5 mL/s)
sites = Dict{String, Any}()
let C = human_circulation()
    println("\ninjection site     PT arrival  PT TTP  PT peak   AA arrival  AA TTP  AA peak (ΔHU, 120 kV)")
    for site in (:antecubital_r, :antecubital_l, :jugular_r, :jugular_l, :svc_top, :femoral_r)
        r = simulate(C, Injection(; site)); pt = curve(r, :rv); aa = curve(r, :lv); thr = 50 / HU[120]
        x = (pt_arr = first_at(r, pt, thr; dt = 0.1f0), pt_ttp = ttp(r, pt), pt_pk = HU[120] * maximum(pt), aa_arr = first_at(r, aa, thr; dt = 0.1f0), aa_ttp = ttp(r, aa), aa_pk = HU[120] * maximum(aa))
        @printf("%-16s %8.1f %8.1f %8.0f %10.1f %8.1f %8.0f\n", site, x...)
        sites[string(site)] = Dict{String, Any}(string(k) => v for (k, v) in pairs(x))
        sites[string(site)]["t"] = r.t[1:5:end]; sites[string(site)]["aa"] = HU[120] .* aa[1:5:end]; sites[string(site)]["pt"] = HU[120] .* pt[1:5:end]
    end
end
np = count(r -> r["pass"] === true, rows); nf = count(r -> r["pass"] === false, rows); nc = count(r -> isequal(r["pass"], "calibrated"), rows)
@printf("\n%d pass, %d calibrated, %d fail, %d info\n", np, nc, nf, length(rows) - np - nf - nc)
if !isempty(ARGS)
    js(x) = x isa AbstractString ? "\"" * replace(x, "\"" => "\\\"") * "\"" : x isa AbstractDict ? "{" * join(["\"$k\":" * js(v) for (k, v) in x], ",") * "}" :
        x isa AbstractVector ? "[" * join(js.(x), ",") * "]" : x === missing ? "null" : x isa Bool ? string(x) : x isa Real ? (isfinite(x) ? @sprintf("%.4g", x) : "null") : string(x)
    path = [
        Dict{String, Any}("name" => string(x.name), "flow" => x.flow, "mtt" => x.mtt, "rd" => x.rd, "v" => s.v)
            for (x, s) in zip(transit_table(human_circulation()), human_circulation().segs) if s.kind !== :bed
    ]
    open(io -> print(io, js(Dict("rows" => rows, "sites" => sites, "path" => path))), ARGS[1], "w")
end
