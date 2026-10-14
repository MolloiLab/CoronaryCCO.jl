# SPDX-License-Identifier: GPL-3.0-or-later
# Figures of HCT-CCO postprocessing/coronary_postprocess.py (commit 9bc141f):
# plot_diameter_vs_order and the figure drawn in compute_hemodynamics. Sizes are in points
# (1/72 in) with matplotlib's default font sizes, saved at 150 dpi like the Python.
module VasculatureMakieExt

using CairoMakie
using Vasculature: Vasculature, LITERATURE_MU, LITERATURE_SD, PERF_REST, PERF_STRESS

const DPI = 150
const INK = colorant"#2c3e50"
const RED = colorant"#c0392b"
const BLUE = colorant"#2c7fb8"

_axis(pos; kw...) = Axis(
    pos; titlesize = 12, titlefont = :regular, xlabelsize = 10, ylabelsize = 10,
    xticklabelsize = 10, yticklabelsize = 10, kw...
)
# matplotlib-like log axis: labels at whole decades, unlabelled ticks at 2..9 x 10^k
_decades(lo, hi) = 10.0 .^ (floor(Int, log10(lo)):ceil(Int, log10(hi)))
_logkw(a) = (
    Symbol(a, :ticks) => _decades, Symbol(a, :minorticksvisible) => true,
    Symbol(a, :minorticks) => IntervalsBetween(9),
    Symbol(a, :tickformat) => vs -> [rich("10", superscript(string(round(Int, log10(v))))) for v in vs],
)
_legend(ax, pos, size) = axislegend(
    ax; position = pos, labelsize = size, framecolor = (:black, 0.2), rowgap = 0,
    padding = (4, 4, 2, 2), patchsize = (20, 10)
)

function Vasculature.plot_postprocess(r::NamedTuple, prefix::AbstractString)
    files = [plot_diameter_vs_order(r.diameter_stats, prefix * "_diameter_vs_order.png")]
    r.hemodynamics === nothing || push!(files, plot_hemodynamics(r.hemodynamics, prefix * "_hemodynamics.png"))
    return files
end

# Ports plot_diameter_vs_order: the tree's per-order median and interquartile band over the
# literature mean +/- SD, log diameter axis.
function plot_diameter_vs_order(st, path)
    fig = Figure(size = (8 * 72, 5.5 * 72), fontsize = 10)
    ax = _axis(
        fig[1, 1]; yscale = log10, xticks = 1:length(LITERATURE_MU),
        xlabel = "Diameter-defined Strahler order", ylabel = "Vessel diameter (mm)",
        title = "Diameter vs order: our tree vs literature", _logkw(:y)...,
        xgridstyle = :dot, ygridstyle = :dot, yminorgridvisible = true, yminorgridstyle = :dot,
        xgridcolor = (:black, 0.12), ygridcolor = (:black, 0.12),
        yminorgridcolor = (:black, 0.12)
    )
    lit = 1:length(LITERATURE_MU)
    scatterlines!(
        ax, lit, LITERATURE_MU; color = INK, linewidth = 2, markersize = 9, marker = :circle,
        label = "Literature (Kassab)"
    )
    band!(ax, lit, LITERATURE_MU .- LITERATURE_SD, LITERATURE_MU .+ LITERATURE_SD; color = (INK, 0.12))
    scatterlines!(
        ax, st.order, st.median; color = RED, linewidth = 2, markersize = 9, marker = :rect,
        linestyle = :dash, label = "Our tree (median)"
    )
    band!(ax, st.order, st.q25, st.q75; color = (RED, 0.15), label = "Our tree IQR")
    _legend(ax, :lt, 10)
    save(path, fig; px_per_unit = DPI / 72)
    return path
end

# Ports the figure of compute_hemodynamics: per-order median WSS and pressure against the
# per-order median diameter, at resting and stress perfusion.
function plot_hemodynamics(h, path)
    fig = Figure(size = (12 * 72, 5 * 72), fontsize = 10)
    states = (
        ("resting", PERF_REST, BLUE, :rect, :dash, h.wss_resting, h.pressure_resting),
        ("stress", PERF_STRESS, RED, :circle, :solid, h.wss_stress, h.pressure_stress),
    )
    ax1 = _axis(
        fig[1, 1]; xscale = log10, yscale = log10, xlabel = "vessel diameter (mm)",
        ylabel = "wall shear stress (dyn/cm2)", title = "Wall shear stress vs diameter",
        xgridvisible = false, ygridvisible = false, _logkw(:x)..., _logkw(:y)...
    )
    hspan!(ax1, 10, 70; color = RGBf(0.88, 0.88, 0.88), label = "physiological arterial (10-70)")
    ax2 = _axis(
        fig[1, 2]; xscale = log10, xlabel = "vessel diameter (mm)", ylabel = "pressure (mmHg)",
        title = "Pressure vs diameter (down the tree)", xgridvisible = false, ygridvisible = false,
        _logkw(:x)...
    )
    for (name, q, col, mk, ls, wss, p) in states
        kw = (color = col, linewidth = 2, markersize = 8, marker = mk, linestyle = ls, label = "$name ($q mL/min/g)")
        scatterlines!(ax1, h.median_diameter, wss; kw...)
        scatterlines!(ax2, h.median_diameter, p; kw...)
    end
    _legend(ax1, :rt, 8)
    _legend(ax2, :lt, 8)
    Label(
        fig[0, :], "Hemodynamics (lumped-parameter solve): inlet 100 mmHg, " *
            "resting $PERF_REST vs stress $PERF_STRESS mL/min/g"; fontsize = 12
    )
    save(path, fig; px_per_unit = DPI / 72)
    return path
end

end
