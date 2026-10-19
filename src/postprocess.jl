# SPDX-License-Identifier: GPL-3.0-or-later
# Post-processing: pure-Julia port of HCT-CCO postprocessing/coronary_postprocess.py (commit
# 9bc141f), entry point run(xml_path, prefix, domain_vol). Writes the data upstream plots, not the plots.

include("postprocess/numpy.jl")
include("postprocess/tree.jl")
include("postprocess/ordering.jl")
include("postprocess/morphometry.jl")
include("postprocess/hemodynamics.jl")
include("postprocess/vtk.jl")

"""
    postprocess(xml_path, prefix; domain_vol = nothing, mass = :mask, elements = :series,
                obj = true, verbose = true) -> NamedTuple

Port of upstream `coronary_postprocess.run(xml_path, prefix, domain_vol)`. Reads a tree XML
written by [`export_gxl`](@ref) and writes

  - `<prefix>_connectivity_matrix.csv`, `<prefix>_diameter_stats.csv`, `<prefix>_diameter_ranges.csv`
  - `<prefix>_tree.vtp` (cell arrays `radius` and `final_order`) and `<prefix>_tree.obj`
    (24-sided capped tube mesh; skipped with `obj = false`)
  - `<prefix>_hemodynamics.csv` when `domain_vol` is given

Returns the computed data: per-segment `diameters` and diameter-defined `orders`, the ordering
internals (`strahler`, `seg2elem`, `elem_order`), the `connectivity` mean/SE matrices
(indexed [daughter order, parent order]), per-order `diameter_stats`, and `hemodynamics`
(`nothing` without a domain).

Two upstream bugs are corrected by default; each has an option that restores the upstream
behaviour, used to verify the port against the Python:

  - `mass = :mask` decodes the domain and counts voxels > 128; `mass = :upstream` counts the
    non-zero bytes of the (possibly compressed) file payload, see [`domain_mass_g`](@ref).
  - `elements = :series` builds elements from same-order segments joined end to end;
    `elements = :upstream` also joins siblings, see [`diameter_defined_orders`](@ref).

`postprocess(...; mass = :upstream, elements = :upstream)` reproduces the Python outputs.
"""
function postprocess(
        xml_path::AbstractString, prefix::AbstractString; domain_vol = nothing,
        mass::Symbol = :mask, elements::Symbol = :series, obj::Bool = true,
        verbose::Bool = true
    )
    mkpath(dirname(abspath(prefix)))
    t = read_tree_xml(xml_path)
    diameters = 2 .* t.radius
    ord = diameter_defined_orders(t, diameters; elements)
    orders = ord.orders

    present, counts = element_connectivity_counts(ord.seg2elem, ord.elem_order, ord.conn)
    cmean, cse = connectivity_mean_se(present, counts)
    stats = diameter_stats_per_order(orders, diameters)
    files = String[]
    push!(files, write_connectivity_csv(prefix * "_connectivity_matrix.csv", present, cmean, cse))
    push!(files, write_diameter_stats_csv(prefix * "_diameter_stats.csv", stats))
    push!(files, write_diameter_ranges_csv(prefix * "_diameter_ranges.csv", stats))
    push!(files, write_tree_vtp(prefix * "_tree.vtp", t, orders))
    obj && push!(files, write_tree_obj(prefix * "_tree.obj", t))

    hemo = nothing
    if domain_vol !== nothing && nsegments(t) > 0
        hemo = hemodynamics(t, orders, domain_vol; mass)
        push!(files, write_hemodynamics_csv(prefix * "_hemodynamics.csv", hemo))
        if verbose
            ws = hemo.wss_stress
            println(
                "[hemo] lumped solve at resting $PERF_REST and stress $PERF_STRESS mL/min/g, ",
                "inlet 100 mmHg; stress per-order WSS ", _pyfmt(minimum(ws), 0), "-",
                _pyfmt(maximum(ws), 0), " dyn/cm2 -> $(prefix)_hemodynamics.csv"
            )
        end
    end

    result = (
        tree = t, diameters = diameters, orders = orders, strahler = ord.strahler,
        seg2elem = ord.seg2elem, elem_order = ord.elem_order,
        connectivity = (order = present, mean = cmean, se = cse),
        diameter_stats = stats, hemodynamics = hemo, files = files,
    )
    verbose && println("[postprocess] root order=", maximum(orders; init = 0), ". Wrote: ", join(files, ", "))
    return result
end

export postprocess
