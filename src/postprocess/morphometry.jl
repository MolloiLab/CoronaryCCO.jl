# SPDX-License-Identifier: GPL-3.0-or-later
# Morphometry tables and their CSV files (coronary_postprocess.py, HCT-CCO commit 9bc141f).

# Ports compute_mean_se: mean and standard error (sample SD / sqrt(n)) of each count list,
# returned as matrices indexed [daughter order, parent order] over `present`.
function connectivity_mean_se(present::Vector{Int}, counts)
    k = length(present)
    mean = zeros(k, k)
    se = zeros(k, k)
    for (i, m) in enumerate(present), (j, p) in enumerate(present)
        c = Float64.(counts[(m, p)])
        isempty(c) && continue
        mean[i, j] = _np_mean(c)
        se[i, j] = _np_std(c; ddof = 1) / sqrt(Float64(length(c)))
    end
    return mean, se
end

# Ports write_matrix_to_csv.
function write_connectivity_csv(path, present, mean, se)
    open(path, "w") do io
        _csvrow(io, ["Order m \\ Order n"; string.(present)])
        for (i, m) in enumerate(present)
            row = [string(m)]
            for j in eachindex(present)
                push!(row, mean[i, j] < 0.005 ? "0" : _pyfmt(mean[i, j], 3) * " ± " * _pyfmt(se[i, j], 3))
            end
            _csvrow(io, row)
        end
    end
    return path
end

# Ports compute_diameter_stats_per_order and compute_diameter_ranges_per_order: per order (in
# increasing order) the mean, sample SD, minimum and maximum diameter.
function diameter_stats_per_order(orders::Vector{Int}, diameters::Vector{Float64})
    ords = sort(unique(orders))
    vals = [diameters[orders .== o] for o in ords]
    return (
        order = ords, mean = _np_mean.(vals), sd = [_np_std(v; ddof = 1) for v in vals],
        min = minimum.(vals), max = maximum.(vals), median = _np_median.(vals),
        q25 = [_np_percentile(v, 25) for v in vals], q75 = [_np_percentile(v, 75) for v in vals],
    )
end

# Ports write_diameter_stats_to_csv and write_diameter_ranges_to_csv.
function write_diameter_stats_csv(path, st)
    open(path, "w") do io
        _csvrow(io, ["Order", "Diameter (units) ± SD"])
        for k in eachindex(st.order)
            _csvrow(io, [st.order[k], _pyfmt(st.mean[k], 2) * " ± " * _pyfmt(st.sd[k], 5)])
        end
    end
    return path
end

function write_diameter_ranges_csv(path, st)
    open(path, "w") do io
        _csvrow(io, ["Order", "Diameter Range (units)"])
        for k in eachindex(st.order)
            _csvrow(io, [st.order[k], _pyfmt(st.min[k], 2) * " - " * _pyfmt(st.max[k], 2)])
        end
    end
    return path
end
