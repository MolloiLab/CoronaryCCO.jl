# SPDX-License-Identifier: GPL-3.0-or-later
# Steady lumped-parameter (Poiseuille) hemodynamics: compute_hemodynamics and
# _read_domain_mass_g of coronary_postprocess.py (HCT-CCO commit 9bc141f). CGS units inside.

const MU_VISC_CGS = 0.036   # blood viscosity 3.6 cP in dyn.s/cm^2
const PERF_REST = 1.0       # resting coronary perfusion, mL/min/g
const PERF_STRESS = 3.5     # stress (vasodilated) coronary perfusion, mL/min/g
const TISSUE_DENS = 1.05    # g/cm^3
const _MMHG = 1333.22

"""
    domain_mass_g(vol_path; mass = :mask) -> (mass_g, volume_cm3)

Myocardial mass and volume of a DGtal `.vol` mask. By default (`mass = :mask`) the volume is
decoded with [`read_vol`](@ref) and the voxels the growth domain treats as inside (value > 128)
are counted, at 1.05 g/cm^3.

`mass = :upstream` reproduces upstream `_read_domain_mass_g` byte for byte, for verification
only. It counts the non-zero bytes that follow the header, which is the voxel count only for an
uncompressed (Version 2) file with a zero background. For the zlib-compressed (Version 3)
sample `case0.vol` it counts bytes of the compressed stream: 41.74 g instead of 267.7 g.
"""
function domain_mass_g(vol_path::AbstractString; mass::Symbol = :mask)
    raw = read(vol_path)
    he = findfirst(b"\n.\n", raw)
    he === nothing && error("domain_mass_g: no header terminator in $vol_path")
    header = String(raw[1:(first(he) - 1)])
    d = Dict{String, Int}()
    for m in eachmatch(r"(X|Y|Z|Voxel-Size):\s*(-?\d+)", header)
        d[m.captures[1]] = parse(Int, m.captures[2])
    end
    vs = get(d, "Voxel-Size", 1)
    if mass === :upstream
        lo = last(he) + 1
        hi = min(length(raw), last(he) + d["X"] * d["Y"] * d["Z"])
        nvox = count(!=(0x00), @view raw[lo:hi])
    elseif mass === :mask
        img, _ = read_vol(vol_path)
        nvox = count(>(MASK_THRESHOLD), img)
    else
        throw(ArgumentError("mass must be :upstream or :mask, got :$mass"))
    end
    v_cm3 = nvox * vs^3 / 1000.0
    return v_cm3 * TISSUE_DENS, v_cm3
end

# Ports compute_hemodynamics' network set-up: tree parent from edge endpoints, Poiseuille
# resistance per segment (dyn.s/cm^5), and the post-order from the first parentless segment.
function _hemo_network(t::PPTree)
    n = nsegments(t)
    node2seg = Dict(t.to[i] => i for i in 1:n)
    parent = [get(node2seg, t.from[i], 0) for i in 1:n]
    children = [Int[] for _ in 1:n]
    for i in 1:n
        parent[i] == 0 || push!(children[parent[i]], i)
    end
    root = findfirst(==(0), parent)
    rseg = map(1:n) do i
        a, b = t.points[t.from[i]], t.points[t.to[i]]
        len = max(norm3(sub3(b, a)), 1.0e-3)
        8.0 * MU_VISC_CGS * (len / 10.0) / (pi * cpow(t.radius[i] / 10.0, 4.0))
    end
    post = _postorder!(Int[], falses(n), root, children)
    return (; n, children, root, rseg, post)
end

# Ports compute_hemodynamics.solve: bisect the distal (terminal) resistance so that a 100 mmHg
# inlet drives the prescribed total flow `qtot` (cm^3/s), then split the flow top-down by
# conductance. Returns per-segment flow (cm^3/s), distal pressure (mmHg) and WSS (dyn/cm^2).
function _hemo_solve(net, radius::Vector{Float64}, qtot::Float64; pin = 100.0, pv = 0.0)
    (; n, children, root, rseg, post) = net
    req = zeros(n)
    function req_of!(rt)
        for i in post
            if isempty(children[i])
                req[i] = rseg[i] + rt
            else
                s = 0.0
                for c in children[i]
                    s += 1.0 / req[c]
                end
                req[i] = rseg[i] + 1.0 / s
            end
        end
        return req
    end
    dp = (pin - pv) * _MMHG
    lo, hi = 0.0, 1.0e14
    for _ in 1:80
        mid = (lo + hi) / 2.0
        if dp / req_of!(mid)[root] > qtot
            lo = mid
        else
            hi = mid
        end
    end
    req_of!((lo + hi) / 2.0)
    q = zeros(n)
    q[root] = qtot
    for i in Iterators.reverse(post)
        isempty(children[i]) && continue
        inv = 0.0
        for c in children[i]
            inv += 1.0 / req[c]
        end
        for c in children[i]
            q[c] = q[i] * (1.0 / req[c]) / inv
        end
    end
    p = zeros(n)
    p[root] = (pin * _MMHG - q[root] * rseg[root]) / _MMHG
    for i in Iterators.reverse(post), c in children[i]
        p[c] = p[i] - q[c] * rseg[c] / _MMHG
    end
    wss = [4.0 * MU_VISC_CGS * q[i] / (pi * cpow(radius[i] / 10.0, 3.0)) for i in 1:n]
    return (flow = q, pressure = p, wss = wss)
end

"""
    hemodynamics(t::PPTree, orders, vol_path; mass = :mask) -> NamedTuple

The lumped resistive-network solve of compute_hemodynamics at resting and stress perfusion,
with per-order medians (orders in increasing order) as written to `<prefix>_hemodynamics.csv`.
"""
function hemodynamics(t::PPTree, orders::Vector{Int}, vol_path::AbstractString; mass::Symbol = :mask)
    mass_g, v_cm3 = domain_mass_g(vol_path; mass)
    net = _hemo_network(t)
    resting = _hemo_solve(net, t.radius, PERF_REST * mass_g / 60.0)
    stress = _hemo_solve(net, t.radius, PERF_STRESS * mass_g / 60.0)
    diam = 2.0 .* t.radius
    ords = sort(unique(filter(>(0), orders)))
    med(x, o) = _np_median(x[orders .== o])
    return (
        mass_g = mass_g, volume_cm3 = v_cm3, inlet_pressure = 100.0, resting = resting,
        stress = stress, order = ords, median_diameter = [med(diam, o) for o in ords],
        wss_resting = [med(resting.wss, o) for o in ords], wss_stress = [med(stress.wss, o) for o in ords],
        pressure_resting = [med(resting.pressure, o) for o in ords],
        pressure_stress = [med(stress.pressure, o) for o in ords],
        flow_stress = [med(stress.flow, o) * 60.0 for o in ords],
    )
end

function write_hemodynamics_csv(path, h)
    open(path, "w") do io
        _csvrow(io, ["myocardial_mass_g", _pyfmt(h.mass_g, 2)])
        _csvrow(io, ["resting_perfusion_mL/min/g", string(PERF_REST)])
        _csvrow(io, ["stress_perfusion_mL/min/g", string(PERF_STRESS)])
        _csvrow(io, ["inlet_pressure_mmHg", _pyfmt(h.inlet_pressure, 0)])
        _csvrow(io, String[])
        _csvrow(
            io, [
                "order", "median_diameter_mm", "WSS_resting_dyn/cm2", "WSS_stress_dyn/cm2",
                "pressure_resting_mmHg", "pressure_stress_mmHg", "flow_stress_mL/min",
            ]
        )
        for k in eachindex(h.order)
            _csvrow(
                io, [
                    h.order[k], _pyfmt(h.median_diameter[k], 3), _pyfmt(h.wss_resting[k], 1),
                    _pyfmt(h.wss_stress[k], 1), _pyfmt(h.pressure_resting[k], 1),
                    _pyfmt(h.pressure_stress[k], 1), _pyfmt(h.flow_stress[k], 4),
                ]
            )
        end
    end
    return path
end
