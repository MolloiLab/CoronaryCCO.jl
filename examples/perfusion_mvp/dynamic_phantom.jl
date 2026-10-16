# The shareable dynamic CT perfusion phantom: ONE compartment volume plus ONE CSV per phantom design (scenario).
#
# Each voxel of the volume (`<phantom>_compartments.u16.zst`, UInt16, column-major, zstd) holds the id of its kinetic
# compartment: a set of voxels that share one base material and one iodine curve (a chamber blood pool, an artery segment
# at one transit delay, a class of myocardium with nearly the same curve, a plaque sub-region, ...). A scenario (healthy,
# a stenosis, ...) is one file `label_<scenario>.csv` on the same volume, holding everything needed for it:
#   columns  `row`, `value`, one column per compartment id, one per summary region (`region:<name>`), and the
#            bolus-tracking ROI (`roi:<name>`)
#   rows     scenario facts (one value each in `value`: volume file, dims, voxel size, origin, scan times, ...);
#            per-column facts (name, base material with density and composition, kinetics, voxels, ground-truth flow);
#            then one row per time point, `rest` and `stress`, with `value` = time in s and each column's iodine in mg per
#            mL of voxel (regions: their mean; ROI: ΔHU at 70 keV).
# The writer half is organ-generic (heart, lung, brain, ...). It holds no anatomy and no licensed phantom data, so it lives
# in the public repo; the phantom-specific step (which voxel belongs to which compartment) is supplied by the caller.
# Mixing rule: iodine is added to a voxel without displacing tissue, ρ = ρ_base + c/1000 g/cm³, so
# μ(E) = μ_base(E) + (c/1000)(μ/ρ)_I(E) exactly; c is mg iodine per mL of voxel.
module DynamicPhantom

using DelimitedFiles, TOML, Unitful, Printf, Random, Statistics
import XrayAttenuation as XA
import AcceleratedKernels as AK

export compartment_classes, class_error, VolumeWriter, write_slice!, write_phantom_csv, scan_times_from_curve, load, iodine,
    scan_times, materials_at, bs_phantom, iodinated

const F32 = Float32

# =============================================================================================== writing

"""
    compartment_classes(X, K; iters = 40, seed = 1) -> (assign, centres)

k-means of voxel curves into `K` kinetic classes. `X` is (features × voxels): typically every state's (and every
scenario's) curve of one voxel stacked into one column, in HU so the error reads in HU. Clustering on the joint curves
keeps one class homogeneous in every state and scenario, so one compartment volume serves them all. k-means++ seeding on
a subsample; `assign[v]` is the class of voxel v, `centres[:, k]` the class mean curve.
"""
function compartment_classes(X::AbstractMatrix{F32}, K::Integer; iters = 40, seed = 1)
    d, n = size(X); rng = MersenneTwister(seed)
    sub = randperm(rng, n)[1:min(n, 30000)]
    C = zeros(F32, d, K); C[:, 1] = X[:, sub[1]]
    D2 = [sum(abs2, X[:, s] .- C[:, 1]) for s in sub]
    for k in 2:K
        w = cumsum(D2); u = rand(rng) * w[end]; j = searchsortedfirst(w, u)
        C[:, k] = X[:, sub[min(j, length(sub))]]
        for (i, s) in enumerate(sub)
            D2[i] = min(D2[i], sum(abs2, X[:, s] .- C[:, k]))
        end
    end
    a = zeros(Int32, n)
    for _ in 1:iters
        c2 = vec(sum(abs2, C; dims = 1)); G = C' * X; anew = similar(a)
        AK.foreachindex(anew; min_elems = 4096) do v
            best = typemax(F32); bk = Int32(1)
            for k in 1:K
                dd = c2[k] - 2 * G[k, v]; dd < best && (best = dd; bk = Int32(k))
            end
            anew[v] = bk
        end
        changed = count(anew .!= a); a .= anew
        S = zeros(F32, d, K); cnt = zeros(Int, K)
        for v in 1:n
            S[:, a[v]] .+= @view X[:, v]; cnt[a[v]] += 1
        end
        for k in 1:K
            cnt[k] > 0 && (C[:, k] = S[:, k] ./ cnt[k])
        end
        changed == 0 && break
    end
    return a, C
end

"""
    class_error(X, assign, centres) -> (; rms_worst, rms_mean, p99, max)

Error of replacing each voxel's curve by its class curve: RMS over voxels at each feature (time point), the worst and
mean of those, and the 99th percentile and maximum of the absolute error, in the units of `X`.
"""
function class_error(X::AbstractMatrix, a::AbstractVector, C::AbstractMatrix)
    E = X .- C[:, a]; r = vec(sqrt.(mean(abs2, E; dims = 2)))
    return (; rms_worst = maximum(r), rms_mean = mean(r), p99 = quantile(vec(abs.(E)), 0.99), max = maximum(abs, E))
end

"""
    VolumeWriter(path, dims; T = UInt16, level = 10)

Streams a compartment volume to `path` one z slice at a time (column-major: x fastest, then y, then z), through zstd if
`path` ends in `.zst` (a compartment volume is mostly uniform regions: 0.2 mm XCAT-size volumes shrink about 100×), raw
otherwise. `write_slice!(w, slice)` appends one slice; `close(w)` finishes. `w.counts[id + 1]` counts voxels per id.
"""
mutable struct VolumeWriter{T}
    io::IO
    dims::NTuple{3, Int}
    k::Int
    counts::Vector{Int}
end
function VolumeWriter(path::AbstractString, dims; T::Type = UInt16, level::Integer = 10)
    io = endswith(path, ".zst") ? open(pipeline(`zstd -q -f -T4 -$level -o $path`), "w") : open(path, "w")
    return VolumeWriter{T}(io, Tuple(dims), 0, zeros(Int, Int(typemax(T)) + 1))
end
function write_slice!(w::VolumeWriter{T}, s::AbstractMatrix) where {T}
    @assert size(s) == w.dims[1:2] && w.k < w.dims[3]
    sl = convert(Matrix{T}, s)
    for v in sl
        w.counts[Int(v) + 1] += 1
    end
    write(w.io, sl); w.k += 1
    return w
end
Base.close(w::VolumeWriter) = (w.k == w.dims[3] || @warn "volume closed after $(w.k) of $(w.dims[3]) slices"; close(w.io))

_cell(x) = x === nothing || (x isa AbstractFloat && isnan(x)) ? "" : x isa AbstractFloat ? @sprintf("%.6g", x) :
    x isa AbstractVector ? join(_cell.(x), " ") : (t = string(x); occursin(r"[,\"\n]", t) ? "\"" * replace(t, "\"" => "\"\"") * "\"" : t)

"""
    write_phantom_csv(path; facts, columns, times, states = ("rest", "stress"), curve, energy_keV = 70)

Writes `label_<scenario>.csv`, the one table of a phantom design (see the file header for the layout).

- `facts`: `key => value` pairs, one row each (scenario, lesion, volume file, dims, voxel size, origin, scan times, ...).
- `columns`: one NamedTuple per column with `key` (the header: a compartment id, `"region:<name>"` or `"roi:<name>"`)
  and any of `name, base_material, kinetics, voxels, mass_g, mbf_rest, mbf_stress, cfr, ffr, notes`. A
  `base_material` that names an `XA.Materials` entry is written out in full (density, Z/A, I, mass fractions, base HU).
- `curve(j, state, f)`: the value of column `j` at `times[f]` in `state` (iodine mg/mL; for an ROI, ΔHU), or `nothing`.
"""
function write_phantom_csv(path::AbstractString; facts, columns, times, states = ("rest", "stress"), curve, energy_keV = 70)
    E = energy_keV * u"keV"; μ(m) = ustrip(XA.linear_attenuation_coeff(m, E)); μw = μ(XA.Materials.water); μa = μ(XA.Materials.air)
    mat(n) = (isempty(n) || !hasproperty(XA.Materials, Symbol(n))) ? nothing : getproperty(XA.Materials, Symbol(n))
    fields = (
        ("name", c -> get(c, :name, "")), ("base_material", c -> get(c, :base_material, "")),
        ("density_g_cm3", c -> (m = mat(get(c, :base_material, "")); m === nothing ? nothing : ustrip(u"g/cm^3", m.density))),
        ("ZA_ratio", c -> (m = mat(get(c, :base_material, "")); m === nothing ? nothing : m.ZA_ratio)),
        ("I_eV", c -> (m = mat(get(c, :base_material, "")); m === nothing ? nothing : ustrip(u"eV", m.I))),
        ("composition_Z_massfrac", c -> (m = mat(get(c, :base_material, "")); m === nothing ? nothing : join(["$(Z):$(round(w; sigdigits = 6))" for (Z, w) in sort(collect(m.composition))], ";"))),
        ("base_hu$(round(Int, energy_keV))", c -> (m = mat(get(c, :base_material, "")); m === nothing ? nothing : 1000 * (μ(m) - μw) / (μw - μa))),
        ("kinetics", c -> get(c, :kinetics, "")), ("voxels", c -> get(c, :voxels, nothing)), ("mass_g", c -> get(c, :mass_g, nothing)),
        ("mbf_rest_ml_min_g", c -> get(c, :mbf_rest, nothing)), ("mbf_stress_ml_min_g", c -> get(c, :mbf_stress, nothing)),
        ("cfr", c -> get(c, :cfr, nothing)), ("ffr", c -> get(c, :ffr, nothing)), ("notes", c -> get(c, :notes, "")),
    )
    open(path, "w") do io
        println(io, join(("row", "value", (string(c.key) for c in columns)...), ","))
        blank = repeat(",", length(columns))
        for (k, v) in facts
            println(io, _cell(k), ",", _cell(v), blank)
        end
        for (n, f) in fields
            println(io, n, ",,", join(_cell.(f.(columns)), ","))
        end     # empty `value`: one cell per column
        for st in states, (fi, t) in enumerate(times)
            println(io, st, ",", _cell(Float64(t)), ",", join((_cell(curve(j, st, fi)) for j in eachindex(columns)), ","))
        end
    end
    return path
end

"""
    scan_times_from_curve(times, dhu; threshold_hu = 80, lead_s = 1) -> (; arrival, peak, prepeak, peak_hu)

Bolus-tracking rule on a ROI curve in ΔHU over baseline: `arrival` is the first sample at or above `threshold_hu`,
`peak` the sample of the maximum and `prepeak` `lead_s` before it (the usual "scan just before the peak").
"""
function scan_times_from_curve(times, dhu; threshold_hu = 80, lead_s = 1)
    k = findfirst(>=(threshold_hu), dhu); k === nothing && error("the ROI never reaches +$threshold_hu HU")
    kp = argmax(dhu)
    return (; arrival = times[k], peak = times[kp], prepeak = times[kp] - lead_s, peak_hu = dhu[kp])
end

# =============================================================================================== reading

struct Dyn
    csv::String
    facts::Dict{String, String}
    labels::Array{UInt16, 3}
    table::Dict{Int, NamedTuple{(:name, :base, :kinetics), Tuple{String, String, String}}}
    materials::Dict{String, XA.Material}
    times::Vector{Float64}
    curves::Dict{Tuple{String, Int}, Vector{Float64}}
    reference::Dict{Tuple{String, String}, Vector{Float64}}
    counts::Dict{Int, Int}
end

function read_volume(path, dims, ::Type{T}) where {T}
    io = endswith(path, ".zst") ? open(`zstd -dc $path`) : open(path)
    try
        # one z slice at a time: a single read of several GB from a pipe fails on macOS (EINVAL)
        A = Array{T}(undef, dims...); nxy = prod(dims[1:(end - 1)]); buf = Vector{T}(undef, nxy)
        for k in 1:dims[end]
            read!(io, buf); copyto!(A, (k - 1) * nxy + 1, buf, 1, nxy)
        end
        return A
    finally
        close(io)
    end
end

"""
    load(csv; volume = true) -> Dyn

Reads one phantom design: `label_<scenario>.csv` and (with `volume = true`) the compartment volume it names, from the
same folder. `ph.facts` holds the scenario facts, `ph.table` each compartment's name / base material / kinetics,
`ph.curves[(state, id)]` its iodine curve, `ph.reference[(state, "roi:…")]` the tracking ROI.
"""
function load(csv::AbstractString; volume::Bool = true)
    A = readdlm(csv, ','; header = false, quotes = true)
    hdr = string.(A[1, :]); rows = string.(A[:, 1])
    facts = Dict{String, String}(); times = Float64[]; tab = Dict{String, Vector{Any}}(); trows = Dict{String, Vector{Int}}()
    for r in 2:size(A, 1)
        key = rows[r]
        if key in ("rest", "stress")
            push!(get!(trows, key, Int[]), r)
        elseif isempty(string(A[r, 2]))             # per-column rows leave `value` empty (even when every cell is blank)
            tab[key] = A[r, :]
        else
            facts[key] = string(A[r, 2])
        end
    end
    st1 = first(keys(trows)); times = Float64.(A[trows[st1], 2])
    ids = Dict(c => parse(Int, hdr[c]) for c in 3:length(hdr) if all(isdigit, hdr[c]))
    nm(c, f) = string(tab[f][c])
    table = Dict(i => (name = nm(c, "name"), base = nm(c, "base_material"), kinetics = nm(c, "kinetics")) for (c, i) in ids)
    mats = Dict{String, XA.Material}()
    for (c, i) in ids
        b = nm(c, "base_material"); (isempty(b) || haskey(mats, b)) && continue
        comp = Dict{Int, Float64}(); for kv in split(nm(c, "composition_Z_massfrac"), ';')
            Z, w = split(kv, ':'); comp[parse(Int, Z)] = parse(Float64, w)
        end
        mats[b] = XA.Material(b, Float64(tab["ZA_ratio"][c]), Float64(tab["I_eV"][c]) * u"eV", Float64(tab["density_g_cm3"][c]) * u"g/cm^3", comp)
    end
    num(x) = x isa Number ? Float64(x) : (isempty(string(x)) ? 0.0 : parse(Float64, string(x)))
    curves = Dict{Tuple{String, Int}, Vector{Float64}}(); ref = Dict{Tuple{String, String}, Vector{Float64}}()
    for (st, rr) in trows
        for (c, i) in ids
            table[i].kinetics == "static" || (curves[(st, i)] = num.(A[rr, c]))
        end
        for c in 3:length(hdr)
            startswith(hdr[c], "roi:") && (ref[(st, hdr[c])] = num.(A[rr, c]))
        end
    end
    lab = Array{UInt16, 3}(undef, 0, 0, 0); cnt = Dict{Int, Int}()
    if volume
        dims = Tuple(parse.(Int, split(facts["dims"])))
        lab = UInt16.(read_volume(joinpath(dirname(abspath(csv)), facts["compartment_volume"]), dims, UInt16))
        cv = zeros(Int, 65536); for v in lab
            cv[v + 1] += 1
        end
        cnt = Dict(i - 1 => n for (i, n) in enumerate(cv) if n > 0)
        for l in keys(cnt)
            haskey(table, l) || error("compartment $l is in the volume but not in $(basename(csv))")
        end
    end
    return Dyn(csv, facts, lab, table, mats, times, curves, ref, cnt)
end

_at(times, y, t) = t <= times[1] ? y[1] : t >= times[end] ? y[end] : (k = searchsortedlast(times, t); w = (t - times[k]) / (times[k + 1] - times[k]); y[k] + w * (y[k + 1] - y[k]))
"""iodine(ph, state, t) -> Dict(label => mg/mL) for every label (0 for static labels); linear between samples."""
iodine(ph::Dyn, state, t) = Dict(l => (haskey(ph.curves, (state, l)) ? _at(ph.times, ph.curves[(state, l)], t) : 0.0) for l in keys(ph.table))

"""
    scan_times(ph; state, roi = first ROI column, threshold_hu = 80, lead_s = 1)

Bolus-tracking times of a loaded phantom from its ROI column (ΔHU at 70 keV): `arrival` (first sample ≥ +`threshold_hu`)
and `peak` (`lead_s` before the ROI maximum, the usual scan time), on the 1 s grid.
"""
function scan_times(ph::Dyn; state, roi = first(k for k in keys(ph.reference) if k[1] == state)[2], threshold_hu = 80.0, lead_s = 1.0)
    s = scan_times_from_curve(ph.times, ph.reference[(state, roi)]; threshold_hu, lead_s)
    return (arrival = s.arrival, peak = s.prepeak, peak_time = s.peak, peak_hu = s.peak_hu)
end

"""
    iodinated(base, c) -> XA.Material

`base` with `c` mg/mL iodine added without displacement (the format's mixing rule).
"""
function iodinated(base::XA.Material, c::Real)
    c <= 0 && return base
    rb = ustrip(u"g/cm^3", base.density); ci = Float64(c) * 1.0e-3; rho = rb + ci; I53 = XA.Elements.Iodine
    comp = Dict{Int, Float64}(Z => rb * w / rho for (Z, w) in base.composition); comp[53] = get(comp, 53, 0.0) + ci / rho
    za = (rb * base.ZA_ratio + ci * I53.ZA_ratio) / rho
    lnI = (rb * base.ZA_ratio * log(ustrip(u"eV", base.I)) + ci * I53.ZA_ratio * log(ustrip(u"eV", I53.I))) / (rho * za)
    return XA.Material(string(base.name, " + I ", round(c; digits = 4), " mg/mL"), za, exp(lnI) * u"eV", rho * u"g/cm^3", comp)
end

function _assign(ph, conc, tol)
    groups = Dict{String, Vector{Int}}(); for l in keys(ph.counts)
        push!(get!(groups, ph.table[l].base, Int[]), l)
    end
    lut = Dict{Int, Int}(); mats = XA.Material[]; err = Dict{Int, Float64}()
    for base in sort(collect(keys(groups)))
        ls = groups[base]; v = [conc[l] for l in ls]; w = [ph.counts[l] for l in ls]; o = sortperm(v); start = -Inf; bins = Vector{Vector{Int}}()
        for i in o
            v[i] - start > tol && (push!(bins, Int[]); start = v[i]); push!(bins[end], i)
        end
        for b in bins
            cbar = sum(v[b] .* w[b]) / sum(w[b]); push!(mats, iodinated(ph.materials[base], cbar))
            for i in b
                lut[ls[i]] = length(mats); err[ls[i]] = v[i] - cbar
            end
        end
    end
    return lut, mats, err
end

"""
    materials_at(ph, state, t; max_materials = 64, tol0 = 1e-4) -> (; lut, materials, tol, err)

Labels that share a base material and whose iodine at time `t` differs by ≤ `tol` mg/mL become one material; `tol`
grows until the table fits `max_materials` (64 keeps BasisSimulator's `:dd_fast` projector). `err` is each label's iodine
error in mg/mL.
"""
function materials_at(ph::Dyn, state, t; max_materials = 64, tol0 = 1.0e-4)
    conc = iodine(ph, state, t); tol = tol0
    nbase = length(unique(ph.table[l].base for l in keys(ph.counts)))
    nbase > max_materials && error("$nbase base materials in the volume, more than max_materials = $max_materials")
    while true
        lut, mats, err = _assign(ph, conc, tol); length(mats) <= max_materials && return (; lut, materials = mats, tol, err); tol *= 1.25
    end
end

"""
    bs_phantom(ph, BS; state, t, max_materials = 64) -> (BS.Phantom, materials_at result)

The BasisSimulator.jl phantom at time `t`: a 0-based UInt8 material mask plus the material table, voxel size in cm. Pass
the BasisSimulator module as `BS` (this file does not depend on it).
"""
function bs_phantom(ph::Dyn, BS; state, t, max_materials = 64)
    m = materials_at(ph, state, t; max_materials); lutv = zeros(UInt8, 65536)
    for (l, k) in m.lut
        lutv[l + 1] = UInt8(k - 1)
    end
    mask = similar(ph.labels, UInt8); @inbounds for i in eachindex(mask)
        mask[i] = lutv[ph.labels[i] + 1]
    end
    return BS.Phantom(mask, Dict{Int, XA.Material}(k - 1 => x for (k, x) in enumerate(m.materials)), Tuple(fill(parse(Float64, ph.facts["voxel_mm"]) / 10, 3))), m
end

end # module
