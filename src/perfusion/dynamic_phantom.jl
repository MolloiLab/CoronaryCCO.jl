# SPDX-License-Identifier: GPL-3.0-or-later
# The shareable dynamic CT perfusion phantom: one folder per scenario (healthy, a stenosis, ...), every file on one grid.
#
#   labels.u16.zst      every voxel's label (UInt16, column-major: x fastest, then y, then z; zstd)
#   rois.u8.zst         every voxel's ROI (UInt8; 0 = none): the compartments ground truth is given for, e.g. LV · LAD
#   mbf_<state>.f32.zst every voxel's ground-truth blood flow (mL/min/g; 0 where the voxel has no perfusion of its own)
#   label_<state>.csv   one row per label: name, ROI, kind, base material, voxels, mass, ground-truth MBF, then its iodine
#                       (mg per mL of voxel) at each time
#   ledger_<state>.csv  the organ's iodine at each time: delivered, held (arteries, beds, veins), outflow, painted
#   phantom.toml        the grid, the files, the states, the ROI names and groups, the tracking ROIs, the scenario's facts
#
# Every label lies in one ROI, so an ROI is a set of whole labels: its iodine and ground truth are sums over them. A
# label's iodine is its members' mass over their painted volume, so painting conserves mass. The writer half is
# organ-generic and holds no anatomy; the phantom-specific step (which voxel is what) is supplied by the caller.
# Mixing rule: iodine is added to a voxel without displacing tissue, ρ = ρ_base + c/1000 g/cm³, so
# μ(E) = μ_base(E) + (c/1000)(μ/ρ)_I(E) exactly; c is mg iodine per mL of voxel.
module DynamicPhantom

using DelimitedFiles, Unitful, Printf, Random, Statistics
import TOML
using CodecZstd: ZstdCompressorStream, ZstdDecompressorStream
import XrayAttenuation as XA
import AcceleratedKernels as AK

export compartment_classes, class_error, VolumeWriter, write_slice!, write_bundle, load, iodine, hu, scan_times, materials_at,
    bs_phantom, iodinated, hu_per_mg_ml, read_volume

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
    io = endswith(path, ".zst") ? ZstdCompressorStream(open(path, "w"); level) : open(path, "w")
    return VolumeWriter{T}(io, Tuple(dims), 0, zeros(Int, T <: Integer ? Int(typemax(T)) + 1 : 0))
end
"""
    write_slice!(w::VolumeWriter, slice)

Appends one z slice to the volume `w` streams.
"""
function write_slice!(w::VolumeWriter{T}, s::AbstractMatrix) where {T}
    @assert size(s) == w.dims[1:2] && w.k < w.dims[3]
    sl = convert(Matrix{T}, s)
    if T <: Integer
        for v in sl
            w.counts[Int(v) + 1] += 1
        end
    end
    write(w.io, sl); w.k += 1
    return w
end
Base.close(w::VolumeWriter) = (w.k == w.dims[3] || @warn "volume closed after $(w.k) of $(w.dims[3]) slices"; close(w.io))

_cell(x) = x === nothing || (x isa AbstractFloat && isnan(x)) ? "" : x isa AbstractFloat ? @sprintf("%.6g", x) :
    (t = string(x); occursin(r"[,\"\n]", t) ? "\"" * replace(t, "\"" => "\"\"") * "\"" : t)

function _material_cells(name, energy_keV)
    hasproperty(XA.Materials, Symbol(name)) || return (nothing, nothing, nothing, "", nothing)
    m = getproperty(XA.Materials, Symbol(name)); E = energy_keV * u"keV"; μ(x) = ustrip(XA.linear_attenuation_coeff(x, E))
    μw = μ(XA.Materials.water); μa = μ(XA.Materials.air)
    comp = join(["$(Z):$(round(w; sigdigits = 6))" for (Z, w) in sort(collect(m.composition))], ";")
    return (ustrip(u"g/cm^3", m.density), m.ZA_ratio, ustrip(u"eV", m.I), comp, 1000 * (μ(m) - μw) / (μw - μa))
end

"""
    write_bundle(out; grid, files, states, times, labels, ledger, toml, energy_keV = 70)

Writes the tables of a phantom folder (the volumes are written by the caller with `VolumeWriter`):

- `labels`: one NamedTuple per label, with `id, name, roi, kind, base_material, voxels, mass_g`, `mbf` (one value per
  state, or `nothing`) and `curve` (state index → vector of iodine, mg/mL, over `times`; `nothing` for a static label).
- `ledger`: state index → NamedTuple of equal-length columns (the first one `t_s`).
- `toml`: a Dict written to `phantom.toml`, with `grid`, `files`, `states` and `energy_keV` added.
"""
function write_bundle(out::AbstractString; grid, files, states, times, labels, ledger, toml, energy_keV = 70)
    tcols = [@sprintf("%g", t) for t in times]
    for (s, st) in enumerate(states)
        open(joinpath(out, "label_$st.csv"), "w") do io
            println(
                io, join(
                    [
                        "id", "name", "roi", "kind", "base_material", "density_g_cm3", "ZA_ratio", "I_eV", "composition_Z_massfrac",
                        "base_hu$(round(Int, energy_keV))", "voxels", "mass_g", "mbf_ml_min_g", tcols...,
                    ], ","
                )
            )
            for l in labels
                mc = _material_cells(l.base_material, energy_keV); c = l.curve === nothing ? fill(0.0, length(times)) : l.curve[s]
                println(
                    io, join(
                        _cell.(
                            [
                                l.id, l.name, l.roi, l.kind, l.base_material, mc..., l.voxels, l.mass_g,
                                l.mbf === nothing ? nothing : l.mbf[s], c...,
                            ]
                        ), ","
                    )
                )
            end
        end
        L = ledger[s]
        open(joinpath(out, "ledger_$st.csv"), "w") do io
            println(io, join(string.(keys(L)), ","))
            for k in eachindex(first(L))
                println(io, join((_cell(Float64(v[k])) for v in values(L)), ","))
            end
        end
    end
    meta = merge(
        toml, Dict{String, Any}(
            "grid" => Dict("dims" => collect(grid.dims), "origin_mm" => collect(Float64.(grid.origin)), "voxel_mm" => Float64(grid.h)),
            "files" => files, "states" => collect(states), "energy_keV" => Float64(energy_keV),
        )
    )
    open(io -> TOML.print(io, meta; sorted = true), joinpath(out, "phantom.toml"), "w")
    return out
end

# =============================================================================================== reading

struct Dyn
    dir::String
    meta::Dict{String, Any}                                      # phantom.toml
    labels::Array{UInt16, 3}
    rois::Array{UInt8, 3}
    origin_mm::NTuple{3, Float64}                               # of the loaded slab
    voxel_mm::Float64
    table::Dict{Int, NamedTuple{(:name, :base, :kind, :roi), Tuple{String, String, String, Int}}}
    materials::Dict{String, XA.Material}
    times::Vector{Float64}
    curves::Dict{Tuple{String, Int}, Vector{Float64}}          # (state, label) -> iodine, mg/mL
    mbf::Dict{Tuple{String, Int}, Float64}                     # (state, label) -> ground-truth MBF, mL/min/g
    label_mass::Dict{Int, Float64}                             # label -> g over the whole volume
    counts::Dict{Int, Int}                                     # label -> voxels in the loaded slab
    roi_table::Dict{Tuple{String, String}, NamedTuple}         # (state, ROI name) -> (; id, kind, voxels, mass_g, mbf, mass)
    ledger::Dict{String, Dict{String, Vector{Float64}}}        # state -> column -> values
    reference::Dict{Tuple{String, String}, Vector{Float64}}    # (state, "roi:<tracking name>") -> ΔHU of the tracking point
end

"""
    read_volume(path, dims, T; zrange = 1:dims[3], xrange = 1:dims[1], yrange = 1:dims[2]) -> Array{T, 3}

Reads a volume of the bundle (zstd or raw), only the box `xrange` × `yrange` × `zrange`.
"""
function read_volume(path, dims, ::Type{T}; zrange = 1:dims[3], xrange = 1:dims[1], yrange = 1:dims[2]) where {T}
    io = endswith(path, ".zst") ? ZstdDecompressorStream(open(path)) : open(path)
    try
        # one z slice at a time (a single read of several GB from a pipe fails on macOS, EINVAL), each cropped to (xrange,
        # yrange) as it is read, so only the crop is ever held
        slice = Matrix{T}(undef, dims[1], dims[2]); A = Array{T}(undef, length(xrange), length(yrange), length(zrange))
        for k in 1:last(zrange)
            read!(io, slice); k >= first(zrange) && (A[:, :, k - first(zrange) + 1] .= view(slice, xrange, yrange))
        end
        return A
    finally
        close(io)
    end
end

_num(x) = x isa Number ? Float64(x) : (isempty(strip(string(x))) ? NaN : parse(Float64, string(x)))
_rows(path) = (A = readdlm(path, ','; header = false, quotes = true); (string.(A[1, :]), A[2:end, :]))

"""
    load(dir; volume = true, rois = false, zrange = nothing, xrange = nothing, yrange = nothing) -> Dyn

Reads one scenario folder. `volume` reads `labels.u16.zst`, `rois` also `rois.u8.zst`; `zrange` (slice numbers, 1 = the
first) reads only those z slices, and `xrange`, `yrange` (voxel numbers) only that box of each (`origin_mm` is then the
box's). `ph.curves[(state, label)]` is a label's iodine (mg/mL),
`ph.mbf[(state, label)]` its ground-truth MBF, `ph.table[label]` its name, base material, kind and ROI,
`ph.roi_table[(state, name)]` a compartment's (or group's) voxels, mass, ground-truth MBF and iodine mass curve (sums over
its labels), `ph.ledger[state]` the organ's iodine ledger,
`ph.reference[(state, "roi:<name>")]` a tracking ROI's ΔHU curve (`scan_times` reads it), `ph.meta` the `phantom.toml`.
"""
function load(
        dir::AbstractString; volume::Bool = true, rois::Bool = false, zrange::Union{Nothing, AbstractUnitRange{<:Integer}} = nothing,
        xrange::Union{Nothing, AbstractUnitRange{<:Integer}} = nothing, yrange::Union{Nothing, AbstractUnitRange{<:Integer}} = nothing
    )
    meta = TOML.parsefile(joinpath(dir, "phantom.toml")); states = String.(meta["states"]); g = meta["grid"]
    dims = Tuple(Int.(g["dims"])); h = Float64(g["voxel_mm"]); o = Float64.(g["origin_mm"])
    table = Dict{Int, NamedTuple{(:name, :base, :kind, :roi), Tuple{String, String, String, Int}}}()
    mats = Dict{String, XA.Material}(); curves = Dict{Tuple{String, Int}, Vector{Float64}}(); mbf = Dict{Tuple{String, Int}, Float64}()
    lmass = Dict{Int, Float64}(); lvox = Dict{Int, Int}(); times = Float64[]
    for st in states
        hdr, A = _rows(joinpath(dir, "label_$st.csv")); col = Dict(h => i for (i, h) in enumerate(hdr))
        it = (col["mbf_ml_min_g"] + 1):length(hdr); times = parse.(Float64, hdr[it])
        for r in axes(A, 1)
            id = Int(A[r, col["id"]]); b = string(A[r, col["base_material"]]); kind = string(A[r, col["kind"]])
            table[id] = (name = string(A[r, col["name"]]), base = b, kind, roi = Int(_num(A[r, col["roi"]])))
            lmass[id] = _num(A[r, col["mass_g"]]); lvox[id] = Int(_num(A[r, col["voxels"]]))
            kind == "static" || (curves[(st, id)] = Float64.(_num.(A[r, it])))
            m = _num(A[r, col["mbf_ml_min_g"]]); isnan(m) || (mbf[(st, id)] = m)
            if !isempty(b) && !haskey(mats, b)
                comp = Dict{Int, Float64}(); for kv in split(string(A[r, col["composition_Z_massfrac"]]), ';')
                    Z, w = split(kv, ':'); comp[parse(Int, Z)] = parse(Float64, w)
                end
                mats[b] = XA.Material(b, _num(A[r, col["ZA_ratio"]]), _num(A[r, col["I_eV"]]) * u"eV", _num(A[r, col["density_g_cm3"]]) * u"g/cm^3", comp)
            end
        end
    end
    led = Dict{String, Dict{String, Vector{Float64}}}()
    for st in states
        hdr, A = _rows(joinpath(dir, "ledger_$st.csv")); led[st] = Dict(hdr[c] => Float64.(_num.(A[:, c])) for c in eachindex(hdr))
    end
    # each compartment and group from its labels: voxels, mass, ground-truth MBF (beds' flow over the mass of everything
    # in it) and iodine mass (mg) at every time
    vml = h^3 / 1000; names = String.(get(get(meta, "rois", Dict()), "names", String[]))
    groups = [String(k) => Int.(v) for (k, v) in get(get(meta, "rois", Dict()), "groups", Dict())]
    rt = Dict{Tuple{String, String}, NamedTuple}()
    for st in states, (name, ids, kind) in [[(names[r], [r], "roi") for r in eachindex(names)]; [(g, r, "group") for (g, r) in groups]]
        ls = [l for (l, row) in table if row.roi in ids]
        mass = sum(lmass[l] for l in ls; init = 0.0); flow = sum(get(mbf, (st, l), 0.0) * lmass[l] for l in ls; init = 0.0)
        io = sum((lvox[l] * vml) .* get(curves, (st, l), zeros(length(times))) for l in ls; init = zeros(length(times)))
        rt[(st, name)] = (
            id = kind == "roi" ? ids[1] : name, kind, voxels = sum(lvox[l] for l in ls; init = 0), mass_g = mass,
            mbf = mass > 0 ? flow / mass : NaN, mass = io,
        )
    end
    ref = Dict{Tuple{String, String}, Vector{Float64}}()
    for (name, tr) in get(meta, "tracking", Dict())
        l = Int(tr["label"]); hpm = hu_per_mg_ml(table[l].base, meta["energy_keV"])
        for st in states
            ref[(st, "roi:$name")] = hpm .* get(curves, (st, l), zeros(length(times)))
        end
    end
    xr, yr, zr = (r === nothing ? (1:dims[d]) : r for (d, r) in enumerate((xrange, yrange, zrange)))
    for (name, r, n) in (("xrange", xr, dims[1]), ("yrange", yr, dims[2]), ("zrange", zr, dims[3]))
        1 <= first(r) <= last(r) <= n || throw(ArgumentError("$name $r is outside the volume's 1:$n"))
    end
    org = (o[1] + (first(xr) - 1) * h, o[2] + (first(yr) - 1) * h, o[3] + (first(zr) - 1) * h)
    lab = Array{UInt16, 3}(undef, 0, 0, 0); roi = Array{UInt8, 3}(undef, 0, 0, 0); cnt = Dict{Int, Int}()
    if volume
        lab = read_volume(joinpath(dir, meta["files"]["labels"]), dims, UInt16; zrange = zr, xrange = xr, yrange = yr)
        cnt = _label_counts(lab)
        for l in keys(cnt)
            haskey(table, l) || error("label $l is in the volume but not in label_$(first(states)).csv")
        end
    end
    rois && (roi = read_volume(joinpath(dir, meta["files"]["rois"]), dims, UInt8; zrange = zr, xrange = xr, yrange = yr))
    return Dyn(abspath(dir), meta, lab, roi, org, h, table, mats, times, curves, mbf, lmass, cnt, rt, led, ref)
end

# a function barrier: counted inline in `load`, the loop ran ~30x slower
function _label_counts(lab)
    cv = zeros(Int, 65536)
    for v in lab
        cv[v + 1] += 1
    end
    return Dict(i - 1 => n for (i, n) in enumerate(cv) if n > 0)
end

_at(times, y, t) = t <= times[1] ? y[1] : t >= times[end] ? y[end] : (k = searchsortedlast(times, t); w = (t - times[k]) / (times[k + 1] - times[k]); y[k] + w * (y[k + 1] - y[k]))
"""iodine(ph, state, t) -> Dict(label => mg/mL) for every label (0 for static labels); linear between samples."""
iodine(ph::Dyn, state, t) = Dict(l => (haskey(ph.curves, (state, l)) ? _at(ph.times, ph.curves[(state, l)], t) : 0.0) for l in keys(ph.table))

"""
    hu(ph, state, t; energy_keV = the phantom's) -> Dict(label => HU)

Every label's CT number at time `t` (s) and `energy_keV`: its base material with its iodine added by the format's mixing
rule. A lookup of the label volume through it is the idealized (noise- and blur-free) image of the phantom.
"""
function hu(ph::Dyn, state, t; energy_keV = ph.meta["energy_keV"])
    E = energy_keV * u"keV"; μ(m) = ustrip(XA.linear_attenuation_coeff(m, E)); μw = μ(XA.Materials.water); μa = μ(XA.Materials.air)
    c = iodine(ph, state, t)
    return Dict(l => 1000 * (μ(iodinated(ph.materials[r.base], c[l])) - μw) / (μw - μa) for (l, r) in ph.table)
end

"""
    scan_times(ph; state, roi = its only tracking ROI, threshold_hu = 80, lead_s = 1, step = 0.5)

Bolus-tracking times of a loaded phantom from its tracking ROI (ΔHU), sampled every `step` s (the tracking scan's
cadence) on the curve interpolated between the stored frames: `arrival` (first sample ≥ +`threshold_hu`), `peak`
(`lead_s` before the ROI maximum, the usual scan time), `peak_time` (the maximum) and `peak_hu`. `roi` is a reference
key, `"roi:<tracking name>"` (e.g. `"roi:descending aorta"`); a phantom with several tracking points must name one.
"""
function scan_times(ph::Dyn; state, roi = nothing, threshold_hu = 80.0, lead_s = 1.0, step = 0.5)
    if roi === nothing
        rois = sort!([k[2] for k in keys(ph.reference) if k[1] == state])
        length(rois) == 1 || throw(ArgumentError("scan_times: the phantom has tracking ROIs $rois; name one with `roi`"))
        roi = only(rois)
    end
    dhu = ph.reference[(state, roi)]; ts = ph.times[1]:step:ph.times[end]; y = [_at(ph.times, dhu, t) for t in ts]
    k = findfirst(>=(threshold_hu), y); k === nothing && error("the ROI never reaches +$threshold_hu HU")
    kp = argmax(y)
    return (arrival = ts[k], peak = ts[kp] - lead_s, peak_time = ts[kp], peak_hu = y[kp])
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

"""
    hu_per_mg_ml(material, energy_keV) -> HU per mg/mL

The CT number change (HU) that 1 mg/mL of iodine adds to `XA.Materials` entry `material` at `energy_keV`.
"""
function hu_per_mg_ml(material::AbstractString, energy_keV::Real)
    E = energy_keV * u"keV"; μ(m) = ustrip(XA.linear_attenuation_coeff(m, E)); b = getproperty(XA.Materials, Symbol(material))
    return 1000 * (μ(iodinated(b, 1.0)) - μ(b)) / (μ(XA.Materials.water) - μ(XA.Materials.air))
end

# a label's merge group: its base material, and for a label that carries iodine also its ROI (so no material spans two
# ROIs and every ROI keeps its iodine mass exactly)
_group(ph, l) = (ph.table[l].base, ph.table[l].kind == "static" ? -2 : ph.table[l].roi)

# optimal 1-D weighted k-means of sorted values `v` (weights `w`) for k = 1..kmax: cost[k] (Σ w (v - bin mean)²) and the
# bin starts of each k, by dynamic programming over prefix sums
function _kmeans1d(v, w, kmax)
    n = length(v); kmax = min(kmax, n)
    W = cumsum([0.0; w]); S = cumsum([0.0; w .* v]); Q = cumsum([0.0; w .* v .^ 2])
    sse(i, j) = (ww = W[j + 1] - W[i]; ww <= 0 ? 0.0 : max(0.0, Q[j + 1] - Q[i] - (S[j + 1] - S[i])^2 / ww))
    D = fill(Inf, kmax, n); B = zeros(Int, kmax, n)
    for j in 1:n
        D[1, j] = sse(1, j); B[1, j] = 1
    end
    for k in 2:kmax, j in k:n, i in k:j
        c = D[k - 1, i - 1] + sse(i, j); c < D[k, j] && (D[k, j] = c; B[k, j] = i)
    end
    starts(k) = (
        st = Int[]; j = n; for kk in k:-1:1
            i = B[kk, j]; pushfirst!(st, i); j = i - 1
        end; st
    )
    return D[:, n], starts
end

# The fewest merges that fit `max_materials`: static labels of one base share a material (no iodine, nothing lost); every
# iodinated group starts with one material, and each further material goes to the group whose iodine error it lowers
# most. A merged material is its members' voxel-weighted mean, so mass is conserved in every group.
function _assign(ph, conc, max_materials)
    groups = Dict{Tuple{String, Int}, Vector{Int}}(); for l in keys(ph.counts)
        push!(get!(groups, _group(ph, l), Int[]), l)
    end
    keys_ = sort(collect(keys(groups)))
    length(keys_) > max_materials && error("$(length(keys_)) base material × ROI groups in the volume, more than max_materials = $max_materials: raise it (BasisSimulator then uses its slower :dd projector)")
    sol = Dict{Tuple{String, Int}, Any}(); nk = Dict(g => 1 for g in keys_)
    for g in keys_
        ls = groups[g]; v = [conc[l] for l in ls]; o = sortperm(v); ls = ls[o]; groups[g] = ls
        sol[g] = g[2] == -2 ? (cost = [0.0], starts = k -> [1]) : (r = _kmeans1d(v[o], Float64[ph.counts[l] for l in ls], max_materials); (cost = r[1], starts = r[2]))
    end
    for _ in 1:(max_materials - length(keys_))
        gain(g) = nk[g] < length(sol[g].cost) ? sol[g].cost[nk[g]] - sol[g].cost[nk[g] + 1] : -Inf
        g = argmax(gain, keys_); gain(g) > 0 || break
        nk[g] += 1
    end
    lut = Dict{Int, Int}(); mats = XA.Material[]; err = Dict{Int, Float64}()
    for g in keys_
        ls = groups[g]; st = [sol[g].starts(nk[g]); length(ls) + 1]
        for b in 1:(length(st) - 1)
            m = ls[st[b]:(st[b + 1] - 1)]; w = [ph.counts[l] for l in m]
            cbar = sum(conc[l] * ph.counts[l] for l in m) / sum(w); push!(mats, iodinated(ph.materials[g[1]], cbar))
            for l in m
                lut[l] = length(mats); err[l] = conc[l] - cbar
            end
        end
    end
    return lut, mats, err
end

"""
    materials_at(ph, state, t; max_materials = 64) -> (; lut, materials, err, max_err)

The material table of the phantom at time `t`, at most `max_materials` long (64 keeps BasisSimulator's `:dd_fast`
projector; more falls back to its slower `:dd`). Different base materials are never merged, and a label that carries
iodine is merged only with labels of its own base material and ROI: the fewest merges that fit, each merged material at
its members' voxel-weighted mean iodine, so the total iodine and every ROI's iodine are kept exactly. `err` is each
label's iodine error (mg/mL), `max_err` the largest.
"""
function materials_at(ph::Dyn, state, t; max_materials = 64)
    conc = iodine(ph, state, t)
    lut, mats, err = _assign(ph, conc, max_materials)
    return (; lut, materials = mats, err, max_err = maximum(abs, values(err); init = 0.0))
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
    return BS.Phantom(mask, Dict{Int, XA.Material}(k - 1 => x for (k, x) in enumerate(m.materials)), Tuple(fill(ph.voxel_mm / 10, 3))), m
end

end # module
