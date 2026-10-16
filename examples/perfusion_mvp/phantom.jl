# SPDX-License-Identifier: GPL-3.0-or-later
# Voxelize and materialize: a label volume plus per-voxel iodine concentration (mg I per mL of
# voxel) at chosen times, from a myocardial mask, a grown tree and the transport masses.
#
# Mass bookkeeping is exact: every milligram from `Flow.frame_masses` lands in exactly one voxel,
# so the phantom's iodine inside the mask equals the network's iodine inside the mask.

module Phantom

import AcceleratedKernels as AK
using Atomix: @atomic
using TOML
import XrayAttenuation as XA
using Unitful: @u_str, ustrip

export PhantomGrid, phantom_grid, lv_cavity, mark_lumen!, deposit!, write_frame, write_meta, materialize, write_materials

const F32 = Float32

"""
    PhantomGrid

Regular grid with `dims` voxels of size `h` mm; voxel (1,1,1) is centered at `origin` (mm).
`labels`: 0 air, 1 soft tissue, 2 myocardium, 3 blood (cavity and resolved vessel lumen).
`myo` marks the myocardial mask, including the vessels inside it. `owner[v]`
is the terminal whose territory contains myocardial voxel `v` (0 if none).
"""
struct PhantomGrid
    dims::NTuple{3, Int}
    h::F32
    origin::NTuple{3, F32}
    labels::Array{UInt8, 3}
    myo::BitArray{3}
    owner::Array{Int32, 3}
    count::Vector{Int32}        # myocardial voxels per terminal
end

@inline _vox(g::PhantomGrid, x, y, z) = (
    unsafe_trunc(Int, floor((x - g.origin[1]) / g.h + 0.5f0)) + 1,
    unsafe_trunc(Int, floor((y - g.origin[2]) / g.h + 0.5f0)) + 1,
    unsafe_trunc(Int, floor((z - g.origin[3]) / g.h + 0.5f0)) + 1,
)

"""
    lv_cavity(mask) -> BitArray

Blood pool enclosed by the myocardium. The case0 mask is open at the valve plane, so a 3D fill
leaks. A voxel is cavity if a 2D hole fill encloses it in at least two of the three axis-aligned
slice stacks.
"""
function lv_cavity(M::AbstractArray{Bool, 3})
    votes = zeros(UInt8, size(M))
    for ax in 1:3
        for k in 1:size(M, ax)
            S = collect(selectdim(M, ax, k))
            H = _holes2(S)
            V = selectdim(votes, ax, k)
            V .+= H
        end
    end
    return votes .>= 2
end

function _holes2(S::AbstractMatrix{Bool})
    sz = size(S); seen = falses(sz); st = CartesianIndex{2}[]
    for I in CartesianIndices(S)
        if !S[I] && (I[1] in (1, sz[1]) || I[2] in (1, sz[2]))
            seen[I] = true; push!(st, I)
        end
    end
    while !isempty(st)
        I = pop!(st)
        for J in (I + CartesianIndex(1, 0), I - CartesianIndex(1, 0), I + CartesianIndex(0, 1), I - CartesianIndex(0, 1))
            checkbounds(Bool, S, J) || continue
            (S[J] || seen[J]) && continue
            seen[J] = true; push!(st, J)
        end
    end
    return .!S .& .!seen
end

"""
    phantom_grid(mask, lower, tpts; h = 0.5, margin = 10, nearest) -> PhantomGrid

`mask` is the 1 mm myocardial mask whose voxel (1,1,1) sits at integer coordinate `lower` (mm).
The phantom covers the mask's bounding box plus `margin` mm at spacing `h`, sampled from the mask
by nearest neighbour. `tpts` (3 × nterm) are terminal points; `nearest(tpts, q)` returns the
nearest terminal of each query point (Territory.nearest_terminal).
"""
function phantom_grid(mask::AbstractArray{Bool, 3}, lower, tpts::Matrix{F32}; h = 0.5f0, margin = 10.0f0, nearest)
    h = F32(h)
    cav = lv_cavity(mask)
    lo = F32.(lower) .- F32(margin)
    hi = F32.(lower) .+ F32.(size(mask)) .- 1.0f0 .+ F32(margin)
    dims = ntuple(d -> Int(floor((hi[d] - lo[d]) / h)) + 1, 3)
    labels = fill(UInt8(1), dims)
    AK.foreachindex(labels; min_elems = 1 << 14) do v
        I = CartesianIndices(labels)[v]
        c = ntuple(d -> unsafe_trunc(Int, floor(lo[d] + (I[d] - 1) * h - lower[d] + 0.5f0)) + 1, 3)
        if all(d -> 1 <= c[d] <= size(mask, d), 1:3)
            mask[c...] && (labels[v] = 0x02)
            cav[c...] && (labels[v] = 0x03)
        end
    end
    myomask = labels .== 0x02
    myo = findall(vec(myomask))
    q = Matrix{F32}(undef, 3, length(myo))
    CI = CartesianIndices(labels)
    for (k, v) in enumerate(myo), d in 1:3
        q[d, k] = lo[d] + (CI[v][d] - 1) * h
    end
    own = nearest(tpts, q)
    owner = zeros(Int32, dims)
    count = zeros(Int32, size(tpts, 2))
    for (k, v) in enumerate(myo)
        owner[v] = own[k]; count[own[k]] += 1
    end
    return PhantomGrid(dims, h, lo, labels, myomask, owner, count)
end

# Local frame of segment j: start point, unit direction, length, and two unit normals.
@inline function _segframe(tree, j)
    ax = tree.a[1, j]; ay = tree.a[2, j]; az = tree.a[3, j]
    dx = tree.b[1, j] - ax; dy = tree.b[2, j] - ay; dz = tree.b[3, j] - az
    L = sqrt(dx * dx + dy * dy + dz * dz)
    L > 0.0f0 || return (ax, ay, az, 0.0f0, 0.0f0, 1.0f0, 0.0f0, 1.0f0, 0.0f0, 0.0f0, 0.0f0, 1.0f0, 0.0f0)
    dx /= L; dy /= L; dz /= L
    ux, uy, uz = abs(dz) < 0.9f0 ? (dy, -dx, 0.0f0) : (0.0f0, dz, -dy)
    nu = sqrt(ux * ux + uy * uy + uz * uz); ux /= nu; uy /= nu; uz /= nu
    wx = dy * uz - dz * uy; wy = dz * ux - dx * uz; wz = dx * uy - dy * ux
    return (ax, ay, az, dx, dy, dz, ux, uy, uz, wx, wy, wz, L)
end

"""
    mark_lumen!(g, tree, venous_scale)

Labels as blood (3) every voxel whose center lies inside a resolved lumen (radius ≥ h/2) of an
artery or its companion vein (offset by r_a + r_v, as in [`deposit!`](@ref)). `myo` is left
unchanged, so intramyocardial vessels stay inside the myocardial mask.
"""
function mark_lumen!(g::PhantomGrid, tree, venous_scale::Real)
    s = F32(venous_scale); labels = g.labels
    AK.foreachindex(tree.r; min_elems = 1024) do j
        r = tree.r[j]
        s * r >= 0.5f0 * g.h || return
        ax, ay, az, dx, dy, dz, ux, uy, uz, _, _, _, L = _segframe(tree, j)
        L > 0.0f0 || return
        for (rr, off) in ((r, 0.0f0), (s * r, r + s * r))
            rr >= 0.5f0 * g.h || continue
            ox = ax + off * ux; oy = ay + off * uy; oz = az + off * uz
            lo = _vox(g, min(ox, ox + L * dx) - rr, min(oy, oy + L * dy) - rr, min(oz, oz + L * dz) - rr)
            hi = _vox(g, max(ox, ox + L * dx) + rr, max(oy, oy + L * dy) + rr, max(oz, oz + L * dz) + rr)
            for l in max(1, lo[3]):min(g.dims[3], hi[3]), k in max(1, lo[2]):min(g.dims[2], hi[2]), i in max(1, lo[1]):min(g.dims[1], hi[1])
                px = g.origin[1] + (i - 1) * g.h - ox; py = g.origin[2] + (k - 1) * g.h - oy; pz = g.origin[3] + (l - 1) * g.h - oz
                t = clamp(px * dx + py * dy + pz * dz, 0.0f0, L)
                qx = px - t * dx; qy = py - t * dy; qz = pz - t * dz
                if qx * qx + qy * qy + qz * qz <= rr * rr && labels[i, k, l] != 0x00
                    labels[i, k, l] = 0x03
                end
            end
        end
    end
    return g
end

"""
    deposit!(conc, g, tree, art, ven, tis, f; venous_offset = true, blood_conc = 0)

Fills `conc` (mg/mL) for frame `f`:

- each terminal's lumped-bed mass is spread evenly over its territory voxels;
- each arterial segment's mass is spread over sample points filling its lumen (rings of
  points every ≤ h/2 along and across the vessel);
- each companion vein is placed alongside its artery, offset by r_a + r_v perpendicular to it;
- LV cavity voxels get `blood_conc`.
"""
function deposit!(
        conc::Array{F32, 3}, g::PhantomGrid, tree, art::AbstractMatrix, ven::AbstractMatrix,
        tis::AbstractMatrix, f::Integer, venous_scale::Real; blood_conc::Real = 0.0f0
    )
    fill!(conc, 0.0f0)
    vv = g.h^3 * 1.0f-3                     # mL per voxel
    # territories
    AK.foreachindex(conc; min_elems = 1 << 14) do v
        o = g.owner[v]
        if o > 0
            conc[v] = tis[f, o] / (g.count[o] * vv)
        elseif g.labels[v] == 0x03 && !g.myo[v]
            conc[v] = F32(blood_conc)
        end
    end
    s = F32(venous_scale)
    AK.foreachindex(tree.r; min_elems = 1024) do j
        ma = art[f, j]; mv = ven[f, j]
        (ma > 0.0f0 || mv > 0.0f0) || return
        ax, ay, az, dx, dy, dz, ux, uy, uz, wx, wy, wz, L = _segframe(tree, j)
        L > 0.0f0 || return
        r = tree.r[j]
        for (m, rr, off) in ((ma, r, 0.0f0), (mv, s * r, r + s * r))
            m > 0.0f0 || continue
            step = g.h * 0.5f0
            nl = max(1, unsafe_trunc(Int, ceil(L / step)))
            nr = rr > step ? unsafe_trunc(Int, ceil(rr / step)) : 0
            # total weight: per axial sample, the center (weight 1 if nr == 0) or rings ∝ radius
            wsum = 0.0f0
            for ir in 0:nr
                wsum += nr == 0 ? 1.0f0 : rr * (ir + 0.5f0) / (nr + 1)
            end
            for il in 1:nl
                t = L * (il - 0.5f0) / nl
                cx = ax + t * dx + off * ux; cy = ay + t * dy + off * uy; cz = az + t * dz + off * uz
                for ir in 0:nr
                    ρ = nr == 0 ? 0.0f0 : rr * (ir + 0.5f0) / (nr + 1)
                    na = nr == 0 ? 1 : max(1, unsafe_trunc(Int, ceil(2.0f0 * F32(π) * ρ / step)))
                    wr = (nr == 0 ? 1.0f0 : ρ) / wsum / nl / na
                    for ia in 1:na
                        θ = 2.0f0 * F32(π) * (ia - 1) / na
                        cθ = cos(θ) * ρ; sθ = sin(θ) * ρ
                        px = cx + cθ * ux + sθ * wx; py = cy + cθ * uy + sθ * wy; pz = cz + cθ * uz + sθ * wz
                        i, k, l = _vox(g, px, py, pz)
                        if 1 <= i <= g.dims[1] && 1 <= k <= g.dims[2] && 1 <= l <= g.dims[3]
                            @atomic conc[i, k, l] += m * wr / vv
                        end
                    end
                end
            end
        end
    end
    return conc
end

write_frame(path, conc::Array{F32, 3}) = (write(path, conc); path)

function write_meta(dir, g::PhantomGrid, times, files; extra = Dict{String, Any}())
    write(joinpath(dir, "labels.u8"), g.labels); write(joinpath(dir, "myo.u8"), UInt8.(g.myo))
    meta = Dict{String, Any}(
        "dims" => collect(g.dims), "voxel_mm" => fill(Float64(g.h), 3), "origin_mm" => collect(Float64.(g.origin)),
        "labels" => Dict("0" => "air", "1" => "soft_tissue", "2" => "myocardium", "3" => "blood"),
        "times_s" => collect(Float64.(times)), "files" => collect(files),
        "index_files" => [replace(f, "iodine_" => "index_", ".f32" => ".u16") for f in files],
        "material_files" => [replace(f, "iodine_" => "materials_", ".f32" => ".toml") for f in files],
        "mixing_rule" => "iodine added without displacement: rho = rho_base + c/1000; mass fractions reweighted (mu = mu_base + c/1000 (mu/rho)_I)",
        "myo_file" => "myo.u8",
    )
    merge!(meta, extra)
    open(io -> TOML.print(io, meta), joinpath(dir, "meta.toml"), "w")
    return dir
end

# ---------------------------------------------------------------------------------------------
# Materialize: every voxel becomes an XrayAttenuation.jl Material
# ---------------------------------------------------------------------------------------------

const BASE = Dict(
    0x00 => XA.Materials.air, 0x01 => XA.Materials.ncat_muscle,
    0x02 => XA.Materials.ncat_heart, 0x03 => XA.Materials.ncat_blood,
)

"""
    iodinated(base::XA.Material, c_mg_ml) -> XA.Material

`base` with `c` mg of iodine added per mL. Iodinated contrast adds mass to a fixed volume, so
ρ = ρ_base + c and the mass fractions are (ρ_base·w_base + c·δ_I)/ρ. The Z/A ratio is the
mass-weighted mean, and the mean excitation energy follows Bragg additivity,
ln I = Σ w (Z/A) ln I / Σ w (Z/A) (ICRU Report 37, 1984). Attenuation is then exactly
μ = ρ_base (μ/ρ)_base + c (μ/ρ)_I, the mixture rule XrayAttenuation applies.
"""
function iodinated(base::XA.Material, c::Real)
    c == 0 && return base
    rb = ustrip(u"g/cm^3", base.density); ci = Float64(c) * 1.0e-3
    rho = rb + ci
    I53 = XA.Elements.Iodine
    comp = Dict{Int, Float64}(Z => rb * w / rho for (Z, w) in base.composition)
    comp[53] = get(comp, 53, 0.0) + ci / rho
    za = (rb * base.ZA_ratio + ci * I53.ZA_ratio) / rho
    lnI = (rb * base.ZA_ratio * log(ustrip(u"eV", base.I)) + ci * I53.ZA_ratio * log(ustrip(u"eV", I53.I))) / (rho * za)
    name = string(base.name, " + I ", round(c; digits = 3), " mg/mL")
    return XA.Material(name, za, exp(lnI) * u"eV", rho * u"g/cm^3", comp)
end

@inline function _hash01(v::UInt32, seed::UInt32)
    x = v * 0x9e3779b1 + seed
    x = (x ⊻ (x >> 16)) * 0x85ebca6b
    x = (x ⊻ (x >> 13)) * 0xc2b2ae35
    x = x ⊻ (x >> 16)
    return Float32(x >> 8) * Float32(1 / 16777216)
end

"""
    materialize(labels, conc; levels, seed) -> (index::Array{UInt16,3}, materials::Vector{XA.Material}, widths)

Maps every voxel to one XrayAttenuation material: its base tissue (air, soft tissue, myocardium,
blood) plus iodine at one of `levels[label]` concentration levels 0, w, 2w, …, with w set by that
label's maximum concentration in the frame. Rounding is dithered (round up with probability equal
to the fractional part, from a fixed hash of the voxel index), so every voxel's expected iodine is
exact and the summed iodine over many voxels is preserved (Roberts 1962). The default budget is
33 materials per frame. BasisSimulator's fast projector allows at most 64, and its cost grows with
the material count. Dither noise (≤ one bin, ~0.1–1 mg/mL) is small against CT quantum noise.
"""
function materialize(
        labels::Array{UInt8, 3}, conc::Array{F32, 3};
        levels = Dict(0x00 => 1, 0x01 => 4, 0x02 => 12, 0x03 => 16), seed::Integer = 1, quantile = 0.999
    )
    labs = sort(collect(keys(levels)))
    # Tissue voxels whose iodine is above their label's `quantile` are mostly unresolved vessel
    # (a sub-voxel vessel's blood deposited in one voxel), so they take the blood base instead.
    # This keeps tissue bins fine without clipping any iodine.
    cut = Dict{UInt8, Float32}()
    for l in (0x01, 0x02)
        vals = conc[labels .== l]
        cut[l] = isempty(vals) ? 0.0f0 : sort(vals)[clamp(ceil(Int, quantile * length(vals)), 1, length(vals))]
    end
    labels = copy(labels)
    for v in eachindex(labels)
        l = labels[v]
        haskey(cut, l) && conc[v] > cut[l] && (labels[v] = 0x03)
    end
    cmax = Dict(l => 0.0f0 for l in labs)
    for v in eachindex(labels)
        l = labels[v]; cmax[l] = max(cmax[l], conc[v])
    end
    widths = Dict(l => (levels[l] > 1 && cmax[l] > 0) ? cmax[l] / (levels[l] - 1) : 0.0f0 for l in labs)
    first_ = Dict{UInt8, Int}(); materials = XA.Material[]
    for l in labs
        first_[l] = length(materials) + 1
        nl = widths[l] > 0 ? levels[l] : 1
        for k in 0:(nl - 1)
            push!(materials, iodinated(BASE[l], k * widths[l]))
        end
    end
    off = zeros(Int, 256); wid = zeros(F32, 256); top = zeros(Int, 256)
    for l in labs
        off[l + 1] = first_[l]; wid[l + 1] = widths[l]; top[l + 1] = (widths[l] > 0 ? levels[l] : 1) - 1
    end
    index = Array{UInt16, 3}(undef, size(labels))
    s32 = UInt32(seed)
    AK.foreachindex(index; min_elems = 1 << 14) do v
        l = labels[v] + 1; w = wid[l]
        k = 0
        if w > 0.0f0
            x = conc[v] / w
            k = unsafe_trunc(Int, floor(x + _hash01(UInt32(v), s32)))
            k = clamp(k, 0, top[l])
        end
        index[v] = UInt16(off[l] + k)
    end
    return index, materials, widths
end

"""
    write_materials(path, materials)

TOML table of materials that rebuild exactly with
`XA.Material(name, ZA_ratio, I_eV*u"eV", density_g_cm3*u"g/cm^3", Dict(Z => w))`.
"""
function write_materials(path, materials)
    tab = [
        Dict(
            "name" => m.name, "density_g_cm3" => ustrip(u"g/cm^3", m.density), "ZA_ratio" => m.ZA_ratio,
            "I_eV" => ustrip(u"eV", m.I), "composition" => Dict(string(Z) => w for (Z, w) in m.composition),
        ) for m in materials
    ]
    open(io -> TOML.print(io, Dict("materials" => tab)), path, "w")
    return path
end

end # module
