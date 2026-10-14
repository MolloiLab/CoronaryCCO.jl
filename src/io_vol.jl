# SPDX-License-Identifier: GPL-3.0-or-later
# Reader for the DGtal .vol format (Version 2 raw, Version 3 zlib-compressed), following
# DGtal::VolReader::importVol (DGtal 1.4, LGPL-3.0). Used for the myocardial domain mask.

import Zlib_jll

"""
    read_vol(path) -> (img::Array{UInt8,3}, lower::NTuple{3,Int32})

Voxel `img[i, j, k]` is the DGtal point `lower .+ (i, j, k) .- 1`. When the header has
`Center-*` fields, the domain is centred as in DGtal: `lower = c - (s - 1) ÷ 2`.
"""
function read_vol(path::AbstractString)
    bytes = read(path)
    header = Dict{String, String}()
    pos = 1
    while true
        nl = findnext(==(UInt8('\n')), bytes, pos)
        nl === nothing && error("read_vol: unterminated header in $path")
        line = String(bytes[pos:(nl - 1)])
        pos = nl + 1
        line == "." && break
        k = findfirst(':', line)
        k === nothing && error("read_vol: invalid header line '$line'")
        header[line[1:(k - 1)]] = strip(line[(k + 1):end])
    end
    sx, sy, sz = (parse(Int, header[k]) for k in ("X", "Y", "Z"))
    version = parse(Int, header["Version"])
    if haskey(header, "Center-X")
        c = (parse(Int, header["Center-X"]), parse(Int, header["Center-Y"]), parse(Int, header["Center-Z"]))
        lower = Int32.((c[1] - (sx - 1) ÷ 2, c[2] - (sy - 1) ÷ 2, c[3] - (sz - 1) ÷ 2))
    else
        lower = (Int32(0), Int32(0), Int32(0))
    end
    total = sx * sy * sz
    payload = @view bytes[pos:end]
    img = Array{UInt8, 3}(undef, sx, sy, sz)
    if version == 3
        destlen = Ref{Culong}(total)
        rc = ccall(
            (:uncompress, Zlib_jll.libz), Cint, (Ptr{UInt8}, Ref{Culong}, Ptr{UInt8}, Culong),
            img, destlen, payload, length(payload)
        )
        # Z_BUF_ERROR (-5) only means trailing data; DGtal reads exactly `total` voxels.
        (rc == 0 || rc == -5) && destlen[] == total || error("read_vol: zlib error $rc")
    elseif version == 2
        length(payload) >= total || error("read_vol: truncated raw data")
        copyto!(img, 1, payload, 1, total)
    else
        error("read_vol: unsupported Version $version")
    end
    return img, lower
end
