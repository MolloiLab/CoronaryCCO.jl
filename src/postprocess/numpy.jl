# SPDX-License-Identifier: GPL-3.0-or-later
# The NumPy and Python primitives coronary_postprocess.py relies on, reproduced so that the
# post-processing prints the same digits: pairwise summation (np.sum / np.mean / np.std),
# np.median, np.percentile (method "linear"), Python float formatting, csv.writer rows, and
# C++17 std::to_chars (VTK's number printing in the OBJ writer).

using Printf

# numpy's DOUBLE_pairwise_sum: plain loop below 8 terms, 8 accumulators up to 128 terms,
# otherwise split in halves rounded down to a multiple of 8.
function _np_pairwise(a::AbstractVector{Float64}, lo::Int, n::Int)
    if n < 8
        r = 0.0
        for i in 0:(n - 1)
            r += a[lo + i]
        end
        return r
    elseif n <= 128
        r1, r2, r3, r4 = a[lo], a[lo + 1], a[lo + 2], a[lo + 3]
        r5, r6, r7, r8 = a[lo + 4], a[lo + 5], a[lo + 6], a[lo + 7]
        i = 8
        while i < n - n % 8
            r1 += a[lo + i]; r2 += a[lo + i + 1]; r3 += a[lo + i + 2]; r4 += a[lo + i + 3]
            r5 += a[lo + i + 4]; r6 += a[lo + i + 5]; r7 += a[lo + i + 6]; r8 += a[lo + i + 7]
            i += 8
        end
        r = ((r1 + r2) + (r3 + r4)) + ((r5 + r6) + (r7 + r8))
        while i < n
            r += a[lo + i]
            i += 1
        end
        return r
    else
        n2 = n ÷ 2
        n2 -= n2 % 8
        return _np_pairwise(a, lo, n2) + _np_pairwise(a, lo + n2, n - n2)
    end
end

# np.add.reduce: the identity 0.0 plus the pairwise sum of the whole contiguous array.
_np_sum(a::AbstractVector{Float64}) = 0.0 + _np_pairwise(a, firstindex(a), length(a))
_np_mean(a::AbstractVector{Float64}) = _np_sum(a) / length(a)

# np.std(a, ddof): sqrt(sum((a - mean)^2) / max(n - ddof, 0)), NaN for n <= ddof.
function _np_std(a::AbstractVector{Float64}; ddof::Int = 0)
    m = _np_mean(a)
    x = [(v - m) * (v - m) for v in a]
    return sqrt(_np_sum(x) / max(length(a) - ddof, 0))
end

# np.median: middle element, or mean of the two middle elements.
function _np_median(a::AbstractVector{Float64})
    any(isnan, a) && return NaN
    s = sort(a)
    n = length(s)
    h = n ÷ 2
    return isodd(n) ? s[h + 1] : (0.0 + s[h] + s[h + 1]) / 2
end

# np.percentile(a, q) with the default "linear" method (Hyndman-Fan 7) and numpy's _lerp.
function _np_percentile(a::AbstractVector{Float64}, q::Real)
    any(isnan, a) && return NaN
    s = sort(a)
    n = length(s)
    v = (n - 1) * (Float64(q) / 100)
    if v >= n - 1
        return s[n]
    end
    k = floor(v)
    t = v - k
    lo, hi = s[Int(k) + 1], s[Int(k) + 2]
    d = hi - lo
    return t >= 0.5 ? hi - d * (1 - t) : lo + d * t
end

# Python f"{x:.{d}f}" (correctly rounded like C printf; "nan"/"inf" spelt the Python way).
function _pyfmt(x::Real, d::Int)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    return Printf.format(Printf.Format("%.$(d)f"), Float64(x))
end

# One csv.writer row (QUOTE_MINIMAL, "\r\n" terminator).
function _csvrow(io::IO, fields)
    for (i, f) in enumerate(fields)
        i > 1 && print(io, ',')
        s = string(f)
        if any(c -> c in (',', '"', '\r', '\n'), s)
            print(io, '"', replace(s, "\"" => "\"\""), '"')
        else
            print(io, s)
        end
    end
    print(io, "\r\n")
    return nothing
end

# Appends the decimal digits of n >= 0 to buf.
function _put_int!(buf::Vector{UInt8}, n::Integer)
    n == 0 && return push!(buf, UInt8('0'))
    start = length(buf) + 1
    while n > 0
        push!(buf, UInt8('0') + UInt8(n % 10))
        n ÷= 10
    end
    return reverse!(buf, start, length(buf))
end

# C++17 std::to_chars(double), appended to buf: the shortest round-trip digits, written in
# fixed or %e style, whichever is shorter (fixed on a tie). This is how VTK's OBJ writer prints
# coordinates. `dg` is scratch space for the digits.
function _put_to_chars!(buf::Vector{UInt8}, dg::Vector{UInt8}, x::Float64)
    signbit(x) && push!(buf, UInt8('-'))
    iszero(x) && return push!(buf, UInt8('0'))
    s = codeunits(string(abs(x)))      # Julia's shortest round-trip spelling
    empty!(dg)
    ndot = -1
    e = 0
    i = 1
    while i <= length(s)
        c = s[i]
        if c == UInt8('.')
            ndot = length(dg)
        elseif c == UInt8('e')
            neg = s[i + 1] == UInt8('-')
            for j in (i + 1 + neg):length(s)
                e = 10e + (s[j] - UInt8('0'))
            end
            neg && (e = -e)
            break
        else
            push!(dg, c)
        end
        i += 1
    end
    ndot < 0 && (ndot = length(dg))
    lead = findfirst(!=(UInt8('0')), dg)
    last = findlast(!=(UInt8('0')), dg)
    xe = ndot - 1 + e - (lead - 1)     # decimal exponent of the first significant digit
    nd = last - lead + 1
    sci = nd + (nd > 1) + 2 + max(2, ndigits(abs(xe)))
    fixed = xe >= nd - 1 ? xe + 1 : xe >= 0 ? nd + 1 : nd + 1 - xe
    if fixed <= sci
        if xe >= nd - 1
            append!(buf, @view dg[lead:last])
            for _ in 1:(xe - nd + 1)
                push!(buf, UInt8('0'))
            end
        elseif xe >= 0
            append!(buf, @view dg[lead:(lead + xe)])
            push!(buf, UInt8('.'))
            append!(buf, @view dg[(lead + xe + 1):last])
        else
            push!(buf, UInt8('0'), UInt8('.'))
            for _ in 1:(-xe - 1)
                push!(buf, UInt8('0'))
            end
            append!(buf, @view dg[lead:last])
        end
    else
        push!(buf, dg[lead])
        nd > 1 && (push!(buf, UInt8('.')); append!(buf, @view dg[(lead + 1):last]))
        push!(buf, UInt8('e'), xe < 0 ? UInt8('-') : UInt8('+'))
        abs(xe) < 10 && push!(buf, UInt8('0'))
        _put_int!(buf, abs(xe))
    end
    return buf
end
_to_chars(x::Float64) = String(_put_to_chars!(UInt8[], UInt8[], x))
