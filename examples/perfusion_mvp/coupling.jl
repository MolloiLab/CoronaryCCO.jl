# SPDX-License-Identifier: GPL-3.0-or-later
import AcceleratedKernels as AK
# Which venule drains each capillary bed. Start from the nearest venous terminal; a capillary
# network connects each terminal arteriole to several venules, and blood leaves through the
# lowest-pressure one, so a bed whose paired venule sits above its arterial pressure (reversed
# flow) is re-paired with the lowest-pressure venous terminal within `reach` mm (widened by 1 mm per pass, at most
# 6 mm: a capillary network spans millimetres, not centimetres). Solved at
# maximal vasodilation, where the beds' own pressure drop is smallest. Stops when the reversed
# beds are at most `tol` of all beds, or when a pass no longer removes at least 10% of them.
function couple_veins(F, terr, VT, terms, vterms; reach = 3.0f0, passes = 24, tol = 5.0f-4, to = identity)
    vleaf = Int32.(vterms[nearest_terminal(VT.b[:, vterms], F.b[:, terms])])
    VB = VT.b[:, vterms]; last = typemax(Int); k0 = 1.0f0
    B = Territory.Buckets(VB, 2.0f0)
    for pass in 1:passes
        t0 = time()
        C = circuit(F, terr; state = :stress, veins = VT, vleaf, kappa_init = k0, to); k0 = C.resist.kappa
        bad = [i for i in eachindex(terms) if C.pa[C.terminals[i]] <= C.pv[vleaf[i]]]
        @printf("coupling pass %d: %d beds with reversed flow (%.1f s)\n", pass, length(bad), time() - t0); flush(stdout)
        (length(bad) <= tol * length(terms) || length(bad) > 0.9 * last) && return vleaf, C
        last = length(bad)
        rr = reach + min(pass - 1, 3); nr = ceil(Int, rr / B.h)
        AK.foreachindex(1:length(bad); min_elems = 1) do ii
            i = bad[ii]; x = F.b[1, terms[i]]; y = F.b[2, terms[i]]; z = F.b[3, terms[i]]
            c0 = (floor(Int, (x - B.lo[1]) / B.h), floor(Int, (y - B.lo[2]) / B.h), floor(Int, (z - B.lo[3]) / B.h))
            bp = Inf32; bk = typemax(Int32)                 # lowest venous pressure, ties to the lowest index
            for kc in (c0[3] - nr):(c0[3] + nr), jc in (c0[2] - nr):(c0[2] + nr), ic in (c0[1] - nr):(c0[1] + nr)
                (0 <= ic < B.dims[1] && 0 <= jc < B.dims[2] && 0 <= kc < B.dims[3]) || continue
                c = 1 + ic + B.dims[1] * (jc + B.dims[2] * kc)
                for s in B.start[c]:(B.start[c + 1] - 1)
                    k = B.items[s]
                    (VB[1, k] - x)^2 + (VB[2, k] - y)^2 + (VB[3, k] - z)^2 <= rr^2 || continue
                    pv = C.pv[vterms[k]]
                    (pv < bp || (pv == bp && k < bk)) && (bp = pv; bk = k)
                end
            end
            bk == typemax(Int32) || (vleaf[i] = Int32(vterms[bk]))
        end
    end
    return vleaf, circuit(F, terr; state = :stress, veins = VT, vleaf, kappa_init = k0, to)
end
