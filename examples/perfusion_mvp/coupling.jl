# SPDX-License-Identifier: GPL-3.0-or-later
# Which venule drains each capillary bed. Start from the nearest venous terminal; a capillary
# network connects each terminal arteriole to several venules, and blood leaves through the
# lowest-pressure one, so a bed whose paired venule sits above its arterial pressure (reversed
# flow) is re-paired with the lowest-pressure venous terminal within `reach` mm (widened by 1 mm per pass, at most
# 6 mm: a capillary network spans millimetres, not centimetres). Solved at
# maximal vasodilation, where the beds' own pressure drop is smallest.
function couple_veins(F, terr, VT, terms, vterms; reach = 3.0f0, passes = 24)
    vleaf = Int32.(vterms[nearest_terminal(VT.b[:, vterms], F.b[:, terms])])
    VB = VT.b[:, vterms]
    for pass in 1:passes
        C = circuit(F, terr; state = :stress, veins = VT, vleaf)
        bad = [i for i in eachindex(terms) if C.pa[C.terminals[i]] <= C.pv[vleaf[i]]]
        @printf("coupling pass %d: %d beds with reversed flow\n", pass, length(bad))
        isempty(bad) && return vleaf, C
        Threads.@threads for i in bad
            p = F.b[:, terms[i]]; pa = C.pa[C.terminals[i]]
            best = (Inf32, vleaf[i])
            for k in eachindex(vterms)
                d2 = (VB[1, k] - p[1])^2 + (VB[2, k] - p[2])^2 + (VB[3, k] - p[3])^2
                d2 <= (reach + min(pass - 1, 3))^2 || continue
                pv = C.pv[vterms[k]]
                pv < best[1] && (best = (pv, Int32(vterms[k])))
            end
            vleaf[i] = best[2]
        end
    end
    return vleaf, circuit(F, terr; state = :stress, veins = VT, vleaf)
end
