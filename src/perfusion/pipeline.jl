# SPDX-License-Identifier: GPL-3.0-or-later
# The one flow entry point: grown arteries + veins + the tissue they supply → blood flow from the pressures.

using .Flow: FlowSpec, VTree, Circuit, circuit, vtree
using .Territory: territory_shares, demand_field, nearest_terminal

# Leaf segments (terminals) of a vessel tree.
_leaves(t::VTree) = (nch = zeros(Int, length(t.r)); foreach(p -> p > 0 && (nch[p] += 1), t.parent); findall(==(0), nch))

"""
    couple_veins(arteries, terr, veins, terms, vterms; phys, reach = 3, passes = 24, tol = 5e-4, to = identity)
        -> (vleaf, circuit)

Which venule drains each capillary bed. Start from the nearest venous terminal; a capillary network connects each
terminal arteriole to several venules, and blood leaves through the lowest-pressure one, so a bed whose paired venule
sits above its arterial pressure (reversed flow) is re-paired with the lowest-pressure venous terminal within `reach` mm
(widened by 1 mm per pass, at most 6 mm: a capillary network spans millimetres, not centimetres). Solved at maximal
vasodilation, where the beds' own pressure drop is smallest. Stops when the reversed beds are at most `tol` of all beds,
or when a pass no longer removes at least 10% of them.
"""
function couple_veins(F, terr, VT, terms, vterms; phys, reach = 3.0f0, passes = 24, tol = 5.0f-4, to = identity)
    vleaf = Int32.(vterms[nearest_terminal(VT.b[:, vterms], F.b[:, terms])])
    VB = VT.b[:, vterms]; last = typemax(Int); k0 = 1.0f0
    B = Territory.Buckets(VB, 2.0f0)
    C = nothing
    for pass in 1:passes
        t0 = time()
        C = nothing; GC.gc(false) # the previous pass's circuit is not needed: at most one full circuit in memory
        C = circuit(F, terr; state = :stress, veins = VT, vleaf, phys, kappa_init = k0, to); k0 = C.resist.kappa
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
    C = nothing; GC.gc(false)
    return vleaf, circuit(F, terr; state = :stress, veins = VT, vleaf, phys, kappa_init = k0, to)
end

"""
    flow(spec::FlowSpec, arteries::VTree, veins::VTree, tissue::AbstractMatrix; to = identity)
        -> (; rest, stress, terr, pairs, vleaf, terms, vterms, demand)

Blood flow in a grown vascular bed from its two pressures alone, the same for every organ; only `spec` differs
(`FlowSpec(organ = "heart")`, `FlowSpec(organ = "lung")`, or your own values). `tissue` is 3 × n points (mm), one per
tissue voxel. Steps: each voxel's tissue is shared among its nearest terminals (`terr`, `pairs`); each capillary bed is
paired with the venule it drains into (`vleaf`); then the circuit is solved. With `spec.regulation = :o2` (arteriolar
tone follows each bed's O₂ balance, heart) the bed is solved at maximal dilation first (`stress`, which sets the tone
constant κ) and then at rest with that κ, under a fractal O₂-demand field (`demand`, Bassingthwaighte et al.); with
`:passive` (lung) there is one state and `rest === stress`.
"""
function flow(spec::FlowSpec, art::VTree, vein::VTree, tissue::AbstractMatrix; to = identity)
    terms = _leaves(art); vterms = _leaves(vein)
    terr, pairs = territory_shares(art.b[:, terms], tissue); terr .= max.(terr, 1.0f-7)
    vleaf, Cc = couple_veins(art, terr, vein, terms, vterms; phys = spec, to)
    # passive vessels have one state (the solver treats rest as maximal dilation), so the coupling's final circuit is the
    # rest circuit; solving it again would only double the memory of a large tree
    spec.regulation === :passive && return (; rest = Cc, stress = Cc, terr, pairs, vleaf, terms, vterms, demand = nothing)
    Cc = nothing
    dem = demand_field(art.b[:, terms], terr)
    stress = circuit(art, terr; state = :stress, veins = vein, vleaf, demand = dem, phys = spec, to)
    rest = circuit(art, terr; state = :rest, veins = vein, vleaf, kappa = stress.resist.kappa, demand = dem, phys = spec, to)
    return (; rest, stress, terr, pairs, vleaf, terms, vterms, demand = dem)
end
