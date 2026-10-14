# Where CoronaryCCO.jl departs from upstream HCT-CCO, and why

CoronaryCCO.jl reproduces HCT-CCO (commit `9bc141f`) bit for bit in its Float64 exact mode:
on the `case0` sample, every segment's topology, coordinates, radius, resistance, flow and
radius ratio match the C++ reference build at 500 and 2,000 terminals. Every behaviour below
is therefore reproduced by default, and each correction is an explicit option.

For each item: what the released code does, what the HCT-CCO paper says (Vajire, Choy, Kassab
and Lee, *Ann. Biomed. Eng.* 2026, doi:10.1007/s10439-026-04383-4), where the behaviour comes
from, why it was easy to miss, what it changes, and what CoronaryCCO offers. These notes are
meant to help interpret upstream results, not to diminish the method; several of the
behaviours are inherited from OpenCCO, and some are invisible at the scales the paper ran.

## 1. Growth stops near 50 µm; orders 1–3 come from post-processing

**Released code.** Terminals are added until the requested count (the paper used 1.2–2.0 × 10⁴,
Table 2). With the code's default γ = 3.0, a 40,000-terminal case0 tree has its smallest
vessels at about 50 µm, i.e. diameter-defined order 4.

**Paper.** Orders 1–11 are obtained after growth: diameter-defined Strahler ordering, then an
order-wise affine transform that sets each order's mean and SD to Kassab's (Eq. 13, "the final
step in the generation pipeline to define the final generated coronary tree"). Table 3 then
lists the "smallest resolved diameter" as ~10 µm. The released post-processing instead assigns
orders from Kassab's diameter boundaries (`literature_order_of`).

**Why it was not addressed.** Cost. The paper reports 1 h for 5 × 10³ terminals, 5 h for 10⁴ and
15–20 h for 2 × 10⁴ ("Computational Performance"). Reaching ~9 µm terminals needs roughly
8 × 10⁵ terminals (item 2), which at that growth rate is out of reach, so the hierarchy below
the grown terminals had to be filled in analytically.

**Effect.** Agreement of mean diameter per order with Kassab is imposed by Eq. 13 (or by the
literature thresholds), so it cannot by itself validate the growth. The hemodynamic check also
lacks the arterioles below ~50 µm, which carry most of the arterial pressure drop.

**CoronaryCCO.** Fast mode grows 2 × 10⁵ terminals in about 100 s on a laptop CPU, so terminal
counts that reach true order-1 arterioles are practical. Orders are computed from the tree
itself (`dds_orders`: topological Strahler, and the Jiang–Kassab–Fung diameter-defined
procedure) with no diameter changes.

**Evidence.** The diameter-defined Strahler system defines order boundaries from the measured
vessels' own diameter statistics (Jiang, Kassab and Fung 1994), which is how Kassab et al.
(1993) classified the porcine data the tree is compared against.

## 2. Bifurcation exponent: the paper uses 2.33, the code defaults to 3.0

**Released code.** `-g/--gamma` defaults to 3.0, inherited from OpenCCO (`my_gamma = 3.0`). The
README example does not set it.

**Paper.** γ = 2.33 throughout (Table 1, "HK coronary scaling law"), and the 0.7 taper floor is
justified against the HK symmetric daughter/parent ratio of 0.74.

**Why it matters.** In a roughly symmetric tree, terminal diameter ≈ root diameter × N^(−1/γ).
From the 3.18 mm root to Kassab's 9 µm order-1 arterioles (a ratio of ~350) takes ~8 × 10⁵
terminals at γ = 2.33 but ~4 × 10⁷ at γ = 3. Measured: at 2 × 10⁵ terminals the median leaf
is 51 µm with γ = 3 and 16.8 µm with γ = 2.33.

**CoronaryCCO.** `gamma` is a parameter; the first-principles runs use the paper's 2.33.

**Evidence.** The paper sets γ = 2.33 and attributes it to the HK coronary scaling law, citing
Huo and Kassab (2009, its reference 20). That paper derives a resistance scaling law; the
flow–diameter law with exponent 7/3 is in Huo and Kassab (2012, Eq. 2.4,
Q_s ∝ D_s^(7/3), "similar to Murray's law but with a different exponent"), validated on the
entire pig coronary trees. γ = 3 is Murray's minimum-work law (Murray 1926). Huo and Kassab
(2012) also report that the measured diameter–branching exponent varies with size, from about
2.1 above 200 µm to about 3 below 70 µm, so a single γ = 7/3 is a whole-tree compromise; a
size-dependent exponent is a possible refinement.

**How many terminals is "true terminal"?** The entire pig LAD and LCx trees reconstructed from
Kassab's morphometry contain 946,937 and 571,383 vessel segments down to the smallest
arterioles, capillaries (< 8 µm) excluded (Huo and Kassab 2009, Fig. 1; 2012, Fig. 2). A
binary tree with S segments has (S + 1)/2 terminals, so an LAD + LCx base tree corresponds to
about 7.6 × 10⁵ terminals, consistent with the ~8 × 10⁵ estimated above for 9 µm leaves at
γ = 2.33.

## 3. The volume objective skips everything below a single-child node

**Released code.** `computeTotalVolume` adds a node's children only at true bifurcations.

**Origin.** Identical in OpenCCO, where it is harmless: CCO grown from a single root never
creates single-child nodes. HCT-CCO introduced imported centerline base trees, which contain
them (case0 has 111, starting at the first segment), and added single-child handling to the
resistance update (`HctccoModel.ih`, "Unary chain handling") and the radius update, but not to
the volume.

**Paper.** Does not discuss it; Section "Local Geometric Optimization" describes the classic
CCO objective.

**Effect.** On case0 the objective counts 1 of 1,242 segments at 500 terminals and 183 of 8,242
at 4,000, i.e. 2.6% of the tree's volume. Candidate bifurcations are then compared on an
almost constant number, so the choice falls to the restriction schedule and to the order in
which nearest segments are tried.

**Why it was easy to miss.** The trees still look plausible: the restriction schedule and
nearest-neighbour attachment produce hierarchy on their own.

**CoronaryCCO.** `volume = :full` counts the whole tree.

**Evidence.** The CCO cost is the total intravascular volume of the tree,
V = π Σᵢ rᵢ² lᵢ over all segments (Schreiner and Buxbaum 1993; Karch et al. 1999; OpenCCO
Eq. 20, Kerautret et al. 2023).

## 4. The wall-distance field is an L3 distance, not Euclidean

**Released code.** `getImageDistance` instantiates
`ExactPredicateLpSeparableMetric<Space, Domain::dimension>` under the name `L2Metric`; the second
template argument is the exponent p, so in 3D the distance is (|dx|³+|dy|³+|dz|³)^(1/3).

**Origin.** Inherited unchanged from OpenCCO (`helpers/GeomHelpers.h`). In 2D, p equals the
dimension and the metric is Euclidean, so the issue only appears in 3D.

**Paper.** Describes "a Euclidean distance transform (EDT)" (citing Danielsson 1980) for the
clearance from the epicardial and endocardial surfaces.

**Effect.** The L3 distance is shorter than the Euclidean one off the axes (1.26 instead of
1.41 at one voxel diagonally from the wall), so the `-m` clearance test is applied with the
wrong geometry near oblique walls. The Euclidean distance and DGtal's output differed for the
first time at step 29 of a 500-terminal run, which is how this was found.

**CoronaryCCO.** `metric = :euclidean`.

**Evidence.** Danielsson (1980); the paper's own Methods.

## 5. Candidate terminals only at integer voxel points

**Released code.** `ImageMaskDomainCtrl::randomPoint` returns integer voxel coordinates.

**Origin.** Inherited from OpenCCO; its analytic sphere domain samples continuously.

**Effect.** The spacing test rejects a second terminal at the same point, so the terminal count
cannot exceed the number of admissible voxels (~2.5 × 10⁵ for case0's 1 mm mask). The paper's
counts (≤ 2 × 10⁴) never came close, so the cap was invisible.

**CoronaryCCO.** `sampling = :continuous` samples uniformly inside admissible voxels and
re-checks the wall clearance at the sampled point.

**Evidence.** CCO places each new terminal at a random location drawn from the perfusion
domain (Schreiner and Buxbaum 1993; Karch et al. 1999), not on a lattice.

## 6. Splitting an imaged vessel moves its node

**Released code.** Adding a bifurcation on segment s moves s's distal node to the optimized
bifurcation point. This is standard CCO for synthetic segments, and HCT-CCO applies it to the
imported base tree as well.

**Paper.** The base tree is "the proximal architecture resolved in CT"; "the algorithm does
not regenerate the visible vasculature … while preserving the extracted macrovascular
geometry".

**Effect.** Imaged vessels bend. On a 20,000-terminal case0 tree, joints of vessels with
r ≥ 0.25 mm on the base tree reached a 124° change of direction, against at most 82° in the
input centerlines.

**Why it was easy to miss.** Kinks are local and the paper's figures show whole hearts; the
reported geometric check is vessel crossings (< 0.1% of segments), not path curvature.

**CoronaryCCO.** `fix_base = true` keeps junctions on imaged centerlines: when the host segment
belongs to the base tree, the bifurcation point is projected onto its centerline (5–95% of
its length). Pieces of a split base vessel stay "base", so later splits keep them in place.

**Evidence.** The paper's own design statement: the segmented base tree is the imaged
anatomy that growth extends. (Jaquet et al. 2019, the hybrid method the paper compares
against, likewise starts growth from the segmented epicardial vessels.)

## 7. Myocardial mass in the hemodynamic check

**Released code.** `_read_domain_mass_g` counts non-zero bytes of the `.vol` payload without
decompressing it. For the shipped case0 (Version 3, zlib) that is 41.74 g; the decoded mask
(voxels > 128) is 254,980 mm³, 267.7 g at 1.05 g/mL. Every flow, pressure drop and wall shear
stress scales with this mass.

**Paper.** Q = q·ρ·V with V the segmented myocardial volume and ρ ≈ 1.05 g/mL.

**Note.** If the paper's own domains were uncompressed (Version 2) masks with a zero
background, the count would have been right there; the error shows up on the public sample.

**CoronaryCCO.** The post-processing defaults to the decoded mask (`mass = :mask`);
`mass = :upstream` reproduces the Python for verification.

## 8. Reproducibility details (no change in behaviour)

- `std::sort` is not stable, so equal distances in the nearest-segment sort are ordered
  differently by libc++ and libstdc++ (C++ standard, [alg.sort]); the reference build and the
  port break ties by segment index.
- `rand()` is seeded from the wall clock; both sides use the same seeded SplitMix64 stream.
- Ceres is pinned (2.2.0 with Eigen 3.4.0, sparse normal Cholesky through Eigen, no FMA
  contraction) so the inner Kamiya solve is reproducible bit for bit.
- Smaller code behaviours reproduced for parity: the constructor's terminal entry for segment
  1 survives the XML import, `myVectParent[0]` is overwritten on leaf splits, XML values are
  read in single precision, and the XML `pTerm` does not update the pressure drop.

## 9. Fast mode implements the paper's incremental update

The paper describes updating resistances "incrementally along the path from the new terminal
to the root" using a radius-independent baseline resistance R̃ = 8μL/π (Eq. 2 and the text
before it); the released code instead recomputes every segment, every radius and the total
volume for every Kamiya iteration of every candidate. Fast mode stores R̃ and a per-node
subtree volume factor so that each candidate and each commit touches only the path to the
root, which is what makes 10⁵–10⁶ terminals practical. It equals exact mode to rounding.

## References

- Danielsson PE. Euclidean distance mapping. *Comput. Graph. Image Process.* 14:227–248, 1980. doi:10.1016/0146-664X(80)90054-4
- Huo Y, Kassab GS. The scaling of blood flow resistance: from a single vessel to the entire distal tree. *Biophys. J.* 96:339–346, 2009. doi:10.1016/j.bpj.2008.09.038
- Huo Y, Kassab GS. Intraspecific scaling laws of vascular trees. *J. R. Soc. Interface* 9:190–200, 2012. doi:10.1098/rsif.2011.0270
- Jaquet C, et al. Generation of patient-specific cardiac vascular networks: a hybrid image-based and synthetic geometric model. *IEEE Trans. Biomed. Eng.* 66:946–955, 2019. doi:10.1109/TBME.2018.2865667
- Jiang ZL, Kassab GS, Fung YC. Diameter-defined Strahler system and connectivity matrix of the pulmonary arterial tree. *J. Appl. Physiol.* 76:882–892, 1994. doi:10.1152/jappl.1994.76.2.882
- Karch R, Neumann F, Neumann M, Schreiner W. A three-dimensional model for arterial tree representation, generated by constrained constructive optimization. *Comput. Biol. Med.* 29:19–38, 1999. doi:10.1016/S0010-4825(98)00045-6
- Kassab GS, Rider CA, Tang NJ, Fung YC. Morphometry of pig coronary arterial trees. *Am. J. Physiol.* 265:H350–H365, 1993. doi:10.1152/ajpheart.1993.265.1.H350
- Kerautret B, Ngo P, Passat N, Talbot H, Jaquet C. OpenCCO: an implementation of constrained constructive optimization for generating 2D and 3D vascular trees. *Image Process. On Line* 13:258–279, 2023. doi:10.5201/ipol.2023.477
- Murray CD. The physiological principle of minimum work. *Proc. Natl. Acad. Sci. USA* 12:207–214, 1926. doi:10.1073/pnas.12.3.207
- Schreiner W, Buxbaum PF. Computer-optimization of vascular trees. *IEEE Trans. Biomed. Eng.* 40:482–491, 1993. doi:10.1109/10.243413
- Vajire SL, Choy JS, Kassab GS, Lee LC. Hierarchy aware CT-driven constrained constructive optimization: a framework for coronary network generation. *Ann. Biomed. Eng.*, 2026. doi:10.1007/s10439-026-04383-4
