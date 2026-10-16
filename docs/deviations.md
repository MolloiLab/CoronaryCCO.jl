# Where Vasculature.jl departs from upstream HCT-CCO, and why

Vasculature.jl reproduces HCT-CCO (commit `9bc141f`) bit for bit in its Float64 exact mode:
on the `case0` sample, every segment's topology, coordinates, radius, resistance, flow and
radius ratio match the C++ reference build at 500 and 2,000 terminals. Every behaviour below
is therefore reproduced by default, and each correction is an explicit option.

For each item: what the released code does, what the HCT-CCO paper says (Vajire, Choy, Kassab
and Lee, *Ann. Biomed. Eng.* 2026, doi:10.1007/s10439-026-04383-4), where the behaviour comes
from, why it was easy to miss, what it changes, and what Vasculature offers. These notes are
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

**Vasculature.** Fast mode grows 2 × 10⁵ terminals in about 100 s on a laptop CPU, so terminal
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

**Vasculature.** `gamma` is a parameter; the first-principles runs use the paper's 2.33.

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

**Vasculature.** `volume = :full` counts the whole tree.

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

**Vasculature.** `metric = :euclidean`.

**Evidence.** Danielsson (1980); the paper's own Methods.

## 5. Candidate terminals only at integer voxel points

**Released code.** `ImageMaskDomainCtrl::randomPoint` returns integer voxel coordinates.

**Origin.** Inherited from OpenCCO; its analytic sphere domain samples continuously.

**Effect.** The spacing test rejects a second terminal at the same point, so the terminal count
cannot exceed the number of admissible voxels (~2.5 × 10⁵ for case0's 1 mm mask). The paper's
counts (≤ 2 × 10⁴) never came close, so the cap was invisible.

**Vasculature.** `sampling = :continuous` samples uniformly inside admissible voxels and
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

**Vasculature.** `fix_base = true` keeps junctions on imaged centerlines: when the host segment
belongs to the base tree, the bifurcation point is projected onto its centerline (5–95% of
its length). Pieces of a split base vessel stay "base", so later splits keep them in place.

**Evidence.** The paper's own design statement: the segmented base tree is the imaged
anatomy that growth extends. (Jaquet et al. 2019, the hybrid method the paper compares
against, likewise starts growth from the segmented epicardial vessels.)

## 7. Collision tests use insertion-time radii, which are far below final radii

**Released code.** A candidate bifurcation is rejected if its three segments intersect an
existing vessel or are shorter than twice their radius, using the radii of the current tree.
The test is never repeated.

**Why it matters.** Each terminal carries q_term = Q_perf / N_final, so early in growth the tree
carries only a fraction of its final flow and every radius is a fraction of its final size.
On case0 the root radius is 0.15 mm after import, 0.70 mm at 1% of growth and 1.42 mm at 5%,
reaching its 1.588 mm cap only at about 10% (2 × 10⁵ terminals, γ = 2.33). The large vessels
form in that first 10% (the unrestricted stage), so they are checked at radii up to ten times
too small and later overlap: at 10⁵ terminals, 69 pairs of large vessels (r ≥ 0.25 mm)
overlapped, by up to 1.66 mm, including branches folded back alongside their grandparents.

**Vasculature.** `collision = :final` multiplies every radius in the intersection, spacing and
degenerate-length tests by (root radius cap / current root radius). For a fixed tree
structure all radii scale with the root radius, so this anticipates each vessel's final
size. The degenerate-length test (2r ≤ l) compares a radius with a model length, which is the
geometric length times a length factor that grows as √(terminals + 1). So that test also
scales lengths to their final value, by √((N_final + 1)/(terminals + 1)), and compares final
radii with final lengths. Scaling the radii alone made every short segment look degenerate
early in growth. On an imaged base tree with few terminals (an XCAT right coronary artery with
3 terminals), no candidate was ever accepted. With both sides scaled, that tree grows normally. At 10⁵ terminals the worst overlap fell to 0.33 mm, and what remains is between
segments two steps apart at closely spaced junctions (cylinders overlapping at a shared node),
which a merged surface mesh absorbs. `tree_qc` reports these checks for any grown tree.

## 8. No limit on branching angles

**Released code.** The bifurcation point comes from the Kamiya optimization and the
collision tests; no angle is checked, so a daughter can leave its parent backwards.

**Evidence.** Murray's minimum-work analysis of branching angles (Murray 1926, *J. Gen.
Physiol.*) predicts that the smaller daughter leaves at most perpendicular to the parent and
the larger one nearly straight on; no daughter turns back.

**Vasculature.** `max_branch_angle = 90` rejects candidate bifurcations in which the
continuing segment, the new branch, or the moved parent segment turns by more than 90°. At
10⁵ terminals this removed every grown-vessel turn above 90° (45 among vessels with
r ≥ 0.1 mm without it).

## 9. Myocardial mass in the hemodynamic check

**Released code.** `_read_domain_mass_g` counts non-zero bytes of the `.vol` payload without
decompressing it. For the shipped case0 (Version 3, zlib) that is 41.74 g; the decoded mask
(voxels > 128) is 254,980 mm³, 267.7 g at 1.05 g/mL. Every flow, pressure drop and wall shear
stress scales with this mass.

**Paper.** Q = q·ρ·V with V the segmented myocardial volume and ρ ≈ 1.05 g/mL.

**Note.** If the paper's own domains were uncompressed (Version 2) masks with a zero
background, the count would have been right there; the error shows up on the public sample.

**Vasculature.** The post-processing defaults to the decoded mask (`mass = :mask`);
`mass = :upstream` reproduces the Python for verification.

## 10. Reproducibility details (no change in behaviour)

- `std::sort` is not stable, so equal distances in the nearest-segment sort are ordered
  differently by libc++ and libstdc++ (C++ standard, [alg.sort]); the reference build and the
  port break ties by segment index.
- `rand()` is seeded from the wall clock; both sides use the same seeded SplitMix64 stream.
- Ceres is pinned (2.2.0 with Eigen 3.4.0, sparse normal Cholesky through Eigen, no FMA
  contraction) so the inner Kamiya solve is reproducible bit for bit.
- Smaller code behaviours reproduced for parity: the constructor's terminal entry for segment
  1 survives the XML import, `myVectParent[0]` is overwritten on leaf splits, XML values are
  read in single precision, and the XML `pTerm` does not update the pressure drop.

## 11. Fast mode implements the paper's incremental update

The paper describes updating resistances "incrementally along the path from the new terminal
to the root" using a radius-independent baseline resistance R̃ = 8μL/π (Eq. 2 and the text
before it); the released code instead recomputes every segment, every radius and the total
volume for every Kamiya iteration of every candidate. Fast mode stores R̃ and a per-node
subtree volume factor so that each candidate and each commit touches only the path to the
root, which is what makes 10⁵–10⁶ terminals practical. It equals exact mode to rounding.

## 12. Batched growth (`growth = :batched`)

**What upstream does.** CCO adds one terminal per step (Schreiner & Buxbaum 1993). Each step
optimizes one bifurcation against the whole current tree, then updates every flow, resistance
and radius before the next step starts. The steps cannot overlap.

**What batched growth changes.** After the first `n_seq` terminals (default 10,000), terminals
are added in rounds. A round draws up to `batch_frac` × (current terminals) locations from the
same random stream (default 1%). Every location's best bifurcation is found in parallel
against the same tree. Then, in draw order, a location is accepted unless it conflicts with
one already accepted in the round. Conflicts are:

- a host segment that is the same as, or the parent, child or sibling of, another accepted host;
- a new terminal closer to another new terminal than the current spacing threshold;
- new vessels that intersect, or pass closer than the usual clearance to, the new vessels of
  another accepted insertion, tested at the same scaled radii (`collision = :final`) as every
  other collision test.

Rejected locations are dropped and new ones are drawn in later rounds, so the final terminal
count is unchanged. The accepted insertions are committed, and flows, resistances, radius ratios
and radii are recomputed once per round over the union of their paths to the root. Every rule
of sequential growth still applies to each insertion: the volume objective, the Kamiya
bifurcation solve, wall distance, the branching-angle limit and the imaged-vessel rule.
Selection is deterministic and does not depend on the number of threads.

**What it approximates.** Insertions in the same round do not see each other's flow when
their bifurcations are optimized. Locations are spread over the domain in proportion to volume,
so a round of 1% of the tree adds on average about 1% to the flow of any subtree. A small
subtree can receive one extra terminal more than its share; adjacent hosts are excluded, so the
segments an insertion is optimized against are never split by another insertion in the round.

**Why this is acceptable.** CCO has no single correct tree. Trees grown from different random
seeds differ in structure but are equally valid and statistically alike (Schreiner et al.
1997). So the test for batched growth is whether its trees fall inside the spread of sequential
trees over seeds, not whether they match one tree bit for bit. Cury et al. (2021) grew CCO
networks concurrently by domain decomposition for the same reason and compared them against
sequentially grown networks up to 200,000 segments. With `batch_frac = 0`, or `n_seq ≥ nterm`,
batched growth is exactly sequential growth.

**Measured agreement.** case0 with the corrected options (γ = 2.33, `volume = :full`,
`metric = :euclidean`, `sampling = :continuous`, `fix_base`, `collision = :final`,
`max_branch_angle = 90`), 10⁵ terminals, seeds 1–5, `n_seq = 2000`, `batch_frac = 0.01`
(`bench/batch_validate.jl`). Ranges are min–max over the five seeds.

| Metric | Sequential | Batched, CPU | Batched, GPU (Metal) |
|---|---|---|---|
| wall time, 8 threads (s) | 35–51 | 6.1–7.8 | 7.4–7.9 |
| total volume (mm³) | 1441–1614 | 1441–1620 | 1442–1619 |
| leaf diameter, median (µm) | 22.56–22.60 | 22.55–22.59 | 22.54–22.58 |
| mean terminal depth (segments) | 42.1–44.7 | 41.9–44.8 | 42.0–44.8 |
| overlaps, r ≥ 0.25 mm | 37–62 | 26–58 | 28–58 |
| grown turns > 90° | 0 | 0 | 0 |
| DDS orders 1–10: count, diameter, length | — | 30 of 30 inside | 30 of 30 inside |
| Strahler orders 1–10: count, diameter, length | — | 26 of 30 inside | 27 of 30 inside |

The rows outside the sequential range miss it by about 1% or less: Strahler order 3–6 mean
diameters and order 2–4 counts. Kassab's per-order SDs are 10–25% of the mean. Batched trees
on the CPU and on a GPU differ in the last bits of exp and log, so their trees differ in
detail but not in statistics. Each backend is deterministic.

Growth time at 8 × 10⁵ terminals (M4, 8 threads): sequential fast mode more than 12 min;
batched 48 s on the CPU and 41 s with the round evaluation on the GPU.

## 13. Junction symmetry floor (`symmetry_floor`, off by default)

**What upstream does.** CCO accepts any bifurcation that lowers the tree volume. Late in growth,
the cheapest place for a new terminal is often the side of a medium arteriole, so capillary-
feeding twigs hang straight off vessels of 30–100 µm. In a 1 M-terminal coronary forest, an
order-5 element (diameter-defined Strahler) carries about 1.7 order-1 twigs; in Kassab's pig
casts it carries none (Kassab et al. 1993, Table 7).

**What the option does.** `symmetry_floor` takes knots `(mother diameter µm, minimum
d_small/d_large)`. A candidate bifurcation is rejected when its two daughters are more
asymmetric than the floor at the mother's final-scale diameter (the radius scale already used
by `collision = :final`). Knots are interpolated in log diameter and held flat beyond the ends.

**Source.** Measured junctions become more symmetric as the mother shrinks: in pig casts
(VanBavel & Spaan 1992; Kalsho & Kassab 2004) and in about 88,000 human coronary bifurcations
(Schwarz et al. 2020). The XCAT growth uses the human lower whisker (Q1 − 1.5 IQR) per mother-
diameter class, 0.28 at 37 µm down to 0.02 at 566 µm
(`data/reference/schwarz2020_human_coronary_symmetry.toml`), so only junctions that would be
outliers in real hearts are excluded. Order-1 twigs per order-6 element fall from 1.42 to 0.48.

**What it does not change.** Per-order diameters. Kassab's own counts and diameters imply a
junction exponent of 1.4–1.9 in orders 3–7, below every direct junction measurement (2.25–3.2),
so no flow-conserving tree with measured junction laws reproduces his orders 4–7.

## 14. Fine forbidden mask (`forbidden`, off by default)

**What upstream does.** Segments are tested against the domain image at its own voxel size
(1 mm here). A straight segment between two allowed voxels can cut through a thin structure the
domain grid cannot resolve, such as the endocardial edge of a chamber.

**What the option does.** `generate(...; forbidden = (mask, origin, h))` (or `set_forbidden!`)
adds a finer mask, nonzero where vessels may not go. When a candidate bifurcation is accepted,
its new branch and the two moved pieces (unless they lie on an imaged centreline) are sampled
every `h/8`; the candidate is rejected if any sample lies in a forbidden cell. Applied to XCAT
with a 0.4 mm blood mask (chamber pools, great arteries and veins, from the 0.2 mm labels),
one growth gives no grown artery or vein touching blood while the 1 mm domain keeps the whole
myocardium, including the thin right-ventricular wall.

## 15. Staged growth only, with the stage schedule deciding where vessels attach

**What upstream does.** HCT-CCO's hierarchy-aware schedule is documented as three stages
(`getConnectionRestrictionLevel`): for the first 10% of terminals a new vessel may attach to any
segment, up to 35% to a terminal or a few generations above it, then only to a terminal or a
terminal's parent. `getN_NearestSegments` applies that schedule to the candidate list. But
`isAddable` also rejects every host that is not a terminal or a terminal's parent, at every
stage, so the first two stages only waste candidate slots: their nearest segments are mostly
interior ones, which are then refused. With an imaged base tree that lists few terminal ends this
stalls growth: XCAT's lung arteries list 333 ends on 3915 nodes, and staged growth accepted no
terminal in 7000 consecutive batched rounds. `--anywhereModel` skips both rules.

**What Vasculature.jl does.** Growth is always staged; there is no unrestricted mode (the
`anywhere_mode` option was removed). The schedule is applied where candidates are chosen, and
`isAddable`'s extra check is dropped (`upstream_attach = true` restores it, for the bitwise
parity test only). XCAT lung, 200k terminals: 89 s (unrestricted: 87 s). XCAT heart, 957k
terminals: 63 s (53 s); Kassab connectivity F1 LCA 0.882, RCA 0.876 (unrestricted 0.876, 0.868),
mean |log10 diameter ratio| 0.054 and 0.075 (0.049, 0.070).

**Connection length.** Upstream refuses any new vessel longer than 25 units
(`getMaxBranchLength`, "same units as post-processing filter"), a coronary value in mm. It is
now `max_connection` (default 25 mm), so another organ can set its own.

## References

- Cury LFM, Maso Talou GD, Younes-Ibrahim M, Blanco PJ. Parallel generation of extensive vascular networks with application to an archetypal human kidney model. *R. Soc. Open Sci.* 8:210973, 2021. doi:10.1098/rsos.210973
- Danielsson PE. Euclidean distance mapping. *Comput. Graph. Image Process.* 14:227–248, 1980. doi:10.1016/0146-664X(80)90054-4
- Huo Y, Kassab GS. The scaling of blood flow resistance: from a single vessel to the entire distal tree. *Biophys. J.* 96:339–346, 2009. doi:10.1016/j.bpj.2008.09.038
- Huo Y, Kassab GS. Intraspecific scaling laws of vascular trees. *J. R. Soc. Interface* 9:190–200, 2012. doi:10.1098/rsif.2011.0270
- Jaquet C, et al. Generation of patient-specific cardiac vascular networks: a hybrid image-based and synthetic geometric model. *IEEE Trans. Biomed. Eng.* 66:946–955, 2019. doi:10.1109/TBME.2018.2865667
- Jiang ZL, Kassab GS, Fung YC. Diameter-defined Strahler system and connectivity matrix of the pulmonary arterial tree. *J. Appl. Physiol.* 76:882–892, 1994. doi:10.1152/jappl.1994.76.2.882
- Karch R, Neumann F, Neumann M, Schreiner W. A three-dimensional model for arterial tree representation, generated by constrained constructive optimization. *Comput. Biol. Med.* 29:19–38, 1999. doi:10.1016/S0010-4825(98)00045-6
- Kalsho G, Kassab GS. Bifurcation asymmetry of the porcine coronary vasculature and its implications on coronary flow heterogeneity. *Am. J. Physiol. Heart Circ. Physiol.* 287:H2493–H2500, 2004. doi:10.1152/ajpheart.00371.2004
- Schwarz JCV, et al. Topologic and hemodynamic characteristics of the human coronary arterial circulation. *Front. Physiol.* 10:1611, 2020. doi:10.3389/fphys.2019.01611
- VanBavel E, Spaan JAE. Branching patterns in the porcine coronary arterial tree: estimation of flow heterogeneity. *Circ. Res.* 71:1200–1212, 1992. doi:10.1161/01.RES.71.5.1200
- Kassab GS, Rider CA, Tang NJ, Fung YC. Morphometry of pig coronary arterial trees. *Am. J. Physiol.* 265:H350–H365, 1993. doi:10.1152/ajpheart.1993.265.1.H350
- Kerautret B, Ngo P, Passat N, Talbot H, Jaquet C. OpenCCO: an implementation of constrained constructive optimization for generating 2D and 3D vascular trees. *Image Process. On Line* 13:258–279, 2023. doi:10.5201/ipol.2023.477
- Murray CD. The physiological principle of minimum work. *Proc. Natl. Acad. Sci. USA* 12:207–214, 1926. doi:10.1073/pnas.12.3.207
- Murray CD. The physiological principle of minimum work applied to the angle of branching of arteries. *J. Gen. Physiol.* 9:835–841, 1926. doi:10.1085/jgp.9.6.835
- Schreiner W, Buxbaum PF. Computer-optimization of vascular trees. *IEEE Trans. Biomed. Eng.* 40:482–491, 1993. doi:10.1109/10.243413
- Schreiner W, Neumann F, Neumann M, End A, Roedler SM. Anatomical variability and functional ability of vascular trees modeled by constrained constructive optimization. *J. Theor. Biol.* 187:147–158, 1997. doi:10.1006/jtbi.1997.0413
- Vajire SL, Choy JS, Kassab GS, Lee LC. Hierarchy aware CT-driven constrained constructive optimization: a framework for coronary network generation. *Ann. Biomed. Eng.*, 2026. doi:10.1007/s10439-026-04383-4
