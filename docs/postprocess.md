# Post-processing

`postprocess` is a pure-Julia port of upstream's `postprocessing/coronary_postprocess.py`
(HCT-CCO commit 9bc141f, entry point `run(xml_path, prefix, domain_vol)`). It reads the tree
XML written by `export_gxl` and runs the same algorithms without Python or VTK.

Two upstream bugs are corrected by default (see [Upstream bugs](#upstream-bugs)). Each has an
option that restores the upstream behaviour. With both options set, the outputs match the
Python's exactly:

```julia
postprocess("tree.xml", "out/case0"; domain_vol = "case0.vol", mass = :upstream, elements = :upstream)
```

```julia
using CoronaryCCO
E = generate("case0.vol", "case0.xml"; nterm = 2000, aperf = 10000, min_dist = 1,
             root = (-26, 41, 13), seed = 1, export_xml = "tree.xml")
using CairoMakie                      # optional: also writes the two figures
r = postprocess("tree.xml", "out/case0"; domain_vol = "case0.vol")
r.orders                              # diameter-defined order of every segment (XML edge order)
```

## Outputs

| File | Content |
|---|---|
| `<prefix>_connectivity_matrix.csv` | mean ± SE number of daughter elements of order m (rows) per parent element of order n (columns) |
| `<prefix>_diameter_stats.csv` | per order: mean diameter ± sample SD |
| `<prefix>_diameter_ranges.csv` | per order: min - max diameter |
| `<prefix>_hemodynamics.csv` | only with `domain_vol`: myocardial mass, then per order the median diameter, WSS, pressure and flow at resting (1.0 mL/min/g) and stress (3.5 mL/min/g) perfusion |
| `<prefix>_tree.vtp` | VTK XML PolyData: one point per node id (Float32), one line per segment, cell arrays `radius` (Float64) and `final_order` (Int32, active scalars) |
| `<prefix>_tree.obj` | 24-sided capped tube mesh, built like `vtkCellDataToPointData` + `vtkTubeFilter` + `vtkOBJWriter` (`obj = false` skips it; it is about 1.3 GB for a 40k-terminal tree) |
| `<prefix>_diameter_vs_order.png` | needs CairoMakie: per-order median and IQR of the tree's diameters over the literature (Kassab) mean ± SD |
| `<prefix>_hemodynamics.png` | needs CairoMakie and `domain_vol`: per-order median WSS and pressure against median diameter |

The CSV files use the same format as the Python, including the `±` signs, `nan` for an
undefined SE, and `\r\n` line endings. The returned NamedTuple holds the computed data:
`orders`, `diameters`, the ordering internals (`strahler`, `seg2elem`, `elem_order`), the
`connectivity` matrices, the `diameter_stats` table, and `hemodynamics` (per-segment flow,
pressure and WSS for both states, plus the per-order medians). The figures are drawn by the
`CoronaryCCOMakieExt` package extension. Without CairoMakie, `postprocess` writes everything
else and logs a note.

## Verification

`oracle/pp_reference.py` runs the upstream script unchanged and also records its per-segment
internals. `test/test_postprocess.jl` compares against these references whenever they exist
under `oracle/out/pp`, running the port with `mass = :upstream, elements = :upstream`. The
trees were grown by this package on `case0` (seed 1) with 200 and 2000 terminals, plus the
40k-terminal benchmark tree (80,242 segments). Results:

- All four CSV files are byte-identical for all three trees.
- Per-segment Strahler orders, element ids and final orders are identical.
- The myocardial mass and the per-segment flow, pressure and WSS of both hemodynamic states
  match bit for bit.
- The VTP files have identical points, lines, `radius` and `final_order` arrays when read
  back with VTK.
- The OBJ files have the same number of lines and the same vertex, normal and face layout.
  Every face line is identical. Between 0.1% and 0.2% of the `v`/`vn` lines differ: 207 of
  182,328 lines (200 terminals), 1,256 of 1,204,728 (2000), and 24,784 of 22,788,728 (40k).
  The largest numeric difference is 1.8e-15. All of these are components close to zero
  (about 1e-9 or smaller). They differ in the last bit because VTK's arm64 build fuses
  multiply-adds into FMAs. Fusing four expressions the same way in Julia removed most of
  them. The port keeps plain IEEE arithmetic so the mesh is the same on every platform.

To match numpy's last digits, the port reproduces numpy's pairwise summation (used by `np.mean`
and `np.std`), its `median` and `percentile` (linear method), and Python's `%f` formatting.
The tube mesh reproduces VTK's sliding normals, bevel vectors and strip layout.

## Upstream bugs

Both are fixed by default. The upstream behaviour stays available for verification.

- **Myocardial mass is read from the compressed file bytes** (`mass`). This bug changes every
  hemodynamics number. Upstream `_read_domain_mass_g` treats the `.vol` payload as raw `uint8`
  voxels and counts the non-zero bytes. The sample `case0.vol` is a Version 3
  (zlib-compressed) file. Upstream therefore counts the non-zero bytes of the compressed
  stream: 39,752 of them, giving **41.74 g**. Decoded, the mask has background 128 and
  foreground 255, so every voxel is non-zero. Even a correct decode followed by upstream's
  "non-zero" rule would count the whole box (2,079.8 g). The default `mass = :mask` decodes
  the file with `read_vol` and counts the voxels the growth domain treats as inside
  (value > 128): 254,980 voxels of 1 mm³ at 1.05 g/cm³ = **267.7 g**. `mass = :upstream`
  reproduces 41.74 g. Total inflow is proportional to the mass, so flows, pressure drops and
  WSS all scale by about 6.4.
- **Elements join sibling segments** (`elements`). Upstream's `build_segment_connections`
  lists every segment that touches a segment's start point as one of its parents, so a
  segment's siblings appear as parents too. Elements ("maximal connected runs of the same
  Strahler order") are grown over these lists. As a result, two daughters of one bifurcation
  with equal Strahler order become one element, e.g. two sibling terminal segments. The
  default `elements = :series` joins only segments that meet end to end, which is what an
  element (a run of segments in series) means. `elements = :upstream` reproduces the
  Python. The fix changes the order of 338 of 4,242 segments in the 2000-terminal tree
  (2,689 elements instead of 2,066). In the 40k tree it changes 4,670 of 80,242 segments
  (53,079 elements instead of 40,066). It also changes the connectivity matrix and the
  per-order statistics.

## Deliberate differences

- **VTP encoding.** The arrays are written inline as uncompressed base64 (`format="binary"`),
  not appended and zlib-compressed as `vtkXMLPolyDataWriter` does. The content is the same,
  and the file opens in ParaView and VTK.
- **OBJ header.** The first comment line names CoronaryCCO.jl instead of VTK.
- **Figures** are redrawn with Makie at the Python figure sizes (8 x 5.5 in and 12 x 5 in,
  150 dpi). They show the same series, bands, scales, labels and titles as the matplotlib
  figures. They are not pixel-identical.
- **Malformed input.** A node without a position is placed at the origin, which is what the
  Python does for the VTK output. In this case the Python's hemodynamics step would instead
  raise a `KeyError`. A missing edge radius reads as 0.
