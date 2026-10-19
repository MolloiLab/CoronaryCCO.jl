# SPDX-License-Identifier: GPL-3.0-or-later
# Reference outputs for test/test_postprocess.jl: runs upstream's
# postprocessing/coronary_postprocess.py run(xml, prefix, vol) unchanged and additionally
# records per-segment internals in <prefix>_internals.txt:
#   line 1: "mass <myocardial mass g>"
#   then per segment: strahler seg2elem(0-based) elem_order final_order
#                     Q_rest P_rest WSS_rest Q_stress P_stress WSS_stress
# (floats as Python float.hex). Needs numpy, vtk and matplotlib; run with `python -I`:
#   python -I oracle/pp_reference.py oracle/build/src/postprocessing/coronary_postprocess.py \
#       oracle/out/pp/tree2000.xml oracle/out/pp/ref2000 oracle/build/src/samples/case0.vol
import sys, types

src_path, xml, prefix, vol = sys.argv[1:5]
src = open(src_path).read()
hooks = [("    diam = 2.0 * rad\n", "    global _SOL, _MASS; _SOL = sol; _MASS = mass_g\n"),
         ("    max_tree_order = int(SOF.max())\n", "    global _RUN; _RUN = (strahler, seg2elem, elem_order, SOF)\n")]
for anchor, hook in hooks:
    assert src.count(anchor) == 1, anchor
    src = src.replace(anchor, hook + anchor)
src = src.replace('if __name__ == "__main__":', "if False:")
M = types.ModuleType("coronary_postprocess")
M.__file__ = src_path
exec(compile(src, src_path, "exec"), M.__dict__)
M.run(xml, prefix, vol)

strahler, seg2elem, elem_order, SOF = M._RUN
with open(prefix + "_internals.txt", "w") as f:
    f.write("mass %s\n" % float(M._MASS).hex())
    for s in range(len(SOF)):
        h = [M._SOL[state][k][s] for state in ("resting", "stress") for k in range(3)]
        f.write("%d %d %d %d %s\n" % (strahler[s], seg2elem[s], elem_order[seg2elem[s]], SOF[s],
                                      " ".join(float(x).hex() for x in h)))
