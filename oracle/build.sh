#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds the upstream HCT-CCO generator, pinned to commit 9bc141f, with the
# instrumentation in cco_oracle.h. The patches below only swap the RNG for a portable
# one, add counters, and dump the final tree; no arithmetic is touched.
#
# Requires: cmake, Ceres Solver, DGtal (set DGtal_DIR if not on the default prefix).
# Usage:    oracle/build.sh            -> oracle/build/bin/HctCco
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
COMMIT=9bc141f
SRC="$HERE/build/src"
mkdir -p "$HERE/build"
if [ ! -d "$SRC/.git" ]; then
  git clone -q https://github.com/Sujal-vajire/HCT-CCO.git "$SRC"
fi
git -C "$SRC" checkout -q -f "$COMMIT"
git -C "$SRC" clean -qfd

# 1. Portable RNG (same RAND_MAX contract) seeded from $CCO_SEED.
sed -i.bak 's/rand()/cco_rand()/g' "$SRC/src/MyocardialDomain.h"
sed -i.bak 's/srand ((unsigned int) time(NULL));/cco_srand_env();/' "$SRC/src/helpers/CoronaryGrowth.h"
sed -i.bak 's/srand ((int) time(NULL));/cco_srand_env();/' "$SRC/bin/HctCco.cpp"

# 2. Work counters and a coarse per-step clock.
G="$SRC/src/helpers/CoronaryGrowth.h"
sed -i.bak 's/DGtal::trace.progressBar(i, nbSeed);/if (i % 250 == 0) std::fprintf(stderr, "[cco] step %u t %.4f segs %zu\\n", i, cco_oracle::now(), aTree.myVectSegments.size());/' "$G"
sed -i.bak 's/auto pt = aTree.generateNewLocation(nbTryCandite);/auto pt = aTree.generateNewLocation(nbTryCandite); CCO_COUNT("location");/' "$G"
sed -i.bak 's/auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it));/auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it)); CCO_COUNT("neighbor");/' "$G"
sed -i.bak 's/HctccoModel< DomCtr, TDim  > cTree1 = aTree;/CCO_COUNT("isAddable"); HctccoModel< DomCtr, TDim  > cTree1 = aTree;/' "$G"
sed -i.bak 's/if (nbSeed != nbSeedFound){/cco_depth_stats(aTree); cco_dump_tree(aTree); cco_oracle::report(); if (nbSeed != nbSeedFound){/' "$G"
M="$SRC/src/HctccoModel.ih"
sed -i.bak 's/res1 = kamyiaOptimization(pCurrent/CCO_COUNT("kamiya_iter"); res1 = kamyiaOptimization(pCurrent/' "$M"
sed -i.bak 's/isDone = true;/isDone = true; CCO_COUNT("accepted");/' "$M"
C="$SRC/src/helpers/CoronaryGeometry.h"
sed -i.bak 's/return summary.IsSolutionUsable();/CCO_COUNT("ceres_solve"); CCO_COUNT_ADD("ceres_iter", summary.iterations.size()); return summary.IsSolutionUsable();/' "$C"

# Every patch must have applied exactly; fail loudly if upstream moved.
for pat in 'cco_rand()' 'cco_srand_env' 'CCO_COUNT("location")' 'CCO_COUNT("neighbor")' \
           'CCO_COUNT("isAddable")' 'cco_dump_tree' 'CCO_COUNT("kamiya_iter")' \
           'CCO_COUNT("accepted")' 'CCO_COUNT("ceres_solve")' '\[cco\] step'; do
  grep -rqF -- "${pat//\\/}" "$SRC/src" "$SRC/bin" || { echo "patch failed: $pat" >&2; exit 1; }
done

# Release build without FMA contraction so float results match Julia's IEEE ops.
cmake -S "$SRC" -B "$HERE/build/cmake" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-ffp-contract=off -include $HERE/cco_oracle.h" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ${DGtal_DIR:+-DDGtal_DIR="$DGtal_DIR"} >/dev/null
cmake --build "$HERE/build/cmake" -j"$(getconf _NPROCESSORS_ONLN)"
mkdir -p "$HERE/build/bin"
cp "$HERE/build/cmake/bin/HctCco" "$HERE/build/bin/HctCco"
echo "built $HERE/build/bin/HctCco"
