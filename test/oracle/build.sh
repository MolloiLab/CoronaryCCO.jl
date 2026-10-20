#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds the upstream HCT-CCO generator, pinned to commit 9bc141f, with the
# instrumentation in cco_oracle.h. The patches below only swap the RNG for a portable
# one, add counters, and dump the final tree; no arithmetic is touched.
#
# Ceres Solver is pinned too: 2.2.0 built from source against Eigen 3.4.0 with
# -ffp-contract=off and only Eigen's sparse Cholesky (no SuiteSparse, no Accelerate, no
# LAPACK), so every Kamiya solve is reproducible bit for bit (see src/kamiya_lm.jl). A
# packaged Ceres (e.g. Homebrew's) uses CHOLMOD and FMA contraction and cannot be matched.
#
# Requires: cmake, git, curl, a C++17 compiler, Boost (for DGtal). Ceres and Eigen are
# fetched and built here once. DGtal 1.4.2 is used from $DGtal_DIR if set, otherwise it
# is fetched and built here once too.
# Usage:    test/oracle/build.sh            -> test/oracle/build/bin/HctCco
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
COMMIT=9bc141f
SRC="$HERE/build/src"
B="$HERE/build"
NPROC="$(getconf _NPROCESSORS_ONLN)"
mkdir -p "$B"

# fetch URL FILE SHA256: download once into build/dl and verify the checksum.
fetch() {
  mkdir -p "$B/dl"
  [ -f "$B/dl/$2" ] || curl -sSLf -o "$B/dl/$2" "$1"
  echo "$3  $B/dl/$2" | shasum -a 256 -c - >/dev/null || { echo "checksum mismatch: $2" >&2; exit 1; }
}

# 0a. Eigen 3.4.0 (header only; installed so Ceres and HCT-CCO find the same copy).
EIGEN_DIR="$B/eigen-install/share/eigen3/cmake"
if [ ! -f "$EIGEN_DIR/Eigen3Config.cmake" ]; then
  fetch https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz eigen-3.4.0.tar.gz \
        8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
  rm -rf "$B/eigen-src" && mkdir -p "$B/eigen-src"
  tar xzf "$B/dl/eigen-3.4.0.tar.gz" -C "$B/eigen-src" --strip-components 1
  cmake -S "$B/eigen-src" -B "$B/eigen-build" -DCMAKE_INSTALL_PREFIX="$B/eigen-install" \
    -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 >/dev/null
  cmake --install "$B/eigen-build" >/dev/null
fi

# 0b. Ceres Solver 2.2.0. With SuiteSparse and Accelerate off, the default
# Solver::Options resolve to SPARSE_NORMAL_CHOLESKY + EIGEN_SPARSE (SimplicialLDLT).
# miniglog replaces glog; log level -1 keeps warnings/errors only, like glog's stderr default.
CERES_DIR="$B/ceres-install/lib/cmake/Ceres"
if [ ! -f "$CERES_DIR/CeresConfig.cmake" ]; then
  fetch http://ceres-solver.org/ceres-solver-2.2.0.tar.gz ceres-solver-2.2.0.tar.gz \
        48b2302a7986ece172898477c3bcd6deb8fb5cf19b3327bc49969aad4cede82d
  rm -rf "$B/ceres-src" && mkdir -p "$B/ceres-src"
  tar xzf "$B/dl/ceres-solver-2.2.0.tar.gz" -C "$B/ceres-src" --strip-components 1
  cmake -S "$B/ceres-src" -B "$B/ceres-build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-ffp-contract=off" -DCMAKE_INSTALL_PREFIX="$B/ceres-install" \
    -DEigen3_DIR="$EIGEN_DIR" -DMINIGLOG=ON -DMINIGLOG_MAX_LOG_LEVEL=-1 -DGFLAGS=OFF -DSUITESPARSE=OFF \
    -DACCELERATESPARSE=OFF -DEIGENSPARSE=ON -DEIGENMETIS=OFF -DUSE_CUDA=OFF -DLAPACK=OFF \
    -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF \
    -DPROVIDE_UNINSTALL_TARGET=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 >/dev/null
  cmake --build "$B/ceres-build" -j"$NPROC"
  cmake --install "$B/ceres-build" >/dev/null
fi
grep -q '^#define CERES_USE_EIGEN_SPARSE' "$B/ceres-install/include/ceres/internal/config.h" &&
grep -q '^#define CERES_NO_SUITESPARSE' "$B/ceres-install/include/ceres/internal/config.h" &&
grep -q '^#define CERES_NO_ACCELERATE_SPARSE' "$B/ceres-install/include/ceres/internal/config.h" ||
  { echo "pinned Ceres has the wrong sparse backend" >&2; exit 1; }

# 0c. DGtal 1.4.2, unless the caller points DGtal_DIR at an existing install.
if [ -z "${DGtal_DIR:-}" ]; then
  DGtal_DIR="$B/dgtal-install/lib/DGtal"
  if [ ! -f "$DGtal_DIR/DGtalConfig.cmake" ]; then
    [ -d "$B/dgtal-src/.git" ] || git clone -q --depth 1 --branch 1.4.2 https://github.com/DGtal-team/DGtal.git "$B/dgtal-src"
    cmake -S "$B/dgtal-src" -B "$B/dgtal-build" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$B/dgtal-install" -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
      -DEIGEN_ROOT="$B/eigen-install/include/eigen3" \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_POLICY_DEFAULT_CMP0167=OLD >/dev/null
    cmake --build "$B/dgtal-build" -j"$NPROC"
    cmake --install "$B/dgtal-build" >/dev/null
  fi
fi
if [ ! -d "$SRC/.git" ]; then
  git clone -q https://github.com/Sujal-vajire/HCT-CCO.git "$SRC"
fi
git -C "$SRC" checkout -q -f "$COMMIT"
git -C "$SRC" clean -qfd

# 1. Portable RNG (same RAND_MAX contract) seeded from $CCO_SEED.
sed -i.bak 's/rand()/cco_rand()/g' "$SRC/src/MyocardialDomain.h"
sed -i.bak 's/srand ((unsigned int) time(NULL));/cco_srand_env();/' "$SRC/src/helpers/CoronaryGrowth.h"
sed -i.bak 's/srand ((int) time(NULL));/cco_srand_env();/' "$SRC/bin/HctCco.cpp"

# Upstream sorts candidate segments with std::sort, whose order for exactly tied distances
# is implementation-defined (libc++ and libstdc++ differ). Use stable_sort so ties break by
# segment index on every platform; Vasculature.jl does the same.
sed -i.bak 's/std::sort(allSegments.begin(), allSegments.end(),/std::stable_sort(allSegments.begin(), allSegments.end(),/' "$SRC/src/HctccoModel.ih"

# 2. Work counters and a coarse per-step clock.
G="$SRC/src/helpers/CoronaryGrowth.h"
sed -i.bak 's/DGtal::trace.progressBar(i, nbSeed);/if (i % 250 == 0) std::fprintf(stderr, "[cco] step %u t %.4f segs %zu\\n", i, cco_oracle::now(), aTree.myVectSegments.size());/' "$G"
sed -i.bak 's/auto pt = aTree.generateNewLocation(nbTryCandite);/auto pt = aTree.generateNewLocation(nbTryCandite); CCO_COUNT("location");/' "$G"
sed -i.bak 's/auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it));/auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it)); CCO_COUNT("neighbor");/' "$G"
sed -i.bak 's/HctccoModel< DomCtr, TDim  > cTree1 = aTree;/CCO_COUNT("isAddable"); HctccoModel< DomCtr, TDim  > cTree1 = aTree;/' "$G"
sed -i.bak 's/        aTree.updateRootRadius();/        aTree.updateRootRadius(); cco_trace_step(aTree, i, nbT, volOpt);/' "$G"
sed -i.bak 's/if (nbSeed != nbSeedFound){/cco_depth_stats(aTree); cco_dump_tree(aTree); cco_oracle::report(); if (nbSeed != nbSeedFound){/' "$G"
M="$SRC/src/HctccoModel.ih"
sed -i.bak 's/res1 = kamyiaOptimization(pCurrent/CCO_COUNT("kamiya_iter"); res1 = kamyiaOptimization(pCurrent/' "$M"
sed -i.bak 's/isDone = true;/isDone = true; CCO_COUNT("accepted");/' "$M"
C="$SRC/src/helpers/CoronaryGeometry.h"
sed -i.bak 's/return summary.IsSolutionUsable();/CCO_COUNT("ceres_solve"); CCO_COUNT_ADD("ceres_iter", summary.iterations.size()); cco_ceres_record(gamma, deltaP1, deltaP2, f0, f1, f2, l0, l1, l2, x1_in, x2_in, xx1, xx2, summary.IsSolutionUsable(), summary.iterations.size(), (int)summary.termination_type); return summary.IsSolutionUsable();/' "$C"
sed -i.bak 's/CostOptPos \*f = new CostOptPos();/const double x1_in = xx1, x2_in = xx2; CostOptPos *f = new CostOptPos();/' "$C"

# Every patch must have applied exactly; fail loudly if upstream moved.
for pat in 'cco_rand()' 'cco_srand_env' 'CCO_COUNT("location")' 'CCO_COUNT("neighbor")' \
           'CCO_COUNT("isAddable")' 'cco_dump_tree' 'CCO_COUNT("kamiya_iter")' \
           'CCO_COUNT("accepted")' 'CCO_COUNT("ceres_solve")' 'cco_ceres_record' 'x1_in = xx1' 'std::stable_sort(allSegments' 'cco_trace_step(aTree' '\[cco\] step'; do
  grep -rqF -- "${pat//\\/}" "$SRC/src" "$SRC/bin" || { echo "patch failed: $pat" >&2; exit 1; }
done

# Release build without FMA contraction so float results match Julia's IEEE ops.
# DGtal puts the package prefix (e.g. /opt/homebrew/include) on the -I path, which would
# shadow the pinned Ceres/miniglog/Eigen headers; upstream includes them with quotes, so
# -iquote makes the pinned copies win.
PIN_INC="-iquote $B/ceres-install/include -iquote $B/ceres-install/include/ceres/internal/miniglog -iquote $B/eigen-install/include/eigen3"
cmake -S "$SRC" -B "$HERE/build/cmake" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-ffp-contract=off $PIN_INC -include $HERE/cco_oracle.h" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCeres_DIR="$CERES_DIR" -DEigen3_DIR="$EIGEN_DIR" \
  -DDGtal_DIR="$DGtal_DIR" >/dev/null
cmake --build "$HERE/build/cmake" -j"$NPROC"
mkdir -p "$HERE/build/bin"
cp "$HERE/build/cmake/bin/HctCco" "$HERE/build/bin/HctCco"
# The generator must have compiled against the pinned headers only.
if grep -rhoE '/[^ ]*/include/(ceres|eigen3)/[^ ]*' "$HERE/build/cmake/bin/CMakeFiles/HctCco.dir" 2>/dev/null |
   grep -v "^$B/" | grep -q .; then
  echo "HctCco picked up a non-pinned Ceres/Eigen header" >&2; exit 1
fi
# Standalone Kamiya/Ceres reference (synthetic records for test/test_kamiya.jl).
c++ -std=c++17 -O3 -DNDEBUG -ffp-contract=off -I"$B/ceres-install/include" \
  -I"$B/ceres-install/include/ceres/internal/miniglog" -I"$B/eigen-install/include/eigen3" \
  "$HERE/kamiya_ref.cpp" "$B/ceres-install/lib/libceres.a" -o "$B/bin/kamiya_ref"
echo "built $HERE/build/bin/HctCco $HERE/build/bin/kamiya_ref"
