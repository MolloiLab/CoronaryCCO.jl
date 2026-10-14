// SPDX-License-Identifier: GPL-3.0-or-later
// Instrumentation force-included (-include) into the upstream HCT-CCO build so it can
// serve as the reference ("oracle") for CoronaryCCO.jl. It changes no arithmetic:
//   1. rand()/srand() are replaced by a portable SplitMix64-based generator with the
//      same RAND_MAX contract, seeded from CCO_SEED, so Julia can replay the exact
//      candidate stream on any platform.
//   2. Lightweight counters record how much work each growth step does.
//   3. The final tree can be dumped with hex floats (CCO_DUMP=path) for bitwise diffs.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <chrono>

namespace cco_oracle {
inline uint64_t &state() { static uint64_t s = 0; return s; }
inline std::map<std::string, unsigned long long> &counters() {
  static std::map<std::string, unsigned long long> c; return c;
}
inline void report() {
  for (auto &kv : counters()) std::fprintf(stderr, "[cco] %s %llu\n", kv.first.c_str(), kv.second);
}
inline double now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}
}  // namespace cco_oracle

// SplitMix64, top 31 bits -> [0, 2^31 - 1]; RAND_MAX is 2^31 - 1 on glibc and macOS.
inline int cco_rand() {
  uint64_t z = (cco_oracle::state() += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  z = z ^ (z >> 31);
  return static_cast<int>(z >> 33);
}
inline void cco_srand_env() {
  const char *s = std::getenv("CCO_SEED");
  cco_oracle::state() = s ? std::strtoull(s, nullptr, 10) : 0;
}
// Records every Kamiya/Ceres solve when CCO_CERES_DUMP=path is set: inputs, outputs,
// usability and iteration count, all as exact hex floats.
inline FILE *ceres_dump_file() {
  static FILE *f = [] {
    const char *p = std::getenv("CCO_CERES_DUMP");
    return p ? std::fopen(p, "w") : static_cast<FILE *>(nullptr);
  }();
  return f;
}
inline void cco_ceres_record(double gamma, double dp1, double dp2, double f0, double f1,
                             double f2, double l0, double l1, double l2, double x1in,
                             double x2in, double x1out, double x2out, bool usable,
                             size_t iters, int termination) {
  FILE *f = ceres_dump_file();
  if (!f) return;
  std::fprintf(f, "%a %a %a %a %a %a %a %a %a %a %a %a %a %d %zu %d\n", gamma, dp1, dp2, f0,
               f1, f2, l0, l1, l2, x1in, x2in, x1out, x2out, usable ? 1 : 0, iters,
               termination);
}

#define CCO_COUNT(name) (cco_oracle::counters()[name] += 1)
#define CCO_COUNT_ADD(name, n) (cco_oracle::counters()[name] += (n))

// Writes every segment with %a (exact hex) floats: index parent left right kterm
// x y z radius resistance flow beta.
template <class TTree>
void cco_dump_tree(const TTree &t) {
  const char *path = std::getenv("CCO_DUMP");
  if (!path) return;
  FILE *f = std::fopen(path, "w");
  if (!f) return;
  for (size_t i = 0; i < t.myVectSegments.size(); i++) {
    const auto &s = t.myVectSegments[i];
    std::fprintf(f, "%zu %u %u %u %u %a %a %a %a %a %a %a\n", i, t.myVectParent[i],
                 t.myVectChildren[i].first, t.myVectChildren[i].second, s.myKTerm,
                 s.myCoordinate[0], s.myCoordinate[1], s.myCoordinate[2], s.myRadius,
                 s.myResistance, s.myFlow, s.myBeta);
  }
  std::fclose(f);
}

// Per-step trace (CCO_TRACE=path): step, locations tried, best volume, segment count, and the
// last two segments' distal points (new left = old distal of the split segment, new terminal).
template <class TTree>
void cco_trace_step(const TTree &t, unsigned step, unsigned nbT, double volOpt) {
  static FILE *f = [] {
    const char *p = std::getenv("CCO_TRACE");
    return p ? std::fopen(p, "w") : static_cast<FILE *>(nullptr);
  }();
  if (!f) return;
  size_t n = t.myVectSegments.size();
  const auto &a = t.myVectSegments[n - 1].myCoordinate;
  unsigned s = t.myVectParent[n - 1];
  const auto &b = t.myVectSegments[s].myCoordinate;
  std::fprintf(f, "%u %u %a %zu %u %a %a %a %a %a %a\n", step, nbT, volOpt, n, s, a[0], a[1], a[2], b[0], b[1], b[2]);
}

// Depth statistics of the final tree (deepest root-to-leaf path, in segments).
template <class TTree>
void cco_depth_stats(const TTree &t) {
  size_t n = t.myVectSegments.size();
  std::vector<unsigned> depth(n, 0);
  unsigned maxd = 0; double sum = 0; size_t nt = 0;
  for (size_t i = 1; i < n; i++) {  // parents always precede children? not guaranteed
    unsigned d = 0; size_t j = i;
    while (j != 0) { j = t.myVectParent[j]; d++; }
    depth[i] = d; if (d > maxd) maxd = d;
  }
  for (auto ti : t.myVectTerminals) { sum += depth[ti]; nt++; }
  std::fprintf(stderr, "[cco] segments %zu terminals %zu max_depth %u mean_terminal_depth %.2f\n",
               n, nt, maxd, nt ? sum / nt : 0.0);
}
