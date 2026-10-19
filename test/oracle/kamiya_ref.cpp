// SPDX-License-Identifier: GPL-3.0-or-later
// Standalone reference for upstream's kamyiaOpt Ceres call (CostOptPos and the solver setup
// are copied from HCT-CCO src/helpers/CoronaryGeometry.h), linked against the pinned Ceres
// 2.2.0 that oracle/build.sh builds. It produces the synthetic reference records that
// test/test_kamiya.jl checks src/kamiya_lm.jl against. Modes:
//   kamiya_ref replay < cases.txt   re-solve cases (first 11 fields of each line used)
//   kamiya_ref fuzz N SEED          N random cases near the recorded operating range
//   kamiya_ref wild N SEED          N extreme cases (signs, tiny/huge magnitudes)
//   kamiya_ref trace < cases.txt    per-iteration summaries (hex) and FullReport per case
// Output lines use the cco_ceres_record format: gamma dP1 dP2 f0 f1 f2 l0 l1 l2 x1_in x2_in
// x1_out x2_out (C99 hex), usable, iterations, termination_type.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include "ceres/ceres.h"

using namespace ceres;

struct CostOptPos {
  double deltap1, deltap2, f0, f1, f2, l0, l1, l2, gamma;
  template <typename T>
  bool operator()(const T* const x1, const T* const x2, T* residual) const {
    residual[0] =  deltap1*x1[0]*x1[0]*(pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0)/f2)),2.0/3.0))
    - (f0*l0*x1[0]*x1[0]) - f1 * l1 * (pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0))/f2),2.0/3.0));
    residual[1] =  deltap2*x2[0]*x2[0]*(pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0)/f2)),2.0/3.0))
    - (f0*l0*x2[0]*x2[0]) - f2 * l2 * (pow(f0*((pow(x1[0],(3.0+gamma)/2.0)/f1) + (pow(x2[0],(3.0+gamma)/2.0))/f2),2.0/3.0));
    return true;
  }
};

struct Tracer : public IterationCallback {
  CallbackReturnType operator()(const IterationSummary& s) override {
    std::printf("it %d succ %d valid %d cost %a dcost %a gmax %a step %a rho %a radius %a\n",
                s.iteration, (int)s.step_is_successful, (int)s.step_is_valid, s.cost,
                s.cost_change, s.gradient_max_norm, s.step_norm, s.relative_decrease,
                s.trust_region_radius);
    return SOLVER_CONTINUE;
  }
};

static void solve(const double* in, bool trace) {
  double xx1 = in[9], xx2 = in[10];
  auto* f = new CostOptPos();
  f->gamma = in[0]; f->deltap1 = in[1]; f->deltap2 = in[2];
  f->f0 = in[3]; f->f1 = in[4]; f->f2 = in[5];
  f->l0 = in[6]; f->l1 = in[7]; f->l2 = in[8];
  Problem problem;
  problem.AddResidualBlock(new AutoDiffCostFunction<CostOptPos, 2, 1, 1>(f), nullptr, &xx1, &xx2);
  Solver::Options options;
  options.minimizer_progress_to_stdout = false;
  Tracer tracer;
  if (trace) options.callbacks.push_back(&tracer);
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  if (trace) std::printf("%s\n", summary.FullReport().c_str());
  std::printf("%a %a %a %a %a %a %a %a %a %a %a %a %a %d %zu %d\n", in[0], in[1], in[2], in[3],
              in[4], in[5], in[6], in[7], in[8], in[9], in[10], xx1, xx2,
              summary.IsSolutionUsable() ? 1 : 0, summary.iterations.size(),
              (int)summary.termination_type);
}

static bool read_case(double* in) {
  char line[4096];
  if (!std::fgets(line, sizeof line, stdin)) return false;
  char* p = line;
  for (int i = 0; i < 11; ++i) in[i] = std::strtod(p, &p);
  return true;
}

int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "replay";
  double in[11];
  if (mode == "replay" || mode == "trace") {
    while (read_case(in)) solve(in, mode == "trace");
  } else if (mode == "fuzz") {
    const long n = std::atol(argv[2]);
    std::mt19937_64 rng(std::strtoull(argv[3], nullptr, 10));
    std::uniform_real_distribution<double> u(0.0, 1.0);
    auto lu = [&](double lo, double hi) { return std::exp(std::log(lo) + u(rng) * (std::log(hi) - std::log(lo))); };
    for (long k = 0; k < n; ++k) {
      // Realistic-ish ranges around the recorded calls, plus wide outliers.
      const bool wild = u(rng) < 0.2;
      in[0] = u(rng) < 0.7 ? 3.0 : 2.0 + 1.5 * u(rng);           // gamma
      in[1] = wild ? lu(1e-3, 1e9) : lu(1e4, 1e7);                // deltaP1
      in[2] = wild ? lu(1e-3, 1e9) : lu(1e4, 1e7);                // deltaP2
      in[3] = lu(0.5, 50.0);                                       // f0
      in[4] = in[3] * (wild ? lu(1e-6, 1.5) : u(rng));             // f1
      in[5] = wild ? lu(1e-6, 1.5) * in[3] : in[3] - in[4];        // f2
      if (in[5] <= 0) in[5] = in[3] * 1e-3;
      in[6] = wild ? lu(1e-4, 1e3) : lu(0.1, 10.0);                // l0
      in[7] = wild ? lu(1e-4, 1e3) : lu(0.1, 10.0);                // l1
      in[8] = wild ? lu(1e-4, 1e3) : lu(0.1, 10.0);                // l2
      in[9] = wild ? lu(1e-6, 10.0) : lu(1e-3, 3e-2);              // x1
      in[10] = wild ? lu(1e-6, 10.0) : (u(rng) < 0.3 ? in[9] : lu(1e-3, 3e-2));  // x2
      if (wild && u(rng) < 0.1) in[9] = -in[9];                    // negative start
      if (wild && u(rng) < 0.05) in[10] = 0.0;
      solve(in, false);
    }
  }
  if (mode == "wild") {
    // Extreme inputs (signs, tiny/huge magnitudes) to drive the failure paths:
    // invalid steps, failed factorizations, failed evaluations mid-run.
    const long n = std::atol(argv[2]);
    std::mt19937_64 rng(std::strtoull(argv[3], nullptr, 10));
    std::uniform_real_distribution<double> u(0.0, 1.0);
    auto lu = [&](double lo, double hi) { return std::exp(std::log(lo) + u(rng) * (std::log(hi) - std::log(lo))); };
    auto sg = [&](double p) { return u(rng) < p ? -1.0 : 1.0; };
    for (long k = 0; k < n; ++k) {
      in[0] = u(rng) < 0.5 ? 3.0 : -2.0 + 8.0 * u(rng);
      in[1] = sg(0.1) * lu(1e-30, 1e30);
      in[2] = sg(0.1) * lu(1e-30, 1e30);
      in[3] = sg(0.05) * lu(1e-20, 1e20);
      in[4] = sg(0.05) * lu(1e-20, 1e20);
      in[5] = sg(0.05) * lu(1e-20, 1e20);
      in[6] = sg(0.05) * lu(1e-20, 1e20);
      in[7] = sg(0.05) * lu(1e-20, 1e20);
      in[8] = sg(0.05) * lu(1e-20, 1e20);
      in[9] = sg(0.05) * lu(1e-12, 1e6);
      in[10] = sg(0.05) * lu(1e-12, 1e6);
      solve(in, false);
    }
  }
  return 0;
}
