// bench_engine.cpp -- standalone timing of the engine (no Stata, no R).
// Times conley::spatial_meat on 100,000 seeded points in a CONUS-like box
// with a 500 km cutoff, Bartlett/haversine and uniform/spherical, at the
// thread counts given on the command line, and prints two checksums.
//
//   make -C stata/plugin bench       (Linux / MSYS2 mingw: builds and runs "1 4")
//   ./bench_engine 1 4 16
//
// The input generator is STL-independent (raw std::mt19937_64 output only;
// std::uniform_real_distribution and std::normal_distribution are
// implementation-defined), so the checksums are comparable across
// libstdc++ and the MSVC STL. Build against another engine header with
// -DCORE_HEADER='"path.h"'. clang-cl / MSVC example (Armadillo headers in
// ARMA_INC):
//   cl /O2 /EHsc /std:c++17 /MT /DNDEBUG /DARMA_64BIT_WORD /DARMA_DONT_USE_WRAPPER
//     /DARMA_DONT_USE_BLAS /DARMA_DONT_USE_LAPACK /DARMA_DONT_USE_SUPERLU
//     /I..\..\src /I%ARMA_INC% bench_engine.cpp /Fe:bench_engine_msvc.exe
// (Contributed by the Windows test session, 2026-09-04.)
#ifndef CORE_HEADER
#define CORE_HEADER "conley_core.h"
#endif
#include CORE_HEADER

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>

struct DetRng {
  std::mt19937_64 e;
  DetRng() : e(7) {}
  double u() { return static_cast<double>(e() >> 11) * (1.0 / 9007199254740992.0); }  // [0,1)
  double u(double lo, double hi) { return lo + (hi - lo) * u(); }
};

int main(int argc, char** argv) {
  DetRng rng;
  const std::size_t n = 100000, k = 4;
  arma::vec lat(n), lon(n), time(n, arma::fill::ones);
  arma::mat S(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    lat[i] = rng.u(25.0, 49.0); lon[i] = rng.u(-124.0, -67.0);
    for (std::size_t j = 0; j < k; ++j) S(i, j) = rng.u(-1.0, 1.0) * 1.7320508075688772;
  }
  const char* kernels[2] = { "bartlett", "uniform" };
  const char* dists[2] = { "haversine", "spherical" };
  for (int a = 1; a < (argc > 1 ? argc : 2); ++a) {
    const int nc = argc > 1 ? std::atoi(argv[a]) : 1;
    for (int kk = 0; kk < 2; ++kk) {
      double best = 1e9, chk = 0.0, chk2 = 0.0;
      for (int rep = 0; rep < 2; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        const arma::mat M = conley::spatial_meat(lat, lon, time, S, 500.0, kernels[kk], dists[kk],
                                                 false, nc, "grid", "double");
        const auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration<double>(t1 - t0).count();
        if (sec < best) best = sec;
        chk = M(0, 0); chk2 = M(1, 2);
      }
      std::printf("threads=%2d %-9s/%-9s  %.3f s   (M[0,0] = %.15g  M[1,2] = %.15g)\n", nc, kernels[kk], dists[kk], best, chk, chk2);
    }
  }
  return 0;
}
