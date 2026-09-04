// bench_engine.cpp -- standalone timing of the engine (no Stata, no R).
// Times conley::spatial_meat on 100,000 seeded points in a CONUS-like box
// with a 500 km cutoff, Bartlett/haversine and uniform/spherical, at the
// thread counts given on the command line. Use it to compare toolchains on
// the same machine (e.g. mingw-w64 vs clang-cl/MSVC on Windows):
//
//   make -C stata/plugin bench            (Linux / MSYS2 mingw: builds bench_engine)
//   ./bench_engine 1 4                     -> seconds per kernel at 1 and 4 threads
//
// clang-cl / MSVC (from a developer prompt, Armadillo headers in ARMA_INC):
//   clang-cl /O2 /EHsc /std:c++17 /DNDEBUG /DARMA_64BIT_WORD /DARMA_DONT_USE_WRAPPER
//     /DARMA_DONT_USE_BLAS /DARMA_DONT_USE_LAPACK /DARMA_DONT_USE_SUPERLU
//     /I..\..\src /I%ARMA_INC% bench_engine.cpp /Fe:bench_engine_clang.exe
#include "conley_core.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>

int main(int argc, char** argv) {
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> ulat(25.0, 49.0), ulon(-124.0, -67.0);
  std::normal_distribution<double> nrm(0.0, 1.0);
  const std::size_t n = 100000, k = 4;
  arma::vec lat(n), lon(n), time(n, arma::fill::ones);
  arma::mat S(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    lat[i] = ulat(rng); lon[i] = ulon(rng);
    for (std::size_t j = 0; j < k; ++j) S(i, j) = nrm(rng);
  }
  const char* kernels[2] = { "bartlett", "uniform" };
  const char* dists[2] = { "haversine", "spherical" };
  for (int a = 1; a < (argc > 1 ? argc : 2); ++a) {
    const int nc = argc > 1 ? std::atoi(argv[a]) : 1;
    for (int kk = 0; kk < 2; ++kk) {
      double best = 1e9, chk = 0.0;
      for (int rep = 0; rep < 2; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        const arma::mat M = conley::spatial_meat(lat, lon, time, S, 500.0, kernels[kk], dists[kk],
                                                 false, nc, "grid", "double");
        const auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration<double>(t1 - t0).count();
        if (sec < best) best = sec;
        chk = M(0, 0);
      }
      std::printf("threads=%2d %-9s/%-9s  %.3f s   (M[0,0] = %.10g)\n", nc, kernels[kk], dists[kk], best, chk);
    }
  }
  return 0;
}
