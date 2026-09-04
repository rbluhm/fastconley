// Standalone (no R) build check for src/conley_core.h -- the compile
// path the Stata plugin uses. Writes a random cross-section + balanced
// panel and the engine's spatial / serial / grid meats to OUT_DIR so
// tests/manual/core_standalone_check.R can compare them with the R
// package bit for bit. Build from the package root with:
//
//   g++ -std=c++14 -O2 -pthread -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS
//       -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU -Isrc
//       tests/manual/core_standalone_check.cpp -o .scratch/core_standalone_check
//   .scratch/core_standalone_check OUT_DIR
//
// (needs the Armadillo headers, e.g. libarmadillo-dev or
// RcppArmadillo/include on the include path).
#include "conley_core.h"

#include <cstdio>
#include <fstream>
#include <random>
#include <string>

static void write_mat(const std::string& path, const arma::mat& M) {
  std::ofstream f(path.c_str());
  f.precision(17);
  for (arma::uword i = 0; i < M.n_rows; ++i) {
    for (arma::uword j = 0; j < M.n_cols; ++j) {
      f << (j ? "," : "") << M(i, j);
    }
    f << "\n";
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: core_standalone_check OUT_DIR\n");
    return 2;
  }
  const std::string out = std::string(argv[1]) + "/";
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> ulat(25.0, 49.0), ulon(-124.0, -67.0);
  std::normal_distribution<double> nrm(0.0, 1.0);

  // Cross-section: n points, k = 3 scores, single time block.
  const std::size_t n = 5000, k = 3;
  arma::vec lat(n), lon(n), time(n, arma::fill::ones);
  arma::mat S(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    lat[i] = ulat(rng); lon[i] = ulon(rng);
    for (std::size_t j = 0; j < k; ++j) S(i, j) = nrm(rng);
  }
  write_mat(out + "cs_input.csv", arma::join_rows(arma::join_rows(lat, lon), S));
  for (int nc = 1; nc <= 4; nc += 3) {
    write_mat(out + "cs_bartlett_haversine_nc" + std::to_string(nc) + ".csv",
              conley::spatial_meat(lat, lon, time, S, 300.0, "bartlett",
                                   "haversine", false, nc, "grid", "double"));
    write_mat(out + "cs_uniform_chord_band_nc" + std::to_string(nc) + ".csv",
              conley::spatial_meat(lat, lon, time, S, 300.0, "uniform",
                                   "chord", false, nc, "band", "double"));
  }

  // Balanced panel: n_u units x T periods, sorted by (time, unit); the
  // serial meat wants (unit, time) order, so build both layouts.
  const std::size_t n_u = 800, T = 3;
  arma::vec plat(n_u * T), plon(n_u * T), ptime(n_u * T);
  arma::mat PS(n_u * T, k);
  std::vector<double> ulat_v(n_u), ulon_v(n_u);
  for (std::size_t u = 0; u < n_u; ++u) { ulat_v[u] = ulat(rng); ulon_v[u] = ulon(rng); }
  for (std::size_t t = 0; t < T; ++t) {
    for (std::size_t u = 0; u < n_u; ++u) {
      const std::size_t r = t * n_u + u;
      plat[r] = ulat_v[u]; plon[r] = ulon_v[u]; ptime[r] = static_cast<double>(t + 1);
      for (std::size_t j = 0; j < k; ++j) PS(r, j) = nrm(rng);
    }
  }
  write_mat(out + "pnl_input.csv",
            arma::join_rows(arma::join_rows(plat, plon), arma::join_rows(ptime, PS)));
  write_mat(out + "pnl_bartlett_spherical_balanced_nc4.csv",
            conley::spatial_meat(plat, plon, ptime, PS, 400.0, "bartlett",
                                 "spherical", true, 4, "grid", "double"));
  // (unit, time) layout for the serial meat.
  arma::vec unit_s(n_u * T), time_s(n_u * T);
  arma::mat S_s(n_u * T, k);
  for (std::size_t u = 0; u < n_u; ++u) {
    for (std::size_t t = 0; t < T; ++t) {
      const std::size_t r = u * T + t;
      unit_s[r] = static_cast<double>(u + 1); time_s[r] = static_cast<double>(t + 1);
      S_s.row(r) = PS.row(t * n_u + u);
    }
  }
  write_mat(out + "pnl_serial_lag2_nc4.csv",
            conley::serial_hac_meat(unit_s, time_s, 2.0, S_s, 4));

  // Raster: 40 rings x 60 cols at 0.5 degrees, single block, no holes.
  const int R = 40, C = 60;
  const std::size_t ng = static_cast<std::size_t>(R * C);
  std::vector<int> ring(ng), col(ng);
  arma::vec gtime(ng, arma::fill::ones);
  arma::mat GS(ng, k);
  for (int r = 0; r < R; ++r) {
    for (int c = 0; c < C; ++c) {
      const std::size_t i = static_cast<std::size_t>(r * C + c);
      ring[i] = r; col[i] = c;
      for (std::size_t j = 0; j < k; ++j) GS(i, j) = nrm(rng);
    }
  }
  write_mat(out + "grid_input.csv", GS);
  write_mat(out + "grid_uniform_spherical_nc4.csv",
            conley::grid_meat(ring.data(), col.data(), gtime, GS, 35.0, 0.5, 0.5,
                              R, C, 0, 250.0, "spherical", "uniform", 4));
  write_mat(out + "grid_bartlett_haversine_nc4.csv",
            conley::grid_meat(ring.data(), col.data(), gtime, GS, 35.0, 0.5, 0.5,
                              R, C, 0, 250.0, "haversine", "bartlett", 4));

  // Error path: the engine must throw conley::Error, never abort.
  try {
    conley::spatial_meat(lat, lon, time, S, 300.0, "gaussian", "haversine",
                         false, 1, "grid", "double");
    std::fprintf(stderr, "expected conley::Error\n");
    return 1;
  } catch (const conley::Error& e) {
    std::printf("error path OK: %s\n", e.what());
  }
  std::printf("wrote outputs to %s\n", out.c_str());
  return 0;
}
