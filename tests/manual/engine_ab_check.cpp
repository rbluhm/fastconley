// A/B bit-identity check of src/conley_core.h against a previous revision of
// the same header, for engine refactors that must not change numbers. Both
// headers are compiled into this one program (the old copy with its
// namespace renamed to conley_old) and spatial_meat / serial_hac_meat are
// compared bit for bit on randomized configurations: every kernel and
// distance, k = 1..20, zero / tiny / capped / antipodal cutoffs, duplicate
// locations, lattice points exactly at the cutoff, balanced (double and
// float weights) and band paths, 1 vs 4 threads, and the row-map entry.
// Run it through tests/manual/engine-ab-check.sh, which prepares the old
// header from a git revision:
//
//   tests/manual/engine-ab-check.sh [REV] [ITERATIONS]
#include "conley_core.h"
#include "conley_core_old.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace {

struct Rng {
  std::mt19937_64 e;
  explicit Rng(std::uint64_t seed) : e(seed) {}
  double u() { return static_cast<double>(e() >> 11) * (1.0 / 9007199254740992.0); }
  double u(double lo, double hi) { return lo + (hi - lo) * u(); }
  std::size_t i(std::size_t n) { return static_cast<std::size_t>(e() % n); }
};

int n_cases = 0, n_failed = 0;

void compare(const std::string& name, const arma::mat& a, const arma::mat& b) {
  ++n_cases;
  const bool same = a.n_rows == b.n_rows && a.n_cols == b.n_cols &&
                    std::memcmp(a.memptr(), b.memptr(), sizeof(double) * a.n_elem) == 0;
  if (!same) {
    ++n_failed;
    std::printf("DIFF  %s  max|d| = %.3g\n", name.c_str(), arma::abs(a - b).max());
  }
}

}  // namespace

int main(int argc, char** argv) {
  const int iterations = argc > 1 ? std::atoi(argv[1]) : 300;
  const char* kernels[2] = {"bartlett", "uniform"};
  const char* dists[3] = {"haversine", "spherical", "chord"};
  const double cutoffs[8] = {0.0, 0.001, 5.0, 60.0, 300.0, 1500.0, 13000.0, 25000.0};
  Rng rng(12345);
  for (int it = 0; it < iterations; ++it) {
    // 0 box, 1 global, 2 lattice, 3 clustered duplicates, 4 tiny n
    const int shape = static_cast<int>(rng.i(5));
    const std::size_t n_per = shape == 4 ? 1 + rng.i(40)
        : (rng.i(4) == 0 ? 60000 + rng.i(90000) : 200 + rng.i(3000));
    const std::size_t k = 1 + rng.i(20);
    const int T = rng.i(3) == 0 ? 1 + static_cast<int>(rng.i(4)) : 1;
    double cutoff = cutoffs[rng.i(8)];
    const std::string kernel = kernels[rng.i(2)], dist = dists[rng.i(3)];
    const int ncores = rng.i(2) ? 4 : 1;
    const bool band = rng.i(6) == 0, wfloat = rng.i(3) == 0, use_map = rng.i(2) == 0;
    // Keep the pair count (and a balanced CSR) bounded for large n.
    if (n_per > 20000) {
      const double small[5] = {0.0, 0.001, 5.0, 20.0, 60.0};
      cutoff = small[rng.i(5)];
      if (shape == 3 && cutoff > 5.0) cutoff = 5.0;
    }

    arma::vec la(n_per), lo(n_per);
    const double step = 0.25;
    for (std::size_t i = 0; i < n_per; ++i) {
      if (shape == 0 || shape == 4) { la[i] = rng.u(25, 49); lo[i] = rng.u(-124, -67); }
      else if (shape == 1) { la[i] = rng.u(-89.9, 89.9); lo[i] = rng.u(-180, 360); }
      else if (shape == 2) { la[i] = 30 + step * static_cast<double>(rng.i(80)); lo[i] = -10 + step * static_cast<double>(rng.i(120)); }
      else {
        const std::size_t c = rng.i(1 + n_per / 7);
        la[i] = 40 + 0.01 * static_cast<double>(c % 97);
        lo[i] = 5 + 0.013 * static_cast<double>(c / 97);
      }
    }
    if (shape == 2 && rng.i(2)) {
      // A cutoff equal to a lattice great-circle distance puts pairs on the
      // accept boundary, where Bartlett weights round to (near) zero.
      const double d_lat = step * static_cast<double>(1 + rng.i(6)) * conley::DE2RA;
      cutoff = 2.0 * conley::AVG_ERAD * std::asin(std::sin(d_lat / 2.0));
    }

    const std::size_t n = n_per * static_cast<std::size_t>(T);
    arma::vec lat(n), lon(n), time(n);
    arma::mat S(n, k);
    for (int t = 0; t < T; ++t) {
      for (std::size_t i = 0; i < n_per; ++i) {
        const std::size_t r = static_cast<std::size_t>(t) * n_per + i;
        lat[r] = la[i]; lon[r] = lo[i]; time[r] = t + 1;
      }
    }
    for (std::size_t j = 0; j < k; ++j) {
      for (std::size_t i = 0; i < n; ++i) S(i, j) = rng.u(-1, 1) * (rng.i(50) == 0 ? 0.0 : 1.0);
    }
    const bool balanced = T > 1 && rng.i(3) != 0;
    char name[256];
    std::snprintf(name, sizeof name, "it%d shape%d n=%zu T=%d k=%zu cutoff=%g %s/%s nc%d%s%s%s",
                  it, shape, n, T, k, cutoff, kernel.c_str(), dist.c_str(), ncores,
                  balanced ? " balanced" : "", band ? " band" : "", wfloat ? " float" : "");
    const std::string neighbor = band ? "band" : "grid", csr = wfloat ? "float" : "double";

    const arma::mat old_meat = conley_old::spatial_meat(lat, lon, time, S, cutoff, kernel, dist,
                                                        balanced, ncores, neighbor, csr);
    compare(name, old_meat, conley::spatial_meat(lat, lon, time, S, cutoff, kernel, dist,
                                                 balanced, ncores, neighbor, csr));
    if (use_map) {
      // Scores stored in shuffled rows and reached through the row map.
      std::vector<std::size_t> map(n);
      std::iota(map.begin(), map.end(), 0);
      for (std::size_t i = n; i > 1; --i) std::swap(map[i - 1], map[rng.i(i)]);
      arma::mat stored(n, k);
      for (std::size_t i = 0; i < n; ++i) stored.row(map[i]) = S.row(i);
      compare(std::string(name) + " [row map]", old_meat,
              conley::spatial_meat(lat, lon, time, stored, cutoff, kernel, dist, balanced,
                                   ncores, neighbor, csr, nullptr, &map));
    }
    if (it % 10 == 0 && T > 1) {
      arma::vec unit(n);
      for (std::size_t i = 0; i < n; ++i) unit[i] = static_cast<double>(i % n_per);
      const arma::uvec o = arma::sort_index(unit * 1000.0 + time);
      const arma::vec u = unit.elem(o), tt = time.elem(o);
      const arma::mat So = S.rows(o);
      compare(std::string(name) + " [serial]", conley_old::serial_hac_meat(u, tt, 2, So, ncores),
              conley::serial_hac_meat(u, tt, 2, So, ncores));
    }
  }
  std::printf("engine_ab_check: %d cases, %d bit-identical, %d differ\n",
              n_cases, n_cases - n_failed, n_failed);
  return n_failed ? 1 : 0;
}
