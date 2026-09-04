// Focused standalone probes for conley_core.h boundary behavior.
// Compile from the package root with the same C++14/header-only Armadillo
// flags used by tests/manual/core_standalone_check.R.
#include "conley_core.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error(message);
}

void require_close(double got, double want, double tol,
                   const std::string& message) {
  if (!std::isfinite(got) || std::fabs(got - want) > tol) {
    throw std::runtime_error(message + ": got " + std::to_string(got) +
                             ", want " + std::to_string(want));
  }
}

bool bit_identical(const arma::mat& a, const arma::mat& b) {
  return a.n_rows == b.n_rows && a.n_cols == b.n_cols &&
         std::memcmp(a.memptr(), b.memptr(),
                     static_cast<std::size_t>(a.n_elem) * sizeof(double)) == 0;
}

template <typename F>
void expect_conley_error(const char* name, const F& fn) {
  try {
    fn();
  } catch (const conley::Error& e) {
    std::printf("throw %-24s %s\n", name, e.what());
    return;
  }
  throw std::runtime_error(std::string(name) + " did not throw conley::Error");
}

int interrupt_polls = 0;
bool interrupt_now() {
  ++interrupt_polls;
  return true;
}

}  // namespace

int main() {
  try {
    const char* distances[] = {"haversine", "spherical", "chord"};
    const char* kernels[] = {"bartlett", "uniform"};

    // Exact duplicate observations at cutoff zero retain their cross term:
    // [1,2]'[1,2] has scalar meat 9.
    arma::vec dup_lat(2), dup_lon(2), dup_time(2, arma::fill::ones);
    arma::mat dup_s(2, 1);
    dup_lat.fill(0.0);
    dup_lon.fill(-124.0);
    dup_s(0, 0) = 1.0;
    dup_s(1, 0) = 2.0;
    for (const char* distance : distances) {
      for (const char* kernel : kernels) {
        for (const char* neighbor : {"grid", "band"}) {
          const arma::mat m = conley::spatial_meat(
              dup_lat, dup_lon, dup_time, dup_s, 0.0, kernel, distance, false,
              4, neighbor, "double");
          require_close(m(0, 0), 9.0, 0.0,
                        std::string("duplicate cutoff zero ") + distance + "/" +
                            kernel + "/" + neighbor);
        }
      }
    }
    std::puts("duplicates cutoff 0: full outer product OK");

    // Exact antipodes are accepted when the cutoff exceeds the maximum
    // distance. Bartlett retains the analytically expected positive weight.
    arma::vec anti_lat(2), anti_lon(2), anti_time(2, arma::fill::ones);
    arma::mat anti_s(2, 1);
    anti_lat.fill(0.0);
    anti_lon[0] = 0.0;
    anti_lon[1] = 180.0;
    anti_s(0, 0) = 1.0;
    anti_s(1, 0) = 2.0;
    const double oversized = 25000.0;
    for (const char* distance : distances) {
      const double max_distance = std::strcmp(distance, "chord") == 0
          ? 2.0 * conley::AVG_ERAD
          : conley::PI * conley::AVG_ERAD;
      for (const char* kernel : kernels) {
        const arma::mat m = conley::spatial_meat(
            anti_lat, anti_lon, anti_time, anti_s, oversized, kernel, distance,
            false, 1, "grid", "double");
        const double expected = std::strcmp(kernel, "uniform") == 0
            ? 9.0
            : 5.0 + 4.0 * (1.0 - max_distance / oversized);
        require_close(m(0, 0), expected, 2e-12,
                      std::string("antipode ") + distance + "/" + kernel);
      }
    }
    std::puts("antipodes oversized cutoff: accepted OK");

    // One occupied lattice cell always contributes its score outer product,
    // even when the accept threshold rounds to the endpoint one.
    std::vector<int> one_ring(1, 0), one_col(1, 0);
    arma::vec one_time(1, arma::fill::ones);
    arma::mat one_s(1, 1);
    one_s(0, 0) = 2.0;
    for (double cutoff : {0.0, 1e-9}) {
      for (const char* distance : distances) {
        for (const char* kernel : kernels) {
          const arma::mat m = conley::grid_meat(
              one_ring.data(), one_col.data(), one_time, one_s, 35.0, 1.0,
              1.0, 1, 1, 0, cutoff, distance, kernel, 4);
          require_close(m(0, 0), 4.0, 0.0,
                        std::string("one-cell grid ") + distance + "/" + kernel);
        }
      }
    }
    std::puts("one-cell grid cutoff 0 / 1e-9: diagonal OK");

    const double qnan = std::numeric_limits<double>::quiet_NaN();
    expect_conley_error("grid dlon = 0", [&]() {
      conley::grid_meat(one_ring.data(), one_col.data(), one_time, one_s,
                        35.0, 1.0, 0.0, 1, 1, 0, 1.0,
                        "haversine", "uniform", 1);
    });
    expect_conley_error("grid dlat = 0", [&]() {
      conley::grid_meat(one_ring.data(), one_col.data(), one_time, one_s,
                        35.0, 0.0, 1.0, 1, 1, 0, 1.0,
                        "haversine", "uniform", 1);
    });
    expect_conley_error("grid dlat = NaN", [&]() {
      conley::grid_meat(one_ring.data(), one_col.data(), one_time, one_s,
                        35.0, qnan, 1.0, 1, 1, 0, 1.0,
                        "haversine", "uniform", 1);
    });
    expect_conley_error("grid dlon = NaN", [&]() {
      conley::grid_meat(one_ring.data(), one_col.data(), one_time, one_s,
                        35.0, 1.0, qnan, 1, 1, 0, 1.0,
                        "haversine", "uniform", 1);
    });

    arma::vec bad_time = dup_time;
    bad_time[0] = qnan;
    expect_conley_error("spatial NaN time", [&]() {
      conley::spatial_meat(dup_lat, dup_lon, bad_time, dup_s, 1.0,
                           "uniform", "haversine", false, 1, "grid", "double");
    });
    expect_conley_error("serial NaN time", [&]() {
      conley::serial_hac_meat(dup_time, bad_time, 1.0, dup_s, 1);
    });
    arma::mat bad_s = dup_s;
    bad_s(0, 0) = qnan;
    expect_conley_error("spatial NaN score", [&]() {
      conley::spatial_meat(dup_lat, dup_lon, dup_time, bad_s, 1.0,
                           "uniform", "haversine", false, 1, "grid", "double");
    });
    arma::vec bad_unit = dup_time;
    bad_unit[0] = qnan;
    expect_conley_error("serial NaN unit", [&]() {
      conley::serial_hac_meat(bad_unit, dup_time, 1.0, dup_s, 1);
    });

    // Capped cell sizing keeps a near-zero-cutoff 100k-row call output
    // sensitive instead of collapsing every point into one bucket.
    const std::size_t big_n = 100000;
    arma::vec big_lat(big_n), big_lon(big_n), big_time(big_n, arma::fill::ones);
    arma::mat big_s(big_n, 1, arma::fill::ones);
    std::mt19937_64 rng(20260904);
    std::uniform_real_distribution<double> lat_dist(25.0, 49.0);
    std::uniform_real_distribution<double> lon_dist(-124.0, -67.0);
    for (std::size_t i = 0; i < big_n; ++i) {
      big_lat[i] = lat_dist(rng);
      big_lon[i] = lon_dist(rng);
    }
    const auto tiny_start = std::chrono::steady_clock::now();
    const arma::mat tiny = conley::spatial_meat(
        big_lat, big_lon, big_time, big_s, 1e-9, "uniform", "haversine",
        false, 4, "grid", "double");
    const double tiny_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - tiny_start).count();
    require_close(tiny(0, 0), static_cast<double>(big_n), 0.0,
                  "100k tiny-cutoff self meat");
    std::printf("100k cutoff 1e-9: %.3f seconds\n", tiny_seconds);

    // Determinism across thread counts for both pairwise neighbor engines.
    const std::size_t n = 3000;
    arma::vec lat(n), lon(n), time(n, arma::fill::ones);
    arma::mat s(n, 2);
    for (std::size_t i = 0; i < n; ++i) {
      lat[i] = lat_dist(rng);
      lon[i] = lon_dist(rng);
      s(i, 0) = std::sin(static_cast<double>(i));
      s(i, 1) = std::cos(static_cast<double>(i) / 7.0);
    }
    for (const char* neighbor : {"grid", "band"}) {
      const arma::mat m1 = conley::spatial_meat(
          lat, lon, time, s, 300.0, "bartlett", "spherical", false, 1,
          neighbor, "double");
      const arma::mat m4 = conley::spatial_meat(
          lat, lon, time, s, 300.0, "bartlett", "spherical", false, 4,
          neighbor, "double");
      require(bit_identical(m1, m4),
              std::string("spatial ncores identity failed: ") + neighbor);
    }

    const std::size_t units = 200, periods = 5, ns = units * periods;
    arma::vec unit(ns), stime(ns);
    arma::mat ss(ns, 2);
    for (std::size_t u = 0; u < units; ++u) {
      for (std::size_t t = 0; t < periods; ++t) {
        const std::size_t i = u * periods + t;
        unit[i] = static_cast<double>(u);
        stime[i] = static_cast<double>(t);
        ss(i, 0) = std::sin(static_cast<double>(i) / 11.0);
        ss(i, 1) = std::cos(static_cast<double>(i) / 13.0);
      }
    }
    require(bit_identical(conley::serial_hac_meat(unit, stime, 2.0, ss, 1),
                          conley::serial_hac_meat(unit, stime, 2.0, ss, 4)),
            "serial ncores identity failed");

    const int nr = 24, nc = 32;
    const std::size_t ng = static_cast<std::size_t>(nr * nc);
    std::vector<int> gr(ng), gc(ng);
    arma::vec gt(ng, arma::fill::ones);
    arma::mat gs(ng, 2);
    for (int r = 0; r < nr; ++r) {
      for (int c = 0; c < nc; ++c) {
        const std::size_t i = static_cast<std::size_t>(r * nc + c);
        gr[i] = r;
        gc[i] = c;
        gs(i, 0) = std::sin(static_cast<double>(i) / 5.0);
        gs(i, 1) = std::cos(static_cast<double>(i) / 9.0);
      }
    }
    for (const char* kernel : kernels) {
      const arma::mat g1 = conley::grid_meat(
          gr.data(), gc.data(), gt, gs, 20.0, 0.5, 0.5, nr, nc, 0,
          250.0, "spherical", kernel, 1);
      const arma::mat g4 = conley::grid_meat(
          gr.data(), gc.data(), gt, gs, 20.0, 0.5, 0.5, nr, nc, 0,
          250.0, "spherical", kernel, 4);
      require(bit_identical(g1, g4),
              std::string("grid-native ncores identity failed: ") + kernel);
    }
    std::puts("ncores 1 vs 4: spatial grid/band, serial, grid-native OK");

    conley::set_interrupt_hook(&interrupt_now);
    bool interrupted = false;
    try {
      conley::parallel_blocks(10000, 4, 1,
                              [](std::size_t, std::size_t) {});
    } catch (const conley::Error& e) {
      interrupted = std::strcmp(e.what(), "interrupted") == 0;
    }
    conley::set_interrupt_hook(nullptr);
    require(interrupted && interrupt_polls > 0,
            "cooperative interrupt hook did not stop parallel_blocks");
    std::printf("interrupt hook: interrupted after %d calling-thread poll(s)\n",
                interrupt_polls);

    std::puts("core_edge_probes: all checks passed");
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "core_edge_probes FAILED: %s\n", e.what());
    return 1;
  }
}
