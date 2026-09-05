// Focused standalone probes for conley_core.h boundary behavior.
// Compile from the package root (create .scratch first):
// ARMA_INC=$(Rscript -e 'cat(system.file("include", package="RcppArmadillo"))')
// g++ -std=c++14 -O2 -pthread -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS \
//   -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU -I"$ARMA_INC" -Isrc \
//   tests/manual/core_edge_probes.cpp -o .scratch/core_edge_probes
// .scratch/core_edge_probes
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

// Check the production cursor ranges against the previous binary searches.
// Synthetic cell centres let us cover arbitrary occupied cell sets, including
// cube faces and IDs near 2^63 - 1, independently of spherical sampling.
void check_cell_ranges(std::vector<std::uint64_t> ids, std::int64_t G,
                       std::mt19937_64& rng) {
  if (!ids.empty()) ids.push_back(ids.front());  // duplicate rows in a cell
  std::shuffle(ids.begin(), ids.end(), rng);
  const std::uint64_t uG = static_cast<std::uint64_t>(G);
  const std::size_t n = ids.size();
  const conley::GridSpec spec{2.0 / static_cast<double>(G), G};
  conley::CoordCache coords;
  coords.x3.resize(2 * n);
  coords.y3.resize(2 * n);
  coords.z3.resize(2 * n);
  for (std::size_t i = 0; i < 2 * n; ++i) {
    const std::uint64_t id = ids[i % n];
    coords.x3[i] = (static_cast<double>(id / (uG * uG)) + 0.5) * spec.cell - 1.0;
    coords.y3[i] = (static_cast<double>((id / uG) % uG) + 0.5) * spec.cell - 1.0;
    coords.z3[i] = (static_cast<double>(id % uG) + 0.5) * spec.cell - 1.0;
  }
  conley::CellGrid grid;
  grid.row_cell.resize(2 * n);
  std::vector<std::size_t> sorted;
  std::vector<std::uint64_t> ordered = ids;
  std::sort(ordered.begin(), ordered.end());
  std::vector<std::uint64_t> ucid;
  std::vector<std::size_t> ustart;
  for (std::size_t i = 0; i < n; ++i) {
    if (i == 0 || ordered[i] != ordered[i - 1]) {
      ucid.push_back(ordered[i]);
      ustart.push_back(i);
    }
  }
  ustart.push_back(n);
  for (std::size_t block = 0; block < 2; ++block) {
    const std::size_t row0 = block * n, cell0 = grid.cell_start.size();
    conley::append_block_grid(coords, row0, row0 + n, spec, sorted, grid,
                             block == 0 ? 1 : 4);
    require(grid.cell_start.size() == cell0 + ucid.size(), "cell count");
    for (std::size_t i = 0; i < n; ++i) {
      require(ids[sorted[row0 + i] - row0] == ordered[i], "cell ordering");
    }
    for (std::size_t cdx = 0; cdx < ucid.size(); ++cdx) {
      require(grid.cell_start[cell0 + cdx] == row0 + ustart[cdx], "cell start");
      const std::uint64_t cu = ucid[cdx];
      const std::int64_t cx = cu / (uG * uG), cy = (cu / uG) % uG, cz = cu % uG;
      const std::int64_t intervals[5][4] = {
        {cx, cy, cz + 1, cz + 1}, {cx, cy + 1, cz - 1, cz + 1},
        {cx + 1, cy - 1, cz - 1, cz + 1}, {cx + 1, cy, cz - 1, cz + 1},
        {cx + 1, cy + 1, cz - 1, cz + 1}
      };
      for (std::size_t run = 0; run < 5; ++run) {
        const auto& in = intervals[run];
        const std::int64_t zlo = std::max<std::int64_t>(0, in[2]);
        const std::int64_t zhi = std::min(G - 1, in[3]);
        std::size_t lo = row0, hi = row0;
        if (in[0] >= 0 && in[0] < G && in[1] >= 0 && in[1] < G && zlo <= zhi) {
          const std::uint64_t base = (static_cast<std::uint64_t>(in[0]) * uG + in[1]) * uG;
          const auto clo = std::lower_bound(ucid.begin(), ucid.end(), base + zlo);
          const auto chi = std::upper_bound(clo, ucid.end(), base + zhi);
          lo += ustart[clo - ucid.begin()];
          hi += ustart[chi - ucid.begin()];
        }
        const std::size_t pos = 10 * (cell0 + cdx) + 2 * run;
        require(grid.nbr[pos] == lo && grid.nbr[pos + 1] == hi,
                "cursor/binary range mismatch G=" + std::to_string(G));
      }
    }
  }
}

void probe_cell_ranges() {
  std::mt19937_64 rng(20260905);
  std::size_t cases = 0;
  const auto capped = conley::make_grid_spec(
      conley::make_screen_params(1e-9, conley::DIST_SPHERICAL));
  require(capped.G == 2097152, "tiny cutoff must exercise capped grid");
  for (std::int64_t G : {1LL, 2LL, 3LL, 4LL, 17LL, 127LL, 2097152LL}) {
    const std::uint64_t uG = static_cast<std::uint64_t>(G);
    std::vector<std::uint64_t> faces;
    const std::vector<std::int64_t> edge = {0, std::min<std::int64_t>(1, G - 1),
                                           G / 2, std::max<std::int64_t>(0, G - 2), G - 1};
    for (auto x : edge) for (auto y : edge) for (auto z : edge) {
      faces.push_back((static_cast<std::uint64_t>(x) * uG + y) * uG + z);
    }
    check_cell_ranges({}, G, rng);
    check_cell_ranges({uG * uG * uG - 1}, G, rng);
    check_cell_ranges(faces, G, rng);
    cases += 3;
    for (int trial = 0; trial < 12; ++trial) {
      std::vector<std::uint64_t> ids = faces;
      for (std::size_t i = 0; i < 10000; ++i) {
        ids.push_back(((rng() % uG) * uG + rng() % uG) * uG + rng() % uG);
      }
      check_cell_ranges(ids, G, rng);
      ++cases;
    }
  }
  // A parallel chunk starts inside cy == 0, where the cy - 1 family is
  // invalid, and reaches cy == 1 later. Its first valid interval must seed
  // the cursor then; the final x face also leaves whole families invalid.
  const std::uint64_t G = static_cast<std::uint64_t>(capped.G);
  std::vector<std::uint64_t> delayed;
  for (std::uint64_t z = 0; z < 10000; ++z) {
    delayed.push_back(z);
    delayed.push_back(G + z);
    delayed.push_back((G - 1) * G * G + z);
  }
  check_cell_ranges(delayed, capped.G, rng);
  ++cases;
  std::printf("neighbour cursor vs binary searches: %zu cell sets, two blocks each, all ranges identical\n", cases);
}

// Compare every weight spectrum on the small benchmark lattice, retaining
// the imaginary component too. Exercise full batches and short tails.
void probe_weight_fft() {
  const std::size_t R = 180, C = 180;
  conley::GridGeom gg;
  gg.n_ring = R; gg.n_col = C;
  gg.dlat_rad = gg.dlam_rad = 0.05 * conley::DE2RA;
  gg.sphi.resize(R); gg.cphi.resize(R);
  gg.cos_dl.resize(C); gg.s2h_dl.resize(C);
  for (std::size_t r = 0; r < R; ++r) {
    const double phi = (35.0 + r * 0.05) * conley::DE2RA;
    gg.sphi[r] = std::sin(phi); gg.cphi[r] = std::cos(phi);
  }
  for (std::size_t d = 0; d < C; ++d) {
    gg.cos_dl[d] = std::cos(d * gg.dlam_rad);
    gg.s2h_dl[d] = conley::sq(std::sin(d * gg.dlam_rad / 2.0));
  }
  std::size_t checked = 0;
  for (int distance : {conley::DIST_HAVERSINE, conley::DIST_SPHERICAL, conley::DIST_CHORD}) {
    const auto screen = conley::make_screen_params(250.0, distance);
    long ds_max = 0;
    for (std::size_t r1 = 0; r1 < R; ++r1) for (std::size_t r2 = r1; r2 < R; ++r2) {
      ds_max = std::max(ds_max, conley::grid_ring_window(
          gg, r1, r2, screen.sin2_half_angular_cutoff, C - 1).ds);
    }
    const std::size_t npad = conley::next_pow2(C + ds_max);
    for (std::size_t batch : {8, 16}) {
      arma::mat weights(npad, batch);
      std::vector<double> w(C);
      for (std::size_t r1 = 0; r1 < R; ++r1) {
        std::size_t nw = 0;
        const auto check = [&]() {
          if (nw == 0) return;
          const arma::mat input(weights.memptr(), npad, nw, false, true);
          const arma::cx_mat got = arma::fft(input);
          for (std::size_t j = 0; j < nw; ++j) {
            const arma::vec old_input = input.col(j);
            const arma::cx_vec want = arma::fft(old_input);
            require(std::memcmp(want.memptr(), got.colptr(j),
                                npad * sizeof(std::complex<double>)) == 0,
                    "batched weight FFT spectrum changed");
            ++checked;
          }
        };
        for (std::size_t r2 = r1; r2 < R; ++r2) {
          const auto win = conley::grid_ring_window(
              gg, r1, r2, screen.sin2_half_angular_cutoff, C - 1);
          if (win.ds < 0) continue;
          conley::grid_bartlett_weights(gg, r1, r2, distance, 250.0, win.ds, w.data());
          if (r2 == r1) for (long d = 0; d <= win.ds; ++d) w[d] *= 0.5;
          weights.col(nw).zeros();
          weights(0, nw) = w[0];
          for (long d = 1; d <= win.ds; ++d) {
            weights(d, nw) = weights(npad - d, nw) = w[d];
          }
          if (++nw == batch) { check(); nw = 0; }
        }
        check();
      }
    }
  }
  std::printf("180x180 Bartlett weight FFTs: %zu spectra bit-identical (batches 8/16, all distances)\n", checked);
}

}  // namespace

int main() {
  try {
    probe_cell_ranges();
    probe_weight_fft();
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
