// Deterministic, header-only numerical check for the shared Conley engine.
// Inputs use only raw mt19937_64 output, so no implementation-defined random
// distributions enter the cross-platform golden comparison.
#ifndef CORE_HEADER
#define CORE_HEADER "conley_core.h"
#endif
#include CORE_HEADER

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct DetRng {
  std::mt19937_64 e;
  explicit DetRng(std::uint64_t seed) : e(seed) {}
  double u() {
    return static_cast<double>(e() >> 11) * (1.0 / 9007199254740992.0);
  }
  double u(double lo, double hi) { return lo + (hi - lo) * u(); }
};

struct CheckCase {
  std::string name;
  arma::mat value;
};

arma::mat random_scores(DetRng& rng, std::size_t n, std::size_t k) {
  arma::mat scores(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < k; ++j) {
      scores(i, j) = rng.u(-1.0, 1.0) * 1.7320508075688772;
    }
  }
  return scores;
}

std::vector<CheckCase> run_cases() {
  const std::size_t n = 200;
  const std::size_t k = 3;
  std::vector<CheckCase> out;

  // Cross-section. Keep this generation order synchronized with load_check.c.
  DetRng cross_rng(101);
  arma::vec lat(n), lon(n), time(n, arma::fill::ones);
  arma::mat scores(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    lat[i] = cross_rng.u(25.0, 49.0);
    lon[i] = cross_rng.u(-124.0, -67.0);
    for (std::size_t j = 0; j < k; ++j) {
      scores(i, j) = cross_rng.u(-1.0, 1.0) * 1.7320508075688772;
    }
  }
  out.push_back({"spatial_cross_bartlett_haversine",
                 conley::spatial_meat(lat, lon, time, scores, 500.0,
                                      "bartlett", "haversine", false, 2,
                                      "grid", "double")});
  out.push_back({"spatial_cross_uniform_spherical",
                 conley::spatial_meat(lat, lon, time, scores, 500.0,
                                      "uniform", "spherical", false, 2,
                                      "grid", "double")});

  // Two equal time blocks with stable unit coordinates: the balanced path.
  DetRng panel_rng(202);
  arma::vec unit_lat(n / 2), unit_lon(n / 2);
  for (std::size_t i = 0; i < n / 2; ++i) {
    unit_lat[i] = panel_rng.u(-60.0, 60.0);
    unit_lon[i] = panel_rng.u(-175.0, 175.0);
  }
  arma::vec panel_lat(n), panel_lon(n), panel_time(n);
  for (std::size_t t = 0; t < 2; ++t) {
    for (std::size_t i = 0; i < n / 2; ++i) {
      const std::size_t row = t * (n / 2) + i;
      panel_lat[row] = unit_lat[i];
      panel_lon[row] = unit_lon[i];
      panel_time[row] = static_cast<double>(t + 1);
    }
  }
  const arma::mat panel_scores = random_scores(panel_rng, n, k);
  out.push_back({"spatial_panel_bartlett_haversine",
                 conley::spatial_meat(panel_lat, panel_lon, panel_time,
                                      panel_scores, 700.0, "bartlett",
                                      "haversine", true, 2, "grid", "double")});
  out.push_back({"spatial_panel_uniform_spherical",
                 conley::spatial_meat(panel_lat, panel_lon, panel_time,
                                      panel_scores, 700.0, "uniform",
                                      "spherical", true, 2, "grid", "double")});

  // Forty units observed for five consecutive periods, sorted by unit/time.
  DetRng serial_rng(303);
  arma::vec unit(n), serial_time(n);
  arma::mat serial_scores(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    unit[i] = static_cast<double>(i / 5 + 1);
    serial_time[i] = static_cast<double>(i % 5 + 1);
    for (std::size_t j = 0; j < k; ++j) {
      serial_scores(i, j) = serial_rng.u(-1.0, 1.0) * 1.7320508075688772;
    }
  }
  out.push_back({"serial_lag2",
                 conley::serial_hac_meat(unit, serial_time, 2.0,
                                         serial_scores, 2)});

  // A complete 10 x 20 globe-spanning lattice. n_col_full makes longitude
  // periodic, so pairs around columns 0 and 19 exercise the wrap path.
  DetRng grid_rng(404);
  std::vector<int> ring(n), col(n);
  arma::vec grid_time(n, arma::fill::ones);
  arma::mat grid_scores(n, k);
  for (std::size_t i = 0; i < n; ++i) {
    ring[i] = static_cast<int>(i / 20);
    col[i] = static_cast<int>(i % 20);
    for (std::size_t j = 0; j < k; ++j) {
      grid_scores(i, j) = grid_rng.u(-1.0, 1.0) * 1.7320508075688772;
    }
  }
  out.push_back({"grid_wrap_uniform_spherical",
                 conley::grid_meat(ring.data(), col.data(), grid_time,
                                   grid_scores, -22.5, 5.0, 18.0, 10, 20, 20,
                                   2500.0, "spherical", "uniform", 2)});
  out.push_back({"grid_wrap_bartlett_haversine",
                 conley::grid_meat(ring.data(), col.data(), grid_time,
                                   grid_scores, -22.5, 5.0, 18.0, 10, 20, 20,
                                   2500.0, "haversine", "bartlett", 2)});
  return out;
}

void write_golden(const char* path, const std::vector<CheckCase>& cases) {
  FILE* fp = std::fopen(path, "wb");
  if (!fp) throw std::runtime_error(std::string("cannot write ") + path);
  for (std::size_t c = 0; c < cases.size(); ++c) {
    std::fprintf(fp, "%s", cases[c].name.c_str());
    for (std::size_t i = 0; i < cases[c].value.n_rows; ++i) {
      for (std::size_t j = 0; j < cases[c].value.n_cols; ++j) {
        std::fprintf(fp, " %.17g", cases[c].value(i, j));
      }
    }
    std::fputc('\n', fp);
  }
  if (std::fclose(fp) != 0) {
    throw std::runtime_error(std::string("cannot close ") + path);
  }
}

std::vector<CheckCase> read_golden(const char* path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("cannot read ") + path);
  std::vector<CheckCase> cases;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream in(line);
    CheckCase c;
    c.value.set_size(3, 3);
    if (!(in >> c.name)) throw std::runtime_error("golden line has no name");
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        if (!(in >> c.value(i, j))) {
          throw std::runtime_error("golden line has fewer than nine values: " + c.name);
        }
      }
    }
    std::string extra;
    if (in >> extra) throw std::runtime_error("golden line has extra values: " + c.name);
    cases.push_back(c);
  }
  return cases;
}

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

void compare_cases(const std::vector<CheckCase>& actual,
                   const std::vector<CheckCase>& expected, bool strict) {
  if (actual.size() != expected.size()) {
    throw std::runtime_error("golden case count mismatch");
  }
  for (std::size_t c = 0; c < actual.size(); ++c) {
    if (actual[c].name != expected[c].name) {
      throw std::runtime_error("golden case order/name mismatch: " + actual[c].name);
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        const double got = actual[c].value(i, j);
        const double want = expected[c].value(i, j);
        const double scale = std::max(1.0, std::fabs(want));
        const double rel = std::fabs(got - want) / scale;
        worst = std::max(worst, rel);
        if ((strict && !same_bits(got, want)) || (!strict && rel > 1e-12)) {
          char msg[512];
          std::snprintf(msg, sizeof(msg),
                        "%s[%lu,%lu]: got %.17g, expected %.17g, rel %.17g%s",
                        actual[c].name.c_str(), static_cast<unsigned long>(i),
                        static_cast<unsigned long>(j), got, want, rel,
                        strict ? " (strict bit mismatch)" : "");
          throw std::runtime_error(msg);
        }
      }
    }
    std::printf("%-39s m00=%.17g m01=%.17g m12=%.17g m22=%.17g rel=%.3g\n",
                actual[c].name.c_str(), actual[c].value(0, 0),
                actual[c].value(0, 1), actual[c].value(1, 2),
                actual[c].value(2, 2), worst);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    bool strict = false;
    bool write = false;
    const char* path = "golden.txt";
    for (int i = 1; i < argc; ++i) {
      if (std::strcmp(argv[i], "--strict") == 0) strict = true;
      else if (std::strcmp(argv[i], "--write") == 0) write = true;
      else path = argv[i];
    }
    if (strict && write) throw std::runtime_error("--strict and --write are mutually exclusive");
    const std::vector<CheckCase> actual = run_cases();
    if (write) {
      write_golden(path, actual);
      std::printf("golden_check: wrote %lu cases to %s\n",
                  static_cast<unsigned long>(actual.size()), path);
      return 0;
    }
    compare_cases(actual, read_golden(path), strict);
    std::printf("golden_check: %lu cases passed (%s)\n",
                static_cast<unsigned long>(actual.size()),
                strict ? "strict bit identity" : "relative tolerance 1e-12");
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "golden_check: %s\n", e.what());
    return 1;
  }
}
