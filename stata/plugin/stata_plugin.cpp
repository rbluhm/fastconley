// stata_plugin.cpp -- Stata plugin front-end for the fastconley engine
// (src/conley_core.h, shared with the R package). Built by the Makefile in
// this directory into stata/src/fastconley_<platform>.plugin.
//
// Protocol (driven by fastconley.ado; every variable is passed in the order
// listed, restricted by `if`/`in` to the prepared rows):
//   plugin call fastconley_plugin, check
//       -> globals FASTCONLEY_ENGINE_VERSION, FASTCONLEY_ENGINE_BUILD
//   plugin call fastconley_plugin lat lon time s1 ... sk in 1/n, ///
//       spatial <cutoff scalar> <kernel> <dist> <balanced 0|1> <threads> <neighbor> <csr_weight> <matname>
//       -> k x k meat stored in Stata matrix <matname> (must exist),
//          local fc_unbalanced_fallback
//   plugin call fastconley_plugin unit time s1 ... sk in 1/n, serial <lag scalar> <threads> <matname>
//   plugin call fastconley_plugin ring col time s1 ... sk in 1/n, ///
//       grid <lat0 scalar> <dlat scalar> <dlon scalar> <n_ring> <n_col> <n_col_full> <cutoff scalar> <dist> <kernel> <threads> <matname>
// Real-valued arguments are passed as names of Stata scalars (negative
// literals do not survive plugin call's argument parsing).
// Errors: conley::Error / std::exception -> SF_error + local fc_plugin_error,
// return code 198. Rows must already be sorted (time blocks contiguous;
// (unit, time) for serial) and free of missing values: the ado prepares
// them in Mata exactly as the R package does.
// The engine header comes first: the SPI header defines SYSTEM as a bare
// number, which clashes with identifiers in some Armadillo dependencies.
#include "conley_core.h"
#include "stplugin.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <vector>

#ifndef FASTCONLEY_BUILD_ID
#define FASTCONLEY_BUILD_ID "local"
#endif

namespace {

const char* const SAMPLE_TOO_LARGE =
    "sample too large for the compiled engine; use engine(mata)";

bool stata_interrupt_requested() {
  SF_poll();
  return SW_stopflag != 0;
}

std::vector<ST_int> sample_rows() {
  std::vector<ST_int> rows;
  const ST_int a = SF_in1();
  const ST_int b = SF_in2();
  const ST_int nobs = SF_nobs();
  if (a < 0 || b < 0 || nobs < 0) throw conley::Error(SAMPLE_TOO_LARGE);
  if (b < a) return rows;
  const std::uint64_t requested =
      static_cast<std::uint64_t>(static_cast<unsigned int>(b)) -
      static_cast<std::uint64_t>(static_cast<unsigned int>(a)) + 1U;
  if (requested > static_cast<std::uint64_t>(std::numeric_limits<ST_int>::max())) {
    throw conley::Error(SAMPLE_TOO_LARGE);
  }
  rows.reserve(static_cast<std::size_t>(requested));
  for (ST_int i = a;; ++i) {
    if (SF_ifobs(i)) rows.push_back(i);
    if (i == b) break;
  }
  return rows;
}

void read_column(ST_int var, const std::vector<ST_int>& rows, double* out) {
  for (std::size_t r = 0; r < rows.size(); ++r) {
    double v = 0.0;
    if (SF_vdata(var, rows[r], &v)) {
      throw conley::Error("could not read plugin input variable " + std::to_string(var));
    }
    if (SF_is_missing(v)) {
      throw conley::Error("missing value in plugin input variable " + std::to_string(var));
    }
    out[r] = v;
  }
}

void store_matrix(const std::string& name, const arma::mat& M) {
  for (arma::uword i = 0; i < M.n_rows; ++i) {
    for (arma::uword j = 0; j < M.n_cols; ++j) {
      if (SF_mat_store(const_cast<char*>(name.c_str()),
                       static_cast<ST_int>(i + 1), static_cast<ST_int>(j + 1), M(i, j))) {
        throw conley::Error("could not store result matrix " + name +
                            " (the ado must create it with the right dimensions first)");
      }
    }
  }
}

void save_local_checked(const char* name, const char* value) {
  if (SF_macro_save(const_cast<char*>(name), const_cast<char*>(value))) {
    throw conley::Error("could not save plugin result local");
  }
}

void report_error_noexcept(const char* what) noexcept {
  if (!what) what = "unknown error";
  SF_macro_save(const_cast<char*>("_fc_plugin_error"), const_cast<char*>(what));
  SF_error(const_cast<char*>("fastconley plugin: "));
  SF_error(const_cast<char*>(what));
  SF_error(const_cast<char*>("\n"));
}

double to_double(const char* s, const char* what) {
  errno = 0;
  char* end = 0;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0' || errno == ERANGE || !std::isfinite(v)) {
    throw conley::Error(std::string("bad numeric argument for ") + what);
  }
  return v;
}

// Real-valued arguments arrive as the names of Stata scalars: `plugin call`
// does not pass negative numbers through argv intact.
double read_scalar(const char* name, const char* what) {
  double v = 0.0;
  if (SF_scal_use(const_cast<char*>(name), &v)) {
    throw conley::Error(std::string("could not read scalar ") + name + " for " + what);
  }
  if (SF_is_missing(v) || !std::isfinite(v)) {
    throw conley::Error(std::string("scalar ") + name + " for " + what +
                        " must be finite and nonmissing");
  }
  return v;
}

int to_int(const char* s, const char* what) {
  const double v = to_double(s, what);
  if (std::trunc(v) != v ||
      v < static_cast<double>(std::numeric_limits<int>::min()) ||
      v > static_cast<double>(std::numeric_limits<int>::max())) {
    throw conley::Error(std::string(what) + " must be an integer in range");
  }
  return static_cast<int>(v);
}

int value_to_int(double v, const char* what) {
  if (!std::isfinite(v) || std::trunc(v) != v ||
      v < static_cast<double>(std::numeric_limits<int>::min()) ||
      v > static_cast<double>(std::numeric_limits<int>::max())) {
    throw conley::Error(std::string(what) + " values must be finite integers in range");
  }
  return static_cast<int>(v);
}

}  // namespace

STDLL stata_call(int argc, char* argv[]) {
  try {
    conley::set_interrupt_hook(&stata_interrupt_requested);
    if (argc < 1) throw conley::Error("missing subcommand");
    const std::string cmd = argv[0];

    if (cmd == "check") {
      // Clear both handshake globals before publishing either value, so a
      // partial failure cannot leave a stale compatible-looking handshake.
      const int clear_version = SF_macro_save(
          const_cast<char*>("FASTCONLEY_ENGINE_VERSION"), const_cast<char*>(""));
      const int clear_build = SF_macro_save(
          const_cast<char*>("FASTCONLEY_ENGINE_BUILD"), const_cast<char*>(""));
      if (clear_version || clear_build) {
        throw conley::Error("could not clear plugin version handshake globals");
      }
      const int save_version = SF_macro_save(
          const_cast<char*>("FASTCONLEY_ENGINE_VERSION"),
          const_cast<char*>(CONLEY_CORE_VERSION));
      const int save_build = SF_macro_save(
          const_cast<char*>("FASTCONLEY_ENGINE_BUILD"),
          const_cast<char*>(FASTCONLEY_BUILD_ID));
      if (save_version || save_build) {
        throw conley::Error("could not save plugin version handshake globals");
      }
      return 0;
    }

    const std::vector<ST_int> rows = sample_rows();
    const std::size_t n = rows.size();
    const int nvars = SF_nvars();

    if (cmd == "spatial") {
      if (argc != 9) throw conley::Error("spatial expects 8 arguments");
      if (nvars < 4) throw conley::Error("spatial expects lat lon time and at least one score column");
      const std::size_t k = static_cast<std::size_t>(nvars - 3);
      arma::vec lat(n), lon(n), time(n);
      arma::mat S(n, k);
      read_column(1, rows, lat.memptr());
      read_column(2, rows, lon.memptr());
      read_column(3, rows, time.memptr());
      for (std::size_t kk = 0; kk < k; ++kk) read_column(static_cast<ST_int>(4 + kk), rows, S.colptr(kk));
      bool fallback = false;
      const arma::mat M = conley::spatial_meat(
          lat, lon, time, S, read_scalar(argv[1], "cutoff"), argv[2], argv[3],
          to_int(argv[4], "balanced") != 0, to_int(argv[5], "threads"), argv[6], argv[7],
          &fallback);
      save_local_checked("_fc_unbalanced_fallback", fallback ? "1" : "0");
      store_matrix(argv[8], M);
      return 0;
    }

    if (cmd == "serial") {
      if (argc != 4) throw conley::Error("serial expects 3 arguments");
      if (nvars < 3) throw conley::Error("serial expects unit time and at least one score column");
      const std::size_t k = static_cast<std::size_t>(nvars - 2);
      arma::vec unit(n), time(n);
      arma::mat S(n, k);
      read_column(1, rows, unit.memptr());
      read_column(2, rows, time.memptr());
      for (std::size_t kk = 0; kk < k; ++kk) read_column(static_cast<ST_int>(3 + kk), rows, S.colptr(kk));
      const arma::mat M = conley::serial_hac_meat(unit, time, read_scalar(argv[1], "lag"), S,
                                                  to_int(argv[2], "threads"));
      store_matrix(argv[3], M);
      return 0;
    }

    if (cmd == "grid") {
      if (argc != 12) throw conley::Error("grid expects 11 arguments");
      if (nvars < 4) throw conley::Error("grid expects ring col time and at least one score column");
      const std::size_t k = static_cast<std::size_t>(nvars - 3);
      arma::vec ringd(n), cold(n), time(n);
      arma::mat S(n, k);
      read_column(1, rows, ringd.memptr());
      read_column(2, rows, cold.memptr());
      read_column(3, rows, time.memptr());
      for (std::size_t kk = 0; kk < k; ++kk) read_column(static_cast<ST_int>(4 + kk), rows, S.colptr(kk));
      std::vector<int> ring(n), col(n);
      for (std::size_t i = 0; i < n; ++i) {
        ring[i] = value_to_int(ringd[i], "ring");
        col[i] = value_to_int(cold[i], "col");
      }
      const arma::mat M = conley::grid_meat(
          ring.data(), col.data(), time, S,
          read_scalar(argv[1], "lat0"), read_scalar(argv[2], "dlat"), read_scalar(argv[3], "dlon"),
          to_int(argv[4], "n_ring"), to_int(argv[5], "n_col"), to_int(argv[6], "n_col_full"),
          read_scalar(argv[7], "cutoff"), argv[8], argv[9], to_int(argv[10], "threads"));
      store_matrix(argv[11], M);
      return 0;
    }

    throw conley::Error("unknown subcommand '" + cmd + "'");
  } catch (const std::bad_alloc&) {
    report_error_noexcept("out of memory");
    return 198;
  } catch (const std::exception& e) {
    report_error_noexcept(e.what());
    return 198;
  } catch (...) {
    report_error_noexcept("unknown error");
    return 198;
  }
}
