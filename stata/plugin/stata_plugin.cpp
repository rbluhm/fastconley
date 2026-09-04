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

#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#ifndef FASTCONLEY_BUILD_ID
#define FASTCONLEY_BUILD_ID "local"
#endif

namespace {

std::vector<ST_int> sample_rows() {
  std::vector<ST_int> rows;
  const ST_int a = SF_in1();
  const ST_int b = SF_in2();
  if (b >= a) rows.reserve(static_cast<std::size_t>(b - a + 1));
  for (ST_int i = a; i <= b; ++i) {
    if (SF_ifobs(i)) rows.push_back(i);
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

void save_local(const char* name, const std::string& value) {
  std::string key = std::string("_") + name;   // leading underscore = local macro
  SF_macro_save(const_cast<char*>(key.c_str()), const_cast<char*>(value.c_str()));
}

void save_global(const char* name, const std::string& value) {
  SF_macro_save(const_cast<char*>(name), const_cast<char*>(value.c_str()));
}

double to_double(const char* s, const char* what) {
  char* end = 0;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0') throw conley::Error(std::string("bad numeric argument for ") + what);
  return v;
}

// Real-valued arguments arrive as the names of Stata scalars: `plugin call`
// does not pass negative numbers through argv intact.
double read_scalar(const char* name, const char* what) {
  double v = 0.0;
  if (SF_scal_use(const_cast<char*>(name), &v)) {
    throw conley::Error(std::string("could not read scalar ") + name + " for " + what);
  }
  return v;
}

int to_int(const char* s, const char* what) {
  return static_cast<int>(to_double(s, what));
}

}  // namespace

STDLL stata_call(int argc, char* argv[]) {
  try {
    if (argc < 1) throw conley::Error("missing subcommand");
    const std::string cmd = argv[0];

    if (cmd == "check") {
      // Globals: the ado runs this from a program and reads them afterwards.
      save_global("FASTCONLEY_ENGINE_VERSION", CONLEY_CORE_VERSION);
      save_global("FASTCONLEY_ENGINE_BUILD", FASTCONLEY_BUILD_ID);
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
      save_local("fc_unbalanced_fallback", fallback ? "1" : "0");
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
        ring[i] = static_cast<int>(ringd[i]);
        col[i] = static_cast<int>(cold[i]);
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
  } catch (const std::exception& e) {
    const std::string what = e.what();
    save_local("fc_plugin_error", what);
    std::string msg = "fastconley plugin: " + what + "\n";
    SF_error(const_cast<char*>(msg.c_str()));
    return 198;
  } catch (...) {
    save_local("fc_plugin_error", "unknown error");
    SF_error(const_cast<char*>("fastconley plugin: unknown error\n"));
    return 198;
  }
}
