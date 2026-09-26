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
//   plugin call fastconley_plugin lat lon [klat klon] time [unit] s1 ... sk in 1/n, ///
//       vce <cutoff scalar> <kernel> <dist> <balanced 0|1> <threads> <neighbor> <csr_weight> ///
//           <method auto|pairwise|grid> <grid_tol scalar> <lag scalar> <haskeys 0|1> <hasunit 0|1> ///
//           <spatial matname> <serial matname>
//       The whole plugin-side VCE step on the raw sample rows (engine 0.11.3+):
//       balanced-panel validation, aggregation of identical (time, klat, klon)
//       keys, lattice detection and the auto rule, the spatial meat (grid or
//       pairwise), and the serial meat after an internal (unit, time) sort.
//       It reproduces fastconley.mata's fastconley_prepare_rows(),
//       fastconley_detect_grid() and fastconley_choose_grid() step for step
//       (see those functions); klat/klon are the pixel-snapped keys computed
//       in Mata (absent when pixel = 0, when the keys are lat/lon). Returns
//       locals fc_n_sp, fc_n_periods, fc_sp_balanced, fc_method_used,
//       fc_agg_fallback, fc_unbalanced_fallback, fc_dateline_fallback,
//       fc_serial_done. User errors (balanced-panel validation, method(grid)
//       without a usable lattice) return 3498 instead of 198.
// Real-valued arguments are passed as names of Stata scalars (negative
// literals do not survive plugin call's argument parsing).
// Errors: conley::Error / std::exception -> SF_error + local fc_plugin_error,
// return code 198. For spatial, serial, and grid the rows must already be
// sorted (time blocks contiguous; (unit, time) for serial); vce takes the
// raw sample rows and sorts them itself. No subcommand accepts missing values.
// The engine header comes first: the SPI header defines SYSTEM as a bare
// number, which clashes with identifiers in some Armadillo dependencies.
#include "conley_core.h"
#include "stplugin.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <new>
#include <numeric>
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

// ---------------------------------------------------------------------------
// vce subcommand helpers
// ---------------------------------------------------------------------------

// Errors in the user's specification (Mata raised these with _error(3498)).
struct UserError : public conley::Error {
  explicit UserError(const std::string& msg) : conley::Error(msg) {}
};

// Order-preserving 64-bit key for a finite double. Both zeros map to the
// key of +0, since Mata's comparisons treat them as equal.
inline std::uint64_t double_sort_key(double v) {
  if (v == 0.0) v = 0.0;
  std::uint64_t b = 0;
  std::memcpy(&b, &v, sizeof b);
  return (b >> 63) ? ~b : (b | (std::uint64_t(1) << 63));
}

// Stable order of rows 0..n-1 by the given keys (most significant first):
// LSD passes of the engine's stable radix sort, so ties keep input order.
// Mata's order() is not stable, so rows with fully tied keys may be visited
// in a different order than in the Mata preparation; with distinct keys the
// order is identical.
std::vector<std::size_t> stable_key_order(const std::vector<const double*>& keys,
                                          std::size_t n, int ncores) {
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), std::size_t(0));
  std::vector<std::uint64_t> kv(n);
  std::vector<std::size_t> q, next(n);
  for (auto it = keys.rbegin(); it != keys.rend(); ++it) {
    const double* x = *it;
    bool constant = true;
    for (std::size_t i = 1; i < n; ++i) {
      if (x[i] != x[0]) {
        constant = false;
        break;
      }
    }
    if (constant) continue;
    for (std::size_t i = 0; i < n; ++i) kv[i] = double_sort_key(x[ord[i]]);
    conley::radix_order(kv, q, ncores);
    for (std::size_t i = 0; i < n; ++i) next[i] = ord[q[i]];
    ord.swap(next);
  }
  return ord;
}

// Mata's round() on the nonnegative arguments used below.
inline double mata_round(double x) { return std::floor(x + 0.5); }

// Sorted distinct values (uniqrows on a column).
std::vector<double> sorted_unique(const arma::vec& x) {
  std::vector<double> u(x.begin(), x.end());
  bool sorted = true;
  for (std::size_t i = 1; i < u.size(); ++i) {
    if (u[i] < u[i - 1]) {
      sorted = false;
      break;
    }
  }
  if (!sorted) std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return u;
}

struct Lattice {
  double lat0 = 0, dlat = 0, lon0 = 0, dlon = 0;
  int n_ring = 0, n_col = 0, n_col_full = 0;
  std::vector<int> ring, col;
};

// Port of fastconley_detect_grid() (fastconley.mata), same arithmetic.
bool detect_lattice(const arma::vec& lat, const arma::vec& lon, double tol, Lattice& g) {
  const std::size_t n = lat.n_elem;
  const std::vector<double> ul = sorted_unique(lat);
  if (ul.size() < 2) return false;
  std::vector<double> dl(ul.size() - 1);
  for (std::size_t i = 0; i + 1 < ul.size(); ++i) dl[i] = ul[i + 1] - ul[i];
  double step_l = *std::min_element(dl.begin(), dl.end());
  if (step_l <= 0) return false;
  double nsteps = mata_round((ul.back() - ul.front()) / step_l);
  if (nsteps >= 1) step_l = (ul.back() - ul.front()) / nsteps;
  for (double d : dl) {
    if (std::abs(d / step_l - mata_round(d / step_l)) > tol) return false;
  }
  const std::vector<double> uo = sorted_unique(lon);
  if (uo.size() < 2) return false;
  std::vector<double> dol(uo.size() - 1);
  for (std::size_t i = 0; i + 1 < uo.size(); ++i) dol[i] = uo[i + 1] - uo[i];
  double step_o = *std::min_element(dol.begin(), dol.end());
  if (step_o <= 0) return false;
  nsteps = mata_round((uo.back() - uo.front()) / step_o);
  if (nsteps >= 1) step_o = (uo.back() - uo.front()) / nsteps;
  for (double d : dol) {
    if (std::abs(d / step_o - mata_round(d / step_o)) > tol) return false;
  }
  std::vector<double> ring_d(n), col_d(n);
  double max_ring = 0, max_col = 0, dev_lat = 0, dev_lon = 0;
  for (std::size_t i = 0; i < n; ++i) {
    ring_d[i] = mata_round((lat[i] - ul.front()) / step_l);
    col_d[i] = mata_round((lon[i] - uo.front()) / step_o);
    max_ring = std::max(max_ring, ring_d[i]);
    max_col = std::max(max_col, col_d[i]);
  }
  const double lim = 2147483648.0 - 2.0;
  if (max_ring > lim || max_col > lim) return false;
  for (std::size_t i = 0; i < n; ++i) {
    dev_lat = std::max(dev_lat, std::abs(lat[i] - (ul.front() + ring_d[i] * step_l)));
    dev_lon = std::max(dev_lon, std::abs(lon[i] - (uo.front() + col_d[i] * step_o)));
  }
  if (dev_lat > tol * step_l) return false;
  if (dev_lon > tol * step_o) return false;
  const double n_col = max_col + 1;
  const double cf = 360 / step_o;
  const double n_col_full = (std::abs(cf - mata_round(cf)) < 1e-6 && mata_round(cf) >= n_col) ? mata_round(cf) : 0;
  g.lat0 = ul.front();
  g.dlat = step_l;
  g.lon0 = uo.front();
  g.dlon = step_o;
  g.n_ring = static_cast<int>(max_ring + 1);
  g.n_col = static_cast<int>(n_col);
  g.n_col_full = static_cast<int>(n_col_full);
  g.ring.resize(n);
  g.col.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    g.ring[i] = static_cast<int>(ring_d[i]);
    g.col[i] = static_cast<int>(col_d[i]);
  }
  return true;
}

// Port of fastconley_choose_grid() (fastconley.mata). `time` is sorted, so
// its runs are the time blocks.
bool choose_lattice(const std::string& method, const std::string& kernel,
                    const arma::vec& lat, const arma::vec& lon, const arma::vec& time,
                    double kk, double cutoff, double tol, Lattice& g) {
  if (method == "pairwise") return false;
  if (!detect_lattice(lat, lon, tol, g)) {
    if (method == "grid") throw UserError("method(grid): no regular lat/lon lattice detected");
    return false;
  }
  const double cells = static_cast<double>(g.n_ring) * static_cast<double>(g.n_col);
  const double n = static_cast<double>(lat.n_elem);
  if (cells > 100 * n || cells > 4e9) {
    if (method == "grid") throw UserError("method(grid): lattice too sparse or degenerate; use method(pairwise)");
    return false;
  }
  if (method == "grid") return true;
  const double pi = 3.14159265358979323846;
  const double ang = cutoff / 6371;
  const double rho_r = std::max(1.0, ang / (g.dlat * pi / 180));
  long double csum = 0;
  for (std::size_t i = 0; i < lat.n_elem; ++i) csum += std::cos(lat[i] * pi / 180);
  const double cosbar = std::max(0.05, static_cast<double>(csum / lat.n_elem));
  const double dbar = std::max(1.0, ang / (g.dlon * pi / 180) / cosbar);
  double sum_nb2 = 0, nblocks = 0;
  for (std::size_t i = 0; i < time.n_elem;) {
    std::size_t j = i + 1;
    while (j < time.n_elem && time[j] == time[i]) ++j;
    const double nb = static_cast<double>(j - i);
    sum_nb2 += nb * nb;
    nblocks += 1;
    i = j;
  }
  const double pairs_est = sum_nb2 * pi * rho_r * dbar / cells / 2;
  double per_pair;
  if (kernel == "uniform") {
    per_pair = g.n_col;
  } else {
    const double npad = std::pow(2.0, std::ceil(std::log(g.n_col + std::min(dbar, static_cast<double>(g.n_col))) / std::log(2.0)));
    per_pair = npad * (5 * std::log(npad) / std::log(2.0) + 2 * kk) / kk;
  }
  const double grid_work = nblocks * g.n_ring * (2 * rho_r + 1) * per_pair;
  return pairs_est > 2 * grid_work;
}

void save_local_num(const char* name, double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  save_local_checked(name, buf);
}

// The vce subcommand (see the protocol comment at the top of the file).
void run_vce(int argc, char* argv[], const std::vector<ST_int>& rows, int nvars) {
  if (argc != 15) throw conley::Error("vce expects 14 arguments");
  const double cutoff = read_scalar(argv[1], "cutoff");
  const std::string kernel = argv[2];
  const std::string dist = argv[3];
  const bool balanced = to_int(argv[4], "balanced") != 0;
  const int threads = to_int(argv[5], "threads");
  const std::string neighbor = argv[6];
  const std::string csr_weight = argv[7];
  const std::string method = argv[8];
  const double grid_tol = read_scalar(argv[9], "grid_tol");
  const double lag = read_scalar(argv[10], "lag");
  const bool haskeys = to_int(argv[11], "haskeys") != 0;
  const bool hasunit = to_int(argv[12], "hasunit") != 0;
  if (method != "auto" && method != "pairwise" && method != "grid") {
    throw conley::Error("vce: unknown method " + method);
  }
  if ((balanced || lag > 0) && !hasunit) throw conley::Error("vce: balanced and lag() need the unit column");
  const int ncoord = 3 + (haskeys ? 2 : 0) + (hasunit ? 1 : 0);
  if (nvars < ncoord + 1) throw conley::Error("vce expects coordinates, time, and at least one score column");
  const std::size_t n = rows.size();
  const std::size_t k = static_cast<std::size_t>(nvars - ncoord);
  if (n == 0) throw conley::Error("vce: no rows");

  ST_int v = 1;
  arma::vec lat(n), lon(n), klat, klon, time(n), unit;
  read_column(v++, rows, lat.memptr());
  read_column(v++, rows, lon.memptr());
  if (haskeys) {
    klat.set_size(n);
    klon.set_size(n);
    read_column(v++, rows, klat.memptr());
    read_column(v++, rows, klon.memptr());
  }
  read_column(v++, rows, time.memptr());
  if (hasunit) {
    unit.set_size(n);
    read_column(v++, rows, unit.memptr());
  }
  arma::mat S(n, k);
  for (std::size_t kk = 0; kk < k; ++kk) read_column(v++, rows, S.colptr(kk));
  const arma::vec& key_lat = haskeys ? klat : lat;
  const arma::vec& key_lon = haskeys ? klon : lon;

  // Aggregation order (time, klat, klon); its time runs also give T, the
  // period count of the full sample (fastconley_prepare_rows computes it
  // before aggregating, which never removes a period).
  const std::vector<std::size_t> agg = stable_key_order(
      {time.memptr(), key_lat.memptr(), key_lon.memptr()}, n, threads);
  std::size_t T = 1;
  for (std::size_t i = 1; i < n; ++i) {
    if (time[agg[i]] != time[agg[i - 1]]) ++T;
  }

  // Balanced-panel validation (same conditions and messages as Mata/R).
  std::vector<std::size_t> tu;
  if (balanced) {
    tu = stable_key_order({time.memptr(), unit.memptr()}, n, threads);
    if (T > 1) {
      if (n % T != 0) throw UserError("balanced requires each period to have the same number of observations");
      const std::size_t n_per = n / T;
      for (std::size_t i = 0; i < n;) {
        std::size_t j = i + 1;
        while (j < n && time[tu[j]] == time[tu[i]]) ++j;
        if (j - i != n_per) throw UserError("balanced requires each period to have the same number of observations");
        i = j;
      }
      for (std::size_t i = 1; i < n_per; ++i) {
        if (unit[tu[i]] == unit[tu[i - 1]]) throw UserError("balanced requires each unit to appear at most once per period");
      }
      for (std::size_t t = 1; t < T; ++t) {
        for (std::size_t i = 0; i < n_per; ++i) {
          if (unit[tu[t * n_per + i]] != unit[tu[i]]) throw UserError("balanced requires every period to contain the same set of units");
        }
      }
      for (std::size_t t = 1; t < T; ++t) {
        for (std::size_t i = 0; i < n_per; ++i) {
          if (lat[tu[t * n_per + i]] != lat[tu[i]] || lon[tu[t * n_per + i]] != lon[tu[i]]) {
            throw UserError("balanced requires time-invariant coordinates per unit");
          }
        }
      }
    }
  }

  // Groups of identical keys in aggregation order.
  std::vector<std::size_t> gstart;
  gstart.reserve(n + 1);
  gstart.push_back(0);
  for (std::size_t i = 1; i < n; ++i) {
    const std::size_t a = agg[i - 1], b = agg[i];
    if (time[a] != time[b] || key_lat[a] != key_lat[b] || key_lon[a] != key_lon[b]) gstart.push_back(i);
  }
  const std::size_t n_groups = gstart.size();
  gstart.push_back(n);

  // Spatial rows: either a row map into S (no merges, or the unaggregated
  // fallback) or a summed score matrix.
  arma::vec sp_lat, sp_lon, sp_time;
  arma::mat S_sum;
  std::vector<std::size_t> sp_map;
  bool use_map = true;
  bool agg_fallback = false;
  bool sp_balanced = false;
  auto gather = [&](const std::vector<std::size_t>& ord, const arma::vec& a, const arma::vec& b) {
    sp_lat.set_size(ord.size());
    sp_lon.set_size(ord.size());
    sp_time.set_size(ord.size());
    for (std::size_t i = 0; i < ord.size(); ++i) {
      sp_lat[i] = a[ord[i]];
      sp_lon[i] = b[ord[i]];
      sp_time[i] = time[ord[i]];
    }
    sp_map = ord;
  };
  if (n_groups == n) {
    gather(agg, key_lat, key_lon);
  } else {
    use_map = false;
    sp_lat.set_size(n_groups);
    sp_lon.set_size(n_groups);
    sp_time.set_size(n_groups);
    S_sum.zeros(n_groups, k);
    for (std::size_t g = 0; g < n_groups; ++g) {
      const std::size_t first = agg[gstart[g]];
      sp_lat[g] = key_lat[first];
      sp_lon[g] = key_lon[first];
      sp_time[g] = time[first];
    }
    for (std::size_t kk = 0; kk < k; ++kk) {
      const double* src = S.colptr(kk);
      double* dst = S_sum.colptr(kk);
      for (std::size_t g = 0; g < n_groups; ++g) {
        double acc = src[agg[gstart[g]]];
        for (std::size_t r = gstart[g] + 1; r < gstart[g + 1]; ++r) acc += src[agg[r]];
        dst[g] = acc;
      }
    }
  }
  if (balanced && T > 1) {
    const std::size_t n_sp = sp_time.n_elem;
    bool ok = (n_sp % T == 0);
    if (ok) {
      const std::size_t n_per = n_sp / T;
      for (std::size_t i = 0; i < n_sp && ok;) {
        std::size_t j = i + 1;
        while (j < n_sp && sp_time[j] == sp_time[i]) ++j;
        if (j - i != n_per) ok = false;
        i = j;
      }
    }
    if (!ok) {
      // Aggregation broke the balance: the unaggregated rows sorted by
      // (time, unit), with the raw coordinates (as fastconley_prepare_rows).
      agg_fallback = true;
      use_map = true;
      S_sum.reset();
      gather(tu, lat, lon);
    }
    sp_balanced = true;
  }

  // Engine choice and the spatial meat.
  Lattice g;
  const bool use_grid = choose_lattice(method, kernel, sp_lat, sp_lon, sp_time,
                                       static_cast<double>(k), cutoff, grid_tol, g);
  arma::mat Msp;
  bool done = false, dateline_fallback = false, unbalanced_fallback = false;
  std::string method_used = "pairwise";
  if (use_grid) {
    arma::mat Sg;
    if (use_map) {
      Sg.set_size(sp_map.size(), k);
      for (std::size_t kk = 0; kk < k; ++kk) {
        const double* src = S.colptr(kk);
        double* dst = Sg.colptr(kk);
        for (std::size_t i = 0; i < sp_map.size(); ++i) dst[i] = src[sp_map[i]];
      }
    }
    try {
      Msp = conley::grid_meat(g.ring.data(), g.col.data(), sp_time, use_map ? Sg : S_sum,
                              g.lat0, g.dlat, g.dlon, g.n_ring, g.n_col, g.n_col_full,
                              cutoff, dist, kernel, threads);
      done = true;
      method_used = "grid";
    } catch (const conley::Error& e) {
      if (method == "auto" && std::string(e.what()).find("dateline") != std::string::npos) {
        dateline_fallback = true;
      } else {
        throw;
      }
    }
  }
  if (!done) {
    Msp = use_map
      ? conley::spatial_meat(sp_lat, sp_lon, sp_time, S, cutoff, kernel, dist, sp_balanced,
                             threads, neighbor, csr_weight, &unbalanced_fallback, &sp_map)
      : conley::spatial_meat(sp_lat, sp_lon, sp_time, S_sum, cutoff, kernel, dist, sp_balanced,
                             threads, neighbor, csr_weight, &unbalanced_fallback);
  }

  // Serial meat on the full sample sorted by (unit, time).
  bool serial_done = false;
  arma::mat Mse;
  if (lag > 0 && T > 1) {
    const std::vector<std::size_t> ut = stable_key_order({unit.memptr(), time.memptr()}, n, threads);
    arma::vec u2(n), t2(n);
    arma::mat S2(n, k);
    for (std::size_t i = 0; i < n; ++i) {
      u2[i] = unit[ut[i]];
      t2[i] = time[ut[i]];
    }
    for (std::size_t kk = 0; kk < k; ++kk) {
      const double* src = S.colptr(kk);
      double* dst = S2.colptr(kk);
      for (std::size_t i = 0; i < n; ++i) dst[i] = src[ut[i]];
    }
    Mse = conley::serial_hac_meat(u2, t2, lag, S2, threads);
    serial_done = true;
  }

  save_local_num("_fc_n_sp", static_cast<double>(sp_time.n_elem));
  save_local_num("_fc_n_periods", static_cast<double>(T));
  save_local_checked("_fc_sp_balanced", sp_balanced ? "1" : "0");
  save_local_checked("_fc_method_used", method_used.c_str());
  save_local_checked("_fc_agg_fallback", agg_fallback ? "1" : "0");
  save_local_checked("_fc_unbalanced_fallback", unbalanced_fallback ? "1" : "0");
  save_local_checked("_fc_dateline_fallback", dateline_fallback ? "1" : "0");
  save_local_checked("_fc_serial_done", serial_done ? "1" : "0");
  store_matrix(argv[13], Msp);
  if (serial_done) store_matrix(argv[14], Mse);
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

    if (cmd == "vce") {
      run_vce(argc, argv, rows, nvars);
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
  } catch (const UserError& e) {
    report_error_noexcept(e.what());
    return 3498;
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
