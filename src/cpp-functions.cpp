#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD 1
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#if defined(RCPP_PARALLEL_USE_TBB) && RCPP_PARALLEL_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr double DE2RA = PI / 180.0;
constexpr double AVG_ERAD = 6371.0;
constexpr double EPS = 1e-14;

enum DistId { DIST_HAVERSINE = 1, DIST_SPHERICAL = 2, DIST_CHORD = 3 };
enum KernelId { KERNEL_BARTLETT = 1, KERNEL_UNIFORM = 2 };

inline double safe_acos(double x) {
  if (x < -1.0) return PI;
  if (x > 1.0) return 0.0;
  return std::acos(x);
}

inline double sq(double x) {
  return x * x;
}

inline double lon_abs_wrapped(double x, double y) {
  double diff = std::fabs(x - y);
  return diff < PI ? diff : TWO_PI - diff;
}

// Deterministic parallel-for over [0, n): body(lo, hi) must write only
// disjoint per-index outputs (pure gather/scatter), so the result is
// independent of the work partition. Small inputs stay serial.
template <typename Body>
struct RangeWorker : public Worker {
  const Body& body;
  explicit RangeWorker(const Body& body) : body(body) {}
  void operator()(std::size_t b, std::size_t e) { body(b, e); }
};

template <typename Body>
void parallel_range(std::size_t n, int ncores, const Body& body) {
  if (ncores > 1 && n > 8192) {
    RangeWorker<Body> w(body);
    parallelFor(0, n, w);
  } else {
    body(0, n);
  }
}

// Sort that uses TBB's parallel_sort when available and worthwhile. The
// comparator must define a strict total order (callers tiebreak on index),
// so the sorted output is unique and identical to std::sort's.
template <typename It, typename Cmp>
void sort_maybe_parallel(It first, It last, Cmp cmp, int ncores) {
#if defined(RCPP_PARALLEL_USE_TBB) && RCPP_PARALLEL_USE_TBB
  if (ncores > 1 && static_cast<std::size_t>(last - first) > 8192) {
    tbb::parallel_sort(first, last, cmp);
    return;
  }
#endif
  (void)ncores;
  std::sort(first, last, cmp);
}

struct CoordCache {
  std::vector<double> lat_rad;
  std::vector<double> lon_rad;
  std::vector<double> cos_lat;
  std::vector<double> sin_lat;
  // 3D unit vectors; always built. The cell grid bins on them and every
  // pair_weight specialization screens with them.
  std::vector<double> x3;
  std::vector<double> y3;
  std::vector<double> z3;
};

int parse_dist_id(const std::string& dist_fn) {
  if (dist_fn == "haversine") return DIST_HAVERSINE;
  if (dist_fn == "spherical") return DIST_SPHERICAL;
  if (dist_fn == "chord") return DIST_CHORD;
  Rcpp::stop("Unknown dist_fn: %s", dist_fn.c_str());
}

int parse_kernel_id(const std::string& kernel) {
  if (kernel == "bartlett") return KERNEL_BARTLETT;
  if (kernel == "uniform") return KERNEL_UNIFORM;
  Rcpp::stop("Unknown kernel: %s", kernel.c_str());
}

struct ScreenParams {
  double lat_cutoff_rad;
  double angular_cutoff_rad;
  double sin2_half_angular_cutoff;
  // cos(cutoff / R): dot(u_i, u_j) >= cos_cutoff iff angular distance <= cutoff/R.
  // Used by SPHERICAL specializations to skip the per-pair acos/cos entirely.
  double cos_cutoff;
  // (cutoff / R)^2, capped at 4. Used by CHORD specializations to compare
  // squared Euclidean distance between unit vectors against the threshold.
  double chord_cutoff_sq;
  // Unit-sphere chord length equivalent to the cutoff: 2*sin(angular/2).
  // All supported distances are monotone in the chord, so a pair is within
  // the cutoff iff |u_i - u_j| <= chord_cell. Used as the cell-grid edge.
  double chord_cell;
};

ScreenParams make_screen_params(double cutoff, int dist_id) {
  ScreenParams p;
  if (cutoff < 0.0) {
    p.lat_cutoff_rad = -1.0;
    p.angular_cutoff_rad = -1.0;
    p.sin2_half_angular_cutoff = 0.0;
    p.cos_cutoff = 1.0;
    p.chord_cutoff_sq = 0.0;
    p.chord_cell = 0.0;
    return p;
  }

  if (dist_id == DIST_CHORD) {
    const double ratio = std::min(1.0, cutoff / (2.0 * AVG_ERAD));
    p.angular_cutoff_rad = 2.0 * std::asin(ratio);
  } else {
    p.angular_cutoff_rad = cutoff / AVG_ERAD;
  }

  if (p.angular_cutoff_rad > PI) p.angular_cutoff_rad = PI;
  p.lat_cutoff_rad = p.angular_cutoff_rad;
  p.sin2_half_angular_cutoff = sq(std::sin(p.angular_cutoff_rad / 2.0));
  p.cos_cutoff = std::cos(p.angular_cutoff_rad);
  // CHORD threshold on unit-vector squared distance: chord_km = R*|u_i - u_j|,
  // so chord_km <= cutoff iff |u_i - u_j|^2 <= (cutoff/R)^2. Cap at 4 (the
  // maximum possible squared distance between unit vectors, antipodal case).
  p.chord_cutoff_sq = std::min(4.0, sq(cutoff / AVG_ERAD));
  p.chord_cell = 2.0 * std::sin(0.5 * p.angular_cutoff_rad);
  return p;
}

// ------------------------------------------------------------------
// pair_weight<D, K>: templated, inlined screen + distance + kernel.
// Returns 0 for rejected pairs and a positive weight otherwise. The
// (D, K) tag is fixed at the top of the spatial routine via an 8-way
// switch in FastSpatialMeat, so each compiled body knows D and K at
// compile time -- no runtime branching, no function-pointer indirection.
// ------------------------------------------------------------------

template<int D, int K>
inline double pair_weight(const CoordCache& c, std::size_t i, std::size_t j,
                          double cutoff, const ScreenParams& screen);

// HAVERSINE × UNIFORM: dot screen + longitude screen, then haversine
// 'a'-threshold (no atan2).
template<>
inline double pair_weight<DIST_HAVERSINE, KERNEL_UNIFORM>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  (void)cutoff;
  // Cheap 3D screen first: in exact arithmetic a = (1 - dot)/2, so
  // dot < cos(cutoff) can never pass the a-test below. The 1e-12 margin
  // keeps the screen conservative; the a-test remains the arbiter, so
  // results are bit-identical with or without this screen.
  {
    const double dot = c.x3[i] * c.x3[j] + c.y3[i] * c.y3[j] + c.z3[i] * c.z3[j];
    if (dot < screen.cos_cutoff - 1e-12) return 0.0;
  }
  const double cos_i = c.cos_lat[i];
  const double cos_j = c.cos_lat[j];
  const double denom = cos_i * cos_j;
  const double dlon = lon_abs_wrapped(c.lon_rad[i], c.lon_rad[j]);
  const double sin_half_dlon = std::sin(dlon / 2.0);
  const double s2_dlon = sin_half_dlon * sin_half_dlon;

  // Multiply-form of the longitude screen: a >= denom * s2_dlon, so
  // denom * s2_dlon > sin2cut implies rejection by the exact a-test below.
  // The absolute 1e-15 margin keeps the screen conservative (denom <= 1, so
  // this admits weakly more candidates than the old divide form); the a-test
  // remains the arbiter, and a multiply is ~4x cheaper than a divide.
  if (screen.angular_cutoff_rad < PI && denom > EPS) {
    if (s2_dlon * denom > screen.sin2_half_angular_cutoff + 1e-15) return 0.0;
  }
  const double dlat = c.lat_rad[j] - c.lat_rad[i];
  const double sin_half_dlat = std::sin(dlat / 2.0);
  const double a = sin_half_dlat * sin_half_dlat + denom * s2_dlon;
  if (a > screen.sin2_half_angular_cutoff) return 0.0;
  return 1.0;
}

// HAVERSINE × BARTLETT: same screens, plus the real distance via atan2.
template<>
inline double pair_weight<DIST_HAVERSINE, KERNEL_BARTLETT>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  // See the UNIFORM specialization: conservative dot screen, a-test arbiter.
  {
    const double dot = c.x3[i] * c.x3[j] + c.y3[i] * c.y3[j] + c.z3[i] * c.z3[j];
    if (dot < screen.cos_cutoff - 1e-12) return 0.0;
  }
  const double cos_i = c.cos_lat[i];
  const double cos_j = c.cos_lat[j];
  const double denom = cos_i * cos_j;
  const double dlon = lon_abs_wrapped(c.lon_rad[i], c.lon_rad[j]);
  const double sin_half_dlon = std::sin(dlon / 2.0);
  const double s2_dlon = sin_half_dlon * sin_half_dlon;

  // Multiply-form of the longitude screen: a >= denom * s2_dlon, so
  // denom * s2_dlon > sin2cut implies rejection by the exact a-test below.
  // The absolute 1e-15 margin keeps the screen conservative (denom <= 1, so
  // this admits weakly more candidates than the old divide form); the a-test
  // remains the arbiter, and a multiply is ~4x cheaper than a divide.
  if (screen.angular_cutoff_rad < PI && denom > EPS) {
    if (s2_dlon * denom > screen.sin2_half_angular_cutoff + 1e-15) return 0.0;
  }
  const double dlat = c.lat_rad[j] - c.lat_rad[i];
  const double sin_half_dlat = std::sin(dlat / 2.0);
  double a = sin_half_dlat * sin_half_dlat + denom * s2_dlon;
  if (a > screen.sin2_half_angular_cutoff) return 0.0;
  if (a < 0.0) a = 0.0;
  if (a > 1.0) a = 1.0;
  const double d = AVG_ERAD * 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  if (cutoff <= 0.0) return d <= 0.0 ? 1.0 : 0.0;
  return 1.0 - d / cutoff;
}

// SPHERICAL × UNIFORM: 3D dot product threshold. No trig in the hot path.
template<>
inline double pair_weight<DIST_SPHERICAL, KERNEL_UNIFORM>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  (void)cutoff;
  const double dot = c.x3[i] * c.x3[j] + c.y3[i] * c.y3[j] + c.z3[i] * c.z3[j];
  return dot >= screen.cos_cutoff ? 1.0 : 0.0;
}

// SPHERICAL × BARTLETT: dot threshold first, then acos only for accepts.
template<>
inline double pair_weight<DIST_SPHERICAL, KERNEL_BARTLETT>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  const double dot = c.x3[i] * c.x3[j] + c.y3[i] * c.y3[j] + c.z3[i] * c.z3[j];
  if (dot < screen.cos_cutoff) return 0.0;
  const double d = AVG_ERAD * safe_acos(dot);
  if (cutoff <= 0.0) return d <= 0.0 ? 1.0 : 0.0;
  return 1.0 - d / cutoff;
}

// CHORD × UNIFORM: squared Euclidean threshold on unit vectors, no sqrt.
template<>
inline double pair_weight<DIST_CHORD, KERNEL_UNIFORM>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  (void)cutoff;
  const double dx = c.x3[i] - c.x3[j];
  const double dy = c.y3[i] - c.y3[j];
  const double dz = c.z3[i] - c.z3[j];
  const double d2 = dx * dx + dy * dy + dz * dz;
  return d2 <= screen.chord_cutoff_sq ? 1.0 : 0.0;
}

// CHORD × BARTLETT: same threshold; sqrt only for accepted pairs.
template<>
inline double pair_weight<DIST_CHORD, KERNEL_BARTLETT>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  const double dx = c.x3[i] - c.x3[j];
  const double dy = c.y3[i] - c.y3[j];
  const double dz = c.z3[i] - c.z3[j];
  const double d2 = dx * dx + dy * dy + dz * dz;
  if (d2 > screen.chord_cutoff_sq) return 0.0;
  const double d = AVG_ERAD * std::sqrt(d2);
  if (cutoff <= 0.0) return d <= 0.0 ? 1.0 : 0.0;
  return 1.0 - d / cutoff;
}

CoordCache make_coord_cache(const arma::vec& lat, const arma::vec& lon,
                            int dist_id, int ncores) {
  const std::size_t n = lat.n_elem;
  // Validate serially first: Rcpp::stop must not be called from worker
  // threads.
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(lat[i]) || !std::isfinite(lon[i])) {
      Rcpp::stop("lat/lon contain non-finite values.");
    }
  }

  CoordCache c;
  c.lat_rad.resize(n);
  c.lon_rad.resize(n);
  c.cos_lat.resize(n);
  c.sin_lat.resize(n);
  // The cell-grid neighbor search bins every distance function by its 3D
  // unit vector, and all pair_weight screens read it — always built.
  (void)dist_id;
  c.x3.resize(n);
  c.y3.resize(n);
  c.z3.resize(n);
  parallel_range(n, ncores, [&](std::size_t lo_i, std::size_t hi_i) {
    for (std::size_t i = lo_i; i < hi_i; ++i) {
      const double la = lat[i] * DE2RA;
      const double lo = lon[i] * DE2RA;
      c.lat_rad[i] = la;
      c.lon_rad[i] = lo;
      c.cos_lat[i] = std::cos(la);
      c.sin_lat[i] = std::sin(la);
      c.x3[i] = c.cos_lat[i] * std::cos(lo);
      c.y3[i] = c.cos_lat[i] * std::sin(lo);
      c.z3[i] = c.sin_lat[i];
    }
  });
  return c;
}

// Row-major flat score buffer, gathered in one pass from the column-major
// score matrix the R side hands over (scores = e * X, possibly aggregated).
// Row `pos` of the buffer is row `perm[pos]` of `Scol`, so the meat workers
// see row-contiguous, permutation-applied reads with no intermediate
// unsorted copy.
struct RowMajorScores {
  std::size_t n;
  std::size_t k;
  std::vector<double> s;

  RowMajorScores(const arma::mat& Scol, const std::vector<std::size_t>& perm,
                 int ncores)
      : n(perm.size()), k(Scol.n_cols),
        s(perm.size() * static_cast<std::size_t>(Scol.n_cols)) {
    const double* base = Scol.memptr();
    const std::size_t ldn = Scol.n_rows;
    const std::size_t kk_n = k;
    parallel_range(n, ncores, [&](std::size_t lo, std::size_t hi) {
      for (std::size_t pos = lo; pos < hi; ++pos) {
        const std::size_t src = perm[pos];
        double* dst = s.data() + pos * kk_n;
        for (std::size_t kk = 0; kk < kk_n; ++kk) dst[kk] = base[kk * ldn + src];
      }
    });
  }

  // Identity order: the caller's rows are already in the right order.
  RowMajorScores(const arma::mat& Scol, int ncores)
      : n(Scol.n_rows), k(Scol.n_cols),
        s(static_cast<std::size_t>(Scol.n_rows) * Scol.n_cols) {
    const double* base = Scol.memptr();
    const std::size_t ldn = Scol.n_rows;
    const std::size_t kk_n = k;
    parallel_range(n, ncores, [&](std::size_t lo, std::size_t hi) {
      for (std::size_t pos = lo; pos < hi; ++pos) {
        double* dst = s.data() + pos * kk_n;
        for (std::size_t kk = 0; kk < kk_n; ++kk) dst[kk] = base[kk * ldn + pos];
      }
    });
  }

  inline const double* row(std::size_t i) const {
    return s.data() + i * k;
  }
};

// ------------------------------------------------------------------
// Deterministic parallel reduction. Rows are processed in fixed-size
// chunks; each chunk accumulates into its own k x k partial, and the
// partials are summed serially in chunk order. The summation order is
// therefore independent of the thread count and of TBB scheduling:
// ncores = 1 and ncores = N give bit-identical results.
// ------------------------------------------------------------------

template <typename Body>
struct ChunkWorker : public Worker {
  const Body& body;
  std::vector<arma::mat>& partials;
  const std::size_t n;
  const std::size_t chunk;

  ChunkWorker(const Body& body, std::vector<arma::mat>& partials,
              std::size_t n, std::size_t chunk)
      : body(body), partials(partials), n(n), chunk(chunk) {}

  void operator()(std::size_t cb, std::size_t ce) {
    for (std::size_t ci = cb; ci < ce; ++ci) {
      const std::size_t lo = ci * chunk;
      const std::size_t hi = std::min(n, lo + chunk);
      body(lo, hi, partials[ci]);
    }
  }
};

// body(lo, hi, meat&) must accumulate rows [lo, hi) into meat and itself be
// deterministic in row order.
template <typename Body>
arma::mat reduce_deterministic(std::size_t n, std::size_t k, std::size_t chunk,
                               int ncores, const Body& body) {
  arma::mat meat(k, k, arma::fill::zeros);
  if (n == 0) return meat;
  const std::size_t nchunks = (n + chunk - 1) / chunk;
  std::vector<arma::mat> partials(nchunks, arma::mat(k, k, arma::fill::zeros));
  ChunkWorker<Body> w(body, partials, n, chunk);
  if (ncores > 1) parallelFor(0, nchunks, w);
  else w(0, nchunks);
  for (std::size_t ci = 0; ci < nchunks; ++ci) meat += partials[ci];
  return meat;
}

// Fixed chunk sizes (rows / unit blocks per reduction task). Constants, so
// results never depend on ncores.
constexpr std::size_t ROW_CHUNK = 1024;
constexpr std::size_t BLOCK_CHUNK = 128;

// Reorder a CoordCache by an index permutation. The output is the same length
// as `perm` and contains `src` values at positions `perm[pos]`. Empty xyz
// vectors in `src` (i.e., dist_id not in {SPHERICAL, CHORD}) are preserved as
// empty in the output.
CoordCache permute_coord_cache(const CoordCache& src,
                               const std::vector<std::size_t>& perm,
                               int ncores) {
  const std::size_t n = perm.size();
  CoordCache out;
  out.lat_rad.resize(n);
  out.lon_rad.resize(n);
  out.cos_lat.resize(n);
  out.sin_lat.resize(n);
  const bool has_xyz = !src.x3.empty();
  if (has_xyz) {
    out.x3.resize(n);
    out.y3.resize(n);
    out.z3.resize(n);
  }
  parallel_range(n, ncores, [&](std::size_t lo, std::size_t hi) {
    for (std::size_t pos = lo; pos < hi; ++pos) {
      const std::size_t i = perm[pos];
      out.lat_rad[pos] = src.lat_rad[i];
      out.lon_rad[pos] = src.lon_rad[i];
      out.cos_lat[pos] = src.cos_lat[i];
      out.sin_lat[pos] = src.sin_lat[i];
      if (has_xyz) {
        out.x3[pos] = src.x3[i];
        out.y3[pos] = src.y3[i];
        out.z3[pos] = src.z3[i];
      }
    }
  });
  return out;
}

struct TimeBlocks {
  std::vector<std::size_t> start;
  std::vector<std::size_t> end;
};

TimeBlocks make_time_blocks(const arma::vec& time) {
  TimeBlocks blocks;
  const std::size_t n = time.n_elem;
  if (n == 0) return blocks;

  blocks.start.push_back(0);
  for (std::size_t i = 1; i < n; ++i) {
    if (time[i] != time[i - 1]) {
      blocks.end.push_back(i);
      blocks.start.push_back(i);
    }
  }
  blocks.end.push_back(n);
  return blocks;
}

bool same_block_size(const TimeBlocks& blocks) {
  if (blocks.start.empty()) return true;
  const std::size_t n0 = blocks.end[0] - blocks.start[0];
  for (std::size_t b = 1; b < blocks.start.size(); ++b) {
    if (blocks.end[b] - blocks.start[b] != n0) return false;
  }
  return true;
}

void sort_block_indices(const std::vector<double>& lat_rad,
                        const std::vector<double>& lon_rad,
                        std::size_t begin, std::size_t end,
                        std::vector<std::size_t>& out) {
  out.resize(end - begin);
  std::iota(out.begin(), out.end(), begin);
  std::sort(out.begin(), out.end(), [&](std::size_t a, std::size_t b) {
    if (lat_rad[a] == lat_rad[b]) return lon_rad[a] < lon_rad[b];
    return lat_rad[a] < lat_rad[b];
  });
}

// col_idx stores within-period sorted positions in [0, n_per); 32-bit halves
// the per-pair index footprint (build_csr guards n_per < 2^32). Weights are
// only populated for KERNEL_BARTLETT — the uniform meat branch never reads
// them — and live in exactly one of `weight` (csr_weight = "double", exact)
// or `weight_f` (csr_weight = "float", 4 bytes/pair, ~6e-8 relative error
// per weight).
struct CsrGraph {
  std::vector<std::size_t> row_ptr;
  std::vector<uint32_t> col_idx;
  std::vector<double> weight;
  std::vector<float> weight_f;
};

template<int D, int K>
struct CsrCountWorkerT : public Worker {
  const std::vector<std::size_t>& sorted_idx;
  const std::vector<std::size_t>& row_end;
  const CoordCache& c;
  std::vector<std::size_t>& counts;
  const double cutoff;
  const ScreenParams screen;

  CsrCountWorkerT(const std::vector<std::size_t>& sorted_idx,
                  const std::vector<std::size_t>& row_end,
                  const CoordCache& c,
                  std::vector<std::size_t>& counts,
                  double cutoff, const ScreenParams& screen)
      : sorted_idx(sorted_idx), row_end(row_end), c(c),
        counts(counts), cutoff(cutoff), screen(screen) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pos = begin; pos < end; ++pos) {
      const std::size_t i = sorted_idx[pos];
      const double lat_i = c.lat_rad[i];
      std::size_t n_found = 0;
      const std::size_t pos_end = row_end[pos];

      for (std::size_t q = pos + 1; q < pos_end; ++q) {
        const std::size_t j = sorted_idx[q];
        const double dlat = c.lat_rad[j] - lat_i;
        if (dlat > screen.lat_cutoff_rad + 1e-15) break;

        const double w = pair_weight<D, K>(c, i, j, cutoff, screen);
        if (w != 0.0) ++n_found;
      }
      counts[pos] = n_found;
    }
  }
};

template<int D, int K>
struct CsrFillWorkerT : public Worker {
  const std::vector<std::size_t>& sorted_idx;
  const std::vector<std::size_t>& row_end;
  const CoordCache& c;
  const std::vector<std::size_t>& row_ptr;
  std::vector<uint32_t>& col_idx;
  std::vector<double>& weight;
  std::vector<float>& weight_f;
  const bool wfloat;
  const double cutoff;
  const ScreenParams screen;

  CsrFillWorkerT(const std::vector<std::size_t>& sorted_idx,
                 const std::vector<std::size_t>& row_end,
                 const CoordCache& c,
                 const std::vector<std::size_t>& row_ptr,
                 std::vector<uint32_t>& col_idx,
                 std::vector<double>& weight,
                 std::vector<float>& weight_f, bool wfloat,
                 double cutoff, const ScreenParams& screen)
      : sorted_idx(sorted_idx), row_end(row_end), c(c),
        row_ptr(row_ptr), col_idx(col_idx),
        weight(weight), weight_f(weight_f), wfloat(wfloat),
        cutoff(cutoff), screen(screen) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pos = begin; pos < end; ++pos) {
      const std::size_t i = sorted_idx[pos];
      const double lat_i = c.lat_rad[i];
      std::size_t out = row_ptr[pos];
      const std::size_t pos_end = row_end[pos];

      for (std::size_t q = pos + 1; q < pos_end; ++q) {
        const std::size_t j = sorted_idx[q];
        const double dlat = c.lat_rad[j] - lat_i;
        if (dlat > screen.lat_cutoff_rad + 1e-15) break;

        const double w = pair_weight<D, K>(c, i, j, cutoff, screen);
        if (w != 0.0) {
          col_idx[out] = static_cast<uint32_t>(j);
          if (K != KERNEL_UNIFORM) {
            if (wfloat) weight_f[out] = static_cast<float>(w);
            else weight[out] = w;
          }
          ++out;
        }
      }
    }
  }
};

template<int D, int K>
CsrGraph build_csr(const std::vector<std::size_t>& sorted_idx,
                   const std::vector<std::size_t>& row_end,
                   const CoordCache& c,
                   double cutoff, bool wfloat, int ncores) {
  const std::size_t n = sorted_idx.size();
  CsrGraph g;
  g.row_ptr.assign(n + 1, 0);

  if (n == 0 || cutoff < 0.0) return g;
  if (n > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
    Rcpp::stop("Balanced CSR path supports at most 2^32 - 1 units per period.");
  }

  std::vector<std::size_t> counts(n, 0);
  const ScreenParams screen = make_screen_params(cutoff, D);

  CsrCountWorkerT<D, K> count_worker(sorted_idx, row_end, c,
                                     counts, cutoff, screen);
  if (ncores > 1) parallelFor(0, n, count_worker);
  else count_worker(0, n);

  for (std::size_t i = 0; i < n; ++i) {
    g.row_ptr[i + 1] = g.row_ptr[i] + counts[i];
  }

  const std::size_t nnz = g.row_ptr[n];
  g.col_idx.assign(nnz, 0);
  if (K != KERNEL_UNIFORM) {
    if (wfloat) g.weight_f.assign(nnz, 0.0f);
    else g.weight.assign(nnz, 0.0);
  }

  CsrFillWorkerT<D, K> fill_worker(sorted_idx, row_end, c,
                                   g.row_ptr, g.col_idx, g.weight,
                                   g.weight_f, wfloat,
                                   cutoff, screen);
  if (ncores > 1) parallelFor(0, n, fill_worker);
  else fill_worker(0, n);

  return g;
}

// ------------------------------------------------------------------
// 3D cell-grid neighbor search (neighbor = "grid", the default).
//
// Points are bucketed by their 3D unit vectors into a cubic grid whose
// edge is >= the unit-sphere chord equivalent of the cutoff. Every
// supported distance is monotone in the chord, so an accepted pair is
// never more than one cell apart per axis. Candidate enumeration is
// therefore output-sensitive (~9/pi candidates per accepted pair,
// independent of geographic extent) instead of the latitude band scan's
// ~2*L_lon/(pi*r). pair_weight<D, K> remains the arbiter for every
// candidate, so accept sets and weights are identical to the band path;
// only summation order differs (FP-level differences in the meat).
//
// Half-pair iteration: rows are sorted by (cell id, original index).
// For a row we scan (a) later rows in its own cell and (b) rows in the
// 13 neighbor cells with lexicographically larger (cx, cy, cz). Those
// 13 cells form five contiguous cell-id intervals (same cx', cy',
// cz' in {cz-1, cz, cz+1}); consecutive occupied cells are contiguous
// in row space, so each interval is a single [begin, end) row range,
// precomputed per cell.
// ------------------------------------------------------------------

struct GridSpec {
  double cell;     // cell edge, unit-sphere chord units
  std::int64_t G;  // cells per axis
};

GridSpec make_grid_spec(const ScreenParams& screen) {
  GridSpec spec;
  // The grid only needs cell >= the true chord cutoff: a slightly larger
  // cell admits a few extra candidates (rejected by pair_weight), never
  // loses a pair. The margin covers FP error in the chord/dot/haversine
  // accept tests relative to the binning coordinates.
  spec.cell = screen.chord_cell * (1.0 + 1e-7) + 1e-13;
  const double gf = 2.0 / spec.cell;
  std::int64_t G = gf >= 1.0 ? static_cast<std::int64_t>(gf) : 1;
  // Cap so cid = ((cx*G)+cy)*G+cz fits 64 bits (G^3 - 1 <= 2^63 - 1). Only
  // binds for sub-metre cutoffs, where cells merely become larger than
  // strictly needed -- still exact, just more candidates.
  const std::int64_t G_MAX = 2097152;  // 2^21
  if (G > G_MAX) G = G_MAX;
  spec.G = G;
  return spec;
}

// All row positions are in the grid-sorted order the caller permutes the
// coord cache / scores into. Cells never span time blocks, and consecutive
// blocks are contiguous, so cell_start[cell + 1] is always the end of
// `cell`'s rows even at block boundaries.
struct CellGrid {
  std::vector<std::size_t> cell_start;   // n_cells + 1 (caller appends total n)
  std::vector<std::uint32_t> row_cell;   // per sorted row: its cell index
  std::vector<std::size_t> nbr;          // 10 * n_cells: five [begin, end) ranges
};

// Bucket one time block, append its rows (original indices, grid order) to
// `sorted_idx` and its cells to `grid`. `grid.row_cell` must be pre-sized to
// the total row count.
void append_block_grid(const CoordCache& c, std::size_t bs, std::size_t be,
                       const GridSpec& spec,
                       std::vector<std::size_t>& sorted_idx, CellGrid& grid,
                       int ncores) {
  const std::size_t nb = be - bs;
  const double cell = spec.cell;
  const std::int64_t G = spec.G;
  const std::uint64_t uG = static_cast<std::uint64_t>(G);

  std::vector<std::uint64_t> cid(nb);
  parallel_range(nb, ncores, [&](std::size_t lo, std::size_t hi) {
    for (std::size_t i = lo; i < hi; ++i) {
      const std::size_t idx = bs + i;
      std::int64_t cx = static_cast<std::int64_t>((c.x3[idx] + 1.0) / cell);
      std::int64_t cy = static_cast<std::int64_t>((c.y3[idx] + 1.0) / cell);
      std::int64_t cz = static_cast<std::int64_t>((c.z3[idx] + 1.0) / cell);
      if (cx < 0) cx = 0; else if (cx >= G) cx = G - 1;
      if (cy < 0) cy = 0; else if (cy >= G) cy = G - 1;
      if (cz < 0) cz = 0; else if (cz >= G) cz = G - 1;
      cid[i] = (static_cast<std::uint64_t>(cx) * uG +
                static_cast<std::uint64_t>(cy)) * uG +
               static_cast<std::uint64_t>(cz);
    }
  });

  // (cell id, original index) order: a strict total order, so the sorted
  // output is unique and deterministic (parallel or not), and rows of a
  // cell are contiguous with consecutive cells contiguous in row space.
  std::vector<std::size_t> ord(nb);
  std::iota(ord.begin(), ord.end(), 0);
  sort_maybe_parallel(ord.begin(), ord.end(),
                      [&](std::size_t a, std::size_t b) {
                        if (cid[a] != cid[b]) return cid[a] < cid[b];
                        return a < b;
                      },
                      ncores);

  const std::size_t row0 = sorted_idx.size();
  std::vector<std::uint64_t> ucid;
  std::vector<std::size_t> ustart;  // block-local row offsets per cell
  for (std::size_t j = 0; j < nb; ++j) {
    if (j == 0 || cid[ord[j]] != cid[ord[j - 1]]) {
      ucid.push_back(cid[ord[j]]);
      ustart.push_back(j);
    }
    sorted_idx.push_back(bs + ord[j]);
  }
  ustart.push_back(nb);

  const std::size_t ncb = ucid.size();
  const std::size_t cell0 = grid.cell_start.size();
  for (std::size_t cdx = 0; cdx < ncb; ++cdx) {
    grid.cell_start.push_back(row0 + ustart[cdx]);
    for (std::size_t j = ustart[cdx]; j < ustart[cdx + 1]; ++j) {
      grid.row_cell[row0 + j] = static_cast<std::uint32_t>(cell0 + cdx);
    }
  }

  // Five forward cell-id intervals per cell; each maps to one contiguous
  // row range because consecutive occupied cells are contiguous rows.
  // Pure per-cell output: safe and deterministic to parallelize.
  grid.nbr.resize(grid.nbr.size() + 10 * ncb);
  parallel_range(ncb, ncores, [&](std::size_t cdx_lo, std::size_t cdx_hi) {
  for (std::size_t cdx = cdx_lo; cdx < cdx_hi; ++cdx) {
    const std::uint64_t cu = ucid[cdx];
    const std::int64_t cz = static_cast<std::int64_t>(cu % uG);
    const std::int64_t cy = static_cast<std::int64_t>((cu / uG) % uG);
    const std::int64_t cx = static_cast<std::int64_t>(cu / (uG * uG));
    std::size_t* out = &grid.nbr[10 * (cell0 + cdx)];
    int run = 0;
    auto add_run = [&](std::int64_t nx, std::int64_t ny,
                       std::int64_t zlo, std::int64_t zhi) {
      if (zlo < 0) zlo = 0;
      if (zhi > G - 1) zhi = G - 1;
      if (nx < 0 || nx >= G || ny < 0 || ny >= G || zlo > zhi) {
        out[2 * run] = out[2 * run + 1] = row0;
        ++run;
        return;
      }
      const std::uint64_t base = (static_cast<std::uint64_t>(nx) * uG +
                                  static_cast<std::uint64_t>(ny)) * uG;
      const std::uint64_t lo = base + static_cast<std::uint64_t>(zlo);
      const std::uint64_t hi = base + static_cast<std::uint64_t>(zhi);
      const std::size_t clo = std::lower_bound(ucid.begin(), ucid.end(), lo) -
                              ucid.begin();
      const std::size_t chi = std::upper_bound(ucid.begin() + clo, ucid.end(), hi) -
                              ucid.begin();
      out[2 * run] = row0 + ustart[clo];
      out[2 * run + 1] = row0 + ustart[chi];
      ++run;
    };
    add_run(cx,     cy,     cz + 1, cz + 1);
    add_run(cx,     cy + 1, cz - 1, cz + 1);
    add_run(cx + 1, cy - 1, cz - 1, cz + 1);
    add_run(cx + 1, cy,     cz - 1, cz + 1);
    add_run(cx + 1, cy + 1, cz - 1, cz + 1);
  }
  });
}

// Invoke fn(q) for every forward candidate row of `pos` (own-cell rows after
// pos, then the five precomputed neighbor row ranges). Each unordered pair
// is visited exactly once across all rows.
template <typename F>
inline void for_each_grid_candidate(const CellGrid& grid, std::size_t pos, F fn) {
  const std::uint32_t cell = grid.row_cell[pos];
  const std::size_t own_end = grid.cell_start[cell + 1];
  for (std::size_t q = pos + 1; q < own_end; ++q) fn(q);
  const std::size_t* r = &grid.nbr[10 * static_cast<std::size_t>(cell)];
  for (int s = 0; s < 5; ++s) {
    for (std::size_t q = r[2 * s]; q < r[2 * s + 1]; ++q) fn(q);
  }
}

// Grid analogue of meat_stream_band: fused candidate scan + score
// accumulation, O(n) memory, no CSR. Inputs are pre-permuted to grid order.
template<int D, int K>
arma::mat meat_stream_grid(const RowMajorScores& S, const CoordCache& coord,
                           const CellGrid& grid, double cutoff,
                           const ScreenParams& screen, int ncores) {
  const std::size_t k = S.k;
  auto body = [&](std::size_t lo, std::size_t hi, arma::mat& meat) {
    std::vector<double> c(k, 0.0);
    for (std::size_t pos = lo; pos < hi; ++pos) {
      const double* si = S.row(pos);
      for (std::size_t kk = 0; kk < k; ++kk) {
        c[kk] = 0.5 * si[kk];
      }

      for_each_grid_candidate(grid, pos, [&](std::size_t q) {
        const double w = pair_weight<D, K>(coord, pos, q, cutoff, screen);
        if (w == 0.0) return;
        const double* sj = S.row(q);
        if (K == KERNEL_UNIFORM) {
          for (std::size_t kk = 0; kk < k; ++kk) c[kk] += sj[kk];
        } else {
          for (std::size_t kk = 0; kk < k; ++kk) c[kk] += w * sj[kk];
        }
      });

      for (std::size_t k1 = 0; k1 < k; ++k1) {
        const double s1 = si[k1];
        for (std::size_t k2 = 0; k2 < k; ++k2) {
          meat(k1, k2) += s1 * c[k2];
        }
      }
    }
  };
  arma::mat meat = reduce_deterministic(S.n, k, ROW_CHUNK, ncores, body);
  return meat + meat.t();
}

// Grid analogues of the CSR count/fill workers (balanced path).
template<int D, int K>
struct GridCsrCountWorkerT : public Worker {
  const CellGrid& grid;
  const CoordCache& c;
  std::vector<std::size_t>& counts;
  const double cutoff;
  const ScreenParams screen;

  GridCsrCountWorkerT(const CellGrid& grid, const CoordCache& c,
                      std::vector<std::size_t>& counts,
                      double cutoff, const ScreenParams& screen)
      : grid(grid), c(c), counts(counts), cutoff(cutoff), screen(screen) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pos = begin; pos < end; ++pos) {
      std::size_t n_found = 0;
      for_each_grid_candidate(grid, pos, [&](std::size_t q) {
        if (pair_weight<D, K>(c, pos, q, cutoff, screen) != 0.0) ++n_found;
      });
      counts[pos] = n_found;
    }
  }
};

template<int D, int K>
struct GridCsrFillWorkerT : public Worker {
  const CellGrid& grid;
  const CoordCache& c;
  const std::vector<std::size_t>& row_ptr;
  std::vector<std::uint32_t>& col_idx;
  std::vector<double>& weight;
  std::vector<float>& weight_f;
  const bool wfloat;
  const double cutoff;
  const ScreenParams screen;

  GridCsrFillWorkerT(const CellGrid& grid, const CoordCache& c,
                     const std::vector<std::size_t>& row_ptr,
                     std::vector<std::uint32_t>& col_idx,
                     std::vector<double>& weight,
                     std::vector<float>& weight_f, bool wfloat,
                     double cutoff, const ScreenParams& screen)
      : grid(grid), c(c), row_ptr(row_ptr), col_idx(col_idx), weight(weight),
        weight_f(weight_f), wfloat(wfloat), cutoff(cutoff), screen(screen) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pos = begin; pos < end; ++pos) {
      std::size_t out = row_ptr[pos];
      for_each_grid_candidate(grid, pos, [&](std::size_t q) {
        const double w = pair_weight<D, K>(c, pos, q, cutoff, screen);
        if (w != 0.0) {
          col_idx[out] = static_cast<std::uint32_t>(q);
          if (K != KERNEL_UNIFORM) {
            if (wfloat) weight_f[out] = static_cast<float>(w);
            else weight[out] = w;
          }
          ++out;
        }
      });
    }
  }
};

// Build the balanced-path CSR by scanning grid candidates. `c` must be the
// block-0 coord cache permuted into grid order; col_idx values are sorted
// positions in [0, n_per), exactly like build_csr.
template<int D, int K>
CsrGraph build_csr_grid(const CellGrid& grid, const CoordCache& c,
                        double cutoff, const ScreenParams& screen,
                        bool wfloat, int ncores) {
  const std::size_t n = grid.row_cell.size();
  CsrGraph g;
  g.row_ptr.assign(n + 1, 0);
  if (n == 0 || cutoff < 0.0) return g;
  if (n > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
    Rcpp::stop("Balanced CSR path supports at most 2^32 - 1 units per period.");
  }

  std::vector<std::size_t> counts(n, 0);
  GridCsrCountWorkerT<D, K> count_worker(grid, c, counts, cutoff, screen);
  if (ncores > 1) parallelFor(0, n, count_worker);
  else count_worker(0, n);

  for (std::size_t i = 0; i < n; ++i) {
    g.row_ptr[i + 1] = g.row_ptr[i] + counts[i];
  }

  const std::size_t nnz = g.row_ptr[n];
  g.col_idx.assign(nnz, 0);
  if (K != KERNEL_UNIFORM) {
    if (wfloat) g.weight_f.assign(nnz, 0.0f);
    else g.weight.assign(nnz, 0.0);
  }

  GridCsrFillWorkerT<D, K> fill_worker(grid, c, g.row_ptr, g.col_idx,
                                       g.weight, g.weight_f, wfloat,
                                       cutoff, screen);
  if (ncores > 1) parallelFor(0, n, fill_worker);
  else fill_worker(0, n);

  return g;
}

// Inputs `S` and `coord` are pre-permuted into lat-sorted order by the
// caller, so rows are indexed by sorted position `pos` directly. This makes
// every read of `S.row(...)` and `coord.lat_rad[...]` sequential across the
// outer loop and within the inner lat window.
template<int D, int K>
arma::mat meat_stream_band(const RowMajorScores& S,
                           const std::vector<std::size_t>& row_end,
                           const CoordCache& coord,
                           double cutoff, const ScreenParams& screen,
                           int ncores) {
  const std::size_t k = S.k;
  auto body = [&](std::size_t lo, std::size_t hi, arma::mat& meat) {
    std::vector<double> c(k, 0.0);
    for (std::size_t pos = lo; pos < hi; ++pos) {
      const double* si = S.row(pos);
      const double lat_i = coord.lat_rad[pos];

      for (std::size_t kk = 0; kk < k; ++kk) {
        c[kk] = 0.5 * si[kk];
      }

      const std::size_t pos_end = row_end[pos];
      for (std::size_t q = pos + 1; q < pos_end; ++q) {
        const double dlat = coord.lat_rad[q] - lat_i;
        if (dlat > screen.lat_cutoff_rad + 1e-15) break;

        const double w = pair_weight<D, K>(coord, pos, q, cutoff, screen);
        if (w == 0.0) continue;

        const double* sj = S.row(q);
        // K == UNIFORM: w is always 1.0 past the w==0 filter, so skip the
        // per-element multiply. Compile-time dead-coded.
        if (K == KERNEL_UNIFORM) {
          for (std::size_t kk = 0; kk < k; ++kk) {
            c[kk] += sj[kk];
          }
        } else {
          for (std::size_t kk = 0; kk < k; ++kk) {
            c[kk] += w * sj[kk];
          }
        }
      }

      for (std::size_t k1 = 0; k1 < k; ++k1) {
        const double s1 = si[k1];
        for (std::size_t k2 = 0; k2 < k; ++k2) {
          meat(k1, k2) += s1 * c[k2];
        }
      }
    }
  };
  arma::mat meat = reduce_deterministic(S.n, k, ROW_CHUNK, ncores, body);
  return meat + meat.t();
}

// Balanced meat over the period-major sorted scores (row `base + pos` is the
// unit at sorted position `pos` in period `b`). The CSR is re-streamed once
// per period -- benchmarking showed this beats a unit-major T*k stacked
// layout: spatial sorting makes the neighbor gathers a sliding window of
// k-wide rows that stays L2-resident per period, whereas any wider stacking
// pushes the window past L2 and thrashes the shared L3 at high thread
// counts. (Tried in the v0.6.0 cycle; see notes/OPTIMIZATION_PLAN.md.)
template<int K, typename WT>
arma::mat meat_from_csr_periodmajor(const RowMajorScores& S,
                                    const std::vector<std::size_t>& block_start,
                                    std::size_t n_per,
                                    const std::vector<std::size_t>& row_ptr,
                                    const std::vector<std::uint32_t>& col_idx,
                                    const std::vector<WT>& weight,
                                    int ncores) {
  const std::size_t k = S.k;
  const std::size_t T = block_start.size();
  auto body = [&](std::size_t lo, std::size_t hi, arma::mat& meat) {
    std::vector<double> c(k, 0.0);
    for (std::size_t b = 0; b < T; ++b) {
      const std::size_t base = block_start[b];
      for (std::size_t pos = lo; pos < hi; ++pos) {
        const double* si = S.row(base + pos);
        for (std::size_t kk = 0; kk < k; ++kk) {
          c[kk] = 0.5 * si[kk];
        }
        for (std::size_t ep = row_ptr[pos]; ep < row_ptr[pos + 1]; ++ep) {
          const double* sj = S.row(base + static_cast<std::size_t>(col_idx[ep]));
          if (K == KERNEL_UNIFORM) {
            for (std::size_t kk = 0; kk < k; ++kk) c[kk] += sj[kk];
          } else {
            const double wgt = static_cast<double>(weight[ep]);
            for (std::size_t kk = 0; kk < k; ++kk) c[kk] += wgt * sj[kk];
          }
        }
        for (std::size_t k1 = 0; k1 < k; ++k1) {
          const double s1 = si[k1];
          for (std::size_t k2 = 0; k2 < k; ++k2) {
            meat(k1, k2) += s1 * c[k2];
          }
        }
      }
    }
  };
  arma::mat meat = reduce_deterministic(n_per, k, ROW_CHUNK, ncores, body);
  return meat + meat.t();
}

// Runtime csr_weight dispatch shared by both balanced variants.
template<int K>
arma::mat meat_from_csr_dispatch(const RowMajorScores& S,
                                 const std::vector<std::size_t>& block_start,
                                 std::size_t n_per, const CsrGraph& graph,
                                 bool wfloat, int ncores) {
  if (K != KERNEL_UNIFORM && wfloat) {
    return meat_from_csr_periodmajor<K, float>(S, block_start, n_per,
                                               graph.row_ptr, graph.col_idx,
                                               graph.weight_f, ncores);
  }
  return meat_from_csr_periodmajor<K, double>(S, block_start, n_per,
                                              graph.row_ptr, graph.col_idx,
                                              graph.weight, ncores);
}

template<int D, int K>
arma::mat fast_spatial_general(const arma::vec& lat, const arma::vec& lon,
                               const arma::vec& time, const arma::mat& S_col,
                               double cutoff, int ncores) {
  const std::size_t n = S_col.n_rows;
  const CoordCache c = make_coord_cache(lat, lon, D, ncores);

  const TimeBlocks blocks = make_time_blocks(time);
  std::vector<std::size_t> sorted_idx;
  sorted_idx.reserve(n);
  std::vector<std::size_t> row_end;
  row_end.reserve(n);

  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    std::vector<std::size_t> block_sorted;
    sort_block_indices(c.lat_rad, c.lon_rad, blocks.start[b], blocks.end[b], block_sorted);
    sorted_idx.insert(sorted_idx.end(), block_sorted.begin(), block_sorted.end());
    row_end.resize(sorted_idx.size(), sorted_idx.size());
  }

  // Reorder coord cache and scores into lat-sorted order. The meat loop then
  // indexes both buffers by sorted position directly, so every read is
  // sequential -- independent of how the caller laid out the input.
  const CoordCache c_sorted = permute_coord_cache(c, sorted_idx, ncores);
  const RowMajorScores S_sorted(S_col, sorted_idx, ncores);

  const ScreenParams screen = make_screen_params(cutoff, D);
  return meat_stream_band<D, K>(S_sorted, row_end, c_sorted, cutoff, screen, ncores);
}

template<int D, int K>
arma::mat fast_spatial_balanced(const arma::vec& lat, const arma::vec& lon,
                                const arma::vec& time, const arma::mat& S_col,
                                double cutoff, bool wfloat, int ncores) {
  const TimeBlocks blocks = make_time_blocks(time);
  const std::size_t n_per = blocks.end[0] - blocks.start[0];
  const std::size_t T = blocks.start.size();
  const CoordCache c = make_coord_cache(lat, lon, D, ncores);

  // Sort block 0 once. Because coordinates are time-invariant in the balanced
  // path, this same permutation applies to every block.
  std::vector<std::size_t> sorted_abs;
  sort_block_indices(c.lat_rad, c.lon_rad, blocks.start[0], blocks.end[0], sorted_abs);
  std::vector<std::size_t> sorted_rel(n_per);
  for (std::size_t i = 0; i < n_per; ++i) {
    sorted_rel[i] = sorted_abs[i] - blocks.start[0];
  }

  // Permute block 0's coord cache into sorted order so build_csr can run on
  // identity indices and store col_idx values in sorted-position space.
  const CoordCache c_block0_sorted = permute_coord_cache(c, sorted_abs, ncores);
  std::vector<std::size_t> identity_perm(n_per);
  std::iota(identity_perm.begin(), identity_perm.end(), 0);
  std::vector<std::size_t> row_end(n_per, n_per);
  CsrGraph graph = build_csr<D, K>(identity_perm, row_end, c_block0_sorted,
                                   cutoff, wfloat, ncores);

  std::vector<std::size_t> global_perm(S_col.n_rows);
  for (std::size_t b = 0; b < T; ++b) {
    const std::size_t base = blocks.start[b];
    for (std::size_t pos = 0; pos < n_per; ++pos) {
      global_perm[base + pos] = base + sorted_rel[pos];
    }
  }
  const RowMajorScores S_sorted(S_col, global_perm, ncores);
  return meat_from_csr_dispatch<K>(S_sorted, blocks.start, n_per, graph,
                                   wfloat, ncores);
}

// Grid version of the general path: per-block cell grids, fused scan +
// accumulate, O(n) memory.
template<int D, int K>
arma::mat fast_spatial_general_grid(const arma::vec& lat, const arma::vec& lon,
                                    const arma::vec& time, const arma::mat& S_col,
                                    double cutoff, int ncores) {
  const std::size_t n = S_col.n_rows;
  if (n > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
    Rcpp::stop("The grid neighbor path supports at most 2^32 - 1 rows; "
               "use neighbor = \"band\".");
  }
  const CoordCache c = make_coord_cache(lat, lon, D, ncores);
  const ScreenParams screen = make_screen_params(cutoff, D);
  const GridSpec spec = make_grid_spec(screen);

  const TimeBlocks blocks = make_time_blocks(time);
  std::vector<std::size_t> sorted_idx;
  sorted_idx.reserve(n);
  CellGrid grid;
  grid.row_cell.resize(n);
  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    append_block_grid(c, blocks.start[b], blocks.end[b], spec, sorted_idx, grid, ncores);
  }
  grid.cell_start.push_back(n);

  const CoordCache c_sorted = permute_coord_cache(c, sorted_idx, ncores);
  const RowMajorScores S_sorted(S_col, sorted_idx, ncores);

  return meat_stream_grid<D, K>(S_sorted, c_sorted, grid, cutoff, screen, ncores);
}

// Grid version of the balanced path: cell grid + CSR from block 0, reused
// across periods; the period-major meat pass is shared with the band version.
template<int D, int K>
arma::mat fast_spatial_balanced_grid(const arma::vec& lat, const arma::vec& lon,
                                     const arma::vec& time, const arma::mat& S_col,
                                     double cutoff, bool wfloat, int ncores) {
  const TimeBlocks blocks = make_time_blocks(time);
  const std::size_t n_per = blocks.end[0] - blocks.start[0];
  const std::size_t T = blocks.start.size();
  const CoordCache c = make_coord_cache(lat, lon, D, ncores);
  const ScreenParams screen = make_screen_params(cutoff, D);
  const GridSpec spec = make_grid_spec(screen);

  std::vector<std::size_t> sorted_abs;
  sorted_abs.reserve(n_per);
  CellGrid g0;
  g0.row_cell.resize(n_per);
  append_block_grid(c, blocks.start[0], blocks.end[0], spec, sorted_abs, g0, ncores);
  g0.cell_start.push_back(n_per);

  const CoordCache c0_sorted = permute_coord_cache(c, sorted_abs, ncores);
  CsrGraph graph = build_csr_grid<D, K>(g0, c0_sorted, cutoff, screen,
                                        wfloat, ncores);

  std::vector<std::size_t> sorted_rel(n_per);
  for (std::size_t i = 0; i < n_per; ++i) {
    sorted_rel[i] = sorted_abs[i] - blocks.start[0];
  }
  std::vector<std::size_t> global_perm(S_col.n_rows);
  for (std::size_t b = 0; b < T; ++b) {
    const std::size_t base = blocks.start[b];
    for (std::size_t pos = 0; pos < n_per; ++pos) {
      global_perm[base + pos] = base + sorted_rel[pos];
    }
  }
  const RowMajorScores S_sorted(S_col, global_perm, ncores);
  return meat_from_csr_dispatch<K>(S_sorted, blocks.start, n_per, graph,
                                   wfloat, ncores);
}

// 8-way dispatch on (dist_id, kernel_id) -> template instantiation. Called
// once per FastSpatialMeat invocation; after this point everything is
// compile-time specialized.
template<int D, int K>
inline arma::mat fast_spatial_dispatch(const arma::vec& lat, const arma::vec& lon,
                                       const arma::vec& time, const arma::mat& S_col,
                                       double cutoff, int ncores, bool balanced,
                                       bool use_grid, bool wfloat) {
  if (balanced) {
    if (use_grid) return fast_spatial_balanced_grid<D, K>(lat, lon, time, S_col, cutoff, wfloat, ncores);
    return fast_spatial_balanced<D, K>(lat, lon, time, S_col, cutoff, wfloat, ncores);
  }
  if (use_grid) return fast_spatial_general_grid<D, K>(lat, lon, time, S_col, cutoff, ncores);
  return fast_spatial_general<D, K>(lat, lon, time, S_col, cutoff, ncores);
}

arma::mat dispatch_spatial(int dist_id, int kernel_id,
                           const arma::vec& lat, const arma::vec& lon,
                           const arma::vec& time, const arma::mat& S_col,
                           double cutoff, int ncores, bool balanced,
                           bool use_grid, bool wfloat) {
  const int tag = (dist_id << 4) | kernel_id;
  switch (tag) {
    case (DIST_HAVERSINE << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_HAVERSINE, KERNEL_UNIFORM>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
    case (DIST_HAVERSINE << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_HAVERSINE, KERNEL_BARTLETT>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
    case (DIST_SPHERICAL << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_SPHERICAL, KERNEL_UNIFORM>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
    case (DIST_SPHERICAL << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_SPHERICAL, KERNEL_BARTLETT>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
    case (DIST_CHORD << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_CHORD, KERNEL_UNIFORM>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
    case (DIST_CHORD << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_CHORD, KERNEL_BARTLETT>(lat, lon, time, S_col, cutoff, ncores, balanced, use_grid, wfloat);
  }
  Rcpp::stop("Unsupported (dist_id, kernel_id) combination.");
}

struct UnitBlocks {
  std::vector<std::size_t> start;
  std::vector<std::size_t> end;
};

UnitBlocks make_unit_blocks(const arma::vec& unit) {
  UnitBlocks b;
  const std::size_t n = unit.n_elem;
  if (n == 0) return b;
  b.start.push_back(0);
  for (std::size_t i = 1; i < n; ++i) {
    if (unit[i] != unit[i - 1]) {
      b.end.push_back(i);
      b.start.push_back(i);
    }
  }
  b.end.push_back(n);
  return b;
}

// Serial (Newey-West in time) meat over contiguous unit blocks. Rows must
// be sorted by time within each unit block (the R layer sorts by
// (unit, time)); a guard below enforces it. O(T_u * k) per unit instead of
// O(T_u^2 * k): the Bartlett weight decomposes as
//   1 - (t_j - t_i)/(L+1) = (1 + t_i'/(L+1)) - t_j'/(L+1)
// (t' = t - block base, an exact shift that keeps the running sums well
// conditioned), so the forward-window contribution is
//   c_i = (1 + t_i'/(L+1)) * A - B/(L+1),  A = sum s_j,  B = sum t_j' s_j
// over the sliding window {j : t_i < t_j <= t_i + L}, maintained with two
// pointers as i advances. Pairs with dt = 0 (including ties) are excluded,
// matching the old per-pair loop; meat = M + M' covers both directions.
arma::mat serial_hac_panel(const arma::vec& times, double cutoff,
                           const RowMajorScores& S, const UnitBlocks& blocks,
                           int ncores) {
  const std::size_t k = S.k;
  if (cutoff < 0.0) return arma::mat(k, k, arma::fill::zeros);

  for (std::size_t bi = 0; bi < blocks.start.size(); ++bi) {
    for (std::size_t i = blocks.start[bi] + 1; i < blocks.end[bi]; ++i) {
      if (times[i] < times[i - 1]) {
        Rcpp::stop("FastSerialHacPanel requires rows sorted by time within "
                   "each unit block.");
      }
    }
  }

  const double Lp1 = cutoff + 1.0;
  const double inv = 1.0 / Lp1;
  auto body = [&](std::size_t blo, std::size_t bhi, arma::mat& meat) {
    std::vector<double> A(k, 0.0), B(k, 0.0), c(k, 0.0);
    for (std::size_t bi = blo; bi < bhi; ++bi) {
      const std::size_t bs = blocks.start[bi];
      const std::size_t be = blocks.end[bi];
      const double tb = times[bs];
      std::fill(A.begin(), A.end(), 0.0);
      std::fill(B.begin(), B.end(), 0.0);
      std::size_t lo = bs, hi = bs;
      for (std::size_t i = bs; i < be; ++i) {
        const double ti = times[i] - tb;
        while (hi < be && times[hi] - tb <= ti + cutoff) {
          const double tj = times[hi] - tb;
          const double* sj = S.row(hi);
          for (std::size_t kk = 0; kk < k; ++kk) {
            A[kk] += sj[kk];
            B[kk] += tj * sj[kk];
          }
          ++hi;
        }
        while (lo < be && times[lo] - tb <= ti) {
          const double tj = times[lo] - tb;
          const double* sj = S.row(lo);
          for (std::size_t kk = 0; kk < k; ++kk) {
            A[kk] -= sj[kk];
            B[kk] -= tj * sj[kk];
          }
          ++lo;
        }
        const double coef = 1.0 + ti * inv;
        for (std::size_t kk = 0; kk < k; ++kk) {
          c[kk] = coef * A[kk] - inv * B[kk];
        }
        const double* si = S.row(i);
        for (std::size_t k1 = 0; k1 < k; ++k1) {
          const double s1 = si[k1];
          for (std::size_t k2 = 0; k2 < k; ++k2) {
            meat(k1, k2) += s1 * c[k2];
          }
        }
      }
    }
  };
  arma::mat meat = reduce_deterministic(blocks.start.size(), k, BLOCK_CHUNK,
                                        ncores, body);
  return meat + meat.t();
}

// ------------------------------------------------------------------
// Grid-native exact meat (workstream C2; see notes/OPTIMIZATION_PLAN.md
// and the C1 prototype tests/manual/proto-conv-c1.R).
//
// For the UNIFORM kernel on a regular lat/lon lattice the accept set
// between two latitude rings is a longitude-index interval, so the inner
// sum over a ring pair is a sliding-window sum over per-ring prefix sums:
// O(n_ring * window * n_col * k) total, independent of the pair count.
// Exact: the dot-product accept threshold is the same constant the
// pairwise engine uses, so the answer agrees to FP summation order.
//
// Layout: one dense prefix tensor per time block, ring-major:
// P[r * (C+1) * k + c * k + kk] = sum of scores over columns < c of ring
// r (k-wide rows). Raw cell scores are recovered as adjacent differences,
// so no separate dense score copy is held. Memory: (C+1) * R * k doubles.
// ------------------------------------------------------------------

struct GridGeom {
  std::size_t n_ring;
  std::size_t n_col;
  double dlam_rad;
  std::vector<double> sphi;  // sin(lat) per ring
  std::vector<double> cphi;  // cos(lat) per ring
};

// Largest lon-index offset accepted between rings r1 and r2 (-1 = none),
// matching the pairwise accept "dot >= coscut" with a boundary fix-up in
// the same arithmetic.
inline long grid_delta_star(const GridGeom& gg, std::size_t r1, std::size_t r2,
                            double coscut) {
  const double base = gg.sphi[r1] * gg.sphi[r2];
  const double cc12 = gg.cphi[r1] * gg.cphi[r2];
  const long dmaxc = static_cast<long>(gg.n_col) - 1;
  if (cc12 < 1e-300) {
    return (base >= coscut) ? dmaxc : -1;
  }
  const double rhs = (coscut - base) / cc12;
  if (rhs > 1.0) return -1;
  long ds;
  if (rhs <= -1.0) {
    ds = dmaxc;
  } else {
    ds = static_cast<long>(std::acos(rhs) / gg.dlam_rad) + 1;
    if (ds > dmaxc) ds = dmaxc;
  }
  while (ds >= 0 && base + cc12 * std::cos(ds * gg.dlam_rad) < coscut) --ds;
  while (ds + 1 <= dmaxc &&
         base + cc12 * std::cos((ds + 1) * gg.dlam_rad) >= coscut) ++ds;
  return ds;
}

} // anonymous namespace

// Entry points take pre-computed scores (scores = e * X, possibly
// pre-aggregated by the R layer) and alias the R memory directly via the
// Armadillo advanced constructors -- no input copies. The R wrappers
// guarantee REALSXP storage.

// [[Rcpp::export]]
arma::mat FastSpatialMeat(Rcpp::NumericVector lat, Rcpp::NumericVector lon,
                          Rcpp::NumericVector time, Rcpp::NumericMatrix scores,
                          double cutoff,
                          std::string kernel = "bartlett",
                          std::string dist_fn = "haversine",
                          bool balanced_pnl = false,
                          int ncores = 1,
                          std::string neighbor = "grid",
                          std::string csr_weight = "double") {
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(lat.size()) != n ||
      static_cast<std::size_t>(lon.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("lat, lon, time, and scores have incompatible lengths.");
  }

  ncores = std::max(1, ncores);
  const int kernel_id = parse_kernel_id(kernel);
  const int dist_id = parse_dist_id(dist_fn);

  bool use_grid;
  if (neighbor == "grid") use_grid = true;
  else if (neighbor == "band") use_grid = false;
  else Rcpp::stop("Unknown neighbor: %s (use \"grid\" or \"band\")", neighbor.c_str());
  // cutoff < 0 is the no-pairs sentinel; the band path's break logic handles
  // it (only the 0.5*S_i diagonal survives), so route it there.
  if (cutoff < 0.0) use_grid = false;

  bool wfloat;
  if (csr_weight == "double") wfloat = false;
  else if (csr_weight == "float") wfloat = true;
  else Rcpp::stop("Unknown csr_weight: %s (use \"double\" or \"float\")",
                  csr_weight.c_str());

  if (n == 0) {
    return arma::mat(k, k, arma::fill::zeros);
  }

  const arma::vec lat_v(lat.begin(), n, false, true);
  const arma::vec lon_v(lon.begin(), n, false, true);
  const arma::vec time_v(time.begin(), n, false, true);
  const arma::mat S_col(scores.begin(), n, k, false, true);

  const TimeBlocks blocks = make_time_blocks(time_v);
  const bool use_balanced = balanced_pnl && blocks.start.size() > 1 && same_block_size(blocks);

  if (balanced_pnl && blocks.start.size() > 1 && !same_block_size(blocks)) {
    Rcpp::warning("balanced_pnl = TRUE but time blocks have unequal sizes; using general CSR path.");
  }

  return dispatch_spatial(dist_id, kernel_id, lat_v, lon_v, time_v, S_col,
                          cutoff, ncores, use_balanced, use_grid, wfloat);
}

// [[Rcpp::export]]
arma::mat FastSerialHacPanel(Rcpp::NumericVector unit, Rcpp::NumericVector time,
                             double cutoff, Rcpp::NumericMatrix scores,
                             int ncores = 1) {
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(unit.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("unit, time, and scores have incompatible lengths.");
  }
  if (n == 0) {
    return arma::mat(k, k, arma::fill::zeros);
  }

  ncores = std::max(1, ncores);
  const arma::vec unit_v(unit.begin(), n, false, true);
  const arma::vec time_v(time.begin(), n, false, true);
  const arma::mat S_col(scores.begin(), n, k, false, true);

  const UnitBlocks blocks = make_unit_blocks(unit_v);
  const RowMajorScores S(S_col, ncores);
  return serial_hac_panel(time_v, cutoff, S, blocks, ncores);
}


// [[Rcpp::export]]
arma::mat FastGridMeat(Rcpp::IntegerVector ring, Rcpp::IntegerVector col,
                       Rcpp::NumericVector time, Rcpp::NumericMatrix scores,
                       double lat0, double dlat, double dlon,
                       int n_ring, int n_col,
                       double cutoff, std::string dist_fn = "spherical",
                       int ncores = 1) {
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(ring.size()) != n ||
      static_cast<std::size_t>(col.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("ring, col, time, and scores have incompatible lengths.");
  }
  if (cutoff < 0.0) Rcpp::stop("FastGridMeat requires a non-negative cutoff.");
  if (n_ring < 1 || n_col < 1) Rcpp::stop("n_ring and n_col must be >= 1.");
  ncores = std::max(1, ncores);

  arma::mat meat_total(k, k, arma::fill::zeros);
  if (n == 0) return meat_total;

  const int dist_id = parse_dist_id(dist_fn);
  const ScreenParams screen = make_screen_params(cutoff, dist_id);
  // Uniform-kernel accept threshold on the unit-vector dot product. All
  // three distances are monotone in the chord, so each maps to a dot
  // threshold: haversine/spherical use cos(angular); chord uses
  // 1 - chord_sq/2.
  const double coscut = (dist_id == DIST_CHORD)
      ? 1.0 - screen.chord_cutoff_sq / 2.0
      : screen.cos_cutoff;

  const std::size_t R = static_cast<std::size_t>(n_ring);
  const std::size_t C = static_cast<std::size_t>(n_col);
  GridGeom gg;
  gg.n_ring = R;
  gg.n_col = C;
  gg.dlam_rad = dlon * DE2RA;
  gg.sphi.resize(R);
  gg.cphi.resize(R);
  for (std::size_t r = 0; r < R; ++r) {
    const double la = (lat0 + static_cast<double>(r) * dlat) * DE2RA;
    gg.sphi[r] = std::sin(la);
    gg.cphi[r] = std::cos(la);
  }
  const double dlat_rad = std::fabs(dlat) * DE2RA;
  const std::size_t rho =
      std::min(R, static_cast<std::size_t>(screen.angular_cutoff_rad /
                                           std::max(dlat_rad, 1e-300)) + 2);

  for (std::size_t i = 0; i < n; ++i) {
    if (ring[i] < 0 || ring[i] >= n_ring || col[i] < 0 || col[i] >= n_col) {
      Rcpp::stop("ring/col indices out of range.");
    }
  }

  const arma::vec time_v(time.begin(), n, false, true);
  const TimeBlocks blocks = make_time_blocks(time_v);
  const std::size_t rowlen = (C + 1) * k;
  std::vector<double> P(R * rowlen);
  std::vector<std::size_t> ring_count(R);

  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    std::fill(P.begin(), P.end(), 0.0);
    std::fill(ring_count.begin(), ring_count.end(), 0);

    // Scatter block scores into the (c+1) slots, then prefix-sum per ring.
    for (std::size_t i = blocks.start[b]; i < blocks.end[b]; ++i) {
      const std::size_t r = static_cast<std::size_t>(ring[i]);
      const std::size_t c = static_cast<std::size_t>(col[i]);
      double* dst = P.data() + r * rowlen + (c + 1) * k;
      for (std::size_t kk = 0; kk < k; ++kk) dst[kk] += scores(i, kk);
      ++ring_count[r];
    }
    parallel_range(R, ncores, [&](std::size_t lo, std::size_t hi) {
      for (std::size_t r = lo; r < hi; ++r) {
        if (ring_count[r] == 0) continue;
        double* Pr = P.data() + r * rowlen;
        for (std::size_t c = 0; c < C; ++c) {
          const double* prev = Pr + c * k;
          double* cur = Pr + (c + 1) * k;
          for (std::size_t kk = 0; kk < k; ++kk) cur[kk] += prev[kk];
        }
      }
    });

    // Deterministic reduce over target rings.
    auto body = [&](std::size_t lo, std::size_t hi, arma::mat& meat) {
      std::vector<double> cf(C * k);
      for (std::size_t r1 = lo; r1 < hi; ++r1) {
        if (ring_count[r1] == 0) continue;
        std::fill(cf.begin(), cf.end(), 0.0);
        bool any = false;
        const std::size_t r2_lo = r1 >= rho ? r1 - rho : 0;
        const std::size_t r2_hi = std::min(R - 1, r1 + rho);
        for (std::size_t r2 = r2_lo; r2 <= r2_hi; ++r2) {
          if (ring_count[r2] == 0) continue;
          const long ds = grid_delta_star(gg, r1, r2, coscut);
          if (ds < 0) continue;
          any = true;
          const double* Pr2 = P.data() + r2 * rowlen;
          for (std::size_t c = 0; c < C; ++c) {
            const std::size_t hi_i =
                std::min(C - 1, c + static_cast<std::size_t>(ds)) + 1;
            const std::size_t lo_i =
                c >= static_cast<std::size_t>(ds)
                    ? c - static_cast<std::size_t>(ds) : 0;
            const double* ph = Pr2 + hi_i * k;
            const double* pl = Pr2 + lo_i * k;
            double* cfc = cf.data() + c * k;
            for (std::size_t kk = 0; kk < k; ++kk) cfc[kk] += ph[kk] - pl[kk];
          }
        }
        if (!any) continue;
        const double* Pr1 = P.data() + r1 * rowlen;
        for (std::size_t c = 0; c < C; ++c) {
          const double* pc = Pr1 + c * k;
          const double* pn = Pr1 + (c + 1) * k;
          const double* cfc = cf.data() + c * k;
          for (std::size_t k1 = 0; k1 < k; ++k1) {
            const double s1 = pn[k1] - pc[k1];
            if (s1 == 0.0) continue;
            for (std::size_t k2 = 0; k2 < k; ++k2) {
              meat(k1, k2) += s1 * cfc[k2];
            }
          }
        }
      }
    };
    meat_total += reduce_deterministic(R, k, 4, ncores, body);
  }

  return meat_total;
}
