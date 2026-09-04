// conley_core.h -- the fastconley spatial / serial HAC engine, front-end
// agnostic. Included by the Rcpp front-end (src/cpp-functions.cpp) and by
// the Stata plugin (stata/plugin/). Depends only on the C++11 standard
// library and on Armadillo used header-only (arma::mat accumulators and
// arma::fft; no BLAS/LAPACK calls), so a standalone build needs just
// -std=c++11 -pthread and the Armadillo headers (with
// -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK).
//
// Errors are reported by throwing conley::Error; each front-end catches
// and translates it (Rcpp's BEGIN_RCPP/END_RCPP turns it into an R error,
// the Stata plugin into SF_error + a return code). Nothing here touches
// R or Stata, and nothing here copies the inputs: callers alias their
// own memory through the Armadillo advanced constructors.
#ifndef FASTCONLEY_CONLEY_CORE_H
#define FASTCONLEY_CONLEY_CORE_H

// Engine version, reported by the Stata plugin's "check" subcommand and
// compared against the ado's expectation. Bump with the package version
// whenever the numerics or the entry-point signatures change.
#define CONLEY_CORE_VERSION "0.11.0"

#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD 1
#endif
// The R front-end includes RcppArmadillo.h first (that header pulls in
// armadillo with R's allocator / RNG hooks); a standalone build gets the
// plain library. ARMA_INCLUDES is armadillo's own include guard.
#ifndef ARMA_INCLUDES
#include <armadillo>
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace conley {

struct Error : public std::runtime_error {
  explicit Error(const std::string& msg) : std::runtime_error(msg) {}
};

[[noreturn]] inline void fail(const std::string& msg) { throw Error(msg); }


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

// ------------------------------------------------------------------
// Threading. A tiny std::thread pool per call: [0, n) is cut into
// fixed blocks handed out through an atomic counter (dynamic load
// balance for irregular per-row work). Every use site is either a pure
// gather/scatter with disjoint outputs or writes chunk-indexed partials
// that are reduced serially in chunk order, so results never depend on
// the thread count or on scheduling. Exceptions thrown by workers are
// rethrown on the calling thread once every worker has joined.
// ------------------------------------------------------------------

template <typename F>
void parallel_blocks(std::size_t n, int ncores, std::size_t grain, F&& fn) {
  if (n == 0) return;
  if (grain == 0) grain = 1;
  const std::size_t nblocks = (n + grain - 1) / grain;
  const std::size_t nthreads = std::min<std::size_t>(
      static_cast<std::size_t>(std::max(1, ncores)), nblocks);
  if (nthreads <= 1) {
    fn(static_cast<std::size_t>(0), n);
    return;
  }
  std::atomic<std::size_t> next(0);
  std::atomic<bool> failed(false);
  std::exception_ptr err;
  std::mutex err_mutex;
  auto work = [&]() {
    try {
      for (;;) {
        if (failed.load(std::memory_order_relaxed)) break;
        const std::size_t b = next.fetch_add(1, std::memory_order_relaxed);
        if (b >= nblocks) break;
        const std::size_t lo = b * grain;
        const std::size_t hi = std::min(n, lo + grain);
        fn(lo, hi);
      }
    } catch (...) {
      std::lock_guard<std::mutex> guard(err_mutex);
      if (!err) err = std::current_exception();
      failed.store(true);
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nthreads - 1);
  for (std::size_t t = 1; t < nthreads; ++t) pool.emplace_back(work);
  work();
  for (std::size_t t = 0; t < pool.size(); ++t) pool[t].join();
  if (err) std::rethrow_exception(err);
}

// Block size for per-row scans: ~16 blocks per thread, never below 1024
// rows, so the atomic hand-off is negligible and irregular rows balance.
inline std::size_t row_grain(std::size_t n, int ncores) {
  const std::size_t t = static_cast<std::size_t>(std::max(1, ncores));
  return std::max<std::size_t>(1024, (n + 16 * t - 1) / (16 * t));
}

// Deterministic parallel-for over [0, n): body(lo, hi) must write only
// disjoint per-index outputs (pure gather/scatter), so the result is
// independent of the work partition. Small inputs stay serial.
template <typename Body>
void parallel_range(std::size_t n, int ncores, const Body& body) {
  if (ncores > 1 && n > 8192) {
    parallel_blocks(n, ncores, row_grain(n, ncores), body);
  } else {
    body(0, n);
  }
}

// Like parallel_range, but for coarse-grained items (e.g. one FFT batch
// per index) where even a handful of items is worth parallelising.
template <typename Body>
void parallel_range_coarse(std::size_t n, int ncores, const Body& body) {
  if (ncores > 1 && n > 1) {
    parallel_blocks(n, ncores, 1, body);
  } else {
    body(0, n);
  }
}

// The comparator must define a strict total order (callers tiebreak on
// index), so the sorted output is unique.
template <typename It, typename Cmp>
void sort_maybe_parallel(It first, It last, Cmp cmp, int ncores) {
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

inline int parse_dist_id(const std::string& dist_fn) {
  if (dist_fn == "haversine") return DIST_HAVERSINE;
  if (dist_fn == "spherical") return DIST_SPHERICAL;
  if (dist_fn == "chord") return DIST_CHORD;
  fail("Unknown dist_fn: " + dist_fn);
}

inline int parse_kernel_id(const std::string& kernel) {
  if (kernel == "bartlett") return KERNEL_BARTLETT;
  if (kernel == "uniform") return KERNEL_UNIFORM;
  fail("Unknown kernel: " + kernel);
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

inline ScreenParams make_screen_params(double cutoff, int dist_id) {
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

inline CoordCache make_coord_cache(const arma::vec& lat, const arma::vec& lon,
                            int dist_id, int ncores) {
  const std::size_t n = lat.n_elem;
  // Validate serially first so the error surfaces before any threads spawn.
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(lat[i]) || !std::isfinite(lon[i])) {
      fail("lat/lon contain non-finite values.");
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
// therefore independent of the thread count and of scheduling:
// ncores = 1 and ncores = N give bit-identical results.
// ------------------------------------------------------------------

// body(lo, hi, meat&) must accumulate rows [lo, hi) into meat and itself be
// deterministic in row order.
template <typename Body>
arma::mat reduce_deterministic(std::size_t n, std::size_t k, std::size_t chunk,
                               int ncores, const Body& body) {
  arma::mat meat(k, k, arma::fill::zeros);
  if (n == 0) return meat;
  const std::size_t nchunks = (n + chunk - 1) / chunk;
  std::vector<arma::mat> partials(nchunks, arma::mat(k, k, arma::fill::zeros));
  parallel_blocks(nchunks, ncores, 1, [&](std::size_t cb, std::size_t ce) {
    for (std::size_t ci = cb; ci < ce; ++ci) {
      const std::size_t lo = ci * chunk;
      const std::size_t hi = std::min(n, lo + chunk);
      body(lo, hi, partials[ci]);
    }
  });
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
inline CoordCache permute_coord_cache(const CoordCache& src,
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

inline TimeBlocks make_time_blocks(const arma::vec& time) {
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

inline bool same_block_size(const TimeBlocks& blocks) {
  if (blocks.start.empty()) return true;
  const std::size_t n0 = blocks.end[0] - blocks.start[0];
  for (std::size_t b = 1; b < blocks.start.size(); ++b) {
    if (blocks.end[b] - blocks.start[b] != n0) return false;
  }
  return true;
}

inline void sort_block_indices(const std::vector<double>& lat_rad,
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
struct CsrCountWorkerT {
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
struct CsrFillWorkerT {
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
    fail("Balanced CSR path supports at most 2^32 - 1 units per period.");
  }

  std::vector<std::size_t> counts(n, 0);
  const ScreenParams screen = make_screen_params(cutoff, D);

  CsrCountWorkerT<D, K> count_worker(sorted_idx, row_end, c,
                                     counts, cutoff, screen);
  parallel_blocks(n, ncores, row_grain(n, ncores), count_worker);

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
  parallel_blocks(n, ncores, row_grain(n, ncores), fill_worker);

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

inline GridSpec make_grid_spec(const ScreenParams& screen) {
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
inline void append_block_grid(const CoordCache& c, std::size_t bs, std::size_t be,
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
struct GridCsrCountWorkerT {
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
struct GridCsrFillWorkerT {
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
    fail("Balanced CSR path supports at most 2^32 - 1 units per period.");
  }

  std::vector<std::size_t> counts(n, 0);
  GridCsrCountWorkerT<D, K> count_worker(grid, c, counts, cutoff, screen);
  parallel_blocks(n, ncores, row_grain(n, ncores), count_worker);

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
  parallel_blocks(n, ncores, row_grain(n, ncores), fill_worker);

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
    fail("The grid neighbor path supports at most 2^32 - 1 rows; "
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

inline arma::mat dispatch_spatial(int dist_id, int kernel_id,
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
  fail("Unsupported (dist_id, kernel_id) combination.");
}

struct UnitBlocks {
  std::vector<std::size_t> start;
  std::vector<std::size_t> end;
};

inline UnitBlocks make_unit_blocks(const arma::vec& unit) {
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
inline arma::mat serial_hac_panel(const arma::vec& times, double cutoff,
                           const RowMajorScores& S, const UnitBlocks& blocks,
                           int ncores) {
  const std::size_t k = S.k;
  if (cutoff < 0.0) return arma::mat(k, k, arma::fill::zeros);

  for (std::size_t bi = 0; bi < blocks.start.size(); ++bi) {
    for (std::size_t i = blocks.start[bi] + 1; i < blocks.end[bi]; ++i) {
      if (times[i] < times[i - 1]) {
        fail("FastSerialHacPanel requires rows sorted by time within "
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
// On a regular lat/lon lattice the accept set between two latitude rings
// is a longitude-index interval. For the UNIFORM kernel the inner sum
// over a ring pair is therefore a sliding-window sum over per-ring prefix
// sums: O(n_ring * window * n_col * k) total, independent of the pair
// count. For the BARTLETT kernel the weight varies with the lon offset,
// so the inner sum is a true 1D convolution per ring pair, computed via
// FFT (arma::fft) with per-ring score spectra cached in a sliding
// latitude band. Both are exact: the dot-product accept threshold is the
// same constant the pairwise engine uses, and the bartlett weights use
// the same per-distance arithmetic as pair_weight, so answers agree to
// FP summation order (plus the inherent acos conditioning for
// spherical x bartlett).
//
// If the lattice wraps the full longitude circle (n_col_full > 0 and the
// accept window reaches across the dateline gap), windows become
// circular: modular prefix-sum arcs for the uniform kernel, and circular
// convolution with period n_col_full for bartlett.
//
// Uniform layout: one dense prefix tensor per time block, ring-major:
// P[r * (C+1) * k + c * k + kk] = sum of scores over columns < c of ring
// r (k-wide rows). Raw cell scores are recovered as adjacent differences,
// so no separate dense score copy is held. Memory: (C+1) * R * k doubles.
// Bartlett layout: dense raw scores DS (R * C * k doubles) plus the
// banded spectrum cache (~(band + 2*rho) * npad * k complex doubles).
// ------------------------------------------------------------------

struct GridGeom {
  std::size_t n_ring;
  std::size_t n_col;
  double dlam_rad;
  double dlat_rad;           // |lat step| (bartlett haversine weights)
  std::vector<double> sphi;  // sin(lat) per ring
  std::vector<double> cphi;  // cos(lat) per ring
  // Per-offset trig tables, delta = 0..cap (shared by every ring pair):
  std::vector<double> cos_dl;  // cos(delta * dlam)
  std::vector<double> s2h_dl;  // sin^2(delta * dlam / 2)
};

// Accept window between rings r1 and r2 against the dot threshold,
// matching the pairwise accept "dot >= coscut" with a boundary fix-up in
// the same arithmetic. ds: largest accepted lon-index offset, capped at
// `cap` (-1 = none). full: every lon offset on the circle is accepted
// (polar / antipodal-in-lon cases; the window is the whole ring).
// dstar_rad: the acos accept boundary in lon radians for interval pairs
// (0 when full/none) -- drives the dateline-wrap feasibility check.
struct RingPairWin {
  long ds;
  bool full;
  double dstar_rad;
};

inline RingPairWin grid_ring_window(const GridGeom& gg, std::size_t r1,
                                    std::size_t r2, double coscut, long cap) {
  const double base = gg.sphi[r1] * gg.sphi[r2];
  const double cc12 = gg.cphi[r1] * gg.cphi[r2];
  RingPairWin out;
  out.dstar_rad = 0.0;
  if (cc12 < 1e-300) {
    out.full = base >= coscut;
    out.ds = out.full ? cap : -1;
    return out;
  }
  const double rhs = (coscut - base) / cc12;
  if (rhs > 1.0) {
    out.full = false;
    out.ds = -1;
    return out;
  }
  if (rhs <= -1.0) {
    out.full = true;
    out.ds = cap;
    return out;
  }
  out.full = false;
  out.dstar_rad = std::acos(rhs);
  long ds = static_cast<long>(out.dstar_rad / gg.dlam_rad) + 1;
  if (ds > cap) ds = cap;
  while (ds >= 0 && base + cc12 * std::cos(ds * gg.dlam_rad) < coscut) --ds;
  while (ds + 1 <= cap &&
         base + cc12 * std::cos((ds + 1) * gg.dlam_rad) >= coscut) ++ds;
  out.ds = ds;
  return out;
}

// Bartlett weights between rings r1 and r2 at lon-index offsets 0..ds,
// using the same per-distance arithmetic as the pair_weight
// specializations. Clamped at 0 so FP overshoot at the window boundary
// cannot produce a negative weight (the pairwise engine would have
// rejected such a pair).
inline void grid_bartlett_weights(const GridGeom& gg, std::size_t r1,
                                  std::size_t r2, int dist_id, double cutoff,
                                  long ds, double* w) {
  const double base = gg.sphi[r1] * gg.sphi[r2];
  const double cc12 = gg.cphi[r1] * gg.cphi[r2];
  // The same-cell distance is exactly 0, but acos/sqrt of the FP dot
  // (sin^2 + cos^2 = 1 +/- eps) would inflate it to ~R*sqrt(2*eps). The
  // pairwise engine weights self-pairs exactly 1 (its 0.5*S_i diagonal),
  // so fix the d = 0 weight after the loop below.
  const bool self_ring = (r1 == r2);
  if (dist_id == DIST_HAVERSINE) {
    const double dphi =
        (static_cast<double>(r2) - static_cast<double>(r1)) * gg.dlat_rad;
    const double s2_dphi = sq(std::sin(dphi / 2.0));
    for (long d = 0; d <= ds; ++d) {
      double a = s2_dphi + cc12 * gg.s2h_dl[d];
      if (a < 0.0) a = 0.0;
      if (a > 1.0) a = 1.0;
      const double dist =
          AVG_ERAD * 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
      const double wt = 1.0 - dist / cutoff;
      w[d] = wt > 0.0 ? wt : 0.0;
    }
  } else if (dist_id == DIST_SPHERICAL) {
    for (long d = 0; d <= ds; ++d) {
      const double dot = base + cc12 * gg.cos_dl[d];
      const double dist = AVG_ERAD * safe_acos(dot);
      const double wt = 1.0 - dist / cutoff;
      w[d] = wt > 0.0 ? wt : 0.0;
    }
  } else {  // DIST_CHORD: chord^2 = 2 - 2*dot in exact arithmetic.
    for (long d = 0; d <= ds; ++d) {
      const double dot = base + cc12 * gg.cos_dl[d];
      double c2 = 2.0 - 2.0 * dot;
      if (c2 < 0.0) c2 = 0.0;
      const double dist = AVG_ERAD * std::sqrt(c2);
      const double wt = 1.0 - dist / cutoff;
      w[d] = wt > 0.0 ? wt : 0.0;
    }
  }
  if (self_ring) w[0] = 1.0;
}

inline std::size_t next_pow2(std::size_t n) {
  std::size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}


// ------------------------------------------------------------------
// Entry points shared by the front-ends. `scores` is the column-major
// score matrix (scores = e * X, possibly pre-aggregated by the caller);
// lat / lon / time / unit have one entry per score row. Callers validate
// lengths; the strings are parsed here so both front-ends accept exactly
// the same vocabulary.
// ------------------------------------------------------------------

// Spatial meat. Rows must be sorted so each time block is contiguous
// (and unit-stable within time when balanced_pnl is true). When
// balanced_pnl is requested but the blocks have unequal sizes, the
// general path is used and *unbalanced_fallback is set so the front-end
// can warn.
inline arma::mat spatial_meat(const arma::vec& lat, const arma::vec& lon,
                              const arma::vec& time, const arma::mat& scores,
                              double cutoff, const std::string& kernel,
                              const std::string& dist_fn, bool balanced_pnl,
                              int ncores, const std::string& neighbor,
                              const std::string& csr_weight,
                              bool* unbalanced_fallback = nullptr) {
  const std::size_t n = scores.n_rows;
  const std::size_t k = scores.n_cols;
  ncores = std::max(1, ncores);
  const int kernel_id = parse_kernel_id(kernel);
  const int dist_id = parse_dist_id(dist_fn);

  bool use_grid;
  if (neighbor == "grid") use_grid = true;
  else if (neighbor == "band") use_grid = false;
  else fail("Unknown neighbor: " + neighbor + " (use \"grid\" or \"band\")");
  // cutoff < 0 is the no-pairs sentinel; the band path's break logic handles
  // it (only the 0.5*S_i diagonal survives), so route it there.
  if (cutoff < 0.0) use_grid = false;

  bool wfloat;
  if (csr_weight == "double") wfloat = false;
  else if (csr_weight == "float") wfloat = true;
  else fail("Unknown csr_weight: " + csr_weight + " (use \"double\" or \"float\")");

  if (unbalanced_fallback) *unbalanced_fallback = false;
  if (n == 0) {
    return arma::mat(k, k, arma::fill::zeros);
  }

  const TimeBlocks blocks = make_time_blocks(time);
  const bool multi = blocks.start.size() > 1;
  const bool use_balanced = balanced_pnl && multi && same_block_size(blocks);
  if (unbalanced_fallback) {
    *unbalanced_fallback = balanced_pnl && multi && !same_block_size(blocks);
  }
  return dispatch_spatial(dist_id, kernel_id, lat, lon, time, scores,
                          cutoff, ncores, use_balanced, use_grid, wfloat);
}

// Serial (within-unit, across-time) HAC meat. Rows must be sorted by
// (unit, time) so unit blocks are contiguous.
inline arma::mat serial_hac_meat(const arma::vec& unit, const arma::vec& time,
                                 double cutoff, const arma::mat& scores,
                                 int ncores) {
  const std::size_t k = scores.n_cols;
  if (scores.n_rows == 0) {
    return arma::mat(k, k, arma::fill::zeros);
  }
  ncores = std::max(1, ncores);
  const UnitBlocks blocks = make_unit_blocks(unit);
  const RowMajorScores S(scores, ncores);
  return serial_hac_panel(time, cutoff, S, blocks, ncores);
}

// Exact meat on a regular lat/lon lattice (ring = latitude index, col =
// longitude index, both 0-based, one entry per score row). See the
// grid-native section above for the algorithm.
inline arma::mat grid_meat(const int* ring, const int* col,
                           const arma::vec& time, const arma::mat& scores,
                           double lat0, double dlat, double dlon,
                           int n_ring, int n_col, int n_col_full,
                           double cutoff, const std::string& dist_fn,
                           const std::string& kernel, int ncores) {
  const std::size_t n = scores.n_rows;
  const std::size_t k = scores.n_cols;
  if (cutoff < 0.0) fail("FastGridMeat requires a non-negative cutoff.");
  if (n_ring < 1 || n_col < 1) fail("n_ring and n_col must be >= 1.");
  if (n_col_full < 0 || (n_col_full > 0 && n_col_full < n_col)) {
    fail("n_col_full must be 0 (no wrap) or >= n_col.");
  }
  ncores = std::max(1, ncores);

  const int kernel_id = parse_kernel_id(kernel);
  if (kernel_id == KERNEL_BARTLETT && cutoff <= 0.0) {
    fail("FastGridMeat requires cutoff > 0 for kernel = 'bartlett'.");
  }

  arma::mat meat_total(k, k, arma::fill::zeros);
  if (n == 0) return meat_total;

  const int dist_id = parse_dist_id(dist_fn);
  const ScreenParams screen = make_screen_params(cutoff, dist_id);
  // Accept threshold on the unit-vector dot product. All three distances
  // are monotone in the chord, so each maps to a dot threshold:
  // haversine/spherical use cos(angular); chord uses 1 - chord_sq/2.
  const double coscut = (dist_id == DIST_CHORD)
      ? 1.0 - screen.chord_cutoff_sq / 2.0
      : screen.cos_cutoff;

  const std::size_t R = static_cast<std::size_t>(n_ring);
  const std::size_t C = static_cast<std::size_t>(n_col);
  GridGeom gg;
  gg.n_ring = R;
  gg.n_col = C;
  gg.dlam_rad = dlon * DE2RA;
  gg.dlat_rad = std::fabs(dlat) * DE2RA;
  gg.sphi.resize(R);
  gg.cphi.resize(R);
  for (std::size_t r = 0; r < R; ++r) {
    const double la = (lat0 + static_cast<double>(r) * dlat) * DE2RA;
    gg.sphi[r] = std::sin(la);
    gg.cphi[r] = std::cos(la);
  }
  const std::size_t rho =
      std::min(R, static_cast<std::size_t>(screen.angular_cutoff_rad /
                                           std::max(gg.dlat_rad, 1e-300)) + 2);

  for (std::size_t i = 0; i < n; ++i) {
    if (ring[i] < 0 || ring[i] >= n_ring || col[i] < 0 || col[i] >= n_col) {
      fail("ring/col indices out of range.");
    }
  }

  // ---- Dateline-wrap feasibility (geometry only, block-independent). ----
  // A pair of columns can be physically within the cutoff "the short way
  // around" yet farther apart than the linear window iff the lattice span
  // plus the widest interval accept window reaches the full circle. Full
  // (whole-ring) windows never miss anything. dstar_max / ds_max are
  // order-independent maxima, so the parallel pre-pass is deterministic.
  const long cap_lin = static_cast<long>(C) - 1;
  double dstar_max = 0.0;
  long ds_max = 0;
  {
    std::vector<double> permax(R, 0.0);
    std::vector<long> perds(R, 0);
    parallel_range(R, ncores, [&](std::size_t lo, std::size_t hi) {
      for (std::size_t r1 = lo; r1 < hi; ++r1) {
        const std::size_t r2_hi = std::min(R - 1, r1 + rho);
        for (std::size_t r2 = r1; r2 <= r2_hi; ++r2) {
          const RingPairWin win = grid_ring_window(gg, r1, r2, coscut, cap_lin);
          if (win.ds < 0) continue;
          if (win.ds > perds[r1]) perds[r1] = win.ds;
          if (!win.full && win.dstar_rad > permax[r1]) {
            permax[r1] = win.dstar_rad;
          }
        }
      }
    });
    for (std::size_t r = 0; r < R; ++r) {
      if (permax[r] > dstar_max) dstar_max = permax[r];
      if (perds[r] > ds_max) ds_max = perds[r];
    }
  }
  const double span_rad = static_cast<double>(C - 1) * gg.dlam_rad;
  const bool wrap_needed = span_rad + dstar_max >= TWO_PI - 1e-9;
  if (wrap_needed && n_col_full == 0) {
    fail("FastGridMeat: the accept window crosses the dateline wrap "
               "but the lattice does not tile the full longitude circle "
               "(360 is not an integer multiple of dlon). "
               "Use method = 'pairwise'.");
  }
  const std::size_t Cf = static_cast<std::size_t>(n_col_full);
  const bool wrap = wrap_needed && n_col_full > 0;
  // Wrap implies span >= pi, so cap <= cap_lin and offsets never exceed
  // the half circle (each unordered column pair is counted exactly once).
  const long cap = wrap ? static_cast<long>(Cf / 2) : cap_lin;

  gg.cos_dl.resize(static_cast<std::size_t>(cap) + 1);
  gg.s2h_dl.resize(static_cast<std::size_t>(cap) + 1);
  for (long d = 0; d <= cap; ++d) {
    gg.cos_dl[d] = std::cos(static_cast<double>(d) * gg.dlam_rad);
    gg.s2h_dl[d] = sq(std::sin(static_cast<double>(d) * gg.dlam_rad / 2.0));
  }

  const TimeBlocks blocks = make_time_blocks(time);
  std::vector<std::size_t> ring_count(R);

  if (kernel_id == KERNEL_UNIFORM) {
    // ---- Uniform kernel: sliding-window prefix sums (boxcar). ----
    const std::size_t rowlen = (C + 1) * k;
    std::vector<double> P(R * rowlen);

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
            const long ds = grid_ring_window(gg, r1, r2, coscut, cap).ds;
            if (ds < 0) continue;
            any = true;
            const double* Pr2 = P.data() + r2 * rowlen;
            if (!wrap) {
              for (std::size_t c = 0; c < C; ++c) {
                const std::size_t hi_i =
                    std::min(C - 1, c + static_cast<std::size_t>(ds)) + 1;
                const std::size_t lo_i =
                    c >= static_cast<std::size_t>(ds)
                        ? c - static_cast<std::size_t>(ds) : 0;
                const double* ph = Pr2 + hi_i * k;
                const double* pl = Pr2 + lo_i * k;
                double* cfc = cf.data() + c * k;
                for (std::size_t kk = 0; kk < k; ++kk) {
                  cfc[kk] += ph[kk] - pl[kk];
                }
              }
            } else {
              // Circular window [c - ds, c + ds] on the Cf-circle,
              // intersected with the occupied columns [0, C). Up to two
              // linear arcs; arcs cannot overlap because 2*ds + 1 < Cf
              // whenever this branch splits.
              const long dsl = ds;
              const long Cl = static_cast<long>(C);
              const long Cfl = static_cast<long>(Cf);
              const bool whole = 2 * dsl + 1 >= Cfl;
              for (std::size_t c = 0; c < C; ++c) {
                double* cfc = cf.data() + c * k;
                const auto add_arc = [&](long a, long barc) {
                  if (a > barc || a >= Cl) return;
                  if (barc >= Cl) barc = Cl - 1;
                  const double* ph = Pr2 + static_cast<std::size_t>(barc + 1) * k;
                  const double* pl = Pr2 + static_cast<std::size_t>(a) * k;
                  for (std::size_t kk = 0; kk < k; ++kk) {
                    cfc[kk] += ph[kk] - pl[kk];
                  }
                };
                if (whole) {
                  add_arc(0, Cl - 1);
                  continue;
                }
                const long a = static_cast<long>(c) - dsl;
                const long bb = static_cast<long>(c) + dsl;
                if (a < 0) {
                  add_arc(0, bb);
                  add_arc(a + Cfl, Cfl - 1);
                } else if (bb >= Cfl) {
                  add_arc(a, Cfl - 1);
                  add_arc(0, bb - Cfl);
                } else {
                  add_arc(a, bb);
                }
              }
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

  // ---- Bartlett kernel: per-ring-pair 1D convolutions via FFT. ----
  // Linear mode zero-pads to the next power of two past C + ds_max so the
  // circular convolution cannot wrap; wrap mode uses period n_col_full so
  // it wraps *correctly*. Score spectra are cached for a sliding latitude
  // band of target rings plus a rho-halo; BAND_H is a fixed constant and
  // the reduction is chunked, so results are bit-identical across ncores.
  //
  // The kernel is even in the lon offset, so T(r2,r1) = T(r1,r2)': only
  // the upper triangle r2 >= r1 is computed, with self-ring weights halved
  // (the pairwise engine's 0.5*S_i trick), and the half-meat is
  // symmetrized at the end. This halves the FFT work and the cache only
  // needs a forward halo.
  const std::size_t npad =
      wrap ? Cf : next_pow2(C + static_cast<std::size_t>(std::max(ds_max, 0L)));
  const std::size_t BAND_H = 128;

  std::vector<double> DS(R * C * k);  // raw dense scores, ring-major
  std::vector<arma::cx_mat> shat(R);  // per-ring score spectra (banded cache)

  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    std::fill(DS.begin(), DS.end(), 0.0);
    std::fill(ring_count.begin(), ring_count.end(), 0);
    for (std::size_t r = 0; r < R; ++r) shat[r].reset();

    for (std::size_t i = blocks.start[b]; i < blocks.end[b]; ++i) {
      const std::size_t r = static_cast<std::size_t>(ring[i]);
      const std::size_t c = static_cast<std::size_t>(col[i]);
      double* dst = DS.data() + (r * C + c) * k;
      for (std::size_t kk = 0; kk < k; ++kk) dst[kk] += scores(i, kk);
      ++ring_count[r];
    }

    for (std::size_t band_lo = 0; band_lo < R; band_lo += BAND_H) {
      const std::size_t band_hi = std::min(R, band_lo + BAND_H);
      const std::size_t need_lo = band_lo;
      const std::size_t need_hi = std::min(R, band_hi + rho);

      // Fill missing spectra for the band + halo (one FFT batch per ring;
      // disjoint writes, so the parallel fill is deterministic).
      std::vector<std::size_t> todo;
      for (std::size_t r = need_lo; r < need_hi; ++r) {
        if (ring_count[r] > 0 && shat[r].n_elem == 0) todo.push_back(r);
      }
      parallel_range_coarse(todo.size(), ncores,
                            [&](std::size_t lo, std::size_t hi) {
        for (std::size_t ti = lo; ti < hi; ++ti) {
          const std::size_t r = todo[ti];
          arma::mat A(npad, k, arma::fill::zeros);
          const double* src = DS.data() + r * C * k;
          for (std::size_t c = 0; c < C; ++c) {
            for (std::size_t kk = 0; kk < k; ++kk) A(c, kk) = src[c * k + kk];
          }
          shat[r] = arma::fft(A);
        }
      });

      // Deterministic reduce over the band's target rings (chunk = 1: each
      // ring is a heavy item -- ~2*rho weight FFTs).
      auto body = [&](std::size_t lo, std::size_t hi, arma::mat& meat) {
        std::vector<double> wbuf(static_cast<std::size_t>(cap) + 1);
        arma::vec wpad(npad);
        std::vector<double> rew(npad);
        arma::cx_mat acc(npad, k);
        std::vector<double> cfr(C * k);
        for (std::size_t ii = lo; ii < hi; ++ii) {
          const std::size_t r1 = band_lo + ii;
          if (ring_count[r1] == 0) continue;
          acc.zeros();
          bool any = false;
          const std::size_t r2_hi = std::min(R - 1, r1 + rho);
          for (std::size_t r2 = r1; r2 <= r2_hi; ++r2) {
            if (ring_count[r2] == 0) continue;
            const RingPairWin win = grid_ring_window(gg, r1, r2, coscut, cap);
            if (win.ds < 0) continue;
            const long ds = win.ds;
            grid_bartlett_weights(gg, r1, r2, dist_id, cutoff, ds,
                                  wbuf.data());
            if (r2 == r1) {
              for (long d = 0; d <= ds; ++d) wbuf[d] *= 0.5;
            }
            wpad.zeros();
            wpad[0] = wbuf[0];
            for (long d = 1; d <= ds; ++d) {
              wpad[static_cast<std::size_t>(d)] = wbuf[d];
              wpad[npad - static_cast<std::size_t>(d)] = wbuf[d];
            }
            // w is even in the offset, so its DFT is real in exact
            // arithmetic; dropping the FP-noise imaginary part halves the
            // accumulate cost.
            const arma::cx_vec what = arma::fft(wpad);
            for (std::size_t m = 0; m < npad; ++m) rew[m] = what[m].real();
            const arma::cx_mat& S2 = shat[r2];
            for (std::size_t kk = 0; kk < k; ++kk) {
              std::complex<double>* a = acc.colptr(kk);
              const std::complex<double>* s = S2.colptr(kk);
              for (std::size_t m = 0; m < npad; ++m) a[m] += rew[m] * s[m];
            }
            any = true;
          }
          if (!any) continue;
          const arma::cx_mat cfm = arma::ifft(acc);
          for (std::size_t kk = 0; kk < k; ++kk) {
            const std::complex<double>* colp = cfm.colptr(kk);
            for (std::size_t c = 0; c < C; ++c) {
              cfr[c * k + kk] = colp[c].real();
            }
          }
          const double* s1base = DS.data() + r1 * C * k;
          for (std::size_t c = 0; c < C; ++c) {
            const double* s1 = s1base + c * k;
            const double* cfc = cfr.data() + c * k;
            for (std::size_t k1 = 0; k1 < k; ++k1) {
              if (s1[k1] == 0.0) continue;
              for (std::size_t k2 = 0; k2 < k; ++k2) {
                meat(k1, k2) += s1[k1] * cfc[k2];
              }
            }
          }
        }
      };
      meat_total += reduce_deterministic(band_hi - band_lo, k, 1, ncores, body);

      // Evict spectra that have left the sliding window (later bands only
      // look forward from their own start).
      for (std::size_t r = need_lo; r < band_hi; ++r) shat[r].reset();
    }
  }

  return meat_total + meat_total.t();
}

} // namespace conley

#endif // FASTCONLEY_CONLEY_CORE_H
