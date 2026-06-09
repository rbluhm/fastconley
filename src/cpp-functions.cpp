#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD 1
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]

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

struct CoordCache {
  std::vector<double> lat_rad;
  std::vector<double> lon_rad;
  std::vector<double> cos_lat;
  std::vector<double> sin_lat;
  // 3D Cartesian coordinates on a sphere of radius AVG_ERAD; only filled when
  // dist_id == DIST_CHORD, otherwise empty.
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

// HAVERSINE × UNIFORM: longitude screen, then haversine 'a'-threshold (no atan2).
template<>
inline double pair_weight<DIST_HAVERSINE, KERNEL_UNIFORM>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
  (void)cutoff;
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

// HAVERSINE × BARTLETT: same screen, plus the real distance via atan2.
template<>
inline double pair_weight<DIST_HAVERSINE, KERNEL_BARTLETT>(
    const CoordCache& c, std::size_t i, std::size_t j,
    double cutoff, const ScreenParams& screen) {
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
                            int dist_id) {
  const std::size_t n = lat.n_elem;
  CoordCache c;
  c.lat_rad.resize(n);
  c.lon_rad.resize(n);
  c.cos_lat.resize(n);
  c.sin_lat.resize(n);
  // SPHERICAL and CHORD consume 3D *unit* vectors per pair (dot-product and
  // squared-distance thresholds), and the cell-grid neighbor search bins
  // every distance function by its unit vector — so xyz is always built.
  (void)dist_id;
  const bool need_xyz = true;
  if (need_xyz) {
    c.x3.resize(n);
    c.y3.resize(n);
    c.z3.resize(n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(lat[i]) || !std::isfinite(lon[i])) {
      Rcpp::stop("lat/lon contain non-finite values.");
    }
    const double la = lat[i] * DE2RA;
    const double lo = lon[i] * DE2RA;
    c.lat_rad[i] = la;
    c.lon_rad[i] = lo;
    c.cos_lat[i] = std::cos(la);
    c.sin_lat[i] = std::sin(la);
    if (need_xyz) {
      c.x3[i] = c.cos_lat[i] * std::cos(lo);
      c.y3[i] = c.cos_lat[i] * std::sin(lo);
      c.z3[i] = c.sin_lat[i];
    }
  }
  return c;
}

// Row-major flat score buffer. `s[i*k + kk] = e[i] * X(i, kk)`. Computing this
// once up-front (a) replaces the per-pair `e * X` recomputation in the hot
// loops, (b) gives the meat workers row-contiguous reads instead of striding
// across Armadillo's column-major storage, and (c) lets us pass a single
// `const RowMajorScores&` through the worker hierarchy instead of `(X, e)`.
struct RowMajorScores {
  std::size_t n;
  std::size_t k;
  std::vector<double> s;

  RowMajorScores(const arma::mat& X, const arma::vec& e)
      : n(X.n_rows), k(X.n_cols), s(static_cast<std::size_t>(X.n_rows) * X.n_cols) {
    for (std::size_t i = 0; i < n; ++i) {
      const double ei = e[i];
      for (std::size_t kk = 0; kk < k; ++kk) {
        s[i * k + kk] = ei * X(i, kk);
      }
    }
  }

  // Build by reordering an existing buffer: row `pos` of the result is row
  // `perm[pos]` of `src`. Used to permute scores into the lat-sort order so
  // the meat workers see sequential access into `s`.
  RowMajorScores(const RowMajorScores& src, const std::vector<std::size_t>& perm)
      : n(perm.size()), k(src.k), s(perm.size() * src.k) {
    for (std::size_t pos = 0; pos < n; ++pos) {
      const double* src_row = src.row(perm[pos]);
      double* dst_row = s.data() + pos * k;
      for (std::size_t kk = 0; kk < k; ++kk) dst_row[kk] = src_row[kk];
    }
  }

  inline const double* row(std::size_t i) const {
    return s.data() + i * k;
  }
};

// Reorder a CoordCache by an index permutation. The output is the same length
// as `perm` and contains `src` values at positions `perm[pos]`. Empty xyz
// vectors in `src` (i.e., dist_id not in {SPHERICAL, CHORD}) are preserved as
// empty in the output.
CoordCache permute_coord_cache(const CoordCache& src,
                               const std::vector<std::size_t>& perm) {
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
  for (std::size_t pos = 0; pos < n; ++pos) {
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
// the per-pair index footprint (build_csr guards n_per < 2^32). weight is
// only populated for KERNEL_BARTLETT — the uniform meat branch never reads
// it, so leaving it empty saves 8 bytes per pair.
struct CsrGraph {
  std::vector<std::size_t> row_ptr;
  std::vector<uint32_t> col_idx;
  std::vector<double> weight;
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
  const double cutoff;
  const ScreenParams screen;

  CsrFillWorkerT(const std::vector<std::size_t>& sorted_idx,
                 const std::vector<std::size_t>& row_end,
                 const CoordCache& c,
                 const std::vector<std::size_t>& row_ptr,
                 std::vector<uint32_t>& col_idx,
                 std::vector<double>& weight,
                 double cutoff, const ScreenParams& screen)
      : sorted_idx(sorted_idx), row_end(row_end), c(c),
        row_ptr(row_ptr), col_idx(col_idx),
        weight(weight), cutoff(cutoff), screen(screen) {}

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
          if (K != KERNEL_UNIFORM) weight[out] = w;
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
                   double cutoff, int ncores) {
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
  if (K != KERNEL_UNIFORM) g.weight.assign(nnz, 0.0);

  CsrFillWorkerT<D, K> fill_worker(sorted_idx, row_end, c,
                                   g.row_ptr, g.col_idx, g.weight,
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
                       std::vector<std::size_t>& sorted_idx, CellGrid& grid) {
  const std::size_t nb = be - bs;
  const double cell = spec.cell;
  const std::int64_t G = spec.G;
  const std::uint64_t uG = static_cast<std::uint64_t>(G);

  std::vector<std::uint64_t> cid(nb);
  for (std::size_t i = 0; i < nb; ++i) {
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

  // (cell id, original index) order: deterministic, and rows of a cell are
  // contiguous with consecutive cells contiguous in row space.
  std::vector<std::size_t> ord(nb);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
    if (cid[a] != cid[b]) return cid[a] < cid[b];
    return a < b;
  });

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
  grid.nbr.resize(grid.nbr.size() + 10 * ncb);
  for (std::size_t cdx = 0; cdx < ncb; ++cdx) {
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

// Grid analogue of StreamMeatGeneralWorkerT: fused candidate scan + score
// accumulation, O(n) memory, no CSR. Inputs are pre-permuted to grid order.
template<int D, int K>
struct StreamMeatGridWorkerT : public Worker {
  const RowMajorScores& S;
  const CoordCache& coord;
  const CellGrid& grid;
  const double cutoff;
  const ScreenParams screen;
  const std::size_t k;
  arma::mat meat;

  StreamMeatGridWorkerT(const RowMajorScores& S, const CoordCache& coord,
                        const CellGrid& grid, double cutoff,
                        const ScreenParams& screen)
      : S(S), coord(coord), grid(grid), cutoff(cutoff), screen(screen),
        k(S.k), meat(k, k, arma::fill::zeros) {}

  StreamMeatGridWorkerT(const StreamMeatGridWorkerT& other, Split)
      : S(other.S), coord(other.coord), grid(other.grid),
        cutoff(other.cutoff), screen(other.screen), k(other.k),
        meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t pos = begin; pos < end; ++pos) {
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
  }

  void join(const StreamMeatGridWorkerT& rhs) {
    meat += rhs.meat;
  }
};

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
  const double cutoff;
  const ScreenParams screen;

  GridCsrFillWorkerT(const CellGrid& grid, const CoordCache& c,
                     const std::vector<std::size_t>& row_ptr,
                     std::vector<std::uint32_t>& col_idx,
                     std::vector<double>& weight,
                     double cutoff, const ScreenParams& screen)
      : grid(grid), c(c), row_ptr(row_ptr), col_idx(col_idx), weight(weight),
        cutoff(cutoff), screen(screen) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pos = begin; pos < end; ++pos) {
      std::size_t out = row_ptr[pos];
      for_each_grid_candidate(grid, pos, [&](std::size_t q) {
        const double w = pair_weight<D, K>(c, pos, q, cutoff, screen);
        if (w != 0.0) {
          col_idx[out] = static_cast<std::uint32_t>(q);
          if (K != KERNEL_UNIFORM) weight[out] = w;
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
                        double cutoff, const ScreenParams& screen, int ncores) {
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
  if (K != KERNEL_UNIFORM) g.weight.assign(nnz, 0.0);

  GridCsrFillWorkerT<D, K> fill_worker(grid, c, g.row_ptr, g.col_idx,
                                       g.weight, cutoff, screen);
  if (ncores > 1) parallelFor(0, n, fill_worker);
  else fill_worker(0, n);

  return g;
}

// Inputs `S` and `coord` are pre-permuted into lat-sorted order by the caller,
// so the worker indexes them by sorted position `pos` directly. This makes
// every read of `S.row(...)` and `coord.lat_rad[...]` sequential across the
// outer loop and within the inner lat window.
template<int D, int K>
struct StreamMeatGeneralWorkerT : public Worker {
  const RowMajorScores& S;
  const std::vector<std::size_t>& row_end;
  const CoordCache& coord;
  const double cutoff;
  const ScreenParams screen;
  const std::size_t k;
  arma::mat meat;

  StreamMeatGeneralWorkerT(const RowMajorScores& S,
                           const std::vector<std::size_t>& row_end,
                           const CoordCache& coord,
                           double cutoff,
                           const ScreenParams& screen)
      : S(S), row_end(row_end), coord(coord),
        cutoff(cutoff), screen(screen),
        k(S.k), meat(k, k, arma::fill::zeros) {}

  StreamMeatGeneralWorkerT(const StreamMeatGeneralWorkerT& other, Split)
      : S(other.S), row_end(other.row_end),
        coord(other.coord), cutoff(other.cutoff), screen(other.screen), k(other.k),
        meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t pos = begin; pos < end; ++pos) {
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
  }

  void join(const StreamMeatGeneralWorkerT& rhs) {
    meat += rhs.meat;
  }
};

// `S` is pre-permuted: within each time block, row at sorted-position `pos`
// lives at absolute index `block_start[b] + pos`. `col_idx` from the CSR
// likewise stores sorted positions in `[0, n_per)`, so the worker indexes `S`
// sequentially without an `sorted_rel`/`col_rel` indirection.
template<int K>
struct MeatBalancedWorkerT : public Worker {
  const RowMajorScores& S;
  const std::vector<std::size_t>& block_start;
  const std::vector<std::size_t>& row_ptr;
  const std::vector<uint32_t>& col_idx;
  const std::vector<double>& weight;
  const std::size_t n_per;
  const std::size_t k;
  arma::mat meat;

  MeatBalancedWorkerT(const RowMajorScores& S,
                      const std::vector<std::size_t>& block_start,
                      std::size_t n_per,
                      const std::vector<std::size_t>& row_ptr,
                      const std::vector<uint32_t>& col_idx,
                      const std::vector<double>& weight)
      : S(S), block_start(block_start),
        row_ptr(row_ptr), col_idx(col_idx), weight(weight),
        n_per(n_per), k(S.k), meat(k, k, arma::fill::zeros) {}

  MeatBalancedWorkerT(const MeatBalancedWorkerT& other, Split)
      : S(other.S), block_start(other.block_start),
        row_ptr(other.row_ptr), col_idx(other.col_idx), weight(other.weight),
        n_per(other.n_per), k(other.k), meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t task = begin; task < end; ++task) {
      const std::size_t b = task / n_per;
      const std::size_t pos = task - b * n_per;
      const std::size_t base = block_start[b];
      const std::size_t i = base + pos;
      const double* si = S.row(i);

      for (std::size_t kk = 0; kk < k; ++kk) {
        c[kk] = 0.5 * si[kk];
      }

      for (std::size_t ep = row_ptr[pos]; ep < row_ptr[pos + 1]; ++ep) {
        const std::size_t j = base + col_idx[ep];
        const double* sj = S.row(j);
        if (K == KERNEL_UNIFORM) {
          for (std::size_t kk = 0; kk < k; ++kk) {
            c[kk] += sj[kk];
          }
        } else {
          const double w = weight[ep];
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
  }

  void join(const MeatBalancedWorkerT& rhs) {
    meat += rhs.meat;
  }
};

template<int K>
arma::mat meat_from_csr_balanced(const RowMajorScores& S,
                                 const std::vector<std::size_t>& block_start,
                                 std::size_t n_per,
                                 const CsrGraph& graph,
                                 int ncores) {
  MeatBalancedWorkerT<K> worker(S, block_start, n_per,
                                graph.row_ptr, graph.col_idx, graph.weight);
  const std::size_t tasks = block_start.size() * n_per;
  if (ncores > 1) parallelReduce(0, tasks, worker);
  else worker(0, tasks);
  return worker.meat + worker.meat.t();
}

template<int D, int K>
arma::mat fast_spatial_general(const arma::vec& lat, const arma::vec& lon,
                               const arma::vec& time, const RowMajorScores& S,
                               double cutoff, int ncores) {
  const std::size_t n = S.n;
  const CoordCache c = make_coord_cache(lat, lon, D);

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

  // Reorder coord cache and scores into lat-sorted order. The worker then
  // indexes both buffers by sorted position directly, so every read is
  // sequential -- independent of how the caller laid out the input.
  const CoordCache c_sorted = permute_coord_cache(c, sorted_idx);
  const RowMajorScores S_sorted(S, sorted_idx);

  const ScreenParams screen = make_screen_params(cutoff, D);
  StreamMeatGeneralWorkerT<D, K> worker(S_sorted, row_end, c_sorted, cutoff, screen);
  if (ncores > 1) parallelReduce(0, n, worker);
  else worker(0, n);
  return worker.meat + worker.meat.t();
}

template<int D, int K>
arma::mat fast_spatial_balanced(const arma::vec& lat, const arma::vec& lon,
                                const arma::vec& time, const RowMajorScores& S,
                                double cutoff, int ncores) {
  const TimeBlocks blocks = make_time_blocks(time);
  const std::size_t n_per = blocks.end[0] - blocks.start[0];
  const CoordCache c = make_coord_cache(lat, lon, D);

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
  const CoordCache c_block0_sorted = permute_coord_cache(c, sorted_abs);
  std::vector<std::size_t> identity_perm(n_per);
  std::iota(identity_perm.begin(), identity_perm.end(), 0);
  std::vector<std::size_t> row_end(n_per, n_per);
  CsrGraph graph = build_csr<D, K>(identity_perm, row_end, c_block0_sorted, cutoff, ncores);

  // Permute scores per block by the same permutation, so the meat worker
  // can index `S_sorted.row(base + pos)` directly.
  std::vector<std::size_t> global_perm(S.n);
  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    const std::size_t base = blocks.start[b];
    for (std::size_t pos = 0; pos < n_per; ++pos) {
      global_perm[base + pos] = base + sorted_rel[pos];
    }
  }
  const RowMajorScores S_sorted(S, global_perm);

  return meat_from_csr_balanced<K>(S_sorted, blocks.start, n_per, graph, ncores);
}

// Grid version of the general path: per-block cell grids, fused scan +
// accumulate, O(n) memory.
template<int D, int K>
arma::mat fast_spatial_general_grid(const arma::vec& lat, const arma::vec& lon,
                                    const arma::vec& time, const RowMajorScores& S,
                                    double cutoff, int ncores) {
  const std::size_t n = S.n;
  if (n > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
    Rcpp::stop("The grid neighbor path supports at most 2^32 - 1 rows; "
               "use neighbor = \"band\".");
  }
  const CoordCache c = make_coord_cache(lat, lon, D);
  const ScreenParams screen = make_screen_params(cutoff, D);
  const GridSpec spec = make_grid_spec(screen);

  const TimeBlocks blocks = make_time_blocks(time);
  std::vector<std::size_t> sorted_idx;
  sorted_idx.reserve(n);
  CellGrid grid;
  grid.row_cell.resize(n);
  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    append_block_grid(c, blocks.start[b], blocks.end[b], spec, sorted_idx, grid);
  }
  grid.cell_start.push_back(n);

  const CoordCache c_sorted = permute_coord_cache(c, sorted_idx);
  const RowMajorScores S_sorted(S, sorted_idx);

  StreamMeatGridWorkerT<D, K> worker(S_sorted, c_sorted, grid, cutoff, screen);
  if (ncores > 1) parallelReduce(0, n, worker);
  else worker(0, n);
  return worker.meat + worker.meat.t();
}

// Grid version of the balanced path: cell grid + CSR from block 0, reused
// across periods; the meat worker is shared with the band version.
template<int D, int K>
arma::mat fast_spatial_balanced_grid(const arma::vec& lat, const arma::vec& lon,
                                     const arma::vec& time, const RowMajorScores& S,
                                     double cutoff, int ncores) {
  const TimeBlocks blocks = make_time_blocks(time);
  const std::size_t n_per = blocks.end[0] - blocks.start[0];
  const CoordCache c = make_coord_cache(lat, lon, D);
  const ScreenParams screen = make_screen_params(cutoff, D);
  const GridSpec spec = make_grid_spec(screen);

  std::vector<std::size_t> sorted_abs;
  sorted_abs.reserve(n_per);
  CellGrid g0;
  g0.row_cell.resize(n_per);
  append_block_grid(c, blocks.start[0], blocks.end[0], spec, sorted_abs, g0);
  g0.cell_start.push_back(n_per);

  const CoordCache c0_sorted = permute_coord_cache(c, sorted_abs);
  CsrGraph graph = build_csr_grid<D, K>(g0, c0_sorted, cutoff, screen, ncores);

  std::vector<std::size_t> sorted_rel(n_per);
  for (std::size_t i = 0; i < n_per; ++i) {
    sorted_rel[i] = sorted_abs[i] - blocks.start[0];
  }
  std::vector<std::size_t> global_perm(S.n);
  for (std::size_t b = 0; b < blocks.start.size(); ++b) {
    const std::size_t base = blocks.start[b];
    for (std::size_t pos = 0; pos < n_per; ++pos) {
      global_perm[base + pos] = base + sorted_rel[pos];
    }
  }
  const RowMajorScores S_sorted(S, global_perm);

  return meat_from_csr_balanced<K>(S_sorted, blocks.start, n_per, graph, ncores);
}

// 8-way dispatch on (dist_id, kernel_id) -> template instantiation. Called
// once per FastSpatialMeat invocation; after this point everything is
// compile-time specialized.
template<int D, int K>
inline arma::mat fast_spatial_dispatch(const arma::vec& lat, const arma::vec& lon,
                                       const arma::vec& time, const RowMajorScores& S,
                                       double cutoff, int ncores, bool balanced,
                                       bool use_grid) {
  if (balanced) {
    if (use_grid) return fast_spatial_balanced_grid<D, K>(lat, lon, time, S, cutoff, ncores);
    return fast_spatial_balanced<D, K>(lat, lon, time, S, cutoff, ncores);
  }
  if (use_grid) return fast_spatial_general_grid<D, K>(lat, lon, time, S, cutoff, ncores);
  return fast_spatial_general<D, K>(lat, lon, time, S, cutoff, ncores);
}

arma::mat dispatch_spatial(int dist_id, int kernel_id,
                           const arma::vec& lat, const arma::vec& lon,
                           const arma::vec& time, const RowMajorScores& S,
                           double cutoff, int ncores, bool balanced,
                           bool use_grid) {
  const int tag = (dist_id << 4) | kernel_id;
  switch (tag) {
    case (DIST_HAVERSINE << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_HAVERSINE, KERNEL_UNIFORM>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
    case (DIST_HAVERSINE << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_HAVERSINE, KERNEL_BARTLETT>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
    case (DIST_SPHERICAL << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_SPHERICAL, KERNEL_UNIFORM>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
    case (DIST_SPHERICAL << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_SPHERICAL, KERNEL_BARTLETT>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
    case (DIST_CHORD << 4) | KERNEL_UNIFORM:
      return fast_spatial_dispatch<DIST_CHORD, KERNEL_UNIFORM>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
    case (DIST_CHORD << 4) | KERNEL_BARTLETT:
      return fast_spatial_dispatch<DIST_CHORD, KERNEL_BARTLETT>(lat, lon, time, S, cutoff, ncores, balanced, use_grid);
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

struct SerialHacPanelWorker : public Worker {
  const arma::vec& times;
  const RowMajorScores& S;
  const UnitBlocks& blocks;
  const double cutoff;
  const std::size_t k;
  arma::mat meat;

  SerialHacPanelWorker(const arma::vec& times, double cutoff,
                       const RowMajorScores& S,
                       const UnitBlocks& blocks)
      : times(times), S(S), blocks(blocks), cutoff(cutoff),
        k(S.k), meat(k, k, arma::fill::zeros) {}

  SerialHacPanelWorker(const SerialHacPanelWorker& other, Split)
      : times(other.times), S(other.S), blocks(other.blocks),
        cutoff(other.cutoff), k(other.k),
        meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t b_begin, std::size_t b_end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t bi = b_begin; bi < b_end; ++bi) {
      const std::size_t bs = blocks.start[bi];
      const std::size_t be = blocks.end[bi];
      for (std::size_t i = bs; i < be; ++i) {
        std::fill(c.begin(), c.end(), 0.0);
        for (std::size_t j = bs; j < be; ++j) {
          const double dt = std::fabs(times[j] - times[i]);
          if (dt <= cutoff && dt != 0.0) {
            const double w = 1.0 - dt / (cutoff + 1.0);
            const double* sj = S.row(j);
            for (std::size_t kk = 0; kk < k; ++kk) {
              c[kk] += w * sj[kk];
            }
          }
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
  }

  void join(const SerialHacPanelWorker& rhs) {
    meat += rhs.meat;
  }
};

} // anonymous namespace

// [[Rcpp::export]]
arma::mat FastSpatialMeat(arma::vec lat, arma::vec lon, arma::vec time,
                          arma::mat X, arma::vec e, double cutoff,
                          std::string kernel = "bartlett",
                          std::string dist_fn = "haversine",
                          bool balanced_pnl = false,
                          int ncores = 1,
                          std::string neighbor = "grid") {
  if (X.n_rows != e.n_elem || X.n_rows != lat.n_elem || X.n_rows != lon.n_elem ||
      X.n_rows != time.n_elem) {
    Rcpp::stop("lat, lon, time, X, and e have incompatible lengths.");
  }

  ncores = std::max(1, ncores);
  const int kernel_id = parse_kernel_id(kernel);
  const int dist_id = parse_dist_id(dist_fn);
  const TimeBlocks blocks = make_time_blocks(time);

  bool use_grid;
  if (neighbor == "grid") use_grid = true;
  else if (neighbor == "band") use_grid = false;
  else Rcpp::stop("Unknown neighbor: %s (use \"grid\" or \"band\")", neighbor.c_str());
  // cutoff < 0 is the no-pairs sentinel; the band path's break logic handles
  // it (only the 0.5*S_i diagonal survives), so route it there.
  if (cutoff < 0.0) use_grid = false;

  if (X.n_rows == 0) {
    return arma::mat(X.n_cols, X.n_cols, arma::fill::zeros);
  }

  const RowMajorScores S(X, e);
  const bool use_balanced = balanced_pnl && blocks.start.size() > 1 && same_block_size(blocks);

  if (balanced_pnl && blocks.start.size() > 1 && !same_block_size(blocks)) {
    Rcpp::warning("balanced_pnl = TRUE but time blocks have unequal sizes; using general CSR path.");
  }

  return dispatch_spatial(dist_id, kernel_id, lat, lon, time, S, cutoff, ncores,
                          use_balanced, use_grid);
}

// [[Rcpp::export]]
arma::mat FastSerialHacPanel(arma::vec unit, arma::vec time, double cutoff,
                             arma::mat X, arma::vec e, int ncores = 1) {
  if (X.n_rows != e.n_elem || X.n_rows != unit.n_elem ||
      X.n_rows != time.n_elem) {
    Rcpp::stop("unit, time, X, and e have incompatible lengths.");
  }
  if (X.n_rows == 0) {
    return arma::mat(X.n_cols, X.n_cols, arma::fill::zeros);
  }

  ncores = std::max(1, ncores);
  const UnitBlocks blocks = make_unit_blocks(unit);
  const RowMajorScores S(X, e);
  SerialHacPanelWorker worker(time, cutoff, S, blocks);
  if (ncores > 1) parallelReduce(0, blocks.start.size(), worker);
  else worker(0, blocks.start.size());
  return worker.meat;
}
