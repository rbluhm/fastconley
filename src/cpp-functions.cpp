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
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double TWO_PI = 2.0 * PI;
constexpr double DE2RA = PI / 180.0;
constexpr double RA2DE = 180.0 / PI;
constexpr double AVG_ERAD = 6371.0;
constexpr double FLAT_KM_PER_DEG = 111.0;
constexpr double EPS = 1e-14;

enum DistId { DIST_HAVERSINE = 1, DIST_SPHERICAL = 2, DIST_CHORD = 3, DIST_FLATEARTH = 4 };
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

inline double dist_haversine_idx(const CoordCache& c, std::size_t i, std::size_t j) {
  const double dlat = c.lat_rad[j] - c.lat_rad[i];
  const double dlon = c.lon_rad[j] - c.lon_rad[i];
  double a = sq(std::sin(dlat / 2.0)) + c.cos_lat[i] * c.cos_lat[j] * sq(std::sin(dlon / 2.0));
  if (a < 0.0) a = 0.0;
  if (a > 1.0) a = 1.0;
  return AVG_ERAD * 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
}

inline double dist_spherical_idx(const CoordCache& c, std::size_t i, std::size_t j) {
  const double d = c.sin_lat[i] * c.sin_lat[j] +
                   c.cos_lat[i] * c.cos_lat[j] * std::cos(c.lon_rad[i] - c.lon_rad[j]);
  return AVG_ERAD * safe_acos(d);
}

inline double dist_chord_idx(const CoordCache& c, std::size_t i, std::size_t j) {
  const double dx = c.x3[i] - c.x3[j];
  const double dy = c.y3[i] - c.y3[j];
  const double dz = c.z3[i] - c.z3[j];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double dist_flatearth_idx(const CoordCache& c, std::size_t i, std::size_t j) {
  const double dlat_km = FLAT_KM_PER_DEG * (c.lat_rad[i] - c.lat_rad[j]) * RA2DE;
  const double dlon_km = c.cos_lat[i] * FLAT_KM_PER_DEG * (c.lon_rad[i] - c.lon_rad[j]) * RA2DE;
  return std::sqrt(dlat_km * dlat_km + dlon_km * dlon_km);
}

inline double distance_idx(int dist_id, const CoordCache& c,
                           std::size_t i, std::size_t j) {
  if (dist_id == DIST_HAVERSINE) return dist_haversine_idx(c, i, j);
  if (dist_id == DIST_SPHERICAL) return dist_spherical_idx(c, i, j);
  if (dist_id == DIST_CHORD)     return dist_chord_idx(c, i, j);
  return dist_flatearth_idx(c, i, j);
}

int parse_dist_id(const std::string& dist_fn) {
  if (dist_fn == "haversine") return DIST_HAVERSINE;
  if (dist_fn == "spherical") return DIST_SPHERICAL;
  if (dist_fn == "chord") return DIST_CHORD;
  if (dist_fn == "flatearth") return DIST_FLATEARTH;
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
};

ScreenParams make_screen_params(double cutoff, int dist_id) {
  ScreenParams p;
  if (cutoff < 0.0) {
    p.lat_cutoff_rad = -1.0;
    p.angular_cutoff_rad = -1.0;
    p.sin2_half_angular_cutoff = 0.0;
    return p;
  }

  if (dist_id == DIST_CHORD) {
    const double ratio = std::min(1.0, cutoff / (2.0 * AVG_ERAD));
    p.angular_cutoff_rad = 2.0 * std::asin(ratio);
  } else if (dist_id == DIST_FLATEARTH) {
    p.angular_cutoff_rad = (cutoff / FLAT_KM_PER_DEG) * DE2RA;
  } else {
    p.angular_cutoff_rad = cutoff / AVG_ERAD;
  }

  if (p.angular_cutoff_rad > PI) p.angular_cutoff_rad = PI;
  p.lat_cutoff_rad = p.angular_cutoff_rad;
  p.sin2_half_angular_cutoff = sq(std::sin(p.angular_cutoff_rad / 2.0));
  return p;
}

inline bool passes_lon_screen(int dist_id, const CoordCache& c,
                              std::size_t i, std::size_t j,
                              double cutoff,
                              const ScreenParams& screen) {
  const double cos_i = c.cos_lat[i];
  const double cos_j = c.cos_lat[j];
  const double lon_i = c.lon_rad[i];
  const double lon_j = c.lon_rad[j];

  if (dist_id == DIST_FLATEARTH) {
    const double cos_abs = std::fabs(cos_i);
    if (cos_abs < EPS) return true;
    const double dlon_deg = std::fabs((lon_j - lon_i) * RA2DE);
    return (cos_abs * FLAT_KM_PER_DEG * dlon_deg) <= cutoff + 1e-12;
  }

  if (screen.angular_cutoff_rad >= PI) return true;
  const double denom = cos_i * cos_j;
  if (denom <= EPS) return true;

  const double dlon = lon_abs_wrapped(lon_i, lon_j);
  const double sin2_half_dlon = sq(std::sin(dlon / 2.0));

  // Conservative necessary condition from the haversine formula, ignoring the
  // latitude term. If this fails, the pair cannot be within the cutoff.
  return sin2_half_dlon <= screen.sin2_half_angular_cutoff / denom + 1e-15;
}

inline double kernel_weight(double dist, double cutoff, int kernel_id) {
  if (dist > cutoff) return 0.0;
  if (kernel_id == KERNEL_UNIFORM) return 1.0;
  if (cutoff <= 0.0) return dist <= 0.0 ? 1.0 : 0.0;
  return 1.0 - dist / cutoff;
}

inline double neighbor_weight(const CoordCache& c,
                              std::size_t i, std::size_t j,
                              double cutoff, int kernel_id, int dist_id,
                              const ScreenParams& screen) {
  if (!passes_lon_screen(dist_id, c, i, j, cutoff, screen)) return 0.0;
  return kernel_weight(distance_idx(dist_id, c, i, j), cutoff, kernel_id);
}

CoordCache make_coord_cache(const arma::vec& lat, const arma::vec& lon,
                            int dist_id) {
  const std::size_t n = lat.n_elem;
  CoordCache c;
  c.lat_rad.resize(n);
  c.lon_rad.resize(n);
  c.cos_lat.resize(n);
  c.sin_lat.resize(n);
  const bool need_xyz = (dist_id == DIST_CHORD);
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
      c.x3[i] = AVG_ERAD * c.cos_lat[i] * std::cos(lo);
      c.y3[i] = AVG_ERAD * c.cos_lat[i] * std::sin(lo);
      c.z3[i] = AVG_ERAD * c.sin_lat[i];
    }
  }
  return c;
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

struct CsrGraph {
  std::vector<std::size_t> row_ptr;
  std::vector<std::size_t> col_idx;
  std::vector<double> weight;
};

struct CsrCountWorker : public Worker {
  const std::vector<std::size_t>& sorted_idx;
  const std::vector<std::size_t>& row_end;
  const CoordCache& c;
  std::vector<std::size_t>& counts;
  const double cutoff;
  const int kernel_id;
  const int dist_id;
  const ScreenParams screen;

  CsrCountWorker(const std::vector<std::size_t>& sorted_idx,
                 const std::vector<std::size_t>& row_end,
                 const CoordCache& c,
                 std::vector<std::size_t>& counts,
                 double cutoff, int kernel_id, int dist_id,
                 const ScreenParams& screen)
      : sorted_idx(sorted_idx), row_end(row_end), c(c),
        counts(counts), cutoff(cutoff),
        kernel_id(kernel_id), dist_id(dist_id), screen(screen) {}

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

        const double w = neighbor_weight(c, i, j, cutoff, kernel_id, dist_id, screen);
        if (w != 0.0) ++n_found;
      }
      counts[pos] = n_found;
    }
  }
};

struct CsrFillWorker : public Worker {
  const std::vector<std::size_t>& sorted_idx;
  const std::vector<std::size_t>& row_end;
  const CoordCache& c;
  const std::vector<std::size_t>& row_ptr;
  std::vector<std::size_t>& col_idx;
  std::vector<double>& weight;
  const double cutoff;
  const int kernel_id;
  const int dist_id;
  const ScreenParams screen;

  CsrFillWorker(const std::vector<std::size_t>& sorted_idx,
                const std::vector<std::size_t>& row_end,
                const CoordCache& c,
                const std::vector<std::size_t>& row_ptr,
                std::vector<std::size_t>& col_idx,
                std::vector<double>& weight,
                double cutoff, int kernel_id, int dist_id,
                const ScreenParams& screen)
      : sorted_idx(sorted_idx), row_end(row_end), c(c),
        row_ptr(row_ptr), col_idx(col_idx),
        weight(weight), cutoff(cutoff), kernel_id(kernel_id), dist_id(dist_id),
        screen(screen) {}

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

        const double w = neighbor_weight(c, i, j, cutoff, kernel_id, dist_id, screen);
        if (w != 0.0) {
          col_idx[out] = j;
          weight[out] = w;
          ++out;
        }
      }
    }
  }
};

CsrGraph build_csr(const std::vector<std::size_t>& sorted_idx,
                   const std::vector<std::size_t>& row_end,
                   const CoordCache& c,
                   double cutoff, int kernel_id, int dist_id,
                   int ncores) {
  const std::size_t n = sorted_idx.size();
  CsrGraph g;
  g.row_ptr.assign(n + 1, 0);

  if (n == 0 || cutoff < 0.0) return g;

  std::vector<std::size_t> counts(n, 0);
  const ScreenParams screen = make_screen_params(cutoff, dist_id);

  CsrCountWorker count_worker(sorted_idx, row_end, c,
                              counts, cutoff, kernel_id, dist_id, screen);
  if (ncores > 1) parallelFor(0, n, count_worker);
  else count_worker(0, n);

  for (std::size_t i = 0; i < n; ++i) {
    g.row_ptr[i + 1] = g.row_ptr[i] + counts[i];
  }

  const std::size_t nnz = g.row_ptr[n];
  g.col_idx.assign(nnz, 0);
  g.weight.assign(nnz, 0.0);

  CsrFillWorker fill_worker(sorted_idx, row_end, c,
                            g.row_ptr, g.col_idx, g.weight,
                            cutoff, kernel_id, dist_id, screen);
  if (ncores > 1) parallelFor(0, n, fill_worker);
  else fill_worker(0, n);

  return g;
}

struct StreamMeatGeneralWorker : public Worker {
  const arma::mat& X;
  const arma::vec& e;
  const std::vector<std::size_t>& sorted_idx;
  const std::vector<std::size_t>& row_end;
  const CoordCache& coord;
  const double cutoff;
  const int kernel_id;
  const int dist_id;
  const ScreenParams screen;
  const std::size_t k;
  arma::mat meat;

  StreamMeatGeneralWorker(const arma::mat& X, const arma::vec& e,
                          const std::vector<std::size_t>& sorted_idx,
                          const std::vector<std::size_t>& row_end,
                          const CoordCache& coord,
                          double cutoff, int kernel_id, int dist_id,
                          const ScreenParams& screen)
      : X(X), e(e), sorted_idx(sorted_idx), row_end(row_end), coord(coord),
        cutoff(cutoff), kernel_id(kernel_id), dist_id(dist_id), screen(screen),
        k(X.n_cols), meat(k, k, arma::fill::zeros) {}

  StreamMeatGeneralWorker(const StreamMeatGeneralWorker& other, Split)
      : X(other.X), e(other.e), sorted_idx(other.sorted_idx), row_end(other.row_end),
        coord(other.coord), cutoff(other.cutoff), kernel_id(other.kernel_id),
        dist_id(other.dist_id), screen(other.screen), k(other.k),
        meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t pos = begin; pos < end; ++pos) {
      const std::size_t i = sorted_idx[pos];
      const double ei = e[i];
      const double lat_i = coord.lat_rad[i];

      for (std::size_t kk = 0; kk < k; ++kk) {
        c[kk] = 0.5 * ei * X(i, kk);
      }

      const std::size_t pos_end = row_end[pos];
      for (std::size_t q = pos + 1; q < pos_end; ++q) {
        const std::size_t j = sorted_idx[q];
        const double dlat = coord.lat_rad[j] - lat_i;
        if (dlat > screen.lat_cutoff_rad + 1e-15) break;

        const double w = neighbor_weight(coord, i, j, cutoff, kernel_id, dist_id, screen);
        if (w == 0.0) continue;

        const double f = w * e[j];
        for (std::size_t kk = 0; kk < k; ++kk) {
          c[kk] += f * X(j, kk);
        }
      }

      for (std::size_t k1 = 0; k1 < k; ++k1) {
        const double s1 = ei * X(i, k1);
        for (std::size_t k2 = 0; k2 < k; ++k2) {
          meat(k1, k2) += s1 * c[k2];
        }
      }
    }
  }

  void join(const StreamMeatGeneralWorker& rhs) {
    meat += rhs.meat;
  }
};

struct MeatBalancedWorker : public Worker {
  const arma::mat& X;
  const arma::vec& e;
  const std::vector<std::size_t>& block_start;
  const std::vector<std::size_t>& sorted_rel;
  const std::vector<std::size_t>& row_ptr;
  const std::vector<std::size_t>& col_rel;
  const std::vector<double>& weight;
  const std::size_t n_per;
  const std::size_t k;
  arma::mat meat;

  MeatBalancedWorker(const arma::mat& X, const arma::vec& e,
                     const std::vector<std::size_t>& block_start,
                     const std::vector<std::size_t>& sorted_rel,
                     const std::vector<std::size_t>& row_ptr,
                     const std::vector<std::size_t>& col_rel,
                     const std::vector<double>& weight)
      : X(X), e(e), block_start(block_start), sorted_rel(sorted_rel),
        row_ptr(row_ptr), col_rel(col_rel), weight(weight), n_per(sorted_rel.size()),
        k(X.n_cols), meat(k, k, arma::fill::zeros) {}

  MeatBalancedWorker(const MeatBalancedWorker& other, Split)
      : X(other.X), e(other.e), block_start(other.block_start), sorted_rel(other.sorted_rel),
        row_ptr(other.row_ptr), col_rel(other.col_rel), weight(other.weight), n_per(other.n_per),
        k(other.k), meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    for (std::size_t task = begin; task < end; ++task) {
      const std::size_t b = task / n_per;
      const std::size_t pos = task - b * n_per;
      const std::size_t base = block_start[b];
      const std::size_t i = base + sorted_rel[pos];
      const double ei = e[i];

      for (std::size_t kk = 0; kk < k; ++kk) {
        c[kk] = 0.5 * ei * X(i, kk);
      }

      for (std::size_t ep = row_ptr[pos]; ep < row_ptr[pos + 1]; ++ep) {
        const std::size_t j = base + col_rel[ep];
        const double f = weight[ep] * e[j];
        for (std::size_t kk = 0; kk < k; ++kk) {
          c[kk] += f * X(j, kk);
        }
      }

      for (std::size_t k1 = 0; k1 < k; ++k1) {
        const double s1 = ei * X(i, k1);
        for (std::size_t k2 = 0; k2 < k; ++k2) {
          meat(k1, k2) += s1 * c[k2];
        }
      }
    }
  }

  void join(const MeatBalancedWorker& rhs) {
    meat += rhs.meat;
  }
};

arma::mat meat_from_csr_balanced(const arma::mat& X, const arma::vec& e,
                                 const std::vector<std::size_t>& block_start,
                                 const std::vector<std::size_t>& sorted_rel,
                                 const CsrGraph& graph,
                                 int ncores) {
  std::vector<std::size_t> col_rel(graph.col_idx.size());
  const std::size_t base0 = block_start[0];
  for (std::size_t i = 0; i < graph.col_idx.size(); ++i) {
    col_rel[i] = graph.col_idx[i] - base0;
  }

  MeatBalancedWorker worker(X, e, block_start, sorted_rel,
                            graph.row_ptr, col_rel, graph.weight);
  const std::size_t tasks = block_start.size() * sorted_rel.size();
  if (ncores > 1) parallelReduce(0, tasks, worker);
  else worker(0, tasks);
  return worker.meat + worker.meat.t();
}

arma::mat fast_spatial_general(const arma::vec& lat, const arma::vec& lon,
                               const arma::vec& time, const arma::mat& X,
                               const arma::vec& e, double cutoff,
                               int kernel_id, int dist_id, int ncores) {
  const std::size_t n = X.n_rows;
  const CoordCache c = make_coord_cache(lat, lon, dist_id);

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

  const ScreenParams screen = make_screen_params(cutoff, dist_id);
  StreamMeatGeneralWorker worker(X, e, sorted_idx, row_end, c,
                                 cutoff, kernel_id, dist_id, screen);
  if (ncores > 1) parallelReduce(0, sorted_idx.size(), worker);
  else worker(0, sorted_idx.size());
  return worker.meat + worker.meat.t();
}

arma::mat fast_spatial_balanced(const arma::vec& lat, const arma::vec& lon,
                                const arma::vec& time, const arma::mat& X,
                                const arma::vec& e, double cutoff,
                                int kernel_id, int dist_id, int ncores) {
  const TimeBlocks blocks = make_time_blocks(time);
  const std::size_t n_per = blocks.end[0] - blocks.start[0];
  const CoordCache c = make_coord_cache(lat, lon, dist_id);

  std::vector<std::size_t> sorted_abs;
  sort_block_indices(c.lat_rad, c.lon_rad, blocks.start[0], blocks.end[0], sorted_abs);
  std::vector<std::size_t> sorted_rel(sorted_abs.size());
  for (std::size_t i = 0; i < sorted_abs.size(); ++i) {
    sorted_rel[i] = sorted_abs[i] - blocks.start[0];
  }

  std::vector<std::size_t> row_end(n_per, n_per);
  CsrGraph graph = build_csr(sorted_abs, row_end, c, cutoff, kernel_id, dist_id, ncores);
  return meat_from_csr_balanced(X, e, blocks.start, sorted_rel, graph, ncores);
}

struct SerialHacWorker : public Worker {
  const arma::vec& times;
  const arma::mat& X;
  const arma::vec& e;
  const double cutoff;
  const std::size_t k;
  arma::mat meat;

  SerialHacWorker(const arma::vec& times, double cutoff,
                  const arma::mat& X, const arma::vec& e)
      : times(times), X(X), e(e), cutoff(cutoff), k(X.n_cols),
        meat(k, k, arma::fill::zeros) {}

  SerialHacWorker(const SerialHacWorker& other, Split)
      : times(other.times), X(other.X), e(other.e), cutoff(other.cutoff),
        k(other.k), meat(k, k, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> c(k, 0.0);
    const std::size_t n = X.n_rows;

    for (std::size_t i = begin; i < end; ++i) {
      std::fill(c.begin(), c.end(), 0.0);

      for (std::size_t j = 0; j < n; ++j) {
        const double dt = std::fabs(times[j] - times[i]);
        if (dt <= cutoff && dt != 0.0) {
          const double w = 1.0 - dt / (cutoff + 1.0);
          const double f = w * e[j];
          for (std::size_t kk = 0; kk < k; ++kk) {
            c[kk] += f * X(j, kk);
          }
        }
      }

      const double ei = e[i];
      for (std::size_t k1 = 0; k1 < k; ++k1) {
        const double s1 = ei * X(i, k1);
        for (std::size_t k2 = 0; k2 < k; ++k2) {
          meat(k1, k2) += s1 * c[k2];
        }
      }
    }
  }

  void join(const SerialHacWorker& rhs) {
    meat += rhs.meat;
  }
};

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
  const arma::mat& X;
  const arma::vec& e;
  const UnitBlocks& blocks;
  const double cutoff;
  const std::size_t k;
  arma::mat meat;

  SerialHacPanelWorker(const arma::vec& times, double cutoff,
                       const arma::mat& X, const arma::vec& e,
                       const UnitBlocks& blocks)
      : times(times), X(X), e(e), blocks(blocks), cutoff(cutoff),
        k(X.n_cols), meat(k, k, arma::fill::zeros) {}

  SerialHacPanelWorker(const SerialHacPanelWorker& other, Split)
      : times(other.times), X(other.X), e(other.e), blocks(other.blocks),
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
            const double f = w * e[j];
            for (std::size_t kk = 0; kk < k; ++kk) {
              c[kk] += f * X(j, kk);
            }
          }
        }
        const double ei = e[i];
        for (std::size_t k1 = 0; k1 < k; ++k1) {
          const double s1 = ei * X(i, k1);
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
                          int ncores = 1) {
  if (X.n_rows != e.n_elem || X.n_rows != lat.n_elem || X.n_rows != lon.n_elem ||
      X.n_rows != time.n_elem) {
    Rcpp::stop("lat, lon, time, X, and e have incompatible lengths.");
  }

  ncores = std::max(1, ncores);
  const int kernel_id = parse_kernel_id(kernel);
  const int dist_id = parse_dist_id(dist_fn);
  const TimeBlocks blocks = make_time_blocks(time);

  if (X.n_rows == 0) {
    return arma::mat(X.n_cols, X.n_cols, arma::fill::zeros);
  }

  if (balanced_pnl && blocks.start.size() > 1 && same_block_size(blocks)) {
    return fast_spatial_balanced(lat, lon, time, X, e, cutoff, kernel_id, dist_id, ncores);
  }

  if (balanced_pnl && blocks.start.size() > 1 && !same_block_size(blocks)) {
    Rcpp::warning("balanced_pnl = TRUE but time blocks have unequal sizes; using general CSR path.");
  }

  return fast_spatial_general(lat, lon, time, X, e, cutoff, kernel_id, dist_id, ncores);
}

// [[Rcpp::export]]
arma::mat TimeDist(arma::vec times, double cutoff,
                   arma::mat X, arma::vec e, int n1, int k, int ncores = 1) {
  if (X.n_rows != e.n_elem || X.n_rows != times.n_elem) {
    Rcpp::stop("times, X, and e have incompatible lengths.");
  }
  (void)n1;
  (void)k;

  ncores = std::max(1, ncores);
  SerialHacWorker worker(times, cutoff, X, e);
  if (ncores > 1) parallelReduce(0, static_cast<std::size_t>(X.n_rows), worker);
  else worker(0, static_cast<std::size_t>(X.n_rows));
  return worker.meat;
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
  SerialHacPanelWorker worker(time, cutoff, X, e, blocks);
  if (ncores > 1) parallelReduce(0, blocks.start.size(), worker);
  else worker(0, blocks.start.size());
  return worker.meat;
}
