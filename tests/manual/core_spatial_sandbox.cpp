// Standalone C++ sandbox check for the fast spatial meat logic.
// This does not use R/Rcpp; it mirrors the CSR + score accumulation algorithm.
// Compile from package root with:
//   g++ -std=c++11 tests/manual/core_spatial_sandbox.cpp -O2 -o /tmp/core_spatial_sandbox && /tmp/core_spatial_sandbox

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

static const double PI = 3.141592653589793238462643383279502884;
static const double DE2RA = PI / 180.0;
static const double R = 6371.0;

struct CSR {
  std::vector<size_t> ptr;
  std::vector<size_t> col;
  std::vector<double> w;
};

double sq(double x) { return x * x; }

double haversine(double lat1, double lon1, double lat2, double lon2) {
  lat1 *= DE2RA; lon1 *= DE2RA; lat2 *= DE2RA; lon2 *= DE2RA;
  double dlat = lat2 - lat1;
  double dlon = lon2 - lon1;
  double a = sq(std::sin(dlat / 2.0)) + std::cos(lat1) * std::cos(lat2) * sq(std::sin(dlon / 2.0));
  if (a < 0) a = 0;
  if (a > 1) a = 1;
  return R * 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
}

std::vector<std::pair<size_t, size_t>> blocks(const std::vector<int>& time) {
  std::vector<std::pair<size_t, size_t>> b;
  size_t start = 0;
  for (size_t i = 1; i < time.size(); ++i) {
    if (time[i] != time[i - 1]) {
      b.push_back({start, i});
      start = i;
    }
  }
  b.push_back({start, time.size()});
  return b;
}

CSR build_csr(const std::vector<size_t>& ord, const std::vector<size_t>& row_end,
              const std::vector<double>& lat, const std::vector<double>& lon,
              double cutoff, bool bartlett) {
  CSR g;
  g.ptr.assign(ord.size() + 1, 0);
  double lat_cut = cutoff / R / DE2RA; // degrees, conservative for this test range

  for (size_t p = 0; p < ord.size(); ++p) {
    size_t i = ord[p];
    for (size_t q = p + 1; q < row_end[p]; ++q) {
      size_t j = ord[q];
      if (lat[j] - lat[i] > lat_cut) break;
      double d = haversine(lat[i], lon[i], lat[j], lon[j]);
      if (d <= cutoff) ++g.ptr[p + 1];
    }
  }
  for (size_t p = 0; p < ord.size(); ++p) g.ptr[p + 1] += g.ptr[p];
  g.col.assign(g.ptr.back(), 0);
  g.w.assign(g.ptr.back(), 0.0);

  for (size_t p = 0; p < ord.size(); ++p) {
    size_t out = g.ptr[p];
    size_t i = ord[p];
    for (size_t q = p + 1; q < row_end[p]; ++q) {
      size_t j = ord[q];
      if (lat[j] - lat[i] > lat_cut) break;
      double d = haversine(lat[i], lon[i], lat[j], lon[j]);
      if (d <= cutoff) {
        g.col[out] = j;
        g.w[out] = bartlett ? 1.0 - d / cutoff : 1.0;
        ++out;
      }
    }
  }
  return g;
}

std::vector<std::vector<double>> fast_general(const std::vector<double>& lat, const std::vector<double>& lon,
                                             const std::vector<int>& time,
                                             const std::vector<std::vector<double>>& X,
                                             const std::vector<double>& e,
                                             double cutoff, bool bartlett) {
  size_t n = X.size(), k = X[0].size();
  std::vector<size_t> ord;
  std::vector<size_t> row_end;
  for (auto be : blocks(time)) {
    std::vector<size_t> tmp(be.second - be.first);
    std::iota(tmp.begin(), tmp.end(), be.first);
    std::sort(tmp.begin(), tmp.end(), [&](size_t a, size_t b) { return lat[a] < lat[b]; });
    size_t end_pos = ord.size() + tmp.size();
    ord.insert(ord.end(), tmp.begin(), tmp.end());
    row_end.resize(end_pos, end_pos);
  }
  CSR g = build_csr(ord, row_end, lat, lon, cutoff, bartlett);

  std::vector<std::vector<double>> M(k, std::vector<double>(k, 0.0));
  std::vector<double> c(k);
  for (size_t p = 0; p < ord.size(); ++p) {
    size_t i = ord[p];
    for (size_t kk = 0; kk < k; ++kk) c[kk] = 0.5 * e[i] * X[i][kk];
    for (size_t ep = g.ptr[p]; ep < g.ptr[p + 1]; ++ep) {
      size_t j = g.col[ep];
      double f = g.w[ep] * e[j];
      for (size_t kk = 0; kk < k; ++kk) c[kk] += f * X[j][kk];
    }
    for (size_t k1 = 0; k1 < k; ++k1) {
      double s1 = e[i] * X[i][k1];
      for (size_t k2 = 0; k2 < k; ++k2) M[k1][k2] += s1 * c[k2];
    }
  }
  auto out = M;
  for (size_t a = 0; a < k; ++a) for (size_t b = 0; b < k; ++b) out[a][b] = M[a][b] + M[b][a];
  return out;
}

std::vector<std::vector<double>> naive(const std::vector<double>& lat, const std::vector<double>& lon,
                                       const std::vector<int>& time,
                                       const std::vector<std::vector<double>>& X,
                                       const std::vector<double>& e,
                                       double cutoff, bool bartlett) {
  size_t n = X.size(), k = X[0].size();
  std::vector<std::vector<double>> out(k, std::vector<double>(k, 0.0));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (time[i] != time[j]) continue;
      double d = haversine(lat[i], lon[i], lat[j], lon[j]);
      if (d <= cutoff) {
        double w = bartlett ? 1.0 - d / cutoff : 1.0;
        for (size_t a = 0; a < k; ++a) {
          double si = e[i] * X[i][a];
          for (size_t b = 0; b < k; ++b) out[a][b] += w * si * e[j] * X[j][b];
        }
      }
    }
  }
  return out;
}

double max_abs_diff(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
  double ans = 0.0;
  for (size_t i = 0; i < A.size(); ++i) for (size_t j = 0; j < A[i].size(); ++j) ans = std::max(ans, std::fabs(A[i][j] - B[i][j]));
  return ans;
}

int main() {
  std::mt19937 gen(123);
  std::normal_distribution<double> nd(0.0, 1.0);
  std::uniform_real_distribution<double> latd(40.0, 43.0), lond(-75.0, -70.0);

  const size_t n = 12, k = 5;
  std::vector<std::vector<double>> X(n, std::vector<double>(k));
  std::vector<double> e(n), lat(n), lon(n);
  for (size_t i = 0; i < n; ++i) {
    X[i][0] = 1.0;
    for (size_t j = 1; j < k; ++j) X[i][j] = nd(gen);
    e[i] = nd(gen);
    lat[i] = latd(gen);
    lon[i] = lond(gen);
  }

  std::vector<int> time_xs(n, 1);
  std::vector<int> time_unbal = {1,1,1,1,1,2,2,2,3,3,3,3};

  auto f1 = fast_general(lat, lon, time_xs, X, e, 250.0, true);
  auto s1 = naive(lat, lon, time_xs, X, e, 250.0, true);
  auto f2 = fast_general(lat, lon, time_unbal, X, e, 250.0, false);
  auto s2 = naive(lat, lon, time_unbal, X, e, 250.0, false);

  double d1 = max_abs_diff(f1, s1);
  double d2 = max_abs_diff(f2, s2);
  std::cout << "cross-section Bartlett max abs diff: " << d1 << "\n";
  std::cout << "unbalanced uniform max abs diff: " << d2 << "\n";
  if (d1 > 1e-10 || d2 > 1e-10) return 1;
  return 0;
}
