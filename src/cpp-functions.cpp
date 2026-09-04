// Rcpp front-end for the fastconley engine. All numerics live in
// conley_core.h; this file only validates lengths, aliases R memory
// through the Armadillo advanced constructors (no input copies -- the R
// wrappers in R/engine.R guarantee REALSXP storage), and lets
// conley::Error propagate to Rcpp's BEGIN_RCPP/END_RCPP (RcppExports.cpp),
// which turns it into an R error carrying the same message.
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD 1
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "conley_core.h"

namespace {

bool r_interrupt_requested() {
  try {
    Rcpp::checkUserInterrupt();
    return false;
  } catch (...) {
    return true;
  }
}

inline void install_interrupt_hook() {
  conley::set_interrupt_hook(&r_interrupt_requested);
}

}  // namespace

// [[Rcpp::export]]
arma::mat FastSpatialMeat_cpp(Rcpp::NumericVector lat, Rcpp::NumericVector lon,
                          Rcpp::NumericVector time, Rcpp::NumericMatrix scores,
                          double cutoff,
                          std::string kernel = "bartlett",
                          std::string dist_fn = "haversine",
                          bool balanced_pnl = false,
                          int ncores = 1,
                          std::string neighbor = "grid",
                          std::string csr_weight = "double") {
  install_interrupt_hook();
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(lat.size()) != n ||
      static_cast<std::size_t>(lon.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("lat, lon, time, and scores have incompatible lengths.");
  }
  if (n == 0) {
    // Still parse the option strings so bad input errors as before.
    return conley::spatial_meat(arma::vec(), arma::vec(), arma::vec(),
                                arma::mat(0, k), cutoff, kernel, dist_fn,
                                balanced_pnl, ncores, neighbor, csr_weight);
  }

  const arma::vec lat_v(lat.begin(), n, false, true);
  const arma::vec lon_v(lon.begin(), n, false, true);
  const arma::vec time_v(time.begin(), n, false, true);
  const arma::mat S_col(scores.begin(), n, k, false, true);

  bool unbalanced_fallback = false;
  arma::mat meat = conley::spatial_meat(lat_v, lon_v, time_v, S_col, cutoff,
                                        kernel, dist_fn, balanced_pnl, ncores,
                                        neighbor, csr_weight,
                                        &unbalanced_fallback);
  if (unbalanced_fallback) {
    Rcpp::warning("balanced_pnl = TRUE but time blocks have unequal sizes; using the streaming path.");
  }
  return meat;
}

// [[Rcpp::export]]
arma::mat FastSerialHacPanel_cpp(Rcpp::NumericVector unit, Rcpp::NumericVector time,
                             double cutoff, Rcpp::NumericMatrix scores,
                             int ncores = 1) {
  install_interrupt_hook();
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(unit.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("unit, time, and scores have incompatible lengths.");
  }
  if (n == 0) {
    return conley::serial_hac_meat(arma::vec(), arma::vec(), cutoff,
                                   arma::mat(0, k), ncores);
  }
  const arma::vec unit_v(unit.begin(), n, false, true);
  const arma::vec time_v(time.begin(), n, false, true);
  const arma::mat S_col(scores.begin(), n, k, false, true);
  return conley::serial_hac_meat(unit_v, time_v, cutoff, S_col, ncores);
}

// [[Rcpp::export]]
arma::mat FastGridMeat_cpp(Rcpp::IntegerVector ring, Rcpp::IntegerVector col,
                       Rcpp::NumericVector time, Rcpp::NumericMatrix scores,
                       double lat0, double dlat, double dlon,
                       int n_ring, int n_col, int n_col_full,
                       double cutoff, std::string dist_fn = "spherical",
                       std::string kernel = "uniform",
                       int ncores = 1) {
  install_interrupt_hook();
  const std::size_t n = static_cast<std::size_t>(scores.nrow());
  const std::size_t k = static_cast<std::size_t>(scores.ncol());
  if (static_cast<std::size_t>(ring.size()) != n ||
      static_cast<std::size_t>(col.size()) != n ||
      static_cast<std::size_t>(time.size()) != n) {
    Rcpp::stop("ring, col, time, and scores have incompatible lengths.");
  }
  const arma::vec time_v(time.begin(), n, false, true);
  const arma::mat S_col(scores.begin(), n, k, false, true);
  return conley::grid_meat(ring.begin(), col.begin(), time_v, S_col,
                           lat0, dlat, dlon, n_ring, n_col, n_col_full,
                           cutoff, dist_fn, kernel, ncores);
}
