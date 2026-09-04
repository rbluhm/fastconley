# R-side wrappers around the C++ entry points (src/cpp-functions.cpp; the
# generated bindings in R/RcppExports.R expose them as FastSpatialMeat_cpp
# etc.). The C++ aliases R memory directly (no input copies) and assumes
# REALSXP / INTSXP storage, so the wrappers coerce here. `ncores` is passed
# straight through: the engine runs its own std::thread pool
# (src/conley_core.h). `scores` is the score matrix e * X (possibly
# pre-aggregated); passing `X` and `e` instead is supported for convenience
# and computes `scores <- X * e` once.

.as_dbl_vector <- function(x) {
  if (is.double(x)) x else as.numeric(x)
}

.as_dbl_matrix <- function(m) {
  if (!is.matrix(m)) m <- as.matrix(m)
  if (!is.double(m)) storage.mode(m) <- "double"
  m
}

FastSpatialMeat <- function(lat, lon, time, X = NULL, e = NULL, cutoff,
                            kernel = "bartlett", dist_fn = "haversine",
                            balanced_pnl = FALSE, ncores = 1L,
                            neighbor = "grid", scores = NULL,
                            csr_weight = "double") {
  if (is.null(scores)) scores <- X * e
  .Call(`_fastconley_FastSpatialMeat_cpp`,
        .as_dbl_vector(lat), .as_dbl_vector(lon), .as_dbl_vector(time),
        .as_dbl_matrix(scores), cutoff, kernel, dist_fn, balanced_pnl,
        ncores, neighbor, csr_weight)
}

FastSerialHacPanel <- function(unit, time, cutoff, X = NULL, e = NULL,
                               ncores = 1L, scores = NULL) {
  if (is.null(scores)) scores <- X * e
  .Call(`_fastconley_FastSerialHacPanel_cpp`,
        .as_dbl_vector(unit), .as_dbl_vector(time), cutoff,
        .as_dbl_matrix(scores), ncores)
}

FastGridMeat <- function(ring, col, time, scores, lat0, dlat, dlon,
                         n_ring, n_col, cutoff, dist_fn = "spherical",
                         kernel = "uniform", n_col_full = 0L, ncores = 1L) {
  if (!is.integer(ring)) ring <- as.integer(ring)
  if (!is.integer(col)) col <- as.integer(col)
  .Call(`_fastconley_FastGridMeat_cpp`, ring, col, .as_dbl_vector(time),
        .as_dbl_matrix(scores), lat0, dlat, dlon,
        as.integer(n_ring), as.integer(n_col), as.integer(n_col_full),
        cutoff, dist_fn, kernel, ncores)
}
