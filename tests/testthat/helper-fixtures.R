# Fixture builders shared across testthat suites. Kept small so the full
# R CMD check test step stays well under a few seconds on a modest machine.

suppressMessages(library(data.table))

make_cross_section <- function(n = 80L, k = 3L, seed = 1L) {
  set.seed(seed)
  X <- cbind(`(Intercept)` = 1,
             matrix(stats::rnorm(n * k), n, k,
                    dimnames = list(NULL, paste0("x", seq_len(k)))))
  data.table(
    y = as.numeric(X[, -1, drop = FALSE] %*% rep(0.1, k)) + stats::rnorm(n),
    X[, -1, drop = FALSE],
    lat = stats::runif(n, 25, 50),
    lon = stats::runif(n, -125, -70)
  )
}

make_balanced_panel <- function(n_unit = 25L, n_time = 4L, k = 3L, seed = 2L,
                                unit_labels = NULL, time_labels = NULL) {
  set.seed(seed)
  ulat <- stats::runif(n_unit, 25, 50); ulon <- stats::runif(n_unit, -125, -70)
  unit_idx <- rep(seq_len(n_unit), times = n_time)
  time_idx <- rep(seq_len(n_time), each = n_unit)
  N <- n_unit * n_time
  X <- matrix(stats::rnorm(N * k), N, k,
              dimnames = list(NULL, paste0("x", seq_len(k))))
  data.table(
    y = as.numeric(X %*% rep(0.1, k)) + stats::rnorm(N),
    X,
    lat = ulat[unit_idx],
    lon = ulon[unit_idx],
    unit = if (is.null(unit_labels)) unit_idx else unit_labels[unit_idx],
    time = if (is.null(time_labels)) time_idx else time_labels[time_idx]
  )
}

review_dense_spatial_meat <- function(scores, lat, lon, time, cutoff,
                                      kernel = c("bartlett", "uniform")) {
  kernel <- match.arg(kernel)
  lat_r <- lat * pi / 180
  lon_r <- lon * pi / 180
  dlat <- outer(lat_r, lat_r, "-")
  dlon <- outer(lon_r, lon_r, "-")
  a <- sin(dlat / 2)^2 + outer(cos(lat_r), cos(lat_r), "*") * sin(dlon / 2)^2
  a <- pmin(pmax(a, 0), 1)
  distance <- 2 * 6371 * asin(sqrt(a))
  diag(distance) <- 0
  weight <- (outer(time, time, "==") & distance <= cutoff) *
    if (kernel == "bartlett") pmax(0, 1 - distance / cutoff) else 1
  crossprod(scores, weight %*% scores)
}

review_naive_serial_meat <- function(scores, unit, time, cutoff) {
  k <- ncol(scores)
  meat <- matrix(0, k, k)
  if (nrow(scores) < 2L) return(meat)
  for (i in seq_len(nrow(scores) - 1L)) {
    for (j in seq.int(i + 1L, nrow(scores))) {
      gap <- abs(time[j] - time[i])
      if (unit[i] != unit[j] || gap <= 0 || gap > cutoff) next
      weight <- 1 - gap / (cutoff + 1)
      meat <- meat + weight * (tcrossprod(scores[i, ], scores[j, ]) +
                                 tcrossprod(scores[j, ], scores[i, ]))
    }
  }
  meat
}

review_dense_vcov <- function(scores, bread, lat, lon, time, cutoff,
                              kernel = "bartlett", unit = NULL,
                              lag_cutoff = 0) {
  meat <- review_dense_spatial_meat(scores, lat, lon, time, cutoff, kernel)
  if (lag_cutoff > 0) {
    meat <- meat + review_naive_serial_meat(scores, unit, time, lag_cutoff)
  }
  out <- bread %*% meat %*% bread
  (out + t(out)) / 2
}
