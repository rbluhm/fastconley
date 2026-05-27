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
