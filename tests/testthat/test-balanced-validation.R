# balanced_pnl = TRUE requires the C++ side to reuse the first period's CSR
# neighbor graph for every subsequent period. That is only safe when each
# period contains exactly the same units with no duplicates and no missing.
# These tests pin the validation that catches inputs that violate that.

test_that("balanced_pnl errors on mismatched unit sets across periods", {
  skip_if_not_installed("lfe")
  set.seed(81)
  # Units {1,2,3,4} in period 1, {2,3,4,5} in period 2 — equal sizes, no
  # repeats, time-invariant coords. Should still error.
  ulat <- stats::runif(5, 25, 50); ulon <- stats::runif(5, -125, -70)
  unit <- c(1, 2, 3, 4, 2, 3, 4, 5)
  time <- c(1, 1, 1, 1, 2, 2, 2, 2)
  N <- length(unit); k <- 2L
  X <- matrix(stats::rnorm(N * k), N, k, dimnames = list(NULL, c("x1", "x2")))
  d <- data.table::data.table(
    y = as.numeric(X %*% c(0.1, 0.1)) + stats::rnorm(N),
    X, lat = ulat[unit], lon = ulon[unit], unit = unit, time = time)
  # 8 rows with two-way FEs leave lfe almost no residual df; felm warns about
  # NaNs in its own summary statistics, which is irrelevant to this test.
  fit <- suppressWarnings(
    lfe::felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE))
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 10000, balanced_pnl = TRUE, ncores = 1L),
    regexp = "same set of units"
  )
})

test_that("balanced_pnl errors on unequal period sizes", {
  skip_if_not_installed("lfe")
  set.seed(82)
  ulat <- stats::runif(4, 25, 50); ulon <- stats::runif(4, -125, -70)
  unit <- c(1, 2, 3, 4, 1, 2, 3)
  time <- c(1, 1, 1, 1, 2, 2, 2)
  N <- length(unit); k <- 2L
  X <- matrix(stats::rnorm(N * k), N, k, dimnames = list(NULL, c("x1", "x2")))
  d <- data.table::data.table(
    y = as.numeric(X %*% c(0.1, 0.1)) + stats::rnorm(N),
    X, lat = ulat[unit], lon = ulon[unit], unit = unit, time = time)
  fit <- suppressWarnings(
    lfe::felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)
  )
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 10000, balanced_pnl = TRUE, ncores = 1L),
    regexp = "same number of observations"
  )
})

test_that("balanced_pnl errors on duplicated unit within a period", {
  skip_if_not_installed("lfe")
  set.seed(83)
  ulat <- stats::runif(3, 25, 50); ulon <- stats::runif(3, -125, -70)
  unit <- c(1, 1, 2, 3, 1, 2, 2, 3)
  time <- c(1, 1, 1, 1, 2, 2, 2, 2)
  N <- length(unit); k <- 2L
  X <- matrix(stats::rnorm(N * k), N, k, dimnames = list(NULL, c("x1", "x2")))
  d <- data.table::data.table(
    y = as.numeric(X %*% c(0.1, 0.1)) + stats::rnorm(N),
    X, lat = ulat[unit], lon = ulon[unit], unit = unit, time = time)
  fit <- lfe::felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 10000, balanced_pnl = TRUE, ncores = 1L),
    regexp = "at most once per period"
  )
})

test_that("balanced_pnl errors on time-varying coordinates", {
  skip_if_not_installed("lfe")
  d <- make_balanced_panel(n_unit = 20L, n_time = 3L, k = 2L, seed = 84L)
  d[unit == 1L & time == 2L, lat := lat + 1]  # break the invariant
  fit <- lfe::felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 500, balanced_pnl = TRUE, ncores = 1L),
    regexp = "time-invariant coordinates"
  )
})


test_that("grid lattice detection rejects degenerate implied grids", {
  lat <- c(0, 1e-12, 1)
  lon <- c(0, 1e-12, 1)
  expect_null(fastconley:::detect_lonlat_grid(lat, lon))
})
