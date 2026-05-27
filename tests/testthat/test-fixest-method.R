# Tests for vcovSpHAC.fixest. Mirrors the highlights of
# tests/manual/test-fixest-method.R; the manual script keeps the larger panel
# workloads for offline runs.

test_that("vcovSpHAC.fixest matches vcovSpHAC.felm to machine precision", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  RcppParallel::setThreadOptions(numThreads = 1L)

  d <- make_balanced_panel(n_unit = 150L, n_time = 4L, k = 3L, seed = 21L)
  rhs <- "x1 + x2 + x3"

  fit_lfe <- lfe::felm(stats::as.formula(paste("y ~", rhs, "| unit + time")),
                       data = d, keepCX = TRUE)
  fit_fx  <- fixest::feols(stats::as.formula(paste("y ~", rhs, "| unit + time")),
                           data = d,
                           ssc = fixest::ssc(adj = FALSE, cluster.adj = FALSE),
                           demeaned = TRUE)

  for (kn in c("uniform", "bartlett")) {
    for (df in c("haversine", "spherical", "chord")) {
      V_lfe <- vcovSpHAC(fit_lfe, unit = "unit", time = "time",
                         lat = "lat", lon = "lon",
                         kernel = kn, dist_fn = df,
                         dist_cutoff = 500, lag_cutoff = 1,
                         balanced_pnl = TRUE, ncores = 1L)
      V_fx  <- vcovSpHAC(fit_fx, unit = "unit", time = "time",
                         lat = "lat", lon = "lon",
                         kernel = kn, dist_fn = df,
                         dist_cutoff = 500, lag_cutoff = 1,
                         balanced_pnl = TRUE, ncores = 1L)
      expect_lt(max(abs(V_lfe - V_fx)), 1e-12,
                label = paste("felm vs fixest", kn, df))
    }
  }
})

test_that("vcovSpHAC.fixest honors subset= in feols", {
  skip_if_not_installed("fixest")
  RcppParallel::setThreadOptions(numThreads = 1L)

  d <- make_cross_section(n = 200L, k = 3L, seed = 51L)
  d$keep <- rep(c(TRUE, FALSE), nrow(d) / 2L)

  fit_sub <- fixest::feols(y ~ x1 + x2 + x3, data = d, subset = d$keep,
                           demeaned = TRUE)
  fit_pre <- fixest::feols(y ~ x1 + x2 + x3, data = d[d$keep],
                           demeaned = TRUE)
  V_sub <- vcovSpHAC(fit_sub, lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, ncores = 1L)
  V_pre <- vcovSpHAC(fit_pre, lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, ncores = 1L)
  expect_lt(max(abs(V_sub - V_pre)), 1e-12)
})

test_that("vcovSpHAC.fixest errors without demeaned=TRUE", {
  skip_if_not_installed("fixest")
  d <- make_cross_section(n = 100L, k = 2L, seed = 61L)
  fit <- fixest::feols(y ~ x1 + x2, data = d)
  expect_error(
    vcovSpHAC(fit, lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 500),
    regexp = "demeaned = TRUE"
  )
})

test_that("vcovSpHAC default method errors on unsupported classes", {
  fit <- stats::lm(mpg ~ wt, data = mtcars)
  expect_error(vcovSpHAC(fit, lat = "x", lon = "y", dist_cutoff = 100),
               regexp = "unsupported model class")
})
