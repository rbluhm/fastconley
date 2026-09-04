# Character and factor unit/time columns must be coerced to integer group IDs
# before the C++ entry points see them — the actual labels are immaterial; only
# equality matters. Regression test for the original bug where the fixest path
# errored on character columns and the felm path silently NA'd them.

test_that("character unit/time produces the same VCOV as integer keys (felm)", {
  skip_if_not_installed("lfe")

  n_unit <- 30L; n_time <- 4L
  d <- make_balanced_panel(n_unit = n_unit, n_time = n_time, k = 2L, seed = 71L)
  d$unit_chr <- paste0("U_", d$unit)
  d$time_chr <- paste0("T_", d$time)

  fit_chr <- lfe::felm(y ~ x1 + x2 | unit_chr + time_chr, data = d, keepCX = TRUE)
  fit_int <- lfe::felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)

  V_chr <- vcovSpHAC(fit_chr, unit = "unit_chr", time = "time_chr",
                     lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, lag_cutoff = 1L,
                     balanced_pnl = TRUE, ncores = 1L)
  V_int <- vcovSpHAC(fit_int, unit = "unit", time = "time",
                     lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, lag_cutoff = 1L,
                     balanced_pnl = TRUE, ncores = 1L)
  expect_lt(max(abs(V_chr - V_int)), 1e-12)
  expect_true(all(is.finite(V_chr)))
})

test_that("character unit/time produces the same VCOV as integer keys (fixest)", {
  skip_if_not_installed("fixest")

  n_unit <- 30L; n_time <- 4L
  d <- make_balanced_panel(n_unit = n_unit, n_time = n_time, k = 2L, seed = 71L)
  d$unit_chr <- paste0("U_", d$unit)
  d$time_chr <- paste0("T_", d$time)

  fit_chr <- fixest::feols(y ~ x1 + x2 | unit_chr + time_chr, data = d,
                           ssc = fixest::ssc(adj = FALSE, cluster.adj = FALSE),
                           demeaned = TRUE)
  fit_int <- fixest::feols(y ~ x1 + x2 | unit + time, data = d,
                           ssc = fixest::ssc(adj = FALSE, cluster.adj = FALSE),
                           demeaned = TRUE)

  V_chr <- vcovSpHAC(fit_chr, unit = "unit_chr", time = "time_chr",
                     lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, lag_cutoff = 1L,
                     balanced_pnl = TRUE, ncores = 1L)
  V_int <- vcovSpHAC(fit_int, unit = "unit", time = "time",
                     lat = "lat", lon = "lon",
                     kernel = "uniform", dist_fn = "spherical",
                     dist_cutoff = 500, lag_cutoff = 1L,
                     balanced_pnl = TRUE, ncores = 1L)
  expect_lt(max(abs(V_chr - V_int)), 1e-12)
  expect_true(all(is.finite(V_chr)))
})
