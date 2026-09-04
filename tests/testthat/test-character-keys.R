# Unit labels are equality keys. Character/factor time labels preserve their
# numeric scale when parsable, because serial-HAC weights depend on time gaps.

test_that("character unit/time produces the same VCOV as integer keys (felm)", {
  skip_if_not_installed("lfe")

  n_unit <- 30L; n_time <- 2L
  d <- make_balanced_panel(n_unit = n_unit, n_time = n_time, k = 2L, seed = 71L)
  d$unit_chr <- paste0("U_", d$unit)
  d$time <- c("2000", "2002")[d$time]
  d$time_chr <- as.character(d$time)
  d$time_fac <- factor(d$time)
  d$time_num <- as.numeric(d$time)

  fit_num <- lfe::felm(y ~ x1 + x2 | unit_chr + time_num, data = d, keepCX = TRUE)
  V_num <- vcovSpHAC(fit_num, unit = "unit_chr", time = "time_num",
                     lat = "lat", lon = "lon", kernel = "uniform",
                     dist_fn = "spherical", dist_cutoff = 500,
                     lag_cutoff = 1L, balanced_pnl = TRUE, ncores = 1L)
  for (time_name in c("time_chr", "time_fac")) {
    fit <- lfe::felm(stats::as.formula(
      paste("y ~ x1 + x2 | unit_chr +", time_name)), data = d, keepCX = TRUE)
    V <- vcovSpHAC(fit, unit = "unit_chr", time = time_name,
                   lat = "lat", lon = "lon", kernel = "uniform",
                   dist_fn = "spherical", dist_cutoff = 500,
                   lag_cutoff = 1L, balanced_pnl = TRUE, ncores = 1L)
    expect_identical(V, V_num)
  }
})

test_that("character unit/time produces the same VCOV as integer keys (fixest)", {
  skip_if_not_installed("fixest")

  n_unit <- 30L; n_time <- 2L
  d <- make_balanced_panel(n_unit = n_unit, n_time = n_time, k = 2L, seed = 71L)
  d$unit_chr <- paste0("U_", d$unit)
  d$time <- c("2000", "2002")[d$time]
  d$time_chr <- as.character(d$time)
  d$time_fac <- factor(d$time)
  d$time_num <- as.numeric(d$time)

  fit_num <- fixest::feols(y ~ x1 + x2 | unit_chr + time_num, data = d,
                           demeaned = TRUE)
  V_num <- vcovSpHAC(fit_num, unit = "unit_chr", time = "time_num",
                     lat = "lat", lon = "lon", kernel = "uniform",
                     dist_fn = "spherical", dist_cutoff = 500,
                     lag_cutoff = 1L, balanced_pnl = TRUE, ncores = 1L)
  for (time_name in c("time_chr", "time_fac")) {
    fit <- fixest::feols(stats::as.formula(
      paste("y ~ x1 + x2 | unit_chr +", time_name)), data = d,
      demeaned = TRUE)
    V <- vcovSpHAC(fit, unit = "unit_chr", time = time_name,
                   lat = "lat", lon = "lon", kernel = "uniform",
                   dist_fn = "spherical", dist_cutoff = 500,
                   lag_cutoff = 1L, balanced_pnl = TRUE, ncores = 1L)
    expect_identical(V, V_num)
  }
})

test_that("non-parsable time is rejected for serial HAC", {
  skip_if_not_installed("fixest")
  d <- make_balanced_panel(n_unit = 20L, n_time = 2L, k = 1L, seed = 72L)
  d$period <- factor(c("before", "after")[d$time])
  fit <- fixest::feols(y ~ x1 | unit + period, data = d, demeaned = TRUE)
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "period", lat = "lat", lon = "lon",
              dist_cutoff = 500, lag_cutoff = 1, ncores = 1),
    "numeric or numeric-parsable"
  )
})

test_that("felm also rejects non-parsable time for serial HAC", {
  skip_if_not_installed("lfe")
  d <- make_balanced_panel(n_unit = 20L, n_time = 2L, k = 1L, seed = 73L)
  d$period <- factor(c("before", "after")[d$time])
  fit <- lfe::felm(y ~ x1 | unit + period, data = d, keepCX = TRUE)
  expect_error(
    vcovSpHAC(fit, unit = "unit", time = "period", lat = "lat", lon = "lon",
              dist_cutoff = 500, lag_cutoff = 1, ncores = 1),
    "numeric or numeric-parsable"
  )
})
