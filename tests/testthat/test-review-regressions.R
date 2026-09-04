test_that("felm resolves requested unit/time names and data-column fallbacks", {
  skip_if_not_installed("lfe")
  d <- make_balanced_panel(n_unit = 24L, n_time = 3L, k = 2L, seed = 91L)
  d$region <- d$unit %% 4L
  d$unit_raw <- d$unit
  d$time_raw <- d$time
  fit <- lfe::felm(y ~ x1 + x2 | region + unit + time,
                   data = d, keepCX = TRUE)
  call_vcov <- function(unit, time) {
    vcovSpHAC(fit, unit = unit, time = time, lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "haversine", dist_cutoff = 600,
              lag_cutoff = 1, balanced_pnl = TRUE, ncores = 1,
              ssc = FALSE, psd_fix = FALSE, data = d)
  }

  V_named <- call_vcov("unit", "time")
  V_raw <- call_vcov("unit_raw", "time_raw")
  expect_identical(V_named, V_raw)
  expect_error(call_vcov("missing_unit", "time"), "missing_unit")
})

test_that("felm uses the first two absorbed FEs only when neither name is passed", {
  skip_if_not_installed("lfe")
  d <- make_balanced_panel(n_unit = 20L, n_time = 3L, k = 1L, seed = 92L)
  d$region <- d$unit %% 4L
  fit <- lfe::felm(y ~ x1 | unit + time + region, data = d, keepCX = TRUE)
  args <- list(lat = "lat", lon = "lon", kernel = "uniform",
               dist_fn = "haversine", dist_cutoff = 500, lag_cutoff = 1,
               balanced_pnl = TRUE, ncores = 1, ssc = FALSE,
               psd_fix = FALSE, data = d)
  V_default <- do.call(vcovSpHAC, c(list(reg = fit), args))
  V_named <- do.call(vcovSpHAC,
                     c(list(reg = fit, unit = "unit", time = "time"), args))
  expect_identical(V_default, V_named)
})

test_that("time without unit defines spatial blocks and lag requires unit", {
  skip_if_not_installed("fixest")
  d <- as.data.frame(make_cross_section(n = 90L, k = 2L, seed = 93L))
  d$period <- rep(1:3, length.out = nrow(d))
  d$row_id <- seq_len(nrow(d))
  fit <- fixest::feols(y ~ x1 + x2, data = d, demeaned = TRUE)
  args <- list(reg = fit, lat = "lat", lon = "lon", kernel = "uniform",
               dist_fn = "haversine", dist_cutoff = 20000, ncores = 1,
               ssc = FALSE, psd_fix = FALSE, data = d)
  V_time <- do.call(vcovSpHAC, c(args, list(time = "period")))
  V_explicit <- do.call(vcovSpHAC,
                        c(args, list(unit = "row_id", time = "period")))
  V_ignored <- do.call(vcovSpHAC, args)
  expect_identical(V_time, V_explicit)
  expect_gt(max(abs(V_time - V_ignored)), 1e-10)
  expect_error(do.call(vcovSpHAC, c(args, list(time = "period", lag_cutoff = 1))),
               "lag_cutoff requires unit")
})

test_that("felm also honours time without unit", {
  skip_if_not_installed("lfe")
  d <- as.data.frame(make_cross_section(n = 80L, k = 1L, seed = 94L))
  d$period <- rep(1:2, each = nrow(d) / 2L)
  d$row_id <- seq_len(nrow(d))
  fit <- lfe::felm(y ~ x1, data = d, keepCX = TRUE)
  args <- list(reg = fit, lat = "lat", lon = "lon", kernel = "uniform",
               dist_fn = "haversine", dist_cutoff = 20000, ncores = 1,
               ssc = FALSE, psd_fix = FALSE, data = d)
  expect_identical(
    do.call(vcovSpHAC, c(args, list(time = "period"))),
    do.call(vcovSpHAC, c(args, list(unit = "row_id", time = "period")))
  )
})

test_that("fixest::obs aligns split, NA, and subset selections", {
  skip_if_not_installed("fixest")
  d <- as.data.frame(make_cross_section(n = 120L, k = 2L, seed = 95L))
  d$group <- rep(c("a", "b"), each = 60L)
  d$y[c(4L, 67L)] <- NA_real_
  multi <- fixest::feols(y ~ x1 + x2, data = d, split = ~group,
                         demeaned = TRUE)
  one <- multi[[1L]]
  idx <- fixest::obs(one)
  pre <- fixest::feols(y ~ x1 + x2, data = d[idx, ], demeaned = TRUE)
  common <- list(lat = "lat", lon = "lon", kernel = "bartlett",
                 dist_fn = "haversine", dist_cutoff = 500, ncores = 1,
                 ssc = FALSE, psd_fix = FALSE)
  expect_identical(
    do.call(vcovSpHAC, c(list(reg = one, data = d), common)),
    do.call(vcovSpHAC, c(list(reg = pre, data = d[idx, ]), common))
  )

  keep <- seq_len(nrow(d)) %% 3L != 0L
  sub <- fixest::feols(y ~ x1 + x2, data = d, subset = keep, demeaned = TRUE)
  idx_sub <- fixest::obs(sub)
  pre_sub <- fixest::feols(y ~ x1 + x2, data = d[idx_sub, ], demeaned = TRUE)
  expect_identical(
    do.call(vcovSpHAC, c(list(reg = sub, data = d), common)),
    do.call(vcovSpHAC, c(list(reg = pre_sub, data = d[idx_sub, ]), common))
  )
})

test_that("unsupported fixest object shapes fail clearly", {
  skip_if_not_installed("fixest")
  d <- as.data.frame(make_cross_section(n = 60L, k = 1L, seed = 96L))
  d$group <- rep(1:2, each = 30L)
  d$g <- rep(1:6, each = 10L)
  multi <- fixest::feols(y ~ x1, data = d, split = ~group, demeaned = TRUE)
  lean <- fixest::feols(y ~ x1, data = d, demeaned = TRUE, lean = TRUE)
  only_fe <- fixest::feols(y ~ 1 | g, data = d, demeaned = TRUE)
  call_vcov <- function(x) vcovSpHAC(x, lat = "lat", lon = "lon",
                                     dist_cutoff = 500, data = d)
  expect_error(call_vcov(multi), "fixest_multi")
  expect_error(call_vcov(lean), "lean = TRUE")
  expect_error(call_vcov(only_fe), "only-fixed-effect")
})

test_that("felm rejects non-2SLS k-class estimates", {
  skip_if_not_installed("lfe")
  set.seed(97)
  n <- 140L
  d <- data.frame(x = rnorm(n), z1 = rnorm(n), z2 = rnorm(n),
                  lat = runif(n, 30, 45), lon = runif(n, -10, 5))
  d$q <- 0.7 * d$z1 + 0.3 * d$z2 + rnorm(n)
  d$y <- 0.4 * d$x + 0.8 * d$q + rnorm(n)
  fit <- lfe::felm(y ~ x | 0 | (q ~ z1 + z2) | 0, data = d,
                   keepCX = TRUE, kclass = "liml")
  expect_error(vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 500,
                         data = d), "k-class")
  fit_exact <- lfe::felm(y ~ x | 0 | (q ~ z1) | 0, data = d,
                         keepCX = TRUE, kclass = "liml")
  expect_error(vcovSpHAC(fit_exact, lat = "lat", lon = "lon",
                         dist_cutoff = 500, data = d),
               "k-class")
})

test_that("ncores options, CRAN cap, rounding, and validation are explicit", {
  old_opt <- options(fastconley.ncores = 7)
  old_env <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit({
    options(old_opt)
    if (is.na(old_env)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_env)
  }, add = TRUE)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  va <- function(ncores = NULL) fastconley:::validate_args(
    "uniform", "haversine", 100, 0, 0, ncores,
    "grid", "double", "pairwise", FALSE, FALSE)
  expect_identical(va()$ncores, 7L)
  expect_identical(va(1.6)$ncores, 2L)
  options(fastconley.ncores = NA_real_)
  expect_error(va(), "ncores")
  options(fastconley.ncores = 7)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_identical(va()$ncores, 2L)
  expect_identical(va(12)$ncores, 2L)
  for (bad in list(NA_real_, Inf, 0, -1, c(1, 2), "2")) {
    expect_error(va(bad), "ncores")
  }
})

test_that("argument typos and deprecated maxobsmem are visible", {
  skip_if_not_installed("lfe")
  d <- as.data.frame(make_cross_section(n = 50L, k = 1L, seed = 98L))
  fit <- lfe::felm(y ~ x1, data = d, keepCX = TRUE)
  expect_error(vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 500,
                         psd_fxi = FALSE, data = d), "psd_fxi")
  expect_warning(vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 500,
                           maxobsmem = 10, data = d),
                 "maxobsmem is deprecated and ignored")
})

test_that("coordinate auto-detection rejects aliases across all candidates", {
  expect_error(fastconley:::detect_coord_names(
    c("lat", "latitude", "lon"), NULL, NULL), "ambiguous")
  expect_error(fastconley:::detect_coord_names(
    c("lat", "lon", "longitude"), NULL, NULL), "ambiguous")
})

test_that("panel keys and coordinates are validated before dispatch", {
  skip_if_not_installed("fixest")
  d <- as.data.frame(make_cross_section(n = 50L, k = 1L, seed = 99L))
  d$unit <- seq_len(nrow(d))
  d$time <- rep(1:2, length.out = nrow(d))
  fit <- fixest::feols(y ~ x1, data = d, demeaned = TRUE)
  call_bad <- function(dd) vcovSpHAC(
    fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
    dist_cutoff = 500, ncores = 1, data = dd)

  cases <- list(
    list(col = "unit", value = NA_real_, msg = "unit contains missing"),
    list(col = "unit", value = Inf, msg = "unit contains non-finite"),
    list(col = "time", value = NA_real_, msg = "time contains missing"),
    list(col = "time", value = Inf, msg = "time contains non-finite"),
    list(col = "lat", value = NA_real_, msg = "latitude contains missing"),
    list(col = "lat", value = Inf, msg = "latitude contains non-finite"),
    list(col = "lat", value = 91, msg = "latitude must lie"),
    list(col = "lon", value = NA_real_, msg = "longitude contains missing"),
    list(col = "lon", value = Inf, msg = "longitude contains non-finite"),
    list(col = "lon", value = 361, msg = "longitude must lie"),
    list(col = "lon", value = -181, msg = "longitude must lie")
  )
  for (case in cases) {
    dd <- d
    dd[[case$col]][1L] <- case$value
    expect_error(call_bad(dd), case$msg, label = paste(case$col, case$value))
  }
})

test_that("felm fits through lfe's k-class path are rejected, default 2SLS accepted", {
  skip_if_not_installed("lfe")
  set.seed(31)
  n <- 150
  d <- data.frame(g = rep(1:5, 30), lat = runif(n, 30, 45), lon = runif(n, -110, -80))
  d$w <- rnorm(n); d$z <- 0.8 * d$w + rnorm(n); d$y <- d$z + rnorm(n)
  f_def <- lfe::felm(y ~ 1 | g | (z ~ w) | 0, d, keepCX = TRUE)
  expect_silent(v_def <- vcovSpHAC(f_def, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1))
  # kclass = 1 reproduces the 2SLS point estimate but lfe stores the raw
  # endogenous regressor in cX, so the object is not usable here.
  kval <- 1
  f_sym <- lfe::felm(y ~ 1 | g | (z ~ w) | 0, d, kclass = kval, keepCX = TRUE)
  expect_equal(unname(coef(f_sym)), unname(coef(f_def)))
  expect_error(vcovSpHAC(f_sym, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1), "k-class")
  f_liml <- lfe::felm(y ~ 1 | g | (z ~ w) | 0, d, kclass = "liml", keepCX = TRUE)
  expect_error(vcovSpHAC(f_liml, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1), "k-class")
})

test_that("fixest path accepts a data frame already aligned with the fit", {
  skip_if_not_installed("fixest")
  set.seed(32)
  n <- 160
  d <- data.frame(x = rnorm(n), g = rep(1:4, 40), lat = runif(n, 30, 45), lon = runif(n, -110, -80))
  d$y <- 0.5 * d$x + rnorm(n)
  d$x[c(5, 17)] <- NA
  fit <- fixest::feols(y ~ x | g, d, subset = ~ lat > 32, demeaned = TRUE)
  v_full <- vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = d)
  d_aligned <- d[fixest::obs(fit), ]
  v_aligned <- vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = d_aligned)
  expect_identical(v_full, v_aligned)
  expect_error(vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = d[1:100, ]),
               "original rows")
  # A subset that only permutes the rows keeps n == nobs_origin: the n-row
  # frame must then be read as ORIGINAL data (rows addressed via obs()), so the
  # reversed fit equals the forward fit's covariance (same rows, same weights).
  d2 <- d[!is.na(d$x), ]
  f_fwd <- fixest::feols(y ~ x | g, d2, demeaned = TRUE)
  f_rev <- fixest::feols(y ~ x | g, d2, subset = nrow(d2):1, demeaned = TRUE)
  v_fwd <- vcovSpHAC(f_fwd, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = d2)
  v_rev <- vcovSpHAC(f_rev, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = d2)
  expect_equal(v_fwd, v_rev, tolerance = 1e-12)
})
