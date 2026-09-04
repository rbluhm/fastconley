test_that("fepois and feglm match an independent dense score sandwich", {
  skip_if_not_installed("fixest")
  set.seed(101)
  n <- 90L
  d <- data.frame(x = rnorm(n), g = rep(1:9, each = 10L),
                  lat = runif(n, 35, 43), lon = runif(n, -8, 4))
  d$count <- rpois(n, exp(0.2 + 0.35 * d$x + rep(rnorm(9, sd = 0.2), each = 10L)))
  d$binary <- rbinom(n, 1, plogis(-0.2 + 0.5 * d$x))
  fits <- list(
    fixest::fepois(count ~ x | g, data = d),
    fixest::feglm(binary ~ x, data = d, family = "binomial")
  )
  for (fit in fits) {
    V <- vcovSpHAC(fit, lat = "lat", lon = "lon", kernel = "bartlett",
                   dist_fn = "haversine", dist_cutoff = 450, method = "pairwise",
                   ncores = 1, ssc = FALSE, psd_fix = FALSE, data = d)
    idx <- fixest::obs(fit)
    V_ref <- review_dense_vcov(
      scores = as.matrix(fit$scores), bread = unname(fit$cov.iid),
      lat = d$lat[idx], lon = d$lon[idx], time = rep(1, length(idx)),
      cutoff = 450, kernel = "bartlett"
    )
    expect_equal(unname(V), unname(V_ref), tolerance = 2e-12)
  }
})

test_that("felm and fixest 2SLS match the same independent dense reference", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  set.seed(102)
  n <- 120L
  d <- data.frame(x = rnorm(n), z1 = rnorm(n), z2 = rnorm(n),
                  lat = runif(n, 28, 48), lon = runif(n, -12, 12))
  d$q <- 0.9 * d$z1 - 0.4 * d$z2 + 0.2 * d$x + rnorm(n)
  d$y <- 0.4 * d$x + 0.8 * d$q + rnorm(n)
  fit_lfe <- lfe::felm(y ~ x | 0 | (q ~ z1 + z2) | 0,
                       data = d, keepCX = TRUE)
  fit_fx <- fixest::feols(y ~ x | q ~ z1 + z2, data = d, demeaned = TRUE)
  common <- list(lat = "lat", lon = "lon", kernel = "uniform",
                 dist_fn = "haversine", dist_cutoff = 600,
                 method = "pairwise", ncores = 1, ssc = FALSE,
                 psd_fix = FALSE, data = d)
  V_lfe <- do.call(vcovSpHAC, c(list(reg = fit_lfe), common))
  V_fx <- do.call(vcovSpHAC, c(list(reg = fit_fx), common))
  scores <- fit_fx$X_demeaned * as.numeric(fit_fx$residuals)
  bread <- solve(crossprod(fit_fx$X_demeaned))
  V_ref <- review_dense_vcov(scores, bread, d$lat, d$lon, rep(1, n),
                             cutoff = 600, kernel = "uniform")
  expect_equal(unname(V_fx), unname(V_ref), tolerance = 2e-12)
  lfe_order_in_fixest <- c(1L, 3L, 2L)
  expect_equal(unname(V_lfe),
               unname(V_ref[lfe_order_in_fixest, lfe_order_in_fixest]),
               tolerance = 2e-12)
})

test_that("serial HAC meat matches an independent naive loop", {
  set.seed(103)
  unit <- rep(seq_len(140L), each = 4L)
  time <- rep(c(1, 2, 4, 7), times = 140L)
  scores <- matrix(rnorm(length(unit) * 3L), ncol = 3L)
  got <- fastconley:::FastSerialHacPanel(unit, time, cutoff = 3,
                                        scores = scores, ncores = 2)
  ref <- review_naive_serial_meat(scores, unit, time, cutoff = 3)
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("pixel aggregation matches manual snapping and dense aggregation", {
  skip_if_not_installed("fixest")
  set.seed(104)
  n <- 100L
  d <- data.frame(x = rnorm(n), lat = runif(n, 39, 42),
                  lon = runif(n, -3, 3))
  d$y <- 0.5 * d$x + rnorm(n)
  fit <- fixest::feols(y ~ x, data = d, demeaned = TRUE)
  pixel <- 20
  V <- vcovSpHAC(fit, lat = "lat", lon = "lon", kernel = "bartlett",
                 dist_fn = "haversine", dist_cutoff = 500, pixel = pixel,
                 method = "pairwise", ncores = 1, ssc = FALSE,
                 psd_fix = FALSE, data = d)

  lat_step <- pixel / 111
  lat_agg <- round(d$lat / lat_step) * lat_step
  cos_lat <- cos(lat_agg * pi / 180)
  cos_lat[abs(cos_lat) < 1e-6] <- 1e-6
  lon_step <- pixel / (111 * cos_lat)
  lon_agg <- round(d$lon / lon_step) * lon_step
  scores <- fit$X_demeaned * as.numeric(fit$residuals)
  key <- interaction(lat_agg, lon_agg, drop = TRUE)
  first <- !duplicated(key)
  scores_agg <- rowsum(scores, key, reorder = FALSE)
  V_ref <- review_dense_vcov(
    scores_agg, solve(crossprod(fit$X_demeaned)),
    lat_agg[first], lon_agg[first], rep(1, sum(first)),
    cutoff = 500, kernel = "bartlett"
  )
  expect_equal(unname(V), unname(V_ref), tolerance = 2e-12)
})
