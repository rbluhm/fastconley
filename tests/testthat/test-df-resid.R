test_that("df_resid overrides the ssc denominator and is validated", {
  skip_if_not_installed("lfe")
  set.seed(3); n <- 400
  d <- data.frame(lat = runif(n, -1, 1), lon = runif(n, 36, 38),
                  g = factor(sample(5, n, TRUE)))
  d$x <- rnorm(n); d$y <- 0.5 * d$x + rnorm(n)
  fit <- lfe::felm(y ~ x | g, data = d, keepCX = TRUE)
  base <- function(...) {
    vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 50, data = d,
              psd_fix = FALSE, ...)
  }
  v0 <- base(ssc = FALSE)
  v1 <- base()
  expect_identical(base(df_resid = NULL), v1)
  expect_equal(v1, v0 * n / (n - fit$p), tolerance = 1e-14)
  expect_equal(base(df_resid = 300), v0 * n / 300, tolerance = 1e-14)
  expect_equal(base(df_resid = n), v0, tolerance = 1e-14)
  expect_error(base(ssc = FALSE, df_resid = 300), "only applies with ssc = TRUE")
  for (bad in list(0, -1, n + 1, NA_real_, Inf, c(100, 200), "300", TRUE)) {
    expect_error(base(df_resid = bad), "df_resid must be a single number")
  }
})

test_that("fixest method honours df_resid the same way", {
  skip_if_not_installed("fixest")
  set.seed(4); n <- 300
  d <- data.frame(lat = runif(n, -1, 1), lon = runif(n, 36, 38),
                  g = factor(sample(5, n, TRUE)))
  d$x <- rnorm(n); d$y <- 0.5 * d$x + rnorm(n)
  fit <- fixest::feols(y ~ x | g, data = d, demeaned = TRUE)
  base <- function(...) {
    vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 50, data = d,
              psd_fix = FALSE, ...)
  }
  v0 <- base(ssc = FALSE)
  expect_equal(base(), v0 * n / (n - fit$nparams), tolerance = 1e-14)
  expect_equal(base(df_resid = 250), v0 * n / 250, tolerance = 1e-14)
  expect_error(base(ssc = FALSE, df_resid = 250), "only applies with ssc = TRUE")
  expect_error(base(df_resid = n + 1), "df_resid must be a single number")
})

test_that("felm and fixest agree under a common df_resid on a disconnected FE graph", {
  skip_if_not_installed("lfe")
  skip_if_not_installed("fixest")
  # Restricting the sample so that ethnic groups 11-12 only occur with ages
  # 9-10 splits the (age, ethnic) level graph into two connected components:
  # lfe counts the exact rank of the first two factors, fixest's default does
  # not, so the two default-ssc results differ by one parameter. The same
  # df_resid removes the difference.
  set.seed(11); n <- 2400
  d <- data.frame(age = sample(10, n, TRUE), ethnic = sample(12, n, TRUE),
                  province = sample(8, n, TRUE))
  d$x <- rnorm(n); d$y <- 0.5 * d$x + rnorm(n) + d$age / 10 + d$ethnic / 10
  d$lat <- runif(n, -2, 2); d$lon <- runif(n, 35, 40)
  sub <- d[!((d$ethnic %in% 11:12) != (d$age %in% 9:10)), ]
  ns <- nrow(sub)
  op <- options(lfe.eps = 1e-11); on.exit(options(op), add = TRUE)
  fl <- lfe::felm(y ~ x | age + ethnic + province, data = sub, keepCX = TRUE)
  fx <- fixest::feols(y ~ x | age + ethnic + province, data = sub,
                      fixef.rm = "none", demeaned = TRUE, fixef.tol = 1e-11)
  K <- 1L + qr(model.matrix(~ factor(age) + factor(ethnic) + factor(province),
                            data = sub))$rank
  df <- ns - K
  cl <- function(fit, ...) {
    vcovSpHAC(fit, lat = "lat", lon = "lon", dist_cutoff = 100, data = sub,
              psd_fix = FALSE, ...)
  }
  vl <- cl(fl, df_resid = df); vx <- cl(fx, df_resid = df)
  expect_equal(unname(vl), unname(vx), tolerance = 1e-6)
  # With the packages' own counts the ratio is exactly the count ratio.
  vl0 <- cl(fl); vx0 <- cl(fx)
  expect_equal(vl0[1, 1] / vx0[1, 1],
               (ns - fx$nparams) / (ns - fl$p) * vl[1, 1] / vx[1, 1],
               tolerance = 1e-12)
})
