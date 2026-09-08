# Extracted from test-review-regressions.R:289

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fastconley", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("lfe")
skip_if_not_installed("fixest")
set.seed(43)
n <- 240
d <- data.frame(x = rnorm(n), g1 = rep(1:6, 40), g2 = rep(1:3, 80),
                  lat = runif(n, 30, 45), lon = runif(n, -110, -80))
d$y <- 0.5 * d$x + rnorm(n)
d$rowid <- seq_len(n)
fl <- lfe::felm(y ~ x | g1 + g2, d, keepCX = TRUE)
v_default <- vcovSpHAC(fl, lat = "lat", lon = "lon", dist_cutoff = 400, ncores = 1, data = d)
v_rowid <- vcovSpHAC(fl, unit = "rowid", lat = "lat", lon = "lon", dist_cutoff = 400, ncores = 1, data = d)
expect_identical(unname(v_default), unname(v_rowid))
fx <- fixest::feols(y ~ x | g1 + g2, d, demeaned = TRUE)
v_fx <- vcovSpHAC(fx, lat = "lat", lon = "lon", dist_cutoff = 400, ncores = 1, data = d)
expect_equal(unname(v_default), unname(v_fx), tolerance = 1e-10)
