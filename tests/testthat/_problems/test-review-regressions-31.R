# Extracted from test-review-regressions.R:31

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fastconley", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("lfe")
d <- make_balanced_panel(n_unit = 20L, n_time = 3L, k = 1L, seed = 92L)
d$region <- d$unit %% 4L
fit <- lfe::felm(y ~ x1 | unit + time + region, data = d, keepCX = TRUE)
args <- list(lat = "lat", lon = "lon", kernel = "uniform",
               dist_fn = "haversine", dist_cutoff = 500, lag_cutoff = 1,
               balanced_pnl = TRUE, ncores = 1, ssc = FALSE,
               psd_fix = FALSE, data = d)
V_default <- do.call(vcovSpHAC, c(list(reg = fit), args))
