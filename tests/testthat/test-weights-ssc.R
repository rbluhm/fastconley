# Weighted (WLS) fits, small-sample correction, psd_fix, lat/lon auto-detect.

make_cs_data <- function(n = 400, seed = 7) {
  set.seed(seed)
  d <- data.frame(latitude  = runif(n, 40, 45),
                  longitude = runif(n, -5, 5),
                  x = rnorm(n),
                  g = sample(letters[1:8], n, TRUE),
                  w = runif(n, 0.5, 2))
  d$y <- 0.5 * d$x + rnorm(n)
  d
}

# Dense brute-force weighted Conley reference. The diagonal of D is set to
# exactly 0: acos of a floating-point self-dot otherwise yields ~9e-5 km of
# phantom distance, which taints the bartlett weights.
dense_ref <- function(X, e, w, lat, lon, cutoff, kernel) {
  n <- length(lat)
  la <- lat * pi / 180; lo <- lon * pi / 180
  D <- outer(seq_len(n), seq_len(n), function(i, j) {
    6371 * acos(pmin(1, sin(la[i]) * sin(la[j]) +
                       cos(la[i]) * cos(la[j]) * cos(lo[i] - lo[j])))
  })
  diag(D) <- 0
  K <- (D <= cutoff) * (if (kernel == "bartlett") 1 - D / cutoff else 1)
  S <- X * (w * e)
  B <- solve(crossprod(X, X * w))
  B %*% (t(S) %*% K %*% S) %*% B
}

test_that("weighted felm fit matches the dense reference", {
  skip_if_not_installed("lfe")
  d <- make_cs_data()
  fit <- lfe::felm(y ~ x | g, data = d, weights = d$w, keepCX = TRUE)
  for (kern in c("uniform", "bartlett")) {
    V <- vcovSpHAC(fit, lat = "latitude", lon = "longitude",
                   kernel = kern, dist_fn = "spherical",
                   dist_cutoff = 200, ncores = 2, ssc = FALSE, data = d)
    # lfe stores sqrt(w); the reference wants w on the original scale.
    Vref <- dense_ref(fit$cX, as.numeric(fit$residuals),
                      as.numeric(fit$weights)^2,
                      d$latitude, d$longitude, 200, kern)
    expect_lt(max(abs(V - Vref)) / max(abs(Vref)), 1e-10)
  }
})

test_that("weighted fixest fit matches the dense reference and felm", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  d <- make_cs_data()
  fx <- fixest::feols(y ~ x | g, data = d, weights = ~w, demeaned = TRUE)
  ff <- lfe::felm(y ~ x | g, data = d, weights = d$w, keepCX = TRUE)
  V_fx <- vcovSpHAC(fx, lat = "latitude", lon = "longitude",
                    kernel = "bartlett", dist_fn = "spherical",
                    dist_cutoff = 200, ncores = 2, ssc = FALSE, data = d)
  V_ff <- vcovSpHAC(ff, lat = "latitude", lon = "longitude",
                    kernel = "bartlett", dist_fn = "spherical",
                    dist_cutoff = 200, ncores = 2, ssc = FALSE, data = d)
  Vref <- dense_ref(fx$X_demeaned, as.numeric(fx$residuals),
                    as.numeric(fx$weights),
                    d$latitude, d$longitude, 200, "bartlett")
  expect_lt(max(abs(V_fx - Vref)) / max(abs(Vref)), 1e-10)
  expect_lt(max(abs(V_fx - V_ff)) / max(abs(V_fx)), 1e-10)
})

test_that("weighted fit matches fixest's own Conley vcov exactly at tiny cutoff", {
  # At a sub-meter cutoff only self-pairs survive, so fixest's approximate
  # spherical distance cannot move any pair across the boundary and the
  # comparison is exact. Both calls use their package defaults: since v0.9.0
  # fastconley's ssc/psd_fix defaults mirror fixest's.
  skip_if_not_installed("fixest")
  d <- make_cs_data()
  fx <- fixest::feols(y ~ x | g, data = d, weights = ~w, demeaned = TRUE)
  V <- vcovSpHAC(fx, lat = "latitude", lon = "longitude",
                 kernel = "uniform", dist_fn = "spherical",
                 dist_cutoff = 1e-4, ncores = 2, data = d)
  V_fx <- vcov(fx, vcov = fixest::conley(1e-4, distance = "spherical"))
  expect_lt(max(abs(V - V_fx)) / max(abs(V_fx)), 1e-12)
})

test_that("ssc (default TRUE) scales by n / df.residual, matching fixest", {
  skip_if_not_installed("fixest")
  d <- make_cs_data()
  fx <- fixest::feols(y ~ x | g, data = d, weights = ~w, demeaned = TRUE)
  n <- nrow(d)
  V0 <- vcovSpHAC(fx, lat = "latitude", lon = "longitude",
                  kernel = "uniform", dist_fn = "spherical",
                  dist_cutoff = 1e-4, ncores = 2, ssc = FALSE, data = d)
  V1 <- vcovSpHAC(fx, lat = "latitude", lon = "longitude",
                  kernel = "uniform", dist_fn = "spherical",
                  dist_cutoff = 1e-4, ncores = 2, data = d)
  expect_equal(V1[1, 1] / V0[1, 1], n / (fx$nobs - fx$nparams),
               tolerance = 1e-14)
  # fixest's default ssc: only the n/(n-K) adjustment binds for Conley.
  V_fx <- vcov(fx, vcov = fixest::conley(1e-4, distance = "spherical"),
               vcov_fix = FALSE)
  expect_lt(max(abs(V1 - V_fx)) / max(abs(V_fx)), 1e-12)
})

test_that("ssc default on felm uses the fit's df.residual", {
  skip_if_not_installed("lfe")
  d <- make_cs_data()
  fit <- lfe::felm(y ~ x | g, data = d, keepCX = TRUE)
  V0 <- vcovSpHAC(fit, lat = "latitude", lon = "longitude",
                  kernel = "bartlett", dist_fn = "haversine",
                  dist_cutoff = 200, ncores = 2, ssc = FALSE, data = d)
  V1 <- vcovSpHAC(fit, lat = "latitude", lon = "longitude",
                  kernel = "bartlett", dist_fn = "haversine",
                  dist_cutoff = 200, ncores = 2, data = d)
  expect_equal(V1[1, 1] / V0[1, 1], nrow(d) / fit$df.residual,
               tolerance = 1e-14)
})

test_that("psd_fix clamps negative eigenvalues with fixest semantics", {
  skip_if_not_installed("fixest")
  # Uniform kernel, generous cutoff, several noise regressors: hunt a seed
  # whose vcov is not PSD, then check the warning and the clamped result.
  found <- FALSE
  for (s in 1:40) {
    set.seed(s)
    m <- 250
    dd <- data.frame(lat = runif(m, 40, 42), lon = runif(m, -2, 2),
                     x1 = rnorm(m), x2 = rnorm(m), x3 = rnorm(m))
    dd$y <- rnorm(m)
    fz <- fixest::feols(y ~ x1 + x2 + x3, data = dd, demeaned = TRUE)
    Vraw <- suppressWarnings(
      vcovSpHAC(fz, lat = "lat", lon = "lon", kernel = "uniform",
                dist_fn = "spherical", dist_cutoff = 150, ncores = 2,
                psd_fix = FALSE, ssc = FALSE, data = dd))
    ev <- eigen(Vraw, symmetric = TRUE)
    if (min(ev$values) <= 0 &&
        max(abs(Vraw - tcrossprod(ev$vectors %*% diag(pmax(ev$values, 1e-16),
                                                      4), ev$vectors))) > 1e-8) {
      found <- TRUE
      expect_warning(
        vcovSpHAC(fz, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 150, ncores = 2,
                  psd_fix = FALSE, ssc = FALSE, data = dd),
        "not positive semi-definite")
      # default psd_fix = TRUE clamps (and warns that it did)
      expect_warning(
        Vfix <- vcovSpHAC(fz, lat = "lat", lon = "lon", kernel = "uniform",
                          dist_fn = "spherical", dist_cutoff = 150,
                          ncores = 2, ssc = FALSE, data = dd),
        "was\\b.*fixed|clamping")
      expect_gte(min(eigen(Vfix, symmetric = TRUE, only.values = TRUE)$values), 0)
      Vref <- tcrossprod(ev$vectors %*% diag(pmax(ev$values, 1e-16), 4),
                         ev$vectors)
      expect_lt(max(abs(Vfix - Vref)), 1e-18)
      break
    }
  }
  skip_if_not(found, "no non-PSD vcov found in 40 seeds")
})

test_that("lat/lon auto-detection finds, reports, and validates columns", {
  skip_if_not_installed("fixest")
  d <- make_cs_data(n = 200)
  fx <- fixest::feols(y ~ x | g, data = d, demeaned = TRUE)
  expect_message(
    V_auto <- vcovSpHAC(fx, kernel = "uniform", dist_fn = "spherical",
                        dist_cutoff = 200, ncores = 2, data = d),
    'lat = "latitude", lon = "longitude"')
  V_named <- vcovSpHAC(fx, lat = "latitude", lon = "longitude",
                       kernel = "uniform", dist_fn = "spherical",
                       dist_cutoff = 200, ncores = 2, data = d)
  expect_identical(V_auto, V_named)
  d2 <- d
  names(d2)[names(d2) == "latitude"] <- "yco"
  fx2 <- fixest::feols(y ~ x | g, data = d2, demeaned = TRUE)
  expect_error(
    vcovSpHAC(fx2, lon = "longitude", kernel = "uniform",
              dist_fn = "spherical", dist_cutoff = 200, data = d2),
    "auto-detect the latitude")
})
