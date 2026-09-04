# Brute-force vs. templated FastSpatialMeat across all six (distance × kernel)
# combinations on a small workload. Mirrors tests/manual/test-equiv-distances.R
# but trimmed for R CMD check.

naive_meat <- function(lat, lon, time, X, e, cutoff, kernel, dist_fn) {
  n <- length(e); k <- ncol(X); AVG_ERAD <- 6371
  d_pair <- function(i, j) {
    lat1 <- lat[i] * pi / 180; lat2 <- lat[j] * pi / 180
    lon1 <- lon[i] * pi / 180; lon2 <- lon[j] * pi / 180
    switch(dist_fn,
      haversine = {
        a <- sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
        AVG_ERAD * 2 * atan2(sqrt(min(max(a, 0), 1)), sqrt(1 - min(max(a, 0), 1)))
      },
      spherical = {
        d <- sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)
        AVG_ERAD * acos(min(max(d, -1), 1))
      },
      chord = {
        ux1 <- cos(lat1) * cos(lon1); uy1 <- cos(lat1) * sin(lon1); uz1 <- sin(lat1)
        ux2 <- cos(lat2) * cos(lon2); uy2 <- cos(lat2) * sin(lon2); uz2 <- sin(lat2)
        AVG_ERAD * sqrt((ux1 - ux2)^2 + (uy1 - uy2)^2 + (uz1 - uz2)^2)
      })
  }
  M <- matrix(0, k, k)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j || time[i] != time[j]) next
    d <- d_pair(i, j); if (d > cutoff) next
    w <- if (kernel == "uniform") 1 else 1 - d / cutoff
    si <- e[i] * X[i, ]; sj <- e[j] * X[j, ]
    M <- M + tcrossprod(si, w * sj)
  }
  for (i in seq_len(n)) {
    si <- e[i] * X[i, ]
    M <- M + tcrossprod(si)
  }
  M
}

test_that("FastSpatialMeat matches brute force across all distance/kernel cells", {

  cases <- list(
    list(label = "cross-section", time = rep(1L, 60L), balanced = FALSE),
    list(label = "balanced panel", time = rep(1:3, each = 20L), balanced = TRUE)
  )

  for (cs in cases) {
    n <- length(cs$time); k <- 3L
    set.seed(11)
    if (cs$balanced) {
      n_unit <- 20L
      ulat <- stats::runif(n_unit, 25, 50); ulon <- stats::runif(n_unit, -125, -70)
      lat <- rep(ulat, times = 3L); lon <- rep(ulon, times = 3L)
    } else {
      lat <- stats::runif(n, 25, 50); lon <- stats::runif(n, -125, -70)
    }
    X <- cbind(1, matrix(stats::rnorm(n * (k - 1L)), n, k - 1L))
    e <- stats::rnorm(n)

    ord <- order(cs$time, lat, lon)
    lat <- lat[ord]; lon <- lon[ord]; time <- cs$time[ord]
    X <- X[ord, , drop = FALSE]; e <- e[ord]

    for (df in c("haversine", "spherical", "chord")) {
      for (kn in c("bartlett", "uniform")) {
        M_fast <- fastconley:::FastSpatialMeat(
          lat = lat, lon = lon, time = time, X = X, e = e,
          cutoff = 400, kernel = kn, dist_fn = df,
          balanced_pnl = cs$balanced, ncores = 1L)
        M_ref <- naive_meat(lat, lon, time, X, e, 400, kn, df)
        expect_lt(max(abs(M_fast - M_ref)), 1e-10,
                  label = paste(cs$label, df, kn))
      }
    }
  }
})
