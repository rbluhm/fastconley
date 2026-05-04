# Manual correctness checks for the fast CSR spatial meat.
# Run after installing the package from this source directory:
#   R CMD INSTALL .
#   Rscript tests/manual/test-fast-spatial.R

library(fastconley)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad; lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad; lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  R * 2 * atan2(sqrt(pmin(1, a)), sqrt(pmax(0, 1 - a)))
}

naive_spatial_meat <- function(lat, lon, time, X, e, cutoff, kernel = "bartlett") {
  S <- X * e
  out <- matrix(0, ncol(X), ncol(X))
  for (tt in unique(time)) {
    ii <- which(time == tt)
    for (a in seq_along(ii)) {
      i <- ii[a]
      for (b in seq_along(ii)) {
        j <- ii[b]
        d <- haversine_km(lat[i], lon[i], lat[j], lon[j])
        if (d <= cutoff) {
          w <- if (kernel == "uniform") 1 else if (cutoff <= 0) 1 else 1 - d / cutoff
          out <- out + w * tcrossprod(S[i, ], S[j, ])
        }
      }
    }
  }
  out
}

set.seed(42)
X <- cbind(1, matrix(rnorm(12 * 4), 12, 4))
e <- rnorm(12)
lat <- runif(12, 40, 43)
lon <- runif(12, -75, -70)
time <- rep(1, 12)

fast <- fastconley:::FastSpatialMeat(lat, lon, time, X, e, cutoff = 250,
                                 kernel = "bartlett", dist_fn = "haversine",
                                 balanced_pnl = FALSE, ncores = 2)
slow <- naive_spatial_meat(lat, lon, time, X, e, cutoff = 250, kernel = "bartlett")
stopifnot(max(abs(fast - slow)) < 1e-8)

# Unbalanced panel: three time blocks with different sizes.
time2 <- c(rep(1, 5), rep(2, 3), rep(3, 4))
fast2 <- fastconley:::FastSpatialMeat(lat, lon, time2, X, e, cutoff = 250,
                                  kernel = "uniform", dist_fn = "haversine",
                                  balanced_pnl = FALSE, ncores = 2)
slow2 <- naive_spatial_meat(lat, lon, time2, X, e, cutoff = 250, kernel = "uniform")
stopifnot(max(abs(fast2 - slow2)) < 1e-8)

# Balanced panel reuse: same coordinates/order in every time block.
n_unit <- 4
n_time <- 3
lat_u <- runif(n_unit, 40, 43)
lon_u <- runif(n_unit, -75, -70)
lat3 <- rep(lat_u, n_time)
lon3 <- rep(lon_u, n_time)
time3 <- rep(seq_len(n_time), each = n_unit)
X3 <- cbind(1, matrix(rnorm(n_unit * n_time * 3), n_unit * n_time, 3))
e3 <- rnorm(n_unit * n_time)
fast3 <- fastconley:::FastSpatialMeat(lat3, lon3, time3, X3, e3, cutoff = 250,
                                  kernel = "bartlett", dist_fn = "haversine",
                                  balanced_pnl = TRUE, ncores = 2)
slow3 <- naive_spatial_meat(lat3, lon3, time3, X3, e3, cutoff = 250, kernel = "bartlett")
stopifnot(max(abs(fast3 - slow3)) < 1e-8)

cat("All fast spatial meat checks passed.\n")
