# Manual smoke test for the fast spatial C++ routines.
# Run after installing the package:
#   R CMD INSTALL .
#   Rscript inst/examples/fast-spatial-smoke.R

library(fastconley)

haversine <- function(lat1, lon1, lat2, lon2) {
  r <- 6371
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad
  lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad
  lon2 <- lon2 * to_rad
  a <- sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

ref_meat <- function(coords, X, e, cutoff, kernel = "bartlett") {
  n <- nrow(X)
  S <- X * e
  W <- diag(1, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      d <- haversine(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2])
      if (d <= cutoff) W[i, j] <- if (kernel == "bartlett") 1 - d / cutoff else 1
    }
  }
  crossprod(S, W %*% S)
}

set.seed(42)
n <- 20
k <- 4
coords <- cbind(lat = runif(n, -5, 5), lon = runif(n, 10, 20))
X <- cbind(1, matrix(rnorm(n * (k - 1)), n, k - 1))
e <- rnorm(n)
cutoff <- 700

for (kernel in c("bartlett", "uniform")) {
  fast <- fastconley:::FastSpatialMeat(coords[, 1], coords[, 2], rep(1, n), X, e, cutoff, kernel, "haversine", FALSE, 1L)
  ref <- ref_meat(coords, X, e, cutoff, kernel)
  cat(kernel, "max abs diff:", max(abs(fast - ref)), "\n")
  stopifnot(isTRUE(all.equal(fast, ref, tolerance = 1e-8)))
}

cat("fast spatial smoke test passed\n")
