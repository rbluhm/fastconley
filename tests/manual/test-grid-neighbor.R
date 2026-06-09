# Correctness checks for the 3D cell-grid neighbor search (neighbor = "grid").
#
# The grid and band paths use identical per-pair accept tests (pair_weight),
# so they must agree on every config up to floating-point summation order.
# Run after installing the package:
#   Rscript tests/manual/test-grid-neighbor.R

suppressMessages(library(fastconley))
RcppParallel::setThreadOptions(numThreads = 1L)

AVG_ERAD <- 6371

dist_km <- function(dist_fn, lat1, lon1, lat2, lon2) {
  to_rad <- pi / 180
  la1 <- lat1 * to_rad; lo1 <- lon1 * to_rad
  la2 <- lat2 * to_rad; lo2 <- lon2 * to_rad
  switch(dist_fn,
    haversine = {
      a <- sin((la2 - la1) / 2)^2 + cos(la1) * cos(la2) * sin((lo2 - lo1) / 2)^2
      a <- min(max(a, 0), 1)
      AVG_ERAD * 2 * atan2(sqrt(a), sqrt(1 - a))
    },
    spherical = {
      d <- sin(la1) * sin(la2) + cos(la1) * cos(la2) * cos(lo1 - lo2)
      d <- min(max(d, -1), 1)
      AVG_ERAD * acos(d)
    },
    chord = {
      dx <- cos(la1) * cos(lo1) - cos(la2) * cos(lo2)
      dy <- cos(la1) * sin(lo1) - cos(la2) * sin(lo2)
      dz <- sin(la1) - sin(la2)
      AVG_ERAD * sqrt(dx^2 + dy^2 + dz^2)
    })
}

naive_meat <- function(lat, lon, time, X, e, cutoff, kernel, dist_fn) {
  S <- X * e
  out <- matrix(0, ncol(X), ncol(X))
  for (tt in unique(time)) {
    ii <- which(time == tt)
    for (a in seq_along(ii)) {
      for (b in seq_along(ii)) {
        i <- ii[a]; j <- ii[b]
        d <- dist_km(dist_fn, lat[i], lon[i], lat[j], lon[j])
        if (d <= cutoff) {
          w <- if (kernel == "uniform") 1 else 1 - d / cutoff
          out <- out + w * tcrossprod(S[i, ], S[j, ])
        }
      }
    }
  }
  out
}

rel_diff <- function(a, b) {
  denom <- max(abs(b), 1e-300)
  max(abs(a - b)) / denom
}

run_both <- function(lat, lon, time, X, e, cutoff, kernel, dist_fn, balanced) {
  g <- fastconley:::FastSpatialMeat(lat, lon, time, X, e, cutoff,
                                    kernel = kernel, dist_fn = dist_fn,
                                    balanced_pnl = balanced, ncores = 1,
                                    neighbor = "grid")
  b <- fastconley:::FastSpatialMeat(lat, lon, time, X, e, cutoff,
                                    kernel = kernel, dist_fn = dist_fn,
                                    balanced_pnl = balanced, ncores = 1,
                                    neighbor = "band")
  list(grid = g, band = b, rd = rel_diff(g, b))
}

kernels <- c("bartlett", "uniform")
dists   <- c("haversine", "spherical", "chord")
worst_gb <- 0; worst_naive <- 0; n_checks <- 0

# ---- 1. Brute-force ground truth, small n, varied geometry -----------------
set.seed(101)
geoms <- list(
  regional = list(lat = runif(40, 40, 44), lon = runif(40, -75, -69), cutoff = 200),
  global   = list(lat = runif(40, -85, 85), lon = runif(40, -180, 180), cutoff = 3000),
  polar    = list(lat = runif(40, 87, 89.99), lon = runif(40, -180, 180), cutoff = 150),
  dateline = list(lat = runif(40, -10, 10),
                  lon = ifelse(runif(40) < 0.5, runif(40, 178, 180),
                                                runif(40, -180, -178)),
                  cutoff = 300),
  dupes    = local({
    base_lat <- runif(8, 30, 35); base_lon <- runif(8, 100, 110)
    idx <- sample(8, 40, replace = TRUE)
    list(lat = base_lat[idx], lon = base_lon[idx], cutoff = 250)
  }),
  allpairs = list(lat = runif(15, -60, 60), lon = runif(15, -180, 180),
                  cutoff = 2.1 * pi * AVG_ERAD),  # > max distance: every pair
  tiny     = list(lat = runif(40, 10, 10.001), lon = runif(40, 10, 10.001),
                  cutoff = 0.05)
)
for (gname in names(geoms)) {
  gg <- geoms[[gname]]
  n <- length(gg$lat)
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  e <- rnorm(n)
  time <- rep(1, n)
  for (kern in kernels) for (df in dists) {
    naive <- naive_meat(gg$lat, gg$lon, time, X, e, gg$cutoff, kern, df)
    res <- run_both(gg$lat, gg$lon, time, X, e, gg$cutoff, kern, df, FALSE)
    rd_naive <- rel_diff(res$grid, naive)
    # spherical+bartlett: acos near dot = 1 is ill-conditioned, so two correct
    # implementations legitimately disagree on short distances by up to
    # ~R*sqrt(2*eps)/cutoff in the Bartlett weight. The uniform kernel (pure
    # accept/reject) and the grid-vs-band check below stay strict.
    tol_naive <- if (kern == "bartlett" && df == "spherical") {
      max(1e-8, 20 * AVG_ERAD * sqrt(2 * .Machine$double.eps) / gg$cutoff)
    } else 1e-10
    stopifnot(rd_naive < tol_naive, res$rd < 1e-12)
    worst_naive <- max(worst_naive, rd_naive)
    worst_gb <- max(worst_gb, res$rd)
    n_checks <- n_checks + 1
  }
}
cat(sprintf("1. brute-force geometries: %d checks, worst vs naive %.2e, worst grid-vs-band %.2e\n",
            n_checks, worst_naive, worst_gb))

# ---- 2. Grid vs band, larger random configs (general + unbalanced) ---------
set.seed(202)
n <- 3000
lat <- runif(n, 25, 49); lon <- runif(n, -125, -67)
X <- cbind(1, matrix(rnorm(n * 4), n, 4)); e <- rnorm(n)
time1 <- rep(1, n)
time3 <- sort(sample(1:3, n, replace = TRUE))   # unbalanced blocks
for (kern in kernels) for (df in dists) for (cut in c(50, 500, 2500)) {
  r1 <- run_both(lat, lon, time1, X, e, cut, kern, df, FALSE)
  r2 <- run_both(lat, lon, time3, X, e, cut, kern, df, FALSE)
  stopifnot(r1$rd < 1e-12, r2$rd < 1e-12)
  worst_gb <- max(worst_gb, r1$rd, r2$rd)
}
cat(sprintf("2. general/unbalanced grid-vs-band: worst rel diff %.2e\n", worst_gb))

# ---- 3. Balanced path: grid CSR vs band CSR vs general grid ----------------
set.seed(303)
n_unit <- 800; n_time <- 4
lat_u <- runif(n_unit, -55, 70); lon_u <- runif(n_unit, -180, 180)
latb <- rep(lat_u, n_time); lonb <- rep(lon_u, n_time)
timeb <- rep(seq_len(n_time), each = n_unit)
nb <- n_unit * n_time
Xb <- cbind(1, matrix(rnorm(nb * 3), nb, 3)); eb <- rnorm(nb)
for (kern in kernels) for (df in dists) for (cut in c(100, 1500)) {
  rb <- run_both(latb, lonb, timeb, Xb, eb, cut, kern, df, TRUE)
  gen <- fastconley:::FastSpatialMeat(latb, lonb, timeb, Xb, eb, cut,
                                      kernel = kern, dist_fn = df,
                                      balanced_pnl = FALSE, ncores = 1,
                                      neighbor = "grid")
  stopifnot(rb$rd < 1e-12, rel_diff(rb$grid, gen) < 1e-12)
  worst_gb <- max(worst_gb, rb$rd)
}
cat(sprintf("3. balanced grid CSR: worst rel diff %.2e\n", worst_gb))

# ---- 4. Degenerate sizes ----------------------------------------------------
for (nn in c(1, 2, 3)) {
  lat <- runif(nn, 0, 1); lon <- runif(nn, 0, 1)
  X <- matrix(rnorm(nn), nn, 1); e <- rnorm(nn)
  r <- run_both(lat, lon, rep(1, nn), X, e, 100, "bartlett", "haversine", FALSE)
  stopifnot(r$rd < 1e-12)
}
cat("4. degenerate sizes (n = 1, 2, 3; k = 1): OK\n")

# ---- 5. Multicore sanity ----------------------------------------------------
set.seed(404)
n <- 5000
lat <- runif(n, -85, 85); lon <- runif(n, -180, 180)
X <- cbind(1, matrix(rnorm(n * 3), n, 3)); e <- rnorm(n)
m1 <- fastconley:::FastSpatialMeat(lat, lon, rep(1, n), X, e, 800,
                                   kernel = "bartlett", dist_fn = "spherical",
                                   balanced_pnl = FALSE, ncores = 1,
                                   neighbor = "grid")
m4 <- fastconley:::FastSpatialMeat(lat, lon, rep(1, n), X, e, 800,
                                   kernel = "bartlett", dist_fn = "spherical",
                                   balanced_pnl = FALSE, ncores = 4,
                                   neighbor = "grid")
stopifnot(rel_diff(m4, m1) < 1e-9)
cat("5. multicore (1 vs 4 threads): OK\n")

# ---- 6. R API: vcovSpHAC neighbor argument ----------------------------------
suppressMessages(library(lfe))
set.seed(505)
n_unit <- 120; n_time <- 3; n <- n_unit * n_time
d <- data.frame(
  unit = rep(seq_len(n_unit), times = n_time),
  time = rep(seq_len(n_time), each = n_unit),
  lat = rep(runif(n_unit, -80, 80), times = n_time),
  lon = rep(runif(n_unit, -180, 180), times = n_time),
  x1 = rnorm(n), x2 = rnorm(n)
)
d$y <- 0.4 * d$x1 - 0.2 * d$x2 + rnorm(n)
fit <- felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)
vg <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                dist_cutoff = 1000, lag_cutoff = 1, balanced_pnl = TRUE,
                ncores = 1, neighbor = "grid")
vb <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                dist_cutoff = 1000, lag_cutoff = 1, balanced_pnl = TRUE,
                ncores = 1, neighbor = "band")
stopifnot(rel_diff(vg, vb) < 1e-12)
bad <- tryCatch({
  vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
            dist_cutoff = 1000, neighbor = "nope")
  FALSE
}, error = function(e) TRUE)
stopifnot(bad)
cat("6. vcovSpHAC neighbor argument (grid == band, bad value errors): OK\n")

cat(sprintf("\nAll grid-neighbor checks passed. Worst grid-vs-band rel diff: %.2e\n",
            worst_gb))
