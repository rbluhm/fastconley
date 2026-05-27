suppressMessages({
  library(lfe)
  library(fastconley)
  library(data.table)
})

RcppParallel::setThreadOptions(numThreads = 1L)

# Naive O(n^2) brute-force reference. Computes the spatial meat from scratch
# using simple per-pair distance formulas, matching what the new templated
# pair_weight specializations should produce.
naive_meat <- function(lat, lon, time, X, e, cutoff, kernel, dist_fn) {
  n <- length(e)
  k <- ncol(X)
  AVG_ERAD <- 6371
  FLAT_KM_PER_DEG <- 111
  d_pair <- function(i, j) {
    lat1 <- lat[i] * pi / 180; lat2 <- lat[j] * pi / 180
    lon1 <- lon[i] * pi / 180; lon2 <- lon[j] * pi / 180
    switch(dist_fn,
      haversine = {
        a <- sin((lat2 - lat1) / 2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)^2
        a <- min(max(a, 0), 1)
        AVG_ERAD * 2 * atan2(sqrt(a), sqrt(1 - a))
      },
      spherical = {
        d <- sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)
        d <- min(max(d, -1), 1)
        AVG_ERAD * acos(d)
      },
      chord = {
        ux1 <- cos(lat1) * cos(lon1); uy1 <- cos(lat1) * sin(lon1); uz1 <- sin(lat1)
        ux2 <- cos(lat2) * cos(lon2); uy2 <- cos(lat2) * sin(lon2); uz2 <- sin(lat2)
        AVG_ERAD * sqrt((ux1 - ux2)^2 + (uy1 - uy2)^2 + (uz1 - uz2)^2)
      },
      flatearth = {
        dlat_km <- FLAT_KM_PER_DEG * (lat[i] - lat[j])
        dlon_km <- cos(lat1) * FLAT_KM_PER_DEG * (lon[i] - lon[j])
        sqrt(dlat_km^2 + dlon_km^2)
      }
    )
  }
  M <- matrix(0, k, k)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (time[i] != time[j] || i == j) next
    d <- d_pair(i, j)
    if (d > cutoff) next
    w <- if (kernel == "uniform") 1 else 1 - d / cutoff
    si <- e[i] * X[i, ]
    sj <- e[j] * X[j, ]
    M <- M + tcrossprod(si, w * sj)
  }
  # diagonal i = j contribution: e_i^2 X_i X_i' (weight 1)
  for (i in seq_len(n)) {
    si <- e[i] * X[i, ]
    M <- M + tcrossprod(si)
  }
  M
}

# Direct call to FastSpatialMeat (the templated path).
fast_meat <- function(lat, lon, time, X, e, cutoff, kernel, dist_fn, balanced) {
  fastconley:::FastSpatialMeat(
    lat = lat, lon = lon, time = time, X = X, e = e, cutoff = cutoff,
    kernel = kernel, dist_fn = dist_fn,
    balanced_pnl = balanced, ncores = 1L
  )
}

run <- function(label, lat, lon, time, X, e, cutoff, balanced) {
  cat(sprintf("\n== %s ==\n", label))
  for (df in c("haversine", "spherical", "chord")) {
    for (kn in c("bartlett", "uniform")) {
      ord <- order(time, lat, lon)
      lat_o <- lat[ord]; lon_o <- lon[ord]; time_o <- time[ord]
      X_o <- X[ord, , drop = FALSE]; e_o <- e[ord]

      M_fast <- fast_meat(lat_o, lon_o, time_o, X_o, e_o, cutoff, kn, df, balanced)
      M_ref  <- naive_meat(lat_o, lon_o, time_o, X_o, e_o, cutoff, kn, df)
      diff <- max(abs(M_fast - M_ref))
      cat(sprintf("  %-10s %-9s  max|fast - naive| = %.3e\n", df, kn, diff))
    }
  }
}

set.seed(1)

# Cross-section
n <- 80
lat <- runif(n, 25, 50); lon <- runif(n, -125, -70)
X <- cbind(1, matrix(rnorm(n * 3), n, 3)); e <- rnorm(n)
time <- rep(1L, n)
run("cross-section n=80, cutoff=400", lat, lon, time, X, e, 400, balanced = FALSE)

# Balanced panel, time-invariant coordinates
G <- 25; TT <- 4
ulat <- runif(G, 25, 50); ulon <- runif(G, -125, -70)
unit <- rep(seq_len(G), times = TT); tt <- rep(seq_len(TT), each = G)
lat2 <- ulat[unit]; lon2 <- ulon[unit]
X2 <- cbind(1, matrix(rnorm(G * TT * 3), G * TT, 3)); e2 <- rnorm(G * TT)
run("balanced panel G=25 T=4, cutoff=400", lat2, lon2, tt, X2, e2, 400, balanced = TRUE)

# Unbalanced (drop ~20%)
keep <- sort(sample.int(G * TT, size = round(0.8 * G * TT)))
run("unbalanced panel ~80 rows, cutoff=400",
    lat2[keep], lon2[keep], tt[keep], X2[keep, , drop = FALSE], e2[keep],
    400, balanced = FALSE)
