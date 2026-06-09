# Checks for the minor-backlog items (v0.6.1):
#   M1 serial HAC O(T*k) sliding window, M2 data= passthrough for felm,
#   M3 deterministic parallel prep.
#   Rscript tests/manual/test-minors.R

suppressMessages({library(lfe); library(fastconley)})
rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)

# ---- M2: data= passthrough ---------------------------------------------------
set.seed(42)
n <- 5000
d <- data.frame(unit = rep(1:1000, 5), time = rep(1:5, each = 1000),
                lat = rep(runif(1000, 25, 49), 5),
                lon = rep(runif(1000, -125, -67), 5),
                x1 = rnorm(n), x2 = rnorm(n))
d$y <- 0.5 * d$x1 - 0.2 * d$x2 + rnorm(n)
fit <- felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)

v_default <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                       dist_cutoff = 500, lag_cutoff = 1, balanced_pnl = TRUE,
                       ncores = 2)
v_data <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                    dist_cutoff = 500, lag_cutoff = 1, balanced_pnl = TRUE,
                    ncores = 2, data = d)
stopifnot(identical(v_default, v_data))

# NA rows dropped at fit time -> model-frame fallback, still consistent.
d2 <- d; d2$x1[c(7, 1234)] <- NA
fit2 <- felm(y ~ x1 + x2 | unit + time, data = d2, keepCX = TRUE)
v2a <- vcovSpHAC(fit2, unit = "unit", time = "time", lat = "lat", lon = "lon",
                 dist_cutoff = 500, ncores = 2)
v2b <- vcovSpHAC(fit2, unit = "unit", time = "time", lat = "lat", lon = "lon",
                 dist_cutoff = 500, ncores = 2, data = d2)
stopifnot(identical(v2a, v2b))

bad <- tryCatch({
  vcovSpHAC(fit, unit = "unit", time = "time", lat = "nope", lon = "lon",
            dist_cutoff = 500, data = d)
  FALSE
}, error = function(e) TRUE)
stopifnot(bad)
cat("M2 data= passthrough (fast path, NA fallback, bad column): OK\n")

# ---- M1: serial HAC sliding window -------------------------------------------
naive_serial <- function(unit, tt, L, S) {
  out <- matrix(0, ncol(S), ncol(S))
  for (u in unique(unit)) {
    ii <- which(unit == u)
    tu <- tt[ii]; Su <- S[ii, , drop = FALSE]
    for (a in seq_along(ii)) {
      w <- 1 - abs(tu - tu[a]) / (L + 1)
      keep <- abs(tu - tu[a]) <= L & tu != tu[a]
      if (any(keep)) {
        out <- out + tcrossprod(Su[a, ], colSums(Su[keep, , drop = FALSE] * w[keep]))
      }
    }
  }
  out
}

set.seed(7)
# Long-T panel (would be O(T^2)-slow before), plus irregular times and ties.
for (cfg in list(list(nu = 50, T = 400, L = 25),
                 list(nu = 20, T = 60,  L = 7))) {
  unit <- rep(seq_len(cfg$nu), each = cfg$T)
  tt <- as.numeric(rep(seq_len(cfg$T), cfg$nu))
  k <- 4
  S <- matrix(rnorm(cfg$nu * cfg$T * k), cfg$nu * cfg$T, k)
  fast <- fastconley:::FastSerialHacPanel(unit, tt, cutoff = cfg$L,
                                          scores = S, ncores = 2)
  stopifnot(rel_diff(fast, naive_serial(unit, tt, cfg$L, S)) < 1e-10)
}
# Gaps + duplicate times within unit (ties must be excluded), large raw times
# (exercises the per-block shift conditioning).
set.seed(8)
unit <- rep(1:30, each = 20)
tt <- as.numeric(2000 + sort(sample(0:40, 20, replace = TRUE)))
tt <- as.numeric(unlist(lapply(1:30, function(u) 2000 + sort(sample(0:40, 20, replace = TRUE)))))
S <- matrix(rnorm(600 * 3), 600, 3)
fast <- fastconley:::FastSerialHacPanel(unit, tt, cutoff = 5, scores = S, ncores = 2)
stopifnot(rel_diff(fast, naive_serial(unit, tt, 5, S)) < 1e-10)

# Unsorted input now errors with a clear message.
err <- tryCatch({
  fastconley:::FastSerialHacPanel(unit, rev(tt), cutoff = 2, scores = S)
  FALSE
}, error = function(e) grepl("sorted by time", conditionMessage(e)))
stopifnot(err)
cat("M1 serial HAC (long T, gaps, ties, big raw times, sorted guard): OK\n")

# ---- M3: parallel prep stays deterministic ------------------------------------
set.seed(9)
n <- 50000
lat <- runif(n, -85, 85); lon <- runif(n, -180, 180)
scores <- matrix(rnorm(n * 3), n, 3)
m1 <- fastconley:::FastSpatialMeat(lat, lon, rep(1, n), scores = scores,
                                   cutoff = 300, kernel = "bartlett",
                                   dist_fn = "haversine", ncores = 1)
m16 <- fastconley:::FastSpatialMeat(lat, lon, rep(1, n), scores = scores,
                                    cutoff = 300, kernel = "bartlett",
                                    dist_fn = "haversine", ncores = 16)
stopifnot(identical(m1, m16))
cat("M3 parallel prep: ncores = 1 vs 16 bit-identical: OK\n")

cat("\nAll minor-backlog checks passed.\n")
