# Phase 2 checks: stacked balanced meat, deterministic reduction,
# csr_weight = "float", and the scores-based entry points.
#   Rscript tests/manual/test-phase2.R

suppressMessages(library(fastconley))

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)

set.seed(77)
n_unit <- 900; n_time <- 6
lat_u <- runif(n_unit, -60, 70); lon_u <- runif(n_unit, -180, 180)
lat <- rep(lat_u, n_time); lon <- rep(lon_u, n_time)
time <- rep(seq_len(n_time), each = n_unit)
n <- n_unit * n_time; k <- 4
X <- cbind(1, matrix(rnorm(n * (k - 1)), n, k - 1))
e <- rnorm(n)
scores <- X * e

# ---- 1. Determinism: ncores = 1 and ncores = 8 are bit-identical -----------
for (bal in c(TRUE, FALSE)) for (kern in c("bartlett", "uniform")) {
  m1 <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores,
                                     cutoff = 800, kernel = kern,
                                     dist_fn = "spherical",
                                     balanced_pnl = bal, ncores = 1)
  m8 <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores,
                                     cutoff = 800, kernel = kern,
                                     dist_fn = "spherical",
                                     balanced_pnl = bal, ncores = 8)
  stopifnot(identical(m1, m8))
}
ord_ut <- order(rep(seq_len(n_unit), n_time), time)
s1 <- fastconley:::FastSerialHacPanel(rep(seq_len(n_unit), n_time)[ord_ut],
                                      time[ord_ut], cutoff = 3,
                                      scores = scores[ord_ut, ], ncores = 1)
s8 <- fastconley:::FastSerialHacPanel(rep(seq_len(n_unit), n_time)[ord_ut],
                                      time[ord_ut], cutoff = 3,
                                      scores = scores[ord_ut, ], ncores = 8)
stopifnot(identical(s1, s8))
cat("1. ncores = 1 vs 8 bit-identical (spatial bal/gen x kernels + serial): OK\n")

# ---- 2. scores= is identical to the X,e convenience path -------------------
a <- fastconley:::FastSpatialMeat(lat, lon, time, X = X, e = e, cutoff = 800,
                                  kernel = "bartlett", dist_fn = "haversine",
                                  balanced_pnl = TRUE, ncores = 2)
b <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores, cutoff = 800,
                                  kernel = "bartlett", dist_fn = "haversine",
                                  balanced_pnl = TRUE, ncores = 2)
stopifnot(identical(a, b))
cat("2. scores= equals X,e wrapper path: OK\n")

# ---- 3. Stacked balanced meat == general path (same math) ------------------
worst <- 0
for (kern in c("bartlett", "uniform")) for (df in c("haversine", "spherical", "chord")) {
  bal <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores,
                                      cutoff = 700, kernel = kern, dist_fn = df,
                                      balanced_pnl = TRUE, ncores = 1)
  gen <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores,
                                      cutoff = 700, kernel = kern, dist_fn = df,
                                      balanced_pnl = FALSE, ncores = 1)
  worst <- max(worst, rel_diff(bal, gen))
}
stopifnot(worst < 1e-12)
cat(sprintf("3. balanced (stacked CSR) vs general: worst rel diff %.2e\n", worst))

# ---- 4. Period chunking: T*k wide enough to force multiple chunks ----------
# PCHUNK = max(1, 1024/k) periods; with k = 40, chunks cover 25 periods, so
# T = 60 forces three CSR passes. Compare against the general path.
nu2 <- 150; nt2 <- 60; k2 <- 40
lat2 <- rep(runif(nu2, 30, 50), nt2); lon2 <- rep(runif(nu2, -120, -70), nt2)
time2 <- rep(seq_len(nt2), each = nu2)
sc2 <- matrix(rnorm(nu2 * nt2 * k2), nu2 * nt2, k2)
bal2 <- fastconley:::FastSpatialMeat(lat2, lon2, time2, scores = sc2, cutoff = 400,
                                     kernel = "bartlett", dist_fn = "spherical",
                                     balanced_pnl = TRUE, ncores = 4)
gen2 <- fastconley:::FastSpatialMeat(lat2, lon2, time2, scores = sc2, cutoff = 400,
                                     kernel = "bartlett", dist_fn = "spherical",
                                     balanced_pnl = FALSE, ncores = 4)
stopifnot(rel_diff(bal2, gen2) < 1e-12)
cat(sprintf("4. period-chunked stacked pass (T=60, k=40): rel diff %.2e\n",
            rel_diff(bal2, gen2)))

# ---- 5. csr_weight = "float": close to double, not equal -------------------
fd <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores, cutoff = 800,
                                   kernel = "bartlett", dist_fn = "spherical",
                                   balanced_pnl = TRUE, ncores = 1,
                                   csr_weight = "double")
ff <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores, cutoff = 800,
                                   kernel = "bartlett", dist_fn = "spherical",
                                   balanced_pnl = TRUE, ncores = 1,
                                   csr_weight = "float")
rdf <- rel_diff(ff, fd)
stopifnot(rdf > 0, rdf < 1e-6)
# uniform ignores csr_weight entirely
fu <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores, cutoff = 800,
                                   kernel = "uniform", dist_fn = "spherical",
                                   balanced_pnl = TRUE, ncores = 1,
                                   csr_weight = "float")
du <- fastconley:::FastSpatialMeat(lat, lon, time, scores = scores, cutoff = 800,
                                   kernel = "uniform", dist_fn = "spherical",
                                   balanced_pnl = TRUE, ncores = 1,
                                   csr_weight = "double")
stopifnot(identical(fu, du))
cat(sprintf("5. csr_weight float vs double: rel diff %.2e (uniform unaffected)\n", rdf))

# ---- 6. Serial HAC against a naive R reference ------------------------------
naive_serial <- function(unit, tt, L, S) {
  out <- matrix(0, ncol(S), ncol(S))
  for (u in unique(unit)) {
    ii <- which(unit == u)
    for (i in ii) for (j in ii) {
      dt <- abs(tt[j] - tt[i])
      if (dt <= L && dt != 0) {
        out <- out + (1 - dt / (L + 1)) * tcrossprod(S[i, ], S[j, ])
      }
    }
  }
  out
}
set.seed(88)
nu3 <- 25; nt3 <- 7
unit3 <- rep(seq_len(nu3), each = nt3); time3 <- rep(seq_len(nt3), nu3)
S3 <- matrix(rnorm(nu3 * nt3 * 3), nu3 * nt3, 3)
fast3 <- fastconley:::FastSerialHacPanel(unit3, time3, cutoff = 2,
                                         scores = S3, ncores = 2)
stopifnot(rel_diff(fast3, naive_serial(unit3, time3, 2, S3)) < 1e-12)
cat("6. serial HAC vs naive reference: OK\n")

# ---- 7. vcovSpHAC end-to-end: csr_weight accepted, determinism --------------
suppressMessages(library(lfe))
set.seed(99)
nu4 <- 100; nt4 <- 4; n4 <- nu4 * nt4
d <- data.frame(unit = rep(seq_len(nu4), nt4), time = rep(seq_len(nt4), each = nu4),
                lat = rep(runif(nu4, -50, 60), nt4), lon = rep(runif(nu4, -150, 150), nt4),
                x1 = rnorm(n4), x2 = rnorm(n4))
d$y <- 0.3 * d$x1 + 0.1 * d$x2 + rnorm(n4)
fit <- felm(y ~ x1 + x2 | unit + time, data = d, keepCX = TRUE)
v1 <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                dist_cutoff = 900, lag_cutoff = 1, balanced_pnl = TRUE, ncores = 1)
v8 <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                dist_cutoff = 900, lag_cutoff = 1, balanced_pnl = TRUE, ncores = 8)
stopifnot(identical(v1, v8))
vf <- vcovSpHAC(fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
                dist_cutoff = 900, lag_cutoff = 1, balanced_pnl = TRUE,
                ncores = 1, csr_weight = "float")
stopifnot(rel_diff(vf, v1) < 1e-6)
cat("7. vcovSpHAC: ncores-invariant results, csr_weight accepted: OK\n")

cat("\nAll Phase 2 checks passed.\n")
