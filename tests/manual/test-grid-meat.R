# Checks for the C2 grid-native exact meat (FastGridMeat + method= dispatch).
#   Rscript tests/manual/test-grid-meat.R

suppressMessages({library(lfe); library(fastconley)})
rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)

# ---- 1. FastGridMeat vs pairwise on a lattice -------------------------------
set.seed(3)
step <- 0.1
lats <- seq(40, 47.9, by = step) + step / 2
lons <- seq(-10, 1.9, by = step) + step / 2
gp <- expand.grid(lat = lats, lon = lons)
n <- nrow(gp); k <- 3
S <- matrix(rnorm(n * k), n, k)
ring <- as.integer(round((gp$lat - lats[1]) / step))
col  <- as.integer(round((gp$lon - lons[1]) / step))
worst <- 0
for (df in c("spherical", "haversine", "chord")) for (cut in c(30, 120, 400, 5000)) {
  pw <- fastconley:::FastSpatialMeat(gp$lat, gp$lon, rep(1, n), scores = S,
                                     cutoff = cut, kernel = "uniform",
                                     dist_fn = df, ncores = 4)
  gm <- fastconley:::FastGridMeat(ring, col, rep(1, n), S, lats[1], step, step,
                                  length(lats), length(lons), cut, df,
                                  ncores = 4)
  worst <- max(worst, rel_diff(gm, pw))
}
stopifnot(worst < 1e-12)
cat(sprintf("1. dense lattice x 3 distances x 4 cutoffs: worst rel diff %.2e\n", worst))

# ---- 2. Sparse occupancy, duplicates in a cell, multi-block, determinism ----
set.seed(4)
keep <- c(sample(n, 1500), sample(n, 600))    # 600 duplicated cells
t2 <- sort(rep(1:3, length.out = length(keep)))
S2 <- matrix(rnorm(length(keep) * k), length(keep), k)
pw <- fastconley:::FastSpatialMeat(gp$lat[keep], gp$lon[keep], t2, scores = S2,
                                   cutoff = 150, kernel = "uniform",
                                   dist_fn = "spherical", ncores = 4)
gm1 <- fastconley:::FastGridMeat(ring[keep], col[keep], t2, S2, lats[1], step,
                                 step, length(lats), length(lons), 150,
                                 "spherical", ncores = 1)
gm16 <- fastconley:::FastGridMeat(ring[keep], col[keep], t2, S2, lats[1], step,
                                  step, length(lats), length(lons), 150,
                                  "spherical", ncores = 16)
stopifnot(rel_diff(gm1, pw) < 1e-12, identical(gm1, gm16))
cat("2. sparse + in-cell duplicates + 3 blocks + ncores-determinism: OK\n")

# ---- 3. Degenerate shapes ----------------------------------------------------
one <- matrix(rnorm(5), 5, 1)
gm <- fastconley:::FastGridMeat(rep(0L, 5), 0:4, rep(1, 5), one,
                                10, 0.5, 0.5, 1, 5, 100, "spherical",
                                ncores = 2)
pw <- fastconley:::FastSpatialMeat(rep(10.0, 5), 0.5 * (0:4), rep(1, 5),
                                   scores = one, cutoff = 100,
                                   kernel = "uniform", dist_fn = "spherical")
stopifnot(rel_diff(gm, pw) < 1e-12)
cat("3. single ring, k = 1: OK\n")

# ---- 4. Lattice detection ----------------------------------------------------
gi <- fastconley:::detect_lonlat_grid(gp$lat[keep], gp$lon[keep])
stopifnot(!is.null(gi), abs(gi$dlat - step) < 1e-9, abs(gi$dlon - step) < 1e-9)
set.seed(5)
stopifnot(is.null(fastconley:::detect_lonlat_grid(runif(500, 0, 10),
                                                  runif(500, 0, 10))))
cat("4. lattice detection (accepts raster subset, rejects jittered): OK\n")

# ---- 5. vcovSpHAC end-to-end: method grid == pairwise, auto picks sanely ----
set.seed(6)
d <- data.frame(lat = gp$lat, lon = gp$lon,
                x1 = rnorm(n), x2 = rnorm(n))
d$y <- 0.3 * d$x1 - 0.1 * d$x2 + rnorm(n)
fit <- felm(y ~ x1 + x2, data = d, keepCX = TRUE)
v_pw <- vcovSpHAC(fit, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 200, ncores = 4,
                  method = "pairwise", data = d)
v_gr <- vcovSpHAC(fit, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 200, ncores = 4,
                  method = "grid", data = d)
v_au <- vcovSpHAC(fit, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 200, ncores = 4,
                  method = "auto", data = d)
stopifnot(rel_diff(v_gr, v_pw) < 1e-12, rel_diff(v_au, v_pw) < 1e-12)
cat("5. vcovSpHAC grid == pairwise == auto: OK\n")

# ---- 6. Guard rails -----------------------------------------------------------
d2 <- d; d2$lat <- d2$lat + runif(n, 0, 1e-3)   # break the lattice
fit2 <- felm(y ~ x1 + x2, data = d2, keepCX = TRUE)
err2 <- tryCatch({
  vcovSpHAC(fit2, lat = "lat", lon = "lon", kernel = "uniform",
            dist_fn = "spherical", dist_cutoff = 200, method = "grid", data = d2)
  FALSE
}, error = function(e) grepl("lattice", conditionMessage(e)))
v2 <- vcovSpHAC(fit2, lat = "lat", lon = "lon", kernel = "uniform",
                dist_fn = "spherical", dist_cutoff = 200, method = "auto",
                data = d2)   # auto falls back to pairwise silently
stopifnot(err2, is.matrix(v2))
cat("6. guard rails (no-lattice error, auto fallback): OK\n")

# ---- 7. Bartlett ring-FFT vs pairwise ----------------------------------------
# spherical/chord tolerances allow the inherent acos/sqrt conditioning of
# near-zero distances (adjacent cells); haversine's atan2 form is exact.
worst_b <- c(spherical = 0, haversine = 0, chord = 0)
for (df in c("spherical", "haversine", "chord")) for (cut in c(30, 120, 400, 5000)) {
  pw <- fastconley:::FastSpatialMeat(gp$lat, gp$lon, rep(1, n), scores = S,
                                     cutoff = cut, kernel = "bartlett",
                                     dist_fn = df, ncores = 4)
  gm <- fastconley:::FastGridMeat(ring, col, rep(1, n), S, lats[1], step, step,
                                  length(lats), length(lons), cut, df,
                                  kernel = "bartlett", ncores = 4)
  worst_b[df] <- max(worst_b[df], rel_diff(gm, pw))
}
stopifnot(worst_b["haversine"] < 1e-12, worst_b["spherical"] < 1e-9,
          worst_b["chord"] < 1e-9)
cat(sprintf(paste0("7. bartlett FFT x 3 distances x 4 cutoffs: worst rel diff ",
                   "sph %.1e / hav %.1e / chd %.1e\n"),
            worst_b["spherical"], worst_b["haversine"], worst_b["chord"]))

# bartlett: sparse + duplicates + multi-block + ncores-determinism
pwb <- fastconley:::FastSpatialMeat(gp$lat[keep], gp$lon[keep], t2, scores = S2,
                                    cutoff = 150, kernel = "bartlett",
                                    dist_fn = "haversine", ncores = 4)
gb1 <- fastconley:::FastGridMeat(ring[keep], col[keep], t2, S2, lats[1], step,
                                 step, length(lats), length(lons), 150,
                                 "haversine", kernel = "bartlett", ncores = 1)
gb16 <- fastconley:::FastGridMeat(ring[keep], col[keep], t2, S2, lats[1], step,
                                  step, length(lats), length(lons), 150,
                                  "haversine", kernel = "bartlett", ncores = 16)
stopifnot(rel_diff(gb1, pwb) < 1e-12, identical(gb1, gb16))
cat("8. bartlett sparse + duplicates + 3 blocks + ncores-determinism: OK\n")

# ---- 9. Dateline wrap ---------------------------------------------------------
# Full-circle lattice: windows must wrap, both kernels, vs pairwise.
set.seed(8)
wlon <- seq(0, 358, by = 2); wlat <- seq(-30, 30, by = 2)
wg <- expand.grid(lon = wlon, lat = wlat)
wn <- nrow(wg)
WS <- matrix(rnorm(wn * k), wn, k)
wgi <- fastconley:::detect_lonlat_grid(wg$lat, wg$lon)
stopifnot(wgi$n_col_full == 180L)
for (kern in c("uniform", "bartlett")) for (df in c("spherical", "haversine")) {
  pw <- fastconley:::FastSpatialMeat(wg$lat, wg$lon, rep(1, wn), scores = WS,
                                     cutoff = 600, kernel = kern,
                                     dist_fn = df, ncores = 4)
  gm <- fastconley:::FastGridMeat(wgi$ring, wgi$col, rep(1, wn), WS,
                                  wgi$lat0, wgi$dlat, wgi$dlon,
                                  wgi$n_ring, wgi$n_col,
                                  n_col_full = wgi$n_col_full,
                                  cutoff = 600, dist_fn = df, kernel = kern,
                                  ncores = 4)
  stopifnot(rel_diff(gm, pw) < if (df == "spherical" && kern == "bartlett")
    1e-9 else 1e-12)
}
# Wrap needed but dlon does not tile 360 (n_col_full = 0): informative stop,
# and method = "auto" falls back to the pairwise engine.
err3 <- tryCatch({
  fastconley:::FastGridMeat(wgi$ring, wgi$col, rep(1, wn), WS,
                            wgi$lat0, wgi$dlat, wgi$dlon,
                            wgi$n_ring, wgi$n_col, n_col_full = 0L,
                            cutoff = 600, dist_fn = "spherical",
                            kernel = "uniform", ncores = 4)
  FALSE
}, error = function(e) grepl("dateline", conditionMessage(e)))
stopifnot(err3)
# End-to-end on the full-circle raster (auto must pick a correct engine).
dw <- data.frame(lat = wg$lat, lon = wg$lon, x1 = rnorm(wn), x2 = rnorm(wn))
dw$y <- 0.3 * dw$x1 - 0.1 * dw$x2 + rnorm(wn)
fitw <- felm(y ~ x1 + x2, data = dw, keepCX = TRUE)
for (kern in c("uniform", "bartlett")) {
  vp <- vcovSpHAC(fitw, lat = "lat", lon = "lon", kernel = kern,
                  dist_fn = "spherical", dist_cutoff = 600, ncores = 4,
                  method = "pairwise", data = dw)
  vg <- vcovSpHAC(fitw, lat = "lat", lon = "lon", kernel = kern,
                  dist_fn = "spherical", dist_cutoff = 600, ncores = 4,
                  method = "grid", data = dw)
  stopifnot(rel_diff(vg, vp) < 1e-9)
}
cat("9. dateline wrap (both kernels vs pairwise, stop + auto fallback, end-to-end): OK\n")

cat("\nAll grid-meat checks passed.\n")
