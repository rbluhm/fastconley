# GLM (feglm/fepois) and IV (2SLS) parity tests for vcovSpHAC.
#
# The "exact yardstick" used throughout is the M-estimation / 2SLS sandwich
# assembled by hand from fixest's own stored pieces — score matrix and bread
# (cov.iid for GLMs, cov.iid/sigma2 for feols) — with the meat computed by
# fastconley's FastSpatialMeat on those scores. Agreement with the yardstick
# verifies that vcovSpHAC extracts exactly the right scores and bread; the
# comparison against fixest::vcov_conley() itself is loose (~1e-2) because
# fixest's spherical distance is approximate while ours is exact.

suppressMessages({
  library(lfe); library(fixest); library(fastconley); library(data.table)
})

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b))

# Align a vcov returned for a felm IV fit (names like `x2(fit)`) with the
# fixest naming (fit_x2).
align <- function(V, ref_names) {
  nm <- gsub("`?([^`(]+)\\(fit\\)`?", "fit_\\1", rownames(V))
  i <- match(ref_names, nm)
  stopifnot(!anyNA(i))
  V[i, i]
}

yardstick <- function(fit, d, cutoff, kernel, dist_fn, lag_cutoff = 0,
                      unit = NULL, time = NULL) {
  sc <- fit$scores
  n  <- nrow(sc)
  tvec <- if (is.null(time)) rep(1, n) else as.numeric(factor(d[[time]]))
  ot <- order(tvec)   # FastSpatialMeat needs time-contiguous rows
  XeeX <- fastconley:::FastSpatialMeat(
    lat = d$lat[ot], lon = d$lon[ot], time = tvec[ot],
    scores = sc[ot, , drop = FALSE], cutoff = cutoff, kernel = kernel,
    dist_fn = dist_fn, ncores = 1)
  if (lag_cutoff > 0) {
    o <- order(d[[unit]], d[[time]])
    XeeX <- XeeX + fastconley:::FastSerialHacPanel(
      unit = as.numeric(factor(d[[unit]]))[o],
      time = as.numeric(factor(d[[time]]))[o],
      cutoff = lag_cutoff, scores = sc[o, , drop = FALSE], ncores = 1)
  }
  B <- if (fit$method_type == "feols") fit$cov.iid / fit$sigma2 else fit$cov.iid
  (B %*% XeeX %*% B) * n / (n - fit$nparams)
}

set.seed(42)

## ------------------------------------------------------------------
cat("=== A. fepois cross-section, no FE ===\n")
n <- 4000
d <- data.table(lat = runif(n, 25, 50), lon = runif(n, -125, -70),
                x = rnorm(n), z = rnorm(n))
d$y <- rpois(n, exp(0.3 * d$x - 0.2 * d$z))

fp <- fepois(y ~ x + z, data = d, glm.tol = 1e-12)

for (kern in c("uniform", "bartlett")) {
  for (df in c("haversine", "spherical", "chord")) {
    V  <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = kern,
                    dist_fn = df, dist_cutoff = 300, ncores = 2, data = d)
    Vy <- yardstick(fp, d, 300, kern, df)
    rd <- rel_diff(V, Vy)
    cat(sprintf("  %-8s %-10s vs yardstick: %.2e\n", kern, df, rd))
    stopifnot(rd < 1e-12)
  }
}

V_u <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = "uniform",
                 dist_fn = "spherical", dist_cutoff = 300, ncores = 2, data = d)
V_fx <- vcov(fp, vcov_conley(lat = "lat", lon = "lon", cutoff = 300,
                             distance = "spherical"))
cat(sprintf("  vs fixest::vcov_conley (approx distance): %.2e\n",
            rel_diff(V_u, V_fx)))
stopifnot(rel_diff(V_u, V_fx) < 2.5e-2)

# ssc = FALSE: exact n/(n-K) ratio off the default
V_nossc <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = "uniform",
                     dist_fn = "spherical", dist_cutoff = 300, ncores = 2,
                     ssc = FALSE, data = d)
stopifnot(rel_diff(V_u, V_nossc * n / (n - fp$nparams)) < 1e-14)
cat("  ssc = FALSE scales by exactly n/(n-K): ok\n")

# bit-identical across core counts
V_c1 <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = "bartlett",
                  dist_fn = "haversine", dist_cutoff = 300, ncores = 1, data = d)
V_c4 <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = "bartlett",
                  dist_fn = "haversine", dist_cutoff = 300, ncores = 4, data = d)
stopifnot(identical(V_c1, V_c4))
cat("  bit-identical across ncores 1 vs 4: ok\n")

## ------------------------------------------------------------------
cat("=== B. feglm families, FE, weights, offset ===\n")
d$id <- factor(sample(1:60, n, TRUE))
d$w  <- runif(n, 0.5, 2)
d$off <- runif(n, 0, 0.5)
d$yb <- as.integer(runif(n) < plogis(0.5 * d$x))
d$yp <- rpois(n, exp(0.3 * d$x + 0.1 * (as.integer(d$id) %% 5)))

fits <- list(
  "poisson + FE"        = fepois(yp ~ x + z | id, data = d, glm.tol = 1e-12),
  "logit + FE"          = feglm(yb ~ x + z | id, data = d, family = binomial(),
                                glm.tol = 1e-12),
  "probit"              = feglm(yb ~ x + z, data = d,
                                family = binomial("probit"), glm.tol = 1e-12),
  "weighted poisson"    = fepois(yp ~ x + z | id, data = d, weights = ~w,
                                 glm.tol = 1e-12),
  "poisson + offset"    = fepois(yp ~ x + z | id, data = d, offset = ~off,
                                 glm.tol = 1e-12)
)
for (nmf in names(fits)) {
  f <- fits[[nmf]]
  V  <- vcovSpHAC(f, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 300, ncores = 2, data = d)
  Vy <- yardstick(f, d, 300, "uniform", "spherical")
  rd <- rel_diff(V, Vy)
  cat(sprintf("  %-18s vs yardstick: %.2e\n", nmf, rd))
  stopifnot(rd < 1e-12)
  V_fx <- vcov(f, vcov_conley(lat = "lat", lon = "lon", cutoff = 300,
                              distance = "spherical"))
  stopifnot(rel_diff(V, V_fx) < 2.5e-2)
}
cat("  all families also within 2.5e-2 of fixest::vcov_conley: ok\n")

## ------------------------------------------------------------------
cat("=== C. fepois panel: spatial + serial HAC ===\n")
n_unit <- 400; n_time <- 6
N <- n_unit * n_time
dp <- data.table(unit = rep(1:n_unit, each = n_time),
                 time = rep(1:n_time, times = n_unit))
ulat <- runif(n_unit, 25, 50); ulon <- runif(n_unit, -125, -70)
dp$lat <- ulat[dp$unit]; dp$lon <- ulon[dp$unit]
dp$x <- rnorm(N)
dp$yp <- rpois(N, exp(0.3 * dp$x + 0.05 * (dp$unit %% 7)))

fpp <- fepois(yp ~ x | unit + time, data = dp, glm.tol = 1e-12)

V_sp  <- vcovSpHAC(fpp, unit = "unit", time = "time", lat = "lat", lon = "lon",
                   kernel = "bartlett", dist_fn = "haversine",
                   dist_cutoff = 400, lag_cutoff = 0, ncores = 2, data = dp)
V_ser <- vcovSpHAC(fpp, unit = "unit", time = "time", lat = "lat", lon = "lon",
                   kernel = "bartlett", dist_fn = "haversine",
                   dist_cutoff = 400, lag_cutoff = 2, ncores = 2, data = dp)
Vy_sp  <- yardstick(fpp, dp, 400, "bartlett", "haversine",
                    unit = "unit", time = "time")
Vy_ser <- yardstick(fpp, dp, 400, "bartlett", "haversine", lag_cutoff = 2,
                    unit = "unit", time = "time")
cat(sprintf("  spatial-only vs yardstick:   %.2e\n", rel_diff(V_sp, Vy_sp)))
cat(sprintf("  spatial+serial vs yardstick: %.2e\n", rel_diff(V_ser, Vy_ser)))
stopifnot(rel_diff(V_sp, Vy_sp) < 1e-12, rel_diff(V_ser, Vy_ser) < 1e-12)

V_bal <- vcovSpHAC(fpp, unit = "unit", time = "time", lat = "lat", lon = "lon",
                   kernel = "bartlett", dist_fn = "haversine",
                   dist_cutoff = 400, lag_cutoff = 2, balanced_pnl = TRUE,
                   ncores = 2, data = dp)
cat(sprintf("  balanced_pnl TRUE vs FALSE:  %.2e\n", rel_diff(V_bal, V_ser)))
stopifnot(rel_diff(V_bal, V_ser) < 1e-10)

## ------------------------------------------------------------------
cat("=== D. IV / 2SLS: felm and feols ===\n")
set.seed(7)
n <- 4000
di <- data.table(lat = runif(n, 25, 50), lon = runif(n, -125, -70),
                 id = factor(sample(1:60, n, TRUE)), x1 = rnorm(n),
                 zi1 = rnorm(n), zi2 = rnorm(n), w = runif(n, 0.5, 2))
u <- rnorm(n)
di$x2 <- 0.8 * di$zi1 + 0.5 * u + rnorm(n)
di$x3 <- 0.6 * di$zi2 - 0.3 * u + rnorm(n)
di$y  <- 0.3 * di$x1 + 0.7 * di$x2 - 0.4 * di$x3 + u + rnorm(n)

cases <- list(
  "one endog"     = list(fl = y ~ x1 | id | (x2 ~ zi1),
                         fx = y ~ x1 | id | x2 ~ zi1, w = FALSE),
  "two endog"     = list(fl = y ~ x1 | id | (x2 | x3 ~ zi1 + zi2),
                         fx = y ~ x1 | id | x2 + x3 ~ zi1 + zi2, w = FALSE),
  "weighted IV"   = list(fl = y ~ x1 | id | (x2 ~ zi1),
                         fx = y ~ x1 | id | x2 ~ zi1, w = TRUE)
)
for (nmc in names(cases)) {
  cs <- cases[[nmc]]
  fl <- if (cs$w) felm(cs$fl, data = di, weights = di$w, keepCX = TRUE)
        else      felm(cs$fl, data = di, keepCX = TRUE)
  fx <- if (cs$w) feols(cs$fx, data = di, weights = ~w, demeaned = TRUE)
        else      feols(cs$fx, data = di, demeaned = TRUE)
  Vl <- vcovSpHAC(fl, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 300, ncores = 2,
                  data = di)
  Vf <- vcovSpHAC(fx, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 300, ncores = 2,
                  data = di)
  Vy <- yardstick(fx, di, 300, "uniform", "spherical")
  Vl_a <- align(Vl, rownames(Vf))
  cat(sprintf("  %-12s felm vs feols: %.2e | feols vs yardstick: %.2e\n",
              nmc, rel_diff(Vl_a, Vf), rel_diff(Vf, Vy)))
  stopifnot(rel_diff(Vl_a, Vf) < 1e-10, rel_diff(Vf, Vy) < 1e-12)
  V_fx <- vcov(fx, vcov_conley(lat = "lat", lon = "lon", cutoff = 300,
                               distance = "spherical"))
  stopifnot(rel_diff(Vf, V_fx) < 2.5e-2)
}
cat("  all IV cases also within 2.5e-2 of fixest::vcov_conley: ok\n")

## IV panel with serial HAC: felm vs feols through our paths
set.seed(11)
n_unit <- 300; n_time <- 5; N <- n_unit * n_time
dip <- data.table(unit = rep(1:n_unit, each = n_time),
                  time = rep(1:n_time, times = n_unit))
ulat <- runif(n_unit, 25, 50); ulon <- runif(n_unit, -125, -70)
dip$lat <- ulat[dip$unit]; dip$lon <- ulon[dip$unit]
dip$x1 <- rnorm(N); dip$zi <- rnorm(N)
u <- rnorm(N)
dip$x2 <- 0.8 * dip$zi + 0.5 * u + rnorm(N)
dip$y  <- 0.3 * dip$x1 + 0.7 * dip$x2 + u + rnorm(N)

fl <- felm(y ~ x1 | unit + time | (x2 ~ zi), data = dip, keepCX = TRUE)
fx <- feols(y ~ x1 | unit + time | x2 ~ zi, data = dip, demeaned = TRUE)
Vl <- vcovSpHAC(fl, unit = "unit", time = "time", lat = "lat", lon = "lon",
                kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 400,
                lag_cutoff = 2, balanced_pnl = TRUE, ncores = 2, data = dip)
Vf <- vcovSpHAC(fx, unit = "unit", time = "time", lat = "lat", lon = "lon",
                kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 400,
                lag_cutoff = 2, balanced_pnl = TRUE, ncores = 2, data = dip)
rd <- rel_diff(align(Vl, rownames(Vf)), Vf)
cat(sprintf("  IV panel + serial HAC, felm vs feols: %.2e\n", rd))
stopifnot(rd < 1e-10)

## ------------------------------------------------------------------
cat("=== E. error paths ===\n")
fp_lean <- fepois(y ~ x + z, data = d, lean = TRUE)
err <- tryCatch(vcovSpHAC(fp_lean, lat = "lat", lon = "lon",
                          kernel = "uniform", dist_fn = "spherical",
                          dist_cutoff = 300, data = d),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("lean", err))
cat("  lean = TRUE fepois rejected with a 'lean' message: ok\n")

fm <- femlm(yp ~ x | id, data = d)
err <- tryCatch(vcovSpHAC(fm, lat = "lat", lon = "lon", kernel = "uniform",
                          dist_fn = "spherical", dist_cutoff = 300, data = d),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("method_type", err))
cat("  femlm rejected with a method_type message: ok\n")

## ------------------------------------------------------------------
cat("=== F. engines on GLM scores: grid vs pairwise, pixel ===\n")
cells <- as.data.table(expand.grid(lat = seq(35, 54.75, by = 0.25),
                                   lon = seq(-20, 19.75, by = 0.25)))
cells$x <- rnorm(nrow(cells))
cells$cnt <- rpois(nrow(cells), exp(0.4 * cells$x))
fg <- fepois(cnt ~ x, data = cells, glm.tol = 1e-12)
for (kern in c("uniform", "bartlett")) {
  Vp <- vcovSpHAC(fg, lat = "lat", lon = "lon", kernel = kern,
                  dist_fn = "spherical", dist_cutoff = 300,
                  method = "pairwise", ncores = 2, data = cells)
  Vg <- vcovSpHAC(fg, lat = "lat", lon = "lon", kernel = kern,
                  dist_fn = "spherical", dist_cutoff = 300,
                  method = "grid", ncores = 2, data = cells)
  rd <- rel_diff(Vg, Vp)
  cat(sprintf("  %-8s grid vs pairwise: %.2e\n", kern, rd))
  stopifnot(rd < 1e-9)
}
V_px <- vcovSpHAC(fp, lat = "lat", lon = "lon", kernel = "uniform",
                  dist_fn = "spherical", dist_cutoff = 300, pixel = 20,
                  ncores = 2, data = d)
stopifnot(all(is.finite(V_px)), all(eigen(V_px)$values > 0))
cat("  pixel = 20 on fepois scores: finite, PSD: ok\n")

## ------------------------------------------------------------------
cat("=== G. obs_selection on the GLM path: NA drops, separated FE groups,\n")
cat("       subset =, and data auto-recovery ===\n")
# feglm drops rows fixest-side: NAs, and FE levels whose outcome is all
# zero (fixef.rm = "perfect"). Coordinates must be re-aligned through
# obs_selection ($subset first, then $obsRemoved); a misalignment would
# show up against a yardstick built on the kept rows.
set.seed(99)
ng <- 3000
dg <- data.table(lat = runif(ng, 25, 50), lon = runif(ng, -125, -70),
                 id = factor(sample(1:40, ng, TRUE)), x = rnorm(ng))
dg$y <- rpois(ng, exp(0.4 * dg$x))
dg$y[dg$id %in% c("3", "7")] <- 0      # separated FE groups -> rows dropped
dg$x[sample(ng, 50)] <- NA             # NA rows -> dropped

kept_rows <- function(fit, n_orig) {
  keep <- seq_len(n_orig)
  if (length(fit$obs_selection$subset))     keep <- keep[fit$obs_selection$subset]
  if (length(fit$obs_selection$obsRemoved)) keep <- keep[fit$obs_selection$obsRemoved]
  keep
}

fg2 <- fepois(y ~ x | id, data = dg, glm.tol = 1e-12)
stopifnot(fg2$nobs < ng)               # something was actually dropped
dk <- dg[kept_rows(fg2, ng)]
stopifnot(nrow(dk) == fg2$nobs)
V  <- vcovSpHAC(fg2, lat = "lat", lon = "lon", kernel = "uniform",
                dist_fn = "spherical", dist_cutoff = 300, ncores = 2, data = dg)
Vy <- yardstick(fg2, dk, 300, "uniform", "spherical")
cat(sprintf("  NA + separated-FE drops (%d -> %d rows) vs yardstick: %.2e\n",
            ng, fg2$nobs, rel_diff(V, Vy)))
stopifnot(rel_diff(V, Vy) < 1e-12)

# subset = stacks on top of the NA/separation drops
fg3 <- fepois(y ~ x | id, data = dg, subset = dg$lon < -80, glm.tol = 1e-12)
dk3 <- dg[kept_rows(fg3, ng)]
stopifnot(nrow(dk3) == fg3$nobs)
V3  <- vcovSpHAC(fg3, lat = "lat", lon = "lon", kernel = "uniform",
                 dist_fn = "spherical", dist_cutoff = 300, ncores = 2, data = dg)
Vy3 <- yardstick(fg3, dk3, 300, "uniform", "spherical")
cat(sprintf("  subset = on top of drops (%d rows) vs yardstick: %.2e\n",
            fg3$nobs, rel_diff(V3, Vy3)))
stopifnot(rel_diff(V3, Vy3) < 1e-12)

# data auto-recovery from the fit's call (no data = passed)
V_auto <- vcovSpHAC(fg2, lat = "lat", lon = "lon", kernel = "uniform",
                    dist_fn = "spherical", dist_cutoff = 300, ncores = 2)
stopifnot(identical(V, V_auto))
cat("  data auto-recovery identical to explicit data =: ok\n")

cat("\nAll GLM/IV parity tests passed.\n")
