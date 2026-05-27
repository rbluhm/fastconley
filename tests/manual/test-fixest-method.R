# Smoke test for vcovSpHAC.fixest(). Compares fastconley's fixest path against
# fixest::vcov_conley() for the no-FE cross-section case (where both tools
# compute the same object), and against the felm path for an FE-with-panel
# case (where the textbook within-period spatial path should produce identical
# numbers regardless of which model class supplied the residuals + cX).

suppressMessages({
  library(lfe); library(fixest); library(fastconley); library(data.table)
})

set.seed(42)
n <- 5000L; k <- 4L
lat <- runif(n, 25, 50); lon <- runif(n, -125, -70)
X <- matrix(rnorm(n * k), n, k); colnames(X) <- paste0("x", seq_len(k))
y <- as.numeric(X %*% rep(0.1, k)) + rnorm(n)
d <- data.table(y = y, X, lat = lat, lon = lon)

cat("=== Cross-section, no FE: vcovSpHAC.fixest vs vcovSpHAC.felm ===\n")

rhs <- paste(paste0("x", seq_len(k)), collapse = " + ")
fit_lfe <- felm(as.formula(paste("y ~", rhs, "| 0")), data = d, keepCX = TRUE)
fit_fx  <- feols(as.formula(paste("y ~", rhs)), data = d,
                 ssc = ssc(adj = FALSE, cluster.adj = FALSE),
                 demeaned = TRUE)

V_lfe <- vcovSpHAC(fit_lfe,
                   lat = "lat", lon = "lon",
                   kernel = "uniform", dist_fn = "spherical",
                   dist_cutoff = 500, lag_cutoff = 0,
                   balanced_pnl = FALSE, ncores = 1L, pixel = 0)

V_fc <- vcovSpHAC(fit_fx,
                  lat = "lat", lon = "lon",
                  kernel = "uniform", dist_fn = "spherical",
                  dist_cutoff = 500, lag_cutoff = 0,
                  balanced_pnl = FALSE, ncores = 1L, pixel = 0)

diff_xs <- max(abs(V_lfe - V_fc))
cat(sprintf("max |felm-path - fixest-path|: %.3e\n", diff_xs))
stopifnot(diff_xs < 1e-12)

cat("\n=== Cross-section sanity: fastconley vs fixest::vcov_conley (close, not identical) ===\n")
V_fx_conley <- vcov(fit_fx,
                    vcov = vcov_conley(lat = "lat", lon = "lon",
                                       cutoff = 500, distance = "spherical",
                                       pixel = 0),
                    ssc = ssc(adj = FALSE, cluster.adj = FALSE),
                    vcov_fix = FALSE)
common <- intersect(rownames(V_fc), rownames(V_fx_conley))
soft_gap <- max(abs(V_fc[common, common] - V_fx_conley[common, common]))
cat(sprintf("max |fastconley - fixest::vcov_conley| on slopes: %.3e\n", soft_gap))
# Tools differ by a small formulation gap; loose bound.
stopifnot(soft_gap < 1e-4)

cat("\n=== Panel, FEs: vcovSpHAC.fixest vs vcovSpHAC.felm (should be identical) ===\n")

n_unit <- 500L; n_time <- 4L
set.seed(7)
ulat <- runif(n_unit, 25, 50); ulon <- runif(n_unit, -125, -70)
unit_id <- rep(seq_len(n_unit), times = n_time)
time_id <- rep(seq_len(n_time), each = n_unit)
N <- n_unit * n_time
X2 <- matrix(rnorm(N * k), N, k); colnames(X2) <- paste0("x", seq_len(k))
y2 <- as.numeric(X2 %*% rep(0.1, k)) + rnorm(N)
dp <- data.table(y = y2, X2,
                 lat = ulat[unit_id], lon = ulon[unit_id],
                 unit = unit_id, time = time_id)

fit_fe_lfe <- felm(as.formula(paste("y ~", rhs, "| unit + time")), data = dp, keepCX = TRUE)
fit_fe_fx  <- feols(as.formula(paste("y ~", rhs, "| unit + time")), data = dp,
                    ssc = ssc(adj = FALSE, cluster.adj = FALSE),
                    demeaned = TRUE)

V_lfe <- vcovSpHAC(fit_fe_lfe,
                   unit = "unit", time = "time",
                   lat = "lat", lon = "lon",
                   kernel = "uniform", dist_fn = "spherical",
                   dist_cutoff = 500, lag_cutoff = 1L,
                   balanced_pnl = TRUE, ncores = 1L, pixel = 0)

V_fx <- vcovSpHAC(fit_fe_fx,
                  unit = "unit", time = "time",
                  lat = "lat", lon = "lon",
                  kernel = "uniform", dist_fn = "spherical",
                  dist_cutoff = 500, lag_cutoff = 1L,
                  balanced_pnl = TRUE, ncores = 1L, pixel = 0)

diff_panel <- max(abs(V_lfe - V_fx))
cat(sprintf("max |felm-path - fixest-path|: %.3e\n", diff_panel))
stopifnot(diff_panel < 1e-9)

cat("\n=== Bartlett kernel, chord distance: fixest path vs felm path ===\n")

V_lfe2 <- vcovSpHAC(fit_fe_lfe, unit = "unit", time = "time",
                    lat = "lat", lon = "lon",
                    kernel = "bartlett", dist_fn = "chord",
                    dist_cutoff = 500, lag_cutoff = 0,
                    balanced_pnl = TRUE, ncores = 1L, pixel = 0)
V_fx2 <- vcovSpHAC(fit_fe_fx, unit = "unit", time = "time",
                   lat = "lat", lon = "lon",
                   kernel = "bartlett", dist_fn = "chord",
                   dist_cutoff = 500, lag_cutoff = 0,
                   balanced_pnl = TRUE, ncores = 1L, pixel = 0)
diff_panel2 <- max(abs(V_lfe2 - V_fx2))
cat(sprintf("max |felm-path - fixest-path|: %.3e\n", diff_panel2))
stopifnot(diff_panel2 < 1e-9)

cat("\n=== Error path: feols without demeaned=TRUE ===\n")
fit_nodm <- feols(as.formula(paste("y ~", rhs)), data = d)
err <- tryCatch(
  vcovSpHAC(fit_nodm, lat = "lat", lon = "lon",
            kernel = "uniform", dist_fn = "spherical",
            dist_cutoff = 500),
  error = function(e) conditionMessage(e))
cat("error message:", err, "\n")
stopifnot(grepl("demeaned = TRUE", err))

cat("\n=== Error path: unsupported model class ===\n")
fit_lm <- lm(y ~ x1 + x2, data = d)
err2 <- tryCatch(vcovSpHAC(fit_lm, lat = "lat", lon = "lon",
                            dist_cutoff = 500),
                 error = function(e) conditionMessage(e))
cat("error message:", err2, "\n")
stopifnot(grepl("unsupported model class", err2))

cat("\n=== fixest subset= handling ===\n")
set.seed(7); n_sub <- 200L
d_sub <- data.table(y = rnorm(n_sub),
                    x1 = rnorm(n_sub), x2 = rnorm(n_sub), x3 = rnorm(n_sub),
                    lat = runif(n_sub, 25, 50), lon = runif(n_sub, -125, -70),
                    keep = rep(c(TRUE, FALSE), n_sub / 2))
fit_sub <- feols(y ~ x1 + x2 + x3, data = d_sub, subset = d_sub$keep,
                 demeaned = TRUE)
V_sub <- vcovSpHAC(fit_sub, lat = "lat", lon = "lon",
                   kernel = "uniform", dist_fn = "spherical",
                   dist_cutoff = 500, ncores = 1L)
# Should equal fitting on the pre-subset data.
fit_pre <- feols(y ~ x1 + x2 + x3, data = d_sub[keep == TRUE],
                 demeaned = TRUE)
V_pre <- vcovSpHAC(fit_pre, lat = "lat", lon = "lon",
                   kernel = "uniform", dist_fn = "spherical",
                   dist_cutoff = 500, ncores = 1L)
diff_subset <- max(abs(V_sub - V_pre))
cat(sprintf("max |subset= vs pre-subset|: %.3e\n", diff_subset))
stopifnot(diff_subset < 1e-12)

# subset= + NAs combined
d_sub$y[c(2, 7)] <- NA
fit_both <- feols(y ~ x1 + x2 + x3, data = d_sub, subset = d_sub$keep,
                  demeaned = TRUE)
V_both <- vcovSpHAC(fit_both, lat = "lat", lon = "lon",
                    kernel = "uniform", dist_fn = "spherical",
                    dist_cutoff = 500, ncores = 1L)
cat(sprintf("subset + NA: V is %d x %d (n_kept = %d)\n",
            nrow(V_both), ncol(V_both), nrow(fit_both$X_demeaned)))
stopifnot(nrow(V_both) == 4L)

cat("\nAll fixest-method smoke tests passed.\n")
