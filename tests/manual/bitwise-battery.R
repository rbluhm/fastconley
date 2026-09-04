# Bitwise regression battery for refactors that must not change numbers.
#
#   Rscript tests/manual/bitwise-battery.R run OUT.rds [LIB_LOC]
#   Rscript tests/manual/bitwise-battery.R compare BASELINE.rds NEW.rds
#
# `run` fits a fixed set of lfe/fixest models on seeded synthetic data and
# stores the vcovSpHAC() matrix for every configuration below (pairwise
# grid/band, balanced CSR incl. float weights, unbalanced streaming, serial
# HAC, pixel aggregation, raster grid engine incl. dateline wrap, weights,
# fepois, ncores 1 vs 4). `compare` reports identical() per config and the
# max abs difference otherwise. Pass LIB_LOC to run a baseline install
# (e.g. the previous release installed into a scratch library).

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

if (mode == "compare") {
  a <- readRDS(args[2]); b <- readRDS(args[3])
  stopifnot(identical(names(a), names(b)))
  bad <- 0L
  for (nm in names(a)) {
    same <- identical(a[[nm]], b[[nm]])
    if (!same) {
      bad <- bad + 1L
      d <- if (is.matrix(a[[nm]]) && is.matrix(b[[nm]]) &&
               all(dim(a[[nm]]) == dim(b[[nm]]))) max(abs(a[[nm]] - b[[nm]])) else NA
      cat(sprintf("DIFF  %-48s max|d| = %s\n", nm, format(d)))
    }
  }
  cat(sprintf("%d configs, %d identical, %d differ\n", length(a), length(a) - bad, bad))
  quit(status = if (bad > 0L) 1L else 0L)
}

stopifnot(mode == "run")
out_file <- args[2]
if (length(args) >= 3 && nzchar(args[3])) {
  library(fastconley, lib.loc = args[3])
} else {
  library(fastconley)
}
cat("fastconley", as.character(packageVersion("fastconley")), "from",
    dirname(system.file(package = "fastconley")), "\n")
suppressPackageStartupMessages({ library(lfe); library(fixest) })

set.seed(20260903)
mk_xy <- function(n) {
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- runif(n)
  y <- 0.5 * x1 - 0.3 * x2 + x3 + rnorm(n)
  data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
}
cs_box <- function(n) {
  cbind(mk_xy(n), lat = runif(n, 25, 49), lon = runif(n, -124, -67),
        region = sample(letters[1:5], n, TRUE))
}
cs  <- cs_box(6000)
csb <- cs_box(20000)

n_u <- 1500; T_ <- 4
bal <- data.frame(unit = rep(seq_len(n_u), each = T_), time = rep(seq_len(T_), n_u),
                  lat = rep(runif(n_u, 35, 60), each = T_),
                  lon = rep(runif(n_u, -10, 30), each = T_))
bal <- cbind(bal, mk_xy(nrow(bal)))
unbal <- bal[-sample(nrow(bal), round(0.2 * nrow(bal))), ]

ras <- expand.grid(lat = 35 + 0.5 * (0:59), lon = -20 + 0.5 * (0:79))
ras <- ras[-sample(nrow(ras), round(0.25 * nrow(ras))), ]
ras <- cbind(ras, mk_xy(nrow(ras)))
wrap <- expand.grid(lat = 60 + (0:11), lon = -180 + (0:359))
wrap <- cbind(wrap, mk_xy(nrow(wrap)))

cs$w <- runif(nrow(cs), 0.5, 2)
cs$cnt <- rpois(nrow(cs), exp(0.3 * cs$x1))

f_cs   <- felm(y ~ x1 + x2 + x3, data = cs, keepCX = TRUE)
f_csb  <- felm(y ~ x1 + x2 + x3, data = csb, keepCX = TRUE)
f_bal  <- felm(y ~ x1 + x2 | unit + time, data = bal, keepCX = TRUE)
f_unb  <- felm(y ~ x1 + x2 | unit + time, data = unbal, keepCX = TRUE)
f_ras  <- felm(y ~ x1 + x2, data = ras, keepCX = TRUE)
f_wrap <- felm(y ~ x1 + x2, data = wrap, keepCX = TRUE)
f_fx_w <- feols(y ~ x1 + x2 | region, data = cs, weights = ~w, demeaned = TRUE)
f_fx_p <- fepois(cnt ~ x1 + x2 | region, data = cs)

res <- list()
add <- function(name, expr) {
  V <- suppressWarnings(suppressMessages(expr))
  res[[name]] <<- unname(V)
}
kd <- expand.grid(kernel = c("bartlett", "uniform"),
                  dist = c("haversine", "spherical", "chord"),
                  stringsAsFactors = FALSE)

for (i in seq_len(nrow(kd))) {
  k <- kd$kernel[i]; d <- kd$dist[i]; tag <- paste(k, d, sep = "/")
  for (nb in c("grid", "band")) {
    add(sprintf("cs %s %s nc1", tag, nb),
        vcovSpHAC(f_cs, lat = "lat", lon = "lon", kernel = k, dist_fn = d,
                  dist_cutoff = 300, ncores = 1, neighbor = nb, data = cs))
  }
  add(sprintf("cs %s grid nc4", tag),
      vcovSpHAC(f_cs, lat = "lat", lon = "lon", kernel = k, dist_fn = d,
                dist_cutoff = 300, ncores = 4, data = cs))
  for (lag in c(0, 2)) {
    add(sprintf("bal %s lag%d nc1", tag, lag),
        vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
                  kernel = k, dist_fn = d, dist_cutoff = 300, lag_cutoff = lag,
                  balanced_pnl = TRUE, ncores = 1, data = bal))
  }
  add(sprintf("unbal %s lag2 nc1", tag),
      vcovSpHAC(f_unb, unit = "unit", time = "time", lat = "lat", lon = "lon",
                kernel = k, dist_fn = d, dist_cutoff = 300, lag_cutoff = 2,
                balanced_pnl = FALSE, ncores = 1, data = unbal))
  add(sprintf("raster-grid %s nc1", tag),
      vcovSpHAC(f_ras, lat = "lat", lon = "lon", kernel = k, dist_fn = d,
                dist_cutoff = 250, ncores = 1, method = "grid", data = ras))
}

add("cs bartlett/haversine pixel25 nc1",
    vcovSpHAC(f_cs, lat = "lat", lon = "lon", dist_cutoff = 300, pixel = 25,
              ncores = 1, data = cs))
add("csb bartlett/haversine grid nc1",
    vcovSpHAC(f_csb, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = csb))
add("csb bartlett/haversine grid nc4",
    vcovSpHAC(f_csb, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 4, data = csb))
add("csb uniform/spherical grid nc4",
    vcovSpHAC(f_csb, lat = "lat", lon = "lon", kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 300, ncores = 4, data = csb))
add("csb bartlett/chord band nc4",
    vcovSpHAC(f_csb, lat = "lat", lon = "lon", dist_fn = "chord", dist_cutoff = 300,
              ncores = 4, neighbor = "band", data = csb))
add("bal bartlett/haversine lag2 float nc1",
    vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
              dist_cutoff = 300, lag_cutoff = 2, balanced_pnl = TRUE,
              csr_weight = "float", ncores = 1, data = bal))
add("bal bartlett/haversine lag2 nc4",
    vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
              dist_cutoff = 300, lag_cutoff = 2, balanced_pnl = TRUE, ncores = 4, data = bal))
add("bal uniform/chord lag0 nc4",
    vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "chord", dist_cutoff = 300,
              balanced_pnl = TRUE, ncores = 4, data = bal))
add("bal bartlett/spherical lag1 band nc1",
    vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
              dist_fn = "spherical", dist_cutoff = 300, lag_cutoff = 1,
              balanced_pnl = TRUE, neighbor = "band", ncores = 1, data = bal))
add("bal-as-general bartlett/haversine lag2 nc1",
    vcovSpHAC(f_bal, unit = "unit", time = "time", lat = "lat", lon = "lon",
              dist_cutoff = 300, lag_cutoff = 2, balanced_pnl = FALSE, ncores = 1, data = bal))
add("unbal bartlett/haversine lag2 nc4",
    vcovSpHAC(f_unb, unit = "unit", time = "time", lat = "lat", lon = "lon",
              dist_cutoff = 300, lag_cutoff = 2, ncores = 4, data = unbal))
add("raster-grid uniform/spherical nc4",
    vcovSpHAC(f_ras, lat = "lat", lon = "lon", kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 250, ncores = 4, method = "grid", data = ras))
add("raster-grid bartlett/spherical nc4",
    vcovSpHAC(f_ras, lat = "lat", lon = "lon", dist_fn = "spherical",
              dist_cutoff = 250, ncores = 4, method = "grid", data = ras))
add("raster-pairwise uniform/spherical nc1",
    vcovSpHAC(f_ras, lat = "lat", lon = "lon", kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 250, ncores = 1, method = "pairwise", data = ras))
add("wrap-grid uniform/spherical nc1",
    vcovSpHAC(f_wrap, lat = "lat", lon = "lon", kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 1500, ncores = 1, method = "grid", data = wrap))
add("wrap-grid bartlett/spherical nc4",
    vcovSpHAC(f_wrap, lat = "lat", lon = "lon", dist_fn = "spherical",
              dist_cutoff = 1500, ncores = 4, method = "grid", data = wrap))
add("fixest weighted FE bartlett/haversine nc1",
    vcovSpHAC(f_fx_w, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = cs))
add("fixest fepois FE bartlett/haversine nc1",
    vcovSpHAC(f_fx_p, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1, data = cs))

saveRDS(res, out_file)
cat(length(res), "configs saved to", out_file, "\n")
