# Compares the outputs of tests/manual/core_standalone_check.cpp (the engine
# header compiled without R) with the installed package's entry points on
# the same inputs. Expect identical() everywhere: same engine, same
# summation order.
#   Rscript tests/manual/core_standalone_check.R OUT_DIR
args <- commandArgs(trailingOnly = TRUE)
d <- if (length(args)) args[1] else stop("pass OUT_DIR")
library(fastconley)
rd <- function(f) unname(as.matrix(read.csv(file.path(d, f), header = FALSE)))
FastSpatialMeat <- fastconley:::FastSpatialMeat
FastSerialHacPanel <- fastconley:::FastSerialHacPanel
FastGridMeat <- fastconley:::FastGridMeat

cs <- rd("cs_input.csv"); lat <- cs[, 1]; lon <- cs[, 2]; S <- cs[, 3:5]
one <- rep(1, nrow(cs))
chk <- list()
for (nc in c(1L, 4L)) {
  chk[[sprintf("cs bartlett/haversine nc%d", nc)]] <- list(
    rd(sprintf("cs_bartlett_haversine_nc%d.csv", nc)),
    FastSpatialMeat(lat, lon, one, scores = S, cutoff = 300, ncores = nc))
  chk[[sprintf("cs uniform/chord band nc%d", nc)]] <- list(
    rd(sprintf("cs_uniform_chord_band_nc%d.csv", nc)),
    FastSpatialMeat(lat, lon, one, scores = S, cutoff = 300, kernel = "uniform",
                    dist_fn = "chord", neighbor = "band", ncores = nc))
}
p <- rd("pnl_input.csv")
chk[["pnl balanced bartlett/spherical nc4"]] <- list(
  rd("pnl_bartlett_spherical_balanced_nc4.csv"),
  FastSpatialMeat(p[, 1], p[, 2], p[, 3], scores = p[, 4:6], cutoff = 400,
                  dist_fn = "spherical", balanced_pnl = TRUE, ncores = 4L))
n_u <- 800L; T_ <- 3L
ord <- order(rep(seq_len(n_u), times = T_), rep(seq_len(T_), each = n_u))
chk[["pnl serial lag2 nc4"]] <- list(
  rd("pnl_serial_lag2_nc4.csv"),
  FastSerialHacPanel(rep(seq_len(n_u), times = T_)[ord], rep(seq_len(T_), each = n_u)[ord],
                     cutoff = 2, scores = p[ord, 4:6], ncores = 4L))
g <- rd("grid_input.csv"); R <- 40L; C <- 60L
ring <- rep(0:(R - 1), each = C); col <- rep(0:(C - 1), times = R)
chk[["grid uniform/spherical nc4"]] <- list(
  rd("grid_uniform_spherical_nc4.csv"),
  FastGridMeat(ring, col, rep(1, R * C), g, 35, 0.5, 0.5, R, C, 250,
               dist_fn = "spherical", kernel = "uniform", ncores = 4L))
chk[["grid bartlett/haversine nc4"]] <- list(
  rd("grid_bartlett_haversine_nc4.csv"),
  FastGridMeat(ring, col, rep(1, R * C), g, 35, 0.5, 0.5, R, C, 250,
               dist_fn = "haversine", kernel = "bartlett", ncores = 4L))
bad <- 0L
for (nm in names(chk)) {
  a <- chk[[nm]][[1]]; b <- unname(chk[[nm]][[2]])
  same <- identical(a, b)
  if (!same) bad <- bad + 1L
  cat(sprintf("%-40s %s  max|d| = %.3g\n", nm, if (same) "identical" else "DIFF",
              max(abs(a - b))))
}
quit(status = if (bad) 1L else 0L)
