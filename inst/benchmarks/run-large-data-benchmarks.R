#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fastconley)
})
for (pkg in c("fixest", "lfe")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this benchmark script.")
  }
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- if (!is.na(script_file)) {
  normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
out_dir <- file.path(root, "inst", "benchmarks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

threads_requested <- as.integer(strsplit(Sys.getenv("FASTCONLEY_LARGE_THREADS", "2,4,8,16"), ",")[[1]])
threads_requested <- threads_requested[!is.na(threads_requested) & threads_requested > 0]
threads_requested <- unique(pmin(threads_requested, parallel::detectCores(logical = TRUE)))
bench_reps <- max(1L, as.integer(Sys.getenv("FASTCONLEY_LARGE_REPS", "1")))
run_fixest_large <- identical(Sys.getenv("FASTCONLEY_LARGE_FIXEST", "0"), "1")

fmt_bytes <- function(bytes) {
  units <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB")
  x <- as.numeric(bytes); i <- 1L
  while (x >= 1024 && i < length(units)) { x <- x / 1024; i <- i + 1L }
  sprintf("%.1f %s", x, units[i])
}

machine_lines <- function() {
  c(
    "# CPU",
    tryCatch(system("lscpu", intern = TRUE), error = function(e) paste("lscpu unavailable:", conditionMessage(e))),
    "", "# Memory",
    tryCatch(system("free -h", intern = TRUE), error = function(e) paste("free unavailable:", conditionMessage(e))),
    "", "# R session",
    capture.output(sessionInfo())
  )
}

time_call <- function(f, reps = bench_reps) {
  times <- numeric(reps)
  value <- NULL
  for (i in seq_len(reps)) {
    gc()
    times[i] <- system.time(value <- f())[["elapsed"]]
  }
  list(seconds = min(times), median_seconds = stats::median(times), value = value)
}

add_result <- function(rows, section, benchmark, method, seconds,
                       median_seconds = seconds, n_obs = NA_real_, n_unit = NA_real_,
                       n_time = NA_real_, k = NA_integer_, cutoff_km = NA_real_,
                       threads = NA_integer_, kernel = NA_character_,
                       dist_fn = NA_character_, implied_dense_weight = NA_character_,
                       implied_pairs = NA_real_, max_abs_diff = NA_real_,
                       rel_diff = NA_real_, notes = NA_character_) {
  rows[[length(rows) + 1L]] <- data.table(
    section = section, benchmark = benchmark, method = method,
    seconds = seconds, median_seconds = median_seconds, n_obs = n_obs,
    n_unit = n_unit, n_time = n_time, k = k, cutoff_km = cutoff_km,
    threads = threads, kernel = kernel, dist_fn = dist_fn,
    implied_dense_weight = implied_dense_weight, implied_pairs = implied_pairs,
    max_abs_diff = max_abs_diff, rel_diff = rel_diff, notes = notes
  )
  rows
}

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)

make_cross_section <- function(n, k, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * k), n, k,
              dimnames = list(NULL, paste0("x", seq_len(k))))
  beta <- seq(0.1, 0.5, length.out = k)
  data.table(
    y = as.numeric(X %*% beta) + stats::rnorm(n),
    X,
    lat = stats::runif(n, -55, 70),
    lon = stats::runif(n, -180, 180)
  )
}

make_panel <- function(n_unit, n_time, k, seed = 2L) {
  set.seed(seed)
  loc <- data.table(
    unit = seq_len(n_unit),
    lat = stats::runif(n_unit, 25, 49),
    lon = stats::runif(n_unit, -125, -67)
  )
  d <- data.table(
    unit = rep(seq_len(n_unit), times = n_time),
    time = rep(seq_len(n_time), each = n_unit)
  )
  d <- merge(d, loc, by = "unit", sort = FALSE)
  setorder(d, time, unit)
  X <- matrix(stats::rnorm(nrow(d) * k), nrow(d), k,
              dimnames = list(NULL, paste0("x", seq_len(k))))
  for (j in seq_len(k)) d[[paste0("x", j)]] <- X[, j]
  beta <- seq(0.1, 0.5, length.out = k)
  unit_fe <- stats::rnorm(n_unit)
  time_fe <- stats::rnorm(n_time)
  d[, y := as.numeric(X %*% beta) + unit_fe[unit] + time_fe[time] + stats::rnorm(.N)]
  d[]
}

rows <- list()
message("Large-data benchmark run")
message("  root: ", root)
message("  threads: ", paste(threads_requested, collapse = ", "))
message("  reps: ", bench_reps)
message("  run fixest on large x-section: ", run_fixest_large)

# 1. Large cross-section: fitted 1M-row cross-section, k = 10, cutoff = 100 km.
# The implied dense matrix is n x n; for n = 1e6 this is about 7.3 TiB.
{
  n <- as.integer(Sys.getenv("FASTCONLEY_LARGE_XSECTION_N", "1000000"))
  k <- 10L
  cutoff <- as.numeric(Sys.getenv("FASTCONLEY_LARGE_XSECTION_CUTOFF", "100"))
  message("1/3 fitted cross-section: n=", n, ", k=", k, ", cutoff=", cutoff)
  d <- make_cross_section(n, k, seed = 9001L)
  rhs <- paste0("x", seq_len(k), collapse = " + ")
  fit_time <- time_call(function() fixest::feols(stats::as.formula(paste("y ~", rhs)),
                                                 data = d, demeaned = TRUE), reps = 1L)
  fit <- fit_time$value
  dense_bytes <- as.numeric(n) * as.numeric(n) * 8
  baseline <- NULL
  for (th in threads_requested) {
    fc <- time_call(function() vcovSpHAC(
      fit, lat = "lat", lon = "lon", kernel = "uniform", dist_fn = "spherical",
      dist_cutoff = cutoff, pixel = 0, method = "pairwise", ncores = th, data = d
    ))
    if (is.null(baseline)) baseline <- fc$value
    rows <- add_result(rows, "large_xsection", "fitted global cross-section",
                       "fastconley", fc$seconds, fc$median_seconds,
                       n_obs = n, k = k, cutoff_km = cutoff, threads = th,
                       kernel = "uniform", dist_fn = "spherical",
                       implied_dense_weight = fmt_bytes(dense_bytes),
                       max_abs_diff = max(abs(fc$value - baseline)),
                       rel_diff = rel_diff(fc$value, baseline),
                       notes = paste0("feols fit time ", sprintf("%.2f", fit_time$seconds), "s; post-estimation VCOV only"))
  }
  if (run_fixest_large) {
    th <- max(threads_requested)
    fixest::setFixest_nthreads(th)
    fx <- time_call(function() fixest::vcov_conley(
      fit, lat = ~lat, lon = ~lon, cutoff = cutoff, distance = "spherical",
      pixel = 0, ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE), vcov_fix = FALSE
    ))
    rows <- add_result(rows, "large_xsection", "fitted global cross-section",
                       "fixest", fx$seconds, fx$median_seconds,
                       n_obs = n, k = k, cutoff_km = cutoff, threads = th,
                       kernel = "uniform", dist_fn = "spherical",
                       implied_dense_weight = fmt_bytes(dense_bytes),
                       notes = "optional large fixest comparison; post-estimation VCOV only")
  }
  rm(d, fit); gc()
}

# 2. Large panel scenario, scaled to fit comfortably on development machines by default.
{
  n_unit <- as.integer(Sys.getenv("FASTCONLEY_LARGE_PANEL_UNITS", "50000"))
  n_time <- as.integer(Sys.getenv("FASTCONLEY_LARGE_PANEL_T", "10"))
  k <- 10L
  cutoff <- as.numeric(Sys.getenv("FASTCONLEY_LARGE_PANEL_CUTOFF", "500"))
  lag_cutoff <- as.numeric(Sys.getenv("FASTCONLEY_LARGE_PANEL_LAG", "1"))
  th <- max(threads_requested)
  message("2/3 balanced panel: n_unit=", n_unit, ", T=", n_time, ", k=", k)
  d <- make_panel(n_unit, n_time, k, seed = 9002L)
  rhs <- paste0("x", seq_len(k), collapse = " + ")
  fit <- fixest::feols(stats::as.formula(paste("y ~", rhs, "| unit + time")),
                       data = d, demeaned = TRUE,
                       ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE))
  fc <- time_call(function() vcovSpHAC(
    fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
    kernel = "uniform", dist_fn = "spherical", dist_cutoff = cutoff,
    lag_cutoff = lag_cutoff, balanced_pnl = TRUE, pixel = 0,
    method = "pairwise", ncores = th, data = d
  ))
  rows <- add_result(rows, "large_panel", "balanced panel SHAC",
                     "fastconley", fc$seconds, fc$median_seconds,
                     n_obs = nrow(d), n_unit = n_unit, n_time = n_time,
                     k = k, cutoff_km = cutoff, threads = th,
                     kernel = "uniform", dist_fn = "spherical",
                     implied_dense_weight = fmt_bytes(as.numeric(n_unit) * as.numeric(n_unit) * 8),
                     notes = paste0("lag_cutoff=", lag_cutoff, "; scaled panel scenario; within-period graph size shown"))
  rm(d, fit); gc()
}

# 3. Large raster scenario: 2.25M-cell lattice, regular-raster calculation only.
{
  n_side <- as.integer(Sys.getenv("FASTCONLEY_LARGE_GRID_SIDE", "1500"))
  k <- 3L
  cutoff <- as.numeric(Sys.getenv("FASTCONLEY_LARGE_GRID_CUTOFF", "250"))
  th <- max(threads_requested)
  message("3/3 regular raster calculation: ", n_side, " x ", n_side, " cells")
  n <- n_side * n_side
  # Approx. 1.1 km latitude step, matching the the 2.25M-cell raster scale.
  step <- 0.01
  lat0 <- 35
  dlat <- step
  dlon <- step
  ring <- rep.int(seq.int(0L, n_side - 1L), times = n_side)
  col <- rep(seq.int(0L, n_side - 1L), each = n_side)
  scores <- matrix(stats::rnorm(n * k), n, k)
  time <- rep(1, n)
  # The raster scenario estimate for this scale/cutoff is ~1.8e11 implied pairs.
  implied_pairs <- 1.8e11 * (as.numeric(n) / 2250000)^2 * (cutoff / 250)^2
  for (kern in c("uniform", "bartlett")) {
    gm <- time_call(function() fastconley:::FastGridMeat(
      ring = ring, col = col, time = time, scores = scores,
      lat0 = lat0, dlat = dlat, dlon = dlon, n_ring = n_side, n_col = n_side,
      n_col_full = 0L, cutoff = cutoff, dist_fn = "spherical", kernel = kern,
      ncores = th
    ))
    rows <- add_result(rows, "large_grid", "regular raster lattice",
                       "fastconley", gm$seconds, gm$median_seconds,
                       n_obs = n, k = k, cutoff_km = cutoff, threads = th,
                       kernel = kern, dist_fn = "spherical",
                       implied_dense_weight = fmt_bytes(as.numeric(n) * as.numeric(n) * 8),
                       implied_pairs = implied_pairs,
                       notes = sprintf("%d x %d lattice; regular-raster calculation only", n_side, n_side))
  }
}

results <- rbindlist(rows, fill = TRUE)
results[, run_date := as.character(Sys.time())]
results[, fastconley_version := as.character(utils::packageVersion("fastconley"))]
results[, fixest_version := as.character(utils::packageVersion("fixest"))]
results[, lfe_version := as.character(utils::packageVersion("lfe"))]
result_file <- file.path(out_dir, "fastconley-large-data-results.csv")
utils::write.csv(results, result_file, row.names = FALSE)

session_file <- file.path(out_dir, "fastconley-large-data-session.txt")
writeLines(c(
  paste("run_date:", as.character(Sys.time())),
  paste("threads:", paste(threads_requested, collapse = ",")),
  paste("reps:", bench_reps),
  paste("fastconley:", as.character(utils::packageVersion("fastconley"))),
  paste("fixest:", as.character(utils::packageVersion("fixest"))),
  paste("lfe:", as.character(utils::packageVersion("lfe"))),
  "",
  machine_lines()
), session_file)
message("wrote ", result_file)
message("wrote ", session_file)
