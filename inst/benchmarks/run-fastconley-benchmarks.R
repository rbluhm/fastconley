#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fastconley)
})

for (pkg in c("fixest", "lfe")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for the benchmark script.")
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

bench_threads <- as.integer(Sys.getenv(
  "FASTCONLEY_BENCH_THREADS",
  as.character(min(8L, parallel::detectCores(logical = TRUE)))
))
bench_reps <- as.integer(Sys.getenv("FASTCONLEY_BENCH_REPS", "2"))
bench_reps <- max(1L, bench_reps)

set.seed(1001)
RcppParallel::setThreadOptions(numThreads = bench_threads)

message("fastconley benchmark run")
message("  root: ", root)
message("  threads: ", bench_threads)
message("  reps per cell: ", bench_reps)

time_call <- function(f, reps = bench_reps) {
  times <- numeric(reps)
  value <- NULL
  for (i in seq_len(reps)) {
    gc()
    elapsed <- system.time(value <- f())[["elapsed"]]
    times[i] <- elapsed
  }
  list(seconds = min(times), median_seconds = stats::median(times), value = value)
}

rel_diff <- function(a, b) {
  max(abs(a - b)) / max(abs(b), 1e-300)
}

add_result <- function(rows, section, benchmark, method, seconds,
                       median_seconds = seconds, n_obs = NA_integer_,
                       n_unit = NA_integer_, n_time = NA_integer_,
                       k = NA_integer_, cutoff_km = NA_real_,
                       threads = NA_integer_, kernel = NA_character_,
                       dist_fn = NA_character_, raw_rows = NA_integer_,
                       cases = NA_integer_, max_abs_diff = NA_real_,
                       rel_diff_value = NA_real_, notes = NA_character_) {
  rows[[length(rows) + 1L]] <- data.table(
    section = section,
    benchmark = benchmark,
    method = method,
    seconds = seconds,
    median_seconds = median_seconds,
    n_obs = n_obs,
    n_unit = n_unit,
    n_time = n_time,
    k = k,
    cutoff_km = cutoff_km,
    threads = threads,
    kernel = kernel,
    dist_fn = dist_fn,
    raw_rows = raw_rows,
    cases = cases,
    max_abs_diff = max_abs_diff,
    rel_diff = rel_diff_value,
    notes = notes
  )
  rows
}

make_cross_section <- function(n, k, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * k), n, k,
              dimnames = list(NULL, paste0("x", seq_len(k))))
  beta <- seq(0.1, 0.5, length.out = k)
  data.table(
    y = as.numeric(X %*% beta) + stats::rnorm(n),
    X,
    lat = stats::runif(n, 25, 49),
    lon = stats::runif(n, -125, -67)
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
  data.table::setorder(d, time, unit)
  X <- matrix(stats::rnorm(nrow(d) * k), nrow(d), k,
              dimnames = list(NULL, paste0("x", seq_len(k))))
  for (j in seq_len(k)) d[[paste0("x", j)]] <- X[, j]
  beta <- seq(0.1, 0.5, length.out = k)
  unit_fe <- stats::rnorm(n_unit)
  time_fe <- stats::rnorm(n_time)
  d[, y := as.numeric(X %*% beta) + unit_fe[unit] + time_fe[time] + stats::rnorm(.N)]
  d[]
}

dense_uniform_spherical <- function(lat, lon, scores, cutoff) {
  rad <- pi / 180
  latr <- lat * rad
  lonr <- lon * rad
  xyz <- cbind(cos(latr) * cos(lonr), cos(latr) * sin(lonr), sin(latr))
  W <- tcrossprod(xyz) >= cos(cutoff / 6371)
  storage.mode(W) <- "double"
  crossprod(scores, W %*% scores)
}

rows <- list()

message("1/5 dense reference scaling")
for (n in c(1000L, 2000L, 4000L)) {
  k <- 5L
  cutoff <- 500
  d <- make_cross_section(n, k, seed = 10L + n)
  X <- as.matrix(d[, paste0("x", seq_len(k)), with = FALSE])
  e <- stats::rnorm(n)
  S <- X * e
  dense <- time_call(function() dense_uniform_spherical(d$lat, d$lon, S, cutoff))
  fast <- time_call(function() fastconley:::FastSpatialMeat(
    lat = d$lat, lon = d$lon, time = rep(1, n), scores = S,
    cutoff = cutoff, kernel = "uniform", dist_fn = "spherical",
    balanced_pnl = FALSE, ncores = bench_threads
  ))
  err <- max(abs(dense$value - fast$value))
  rows <- add_result(rows, "dense", "dense uniform spherical", "dense",
                     dense$seconds, dense$median_seconds, n_obs = n, k = k,
                     cutoff_km = cutoff, threads = 1L, kernel = "uniform",
                     dist_fn = "spherical", max_abs_diff = err,
                     notes = sprintf("dense dot matrix is about %.1f MB", n * n * 8 / 1024^2))
  rows <- add_result(rows, "dense", "dense uniform spherical", "fastconley",
                     fast$seconds, fast$median_seconds, n_obs = n, k = k,
                     cutoff_km = cutoff, threads = bench_threads, kernel = "uniform",
                     dist_fn = "spherical", max_abs_diff = err)
}

message("2/5 scattered cross-section against fixest")
xs_cases <- data.table(
  n = c(50000L, 50000L, 100000L, 100000L),
  cutoff = c(100, 500, 100, 500),
  threads = bench_threads
)
if (bench_threads > 1L) {
  xs_cases <- rbind(
    data.table(n = 50000L, cutoff = 500, threads = 1L),
    xs_cases
  )
}
for (ii in seq_len(nrow(xs_cases))) {
  cfg <- xs_cases[ii]
  k <- 10L
  d <- make_cross_section(cfg$n, k, seed = 200L + ii)
  rhs <- paste0("x", seq_len(k), collapse = " + ")
  fit <- fixest::feols(stats::as.formula(paste("y ~", rhs)),
                       data = d, demeaned = TRUE)
  fixest::setFixest_nthreads(cfg$threads)
  fc <- time_call(function() vcovSpHAC(
    fit, lat = "lat", lon = "lon", kernel = "uniform",
    dist_fn = "spherical", dist_cutoff = cfg$cutoff,
    pixel = 0, method = "pairwise", ncores = cfg$threads, data = d
  ))
  fx <- time_call(function() fixest::vcov_conley(
    fit, lat = ~lat, lon = ~lon, cutoff = cfg$cutoff,
    distance = "spherical", pixel = 0,
    ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE),
    vcov_fix = FALSE
  ))
  err <- max(abs(fc$value - fx$value))
  rd <- rel_diff(fc$value, fx$value)
  label <- sprintf("n=%s cutoff=%s", format(cfg$n, scientific = FALSE), cfg$cutoff)
  rows <- add_result(rows, "cross_section", label, "fastconley",
                     fc$seconds, fc$median_seconds, n_obs = cfg$n, k = k,
                     cutoff_km = cfg$cutoff, threads = cfg$threads,
                     kernel = "uniform", dist_fn = "spherical",
                     max_abs_diff = err, rel_diff_value = rd,
                     notes = "post-estimation only")
  rows <- add_result(rows, "cross_section", label, "fixest",
                     fx$seconds, fx$median_seconds, n_obs = cfg$n, k = k,
                     cutoff_km = cfg$cutoff, threads = cfg$threads,
                     kernel = "uniform", dist_fn = "spherical",
                     max_abs_diff = err, rel_diff_value = rd,
                     notes = "post-estimation only")
}

message("3/5 balanced panel SHAC against fixest composition")
{
  n_unit <- 10000L
  n_time <- 4L
  k <- 5L
  cutoff <- 500
  lag_cutoff <- 1
  d <- make_panel(n_unit, n_time, k, seed = 300L)
  rhs <- paste0("x", seq_len(k), collapse = " + ")
  fit <- fixest::feols(stats::as.formula(paste("y ~", rhs, "| unit + time")),
                       data = d, demeaned = TRUE,
                       ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE))
  fixest::setFixest_nthreads(bench_threads)
  fc <- time_call(function() vcovSpHAC(
    fit, unit = "unit", time = "time", lat = "lat", lon = "lon",
    kernel = "uniform", dist_fn = "spherical", dist_cutoff = cutoff,
    lag_cutoff = lag_cutoff, balanced_pnl = TRUE, pixel = 0,
    ncores = bench_threads, data = d
  ))
  fx <- time_call(function() {
    v_spatial <- fixest::vcov_conley(
      fit, lat = ~lat, lon = ~lon, cutoff = cutoff,
      distance = "spherical", pixel = 0,
      ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE),
      vcov_fix = FALSE
    )
    v_serial <- fixest::vcov_NW(
      fit, unit = ~unit, time = ~time, lag = lag_cutoff,
      ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE),
      vcov_fix = FALSE
    )
    v_robust <- fixest::vcov_hetero(
      fit, ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE),
      vcov_fix = FALSE
    )
    v_spatial + v_serial - v_robust
  })
  rd <- rel_diff(fc$value, fx$value)
  err <- max(abs(fc$value - fx$value))
  rows <- add_result(rows, "panel", "balanced panel SHAC", "fastconley",
                     fc$seconds, fc$median_seconds, n_obs = nrow(d),
                     n_unit = n_unit, n_time = n_time, k = k,
                     cutoff_km = cutoff, threads = bench_threads,
                     kernel = "uniform", dist_fn = "spherical",
                     max_abs_diff = err, rel_diff_value = rd,
                     notes = paste0("lag_cutoff=", lag_cutoff))
  rows <- add_result(rows, "panel", "balanced panel SHAC", "fixest composition",
                     fx$seconds, fx$median_seconds, n_obs = nrow(d),
                     n_unit = n_unit, n_time = n_time, k = k,
                     cutoff_km = cutoff, threads = bench_threads,
                     kernel = "uniform", dist_fn = "spherical",
                     max_abs_diff = err, rel_diff_value = rd,
                     notes = "vcov_conley + vcov_NW - vcov_hetero")
}

message("4/5 raster grid engine")
{
  n_side <- 180L
  step <- 0.05
  lats <- seq(35, by = step, length.out = n_side)
  lons <- seq(-10, by = step, length.out = n_side)
  d <- as.data.table(expand.grid(lat = lats, lon = lons))
  n <- nrow(d)
  k <- 3L
  for (j in seq_len(k)) d[[paste0("x", j)]] <- stats::rnorm(n)
  d[, y := 0.3 * x1 - 0.2 * x2 + 0.1 * x3 + stats::rnorm(.N)]
  fit <- lfe::felm(y ~ x1 + x2 + x3, data = d, keepCX = TRUE)
  cutoff <- 250
  for (kern in c("uniform", "bartlett")) {
    grid <- time_call(function() vcovSpHAC(
      fit, lat = "lat", lon = "lon", kernel = kern, dist_fn = "spherical",
      dist_cutoff = cutoff, method = "grid", ncores = bench_threads, data = d
    ))
    pair <- time_call(function() vcovSpHAC(
      fit, lat = "lat", lon = "lon", kernel = kern, dist_fn = "spherical",
      dist_cutoff = cutoff, method = "pairwise", ncores = bench_threads, data = d
    ))
    err <- max(abs(grid$value - pair$value))
    rd <- rel_diff(grid$value, pair$value)
    rows <- add_result(rows, "grid", "regular lat/lon lattice", "grid",
                       grid$seconds, grid$median_seconds, n_obs = n, k = k,
                       cutoff_km = cutoff, threads = bench_threads,
                       kernel = kern, dist_fn = "spherical",
                       max_abs_diff = err, rel_diff_value = rd,
                       notes = sprintf("%d x %d lattice", n_side, n_side))
    rows <- add_result(rows, "grid", "regular lat/lon lattice", "pairwise",
                       pair$seconds, pair$median_seconds, n_obs = n, k = k,
                       cutoff_km = cutoff, threads = bench_threads,
                       kernel = kern, dist_fn = "spherical",
                       max_abs_diff = err, rel_diff_value = rd,
                       notes = "3D cell-grid pair enumeration")
  }
}

message("5/5 coordinate aggregation and pixel")
{
  n_unique <- 20000L
  repeats <- 5L
  k <- 5L
  cutoff <- 250
  set.seed(500L)
  base <- data.table(
    lat = stats::runif(n_unique, 25, 49),
    lon = stats::runif(n_unique, -125, -67)
  )
  d <- base[rep(seq_len(.N), each = repeats)]
  for (j in seq_len(k)) d[[paste0("x", j)]] <- stats::rnorm(nrow(d))
  beta <- seq(0.1, 0.5, length.out = k)
  X <- as.matrix(d[, paste0("x", seq_len(k)), with = FALSE])
  d[, y := as.numeric(X %*% beta) + stats::rnorm(.N)]
  rhs <- paste0("x", seq_len(k), collapse = " + ")
  fit <- fixest::feols(stats::as.formula(paste("y ~", rhs)),
                       data = d, demeaned = TRUE)
  count_cases <- function(pixel) {
    dt <- data.table(fit$X_demeaned, unit = seq_len(nrow(d)), time = 1,
                     lat = d$lat, lon = d$lon, e = as.numeric(fit$residuals))
    xvars <- names(fit$coefficients)
    length(fastconley:::aggregate_scores(dt, xvars, pixel,
                                          balanced_pnl = FALSE,
                                          verbose = FALSE)$lat)
  }
  baseline <- NULL
  for (px in c(0, 10, 25)) {
    tm <- time_call(function() vcovSpHAC(
      fit, lat = "lat", lon = "lon", kernel = "bartlett",
      dist_fn = "spherical", dist_cutoff = cutoff, pixel = px,
      method = "pairwise", ncores = bench_threads, data = d
    ))
    if (is.null(baseline)) baseline <- tm$value
    rows <- add_result(rows, "pixel", "repeated-coordinate cross-section",
                       paste0("pixel=", px), tm$seconds, tm$median_seconds,
                       n_obs = nrow(d), k = k, cutoff_km = cutoff,
                       threads = bench_threads, kernel = "bartlett",
                       dist_fn = "spherical", raw_rows = nrow(d),
                       cases = count_cases(px),
                       max_abs_diff = max(abs(tm$value - baseline)),
                       rel_diff_value = rel_diff(tm$value, baseline),
                       notes = if (px == 0) "exact coordinate dedupe" else "approximate snapping")
  }
}

results <- rbindlist(rows, fill = TRUE)
results[, run_date := as.character(Sys.time())]
results[, fastconley_version := as.character(utils::packageVersion("fastconley"))]
results[, fixest_version := as.character(utils::packageVersion("fixest"))]
results[, lfe_version := as.character(utils::packageVersion("lfe"))]
result_file <- file.path(out_dir, "fastconley-benchmark-results.csv")
utils::write.csv(results, result_file, row.names = FALSE)

session_file <- file.path(out_dir, "fastconley-benchmark-session.txt")
writeLines(c(
  paste("run_date:", as.character(Sys.time())),
  paste("threads:", bench_threads),
  paste("reps:", bench_reps),
  paste("fastconley:", as.character(utils::packageVersion("fastconley"))),
  paste("fixest:", as.character(utils::packageVersion("fixest"))),
  paste("lfe:", as.character(utils::packageVersion("lfe"))),
  "",
  capture.output(sessionInfo())
), session_file)

message("wrote ", result_file)
message("wrote ", session_file)
