test_that("balanced coordinate comparisons preserve exact and rounded equality", {
  skip_if_not_installed("fixest")
  old_rounding <- data.table::getNumericRounding()
  on.exit(data.table::setNumericRounding(old_rounding), add = TRUE)
  d <- make_balanced_panel(n_unit = 12L, n_time = 3L, k = 1L, seed = 501L)
  d[unit == 1L, `:=`(lat = 32, lon = -96)]
  d <- d[sample(.N)]
  fit <- fixest::feols(y ~ x1 | unit + time, d, demeaned = TRUE)
  run <- function(data) suppressWarnings(vcovSpHAC(
    fit, data = data, lat = "lat", lon = "lon", unit = "unit", time = "time",
    dist_cutoff = 100, balanced_pnl = TRUE, method = "pairwise", ncores = 1L))
  message <- paste0("balanced_pnl = TRUE requires time-invariant coordinates per unit, ",
                    "but some units have varying lat or lon across periods.")
  for (rounding in 0:2) {
    data.table::setNumericRounding(rounding)
    expect_error(run(d), NA)
    for (coord in c("lat", "lon")) {
      changed <- data.table::copy(d)
      idx <- which(changed$unit == 1L & changed$time == 3L)
      data.table::set(changed, i = idx, j = coord, value = changed[[coord]][idx] + 1e-13)
      expect_false(changed[[coord]][idx] == d[[coord]][idx])
      if (rounding == 0L) {
        expect_error(run(changed), message, fixed = TRUE)
      } else {
        expect_identical(data.table::uniqueN(changed[[coord]][changed$unit == 1L]), 1L)
        expect_error(run(changed), NA)
      }
      data.table::set(changed, i = idx, j = coord, value = changed[[coord]][idx] + 1)
      expect_error(run(changed), message, fixed = TRUE)
    }
  }
  data.table::setNumericRounding(0L)
  d[unit == 1L, `:=`(lat = 0, lon = 0)]
  d[unit == 1L & time == 3L, `:=`(lat = -0, lon = -0)]
  expect_error(run(d), NA)
})

test_that("matrix extraction selects coefficient order without mutating fits", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  d <- make_cross_section(n = 80L, k = 2L, seed = 502L)
  fits <- list(fixest::feols(y ~ x1 + x2, d, demeaned = TRUE),
               lfe::felm(y ~ x1 + x2, d, keepCX = TRUE))
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    field <- if (i == 1L) "X_demeaned" else "cX"
    run <- function(f) suppressWarnings(vcovSpHAC(
      f, data = d, lat = "lat", lon = "lon", dist_cutoff = 300, ncores = 1L))
    design_before <- serialize(fit[[field]], NULL)
    data_before <- serialize(d, NULL)
    ref <- run(fit)
    expect_identical(serialize(fit[[field]], NULL), design_before)
    expect_identical(serialize(d, NULL), data_before)
    fit[[field]] <- cbind(unused = seq_len(nrow(d)), fit[[field]][, 3:1])
    reordered_before <- serialize(fit[[field]], NULL)
    expect_identical(run(fit), ref)
    expect_identical(serialize(fit[[field]], NULL), reordered_before)
  }
})

test_that("serial scores retain raw rows, weighted association and stable tie order", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  # The package supports older testthat releases without binding mocks.
  skip_if_not_installed("testthat", minimum_version = "3.1.7")
  d <- make_balanced_panel(n_unit = 12L, n_time = 4L, k = 2L, seed = 503L)
  # Every unit shares coordinates: spatial aggregation must not merge the
  # units' serial histories. Include time ties and a nontrivial row order.
  d[, `:=`(lat = 35, lon = -80, w = runif(.N, 0.3, 2))]
  d <- rbind(d, d[1:4])[sample(.N)]
  d[, time := c(1, 2, 4, 7)[time]]
  fits <- list(fixest::feols(y ~ x1 + x2, d, weights = ~w, demeaned = TRUE),
               lfe::felm(y ~ x1 + x2, d, weights = d$w, keepCX = TRUE),
               fixest::fepois(I(round(abs(y) * 3)) ~ x1 + x2, d, weights = ~w))
  serial_engine <- fastconley:::FastSerialHacPanel
  seen <- NULL
  testthat::local_mocked_bindings(FastSerialHacPanel = function(...) {
    seen <<- list(...)
    serial_engine(...)
  }, .package = "fastconley")
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    scores <- if (i == 3L) {
      fit$scores * rep(1, nrow(d))
    } else {
      X <- if (i == 1L) fit$X_demeaned else fit$cX
      w <- if (i == 1L) as.numeric(fit$weights) else as.numeric(fit$weights)^2
      X * (as.numeric(fit$residuals) * w)
    }
    # This is the original two-stage metadata ordering, including tied
    # unit/time rows; verify exact scores at the actual engine boundary.
    keys <- data.table::data.table(unit = d$unit, time = d$time, row = seq_len(nrow(d)))
    data.table::setorderv(keys, "time")
    data.table::setorderv(keys, c("unit", "time"))
    expected <- scores[keys$row, , drop = FALSE]
    dimnames(expected) <- list(NULL, names(fit$coefficients))
    if (i == 2L) colnames(expected) <- rownames(fit$coefficients)
    for (pixel in c(0, 20)) {
      suppressWarnings(vcovSpHAC(
        fit, data = d, lat = "lat", lon = "lon", unit = "unit", time = "time",
        dist_cutoff = 300, lag_cutoff = 3, pixel = pixel, ncores = 1L))
      expect_identical(seen$scores, expected)
      expect_identical(seen$unit, keys$unit)
      expect_identical(seen$time, keys$time)
      expect_equal(serial_engine(unit = seen$unit, time = seen$time,
                                cutoff = 3, scores = seen$scores),
                   review_naive_serial_meat(scores, d$unit, d$time, 3),
                   tolerance = 1e-12)
    }
  }
})

test_that("lattice detection rejects latitude before evaluating longitude", {
  detect <- fastconley:::detect_lonlat_grid
  expect_null(detect(c(0, 1, 2.1), stop("longitude should remain lazy")))
  expect_null(detect(rep(1, 3), stop("longitude should remain lazy")))
  expect_null(detect(c(0, 1e-12, 1), stop("longitude should remain lazy")))
  expect_null(detect(c(0, 1, 2), c(0, 1, 2.1)))
  expect_null(detect(c(0, 1, 2), rep(1, 3)))
  # Missing cells and small coordinate noise still produce the same lattice.
  gi <- detect(c(0, 1, 3, 0), c(2, 4, 6, 6))
  expect_identical(gi, list(lat0 = 0, dlat = 1, lon0 = 2, dlon = 2,
    ring = c(0L, 1L, 3L, 0L), col = c(0L, 1L, 2L, 2L),
    n_ring = 4L, n_col = 3L, n_col_full = 180L))
  expect_equal(detect(c(0, 1, 2 + 1e-8), c(0, 1, 2))$ring, 0:2)
})

# Reference for aggregate_scores: the table aggregation (group sums by
# time, lat, lon, then sorted by the keys), on scores in the row order of dt.
reference_aggregate <- function(dt, scores) {
  tab <- data.table::data.table(time = dt$time, lat = dt$lat, lon = dt$lon,
                                scores[dt$score_row, , drop = FALSE])
  vars <- colnames(scores)
  agg <- tab[, lapply(.SD, sum), by = c("time", "lat", "lon"), .SDcols = vars]
  data.table::setorderv(agg, c("time", "lat", "lon"))
  list(lat = agg$lat, lon = agg$lon, time = agg$time,
       scores = unname(as.matrix(agg[, vars, with = FALSE])))
}
engine_scores <- function(agg) {
  unname(if (is.null(agg$rows)) agg$scores else agg$scores[agg$rows, , drop = FALSE])
}

test_that("unique locations skip aggregation with identical engine inputs", {
  set.seed(503L)
  n_unit <- 40L; n_time <- 3L
  lat <- stats::runif(n_unit, 25, 50); lon <- stats::runif(n_unit, -125, -70)
  panel <- data.table::data.table(
    unit = rep(seq_len(n_unit), times = n_time),
    time = rep(seq_len(n_time), each = n_unit),
    lat = rep(lat, times = n_time), lon = rep(lon, times = n_time))
  cross <- data.table::data.table(unit = seq_len(60L), time = 1L,
                                  lat = stats::runif(60L, 25, 50),
                                  lon = stats::runif(60L, -125, -70))
  for (dt in list(cross, panel)) {
    n <- nrow(dt)
    scores <- matrix(stats::rnorm(n * 3L), n, 3L,
                     dimnames = list(NULL, c("a", "b", "c")))
    # score_row as left by vcovSpHAC_core's sort: a permutation of the rows.
    data.table::set(dt, j = "score_row", value = sample.int(n))
    for (balanced in c(FALSE, TRUE)) {
      agg <- fastconley:::aggregate_scores(dt, scores, colnames(scores), 0,
                                           balanced, FALSE)
      ref <- reference_aggregate(dt, scores)
      expect_false(is.null(agg$rows))
      expect_identical(agg$lat, ref$lat)
      expect_identical(agg$lon, ref$lon)
      expect_identical(agg$time, ref$time)
      expect_identical(engine_scores(agg), ref$scores)
    }
  }
})

test_that("shared locations, including negative-zero keys, still aggregate", {
  set.seed(504L)
  dt <- data.table::data.table(unit = 1:8, time = 1L,
                               lat = c(10, 10, 20, 0, -0, 30, 30, 40),
                               lon = c(5, 5, 6, 7, 7, 8, 9, 10),
                               score_row = 8:1)
  scores <- matrix(stats::rnorm(16L), 8L, 2L, dimnames = list(NULL, c("a", "b")))
  agg <- fastconley:::aggregate_scores(dt, scores, colnames(scores), 0, FALSE, FALSE)
  ref <- reference_aggregate(dt, scores)
  expect_null(agg$rows)
  expect_identical(length(agg$lat), 6L)
  expect_identical(agg$lat, ref$lat)
  expect_identical(agg$lon, ref$lon)
  expect_identical(engine_scores(agg), ref$scores)
  # Only the negative zero is shared: grouping treats -0 as 0, so the table
  # path must run and merge the two rows.
  dt2 <- data.table::copy(dt)[c(4L, 5L, 3L)]
  data.table::set(dt2, j = "score_row", value = 1:3)
  agg2 <- fastconley:::aggregate_scores(dt2, scores[1:3, ], colnames(scores),
                                        0, FALSE, FALSE)
  expect_identical(length(agg2$lat), 2L)
  expect_null(agg2$rows)
})

test_that("FastSpatialMeat row maps match a gathered score matrix", {
  set.seed(505L)
  n_unit <- 150L; n_time <- 3L
  lat <- rep(stats::runif(n_unit, 25, 50), times = n_time)
  lon <- rep(stats::runif(n_unit, -125, -70), times = n_time)
  time <- rep(seq_len(n_time), each = n_unit)
  n <- length(lat)
  scores <- matrix(stats::rnorm(n * 11L), n, 11L)
  perm <- sample.int(n)
  stored <- matrix(0, n, 11L)
  stored[perm, ] <- scores
  for (neighbor in c("grid", "band")) for (balanced in c(FALSE, TRUE)) {
    for (kernel in c("bartlett", "uniform")) {
      direct <- fastconley:::FastSpatialMeat(
        lat, lon, time, scores = scores, cutoff = 400, kernel = kernel,
        balanced_pnl = balanced, neighbor = neighbor, ncores = 2L)
      mapped <- fastconley:::FastSpatialMeat(
        lat, lon, time, scores = stored, cutoff = 400, kernel = kernel,
        balanced_pnl = balanced, neighbor = neighbor, ncores = 2L, rows = perm)
      expect_identical(mapped, direct)
    }
  }
  bad <- perm; bad[1L] <- n + 1L
  expect_error(fastconley:::FastSpatialMeat(lat, lon, time, scores = stored,
                                            cutoff = 400, rows = bad),
               "rows must index rows of scores")
  bad[1L] <- NA_integer_
  expect_error(fastconley:::FastSpatialMeat(lat, lon, time, scores = stored,
                                            cutoff = 400, rows = bad),
               "rows must index rows of scores")
  expect_error(fastconley:::FastSpatialMeat(lat, lon, time, scores = stored,
                                            cutoff = 400, rows = perm[-1L]),
               "incompatible lengths")
})

test_that("the unique-location path and the table path give identical covariances", {
  skip_if_not_installed("fixest")
  old_rounding <- data.table::getNumericRounding()
  on.exit(data.table::setNumericRounding(old_rounding), add = TRUE)
  # Unit-major rows make vcovSpHAC_core's time sort permute them, so the
  # unique-location path hands the engine a composed row map. Rounded
  # numeric equality turns that path off, so the same call then runs the
  # table aggregation (random coordinates share no rounded key either).
  panel <- make_balanced_panel(n_unit = 30L, n_time = 4L, k = 2L, seed = 506L)
  panel <- panel[order(unit, time)]
  fit_panel <- fixest::feols(y ~ x1 + x2 | unit + time, panel, demeaned = TRUE)
  cross <- make_cross_section(n = 120L, k = 2L, seed = 507L)
  fit_cross <- fixest::feols(y ~ x1 + x2, cross, demeaned = TRUE)
  runs <- list(
    function() vcovSpHAC(fit_panel, data = panel, lat = "lat", lon = "lon",
                         unit = "unit", time = "time", dist_cutoff = 500,
                         lag_cutoff = 1, balanced_pnl = TRUE, ncores = 1L),
    function() vcovSpHAC(fit_panel, data = panel, lat = "lat", lon = "lon",
                         unit = "unit", time = "time", dist_cutoff = 500,
                         kernel = "uniform", lag_cutoff = 1, ncores = 2L),
    function() vcovSpHAC(fit_cross, data = cross, lat = "lat", lon = "lon",
                         dist_cutoff = 800, method = "pairwise", ncores = 2L))
  for (run in runs) {
    data.table::setNumericRounding(0L)
    fast <- suppressWarnings(run())
    data.table::setNumericRounding(2L)
    table <- suppressWarnings(run())
    expect_identical(fast, table)
  }
})
