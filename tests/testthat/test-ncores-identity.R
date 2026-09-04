test_that("pairwise grid and band streaming are bit-identical across ncores", {
  set.seed(111)
  n <- 2200L
  lat <- runif(n, -55, 55)
  lon <- runif(n, -170, 170)
  time <- rep(1:2, each = n / 2L)
  scores <- matrix(rnorm(n * 2L), ncol = 2L)
  for (neighbor in c("grid", "band")) {
    one <- fastconley:::FastSpatialMeat(
      lat, lon, time, scores = scores, cutoff = 150,
      kernel = "bartlett", dist_fn = "haversine", balanced_pnl = FALSE,
      ncores = 1, neighbor = neighbor)
    two <- fastconley:::FastSpatialMeat(
      lat, lon, time, scores = scores, cutoff = 150,
      kernel = "bartlett", dist_fn = "haversine", balanced_pnl = FALSE,
      ncores = 2, neighbor = neighbor)
    expect_identical(two, one, label = neighbor)
  }
})

test_that("balanced CSR is bit-identical across ncores", {
  set.seed(112)
  n_unit <- 1200L
  n_time <- 2L
  lat_unit <- runif(n_unit, 25, 50)
  lon_unit <- runif(n_unit, -125, -70)
  lat <- rep(lat_unit, times = n_time)
  lon <- rep(lon_unit, times = n_time)
  time <- rep(seq_len(n_time), each = n_unit)
  scores <- matrix(rnorm(length(time) * 2L), ncol = 2L)
  run <- function(ncores) fastconley:::FastSpatialMeat(
    lat, lon, time, scores = scores, cutoff = 250,
    kernel = "uniform", dist_fn = "spherical", balanced_pnl = TRUE,
    ncores = ncores, neighbor = "grid")
  expect_identical(run(2), run(1))
})

test_that("serial HAC is bit-identical across ncores", {
  set.seed(113)
  unit <- rep(seq_len(300L), each = 4L)
  time <- rep(seq_len(4L), times = 300L)
  scores <- matrix(rnorm(length(unit) * 3L), ncol = 3L)
  run <- function(ncores) fastconley:::FastSerialHacPanel(
    unit, time, cutoff = 2, scores = scores, ncores = ncores)
  expect_identical(run(2), run(1))
})

test_that("native grid uniform and Bartlett are bit-identical across ncores", {
  set.seed(114)
  n_side <- 28L
  gp <- expand.grid(ring = 0:(n_side - 1L), col = 0:(n_side - 1L))
  scores <- matrix(rnorm(nrow(gp) * 2L), ncol = 2L)
  time <- rep(1, nrow(gp))
  run <- function(kernel, ncores) fastconley:::FastGridMeat(
    ring = gp$ring, col = gp$col, time = time, scores = scores,
    lat0 = 30, dlat = 0.25, dlon = 0.25,
    n_ring = n_side, n_col = n_side, n_col_full = 0L,
    cutoff = 300, dist_fn = "haversine", kernel = kernel,
    ncores = ncores)
  for (kernel in c("uniform", "bartlett")) {
    expect_identical(run(kernel, 2), run(kernel, 1), label = kernel)
  }
})
