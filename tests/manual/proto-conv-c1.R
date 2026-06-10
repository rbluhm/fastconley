# C1 prototype + error study (workstream C, notes/OPTIMIZATION_PLAN.md).
#
# The spatial meat on a lat/lon grid is sum_{cells i,j} w(d_ij) s_i' s_j with
# a radially symmetric, compactly supported kernel — a convolution. This
# prototype implements three grid-native ways to compute it and measures
# accuracy + speed against the exact pairwise engine (FastSpatialMeat on the
# cell centers, kernel = spherical great-circle, same dot-threshold accept
# logic):
#
#   ring_fft : exact for any kernel. For a fixed pair of latitude rings the
#              kernel depends only on the lon offset, so each ring pair is a
#              1D circular convolution; one inverse FFT per (target ring, k)
#              after accumulating ring-pair products in the Fourier domain.
#   boxcar   : exact, uniform kernel only. The accept set per ring pair is a
#              lon-index interval, so the inner sum is a sliding-window sum
#              over prefix sums — no FFT, no transcendentals in the loop.
#   band_fft : approximate. Latitude bands of B rings share one 2D kernel
#              evaluated at the band-center latitude (the cos(dphi) term is
#              exact; only the cos(phi1)cos(phi2) longitude scaling is
#              frozen). One padded 2D FFT convolution per (band, k). Error
#              shrinks with B; cost overhead grows as ~(1 + 2*rho/B).
#
# Run after installing the package:  Rscript tests/manual/proto-conv-c1.R
# Takes ~10-15 minutes (dominated by the pairwise reference at the dense
# scale config).

suppressMessages(library(fastconley))

ERAD <- 6371

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)

# Kernel weight from the unit-vector dot product, matching the package's
# spherical pair logic (accept iff dot >= cos(cutoff/R)).
kernel_weight <- function(dot, kernel, cutoff) {
  coscut <- cos(cutoff / ERAD)
  if (kernel == "uniform") {
    as.numeric(dot >= coscut)
  } else {
    ifelse(dot >= coscut,
           1 - ERAD * acos(pmin(1, pmax(-1, dot))) / cutoff, 0)
  }
}

# ---------------------------------------------------------------------------
# Grid container: lats[r], lons[c] in degrees (fixed steps), scores S[r, c, k]
# (dense array; unoccupied cells hold zero scores, which is exact).
# ---------------------------------------------------------------------------

make_grid <- function(lat0, lat1, lon0, lon1, step, k, occupancy = 1, seed = 1) {
  set.seed(seed)
  lats <- seq(lat0, lat1 - step, by = step) + step / 2
  lons <- seq(lon0, lon1 - step, by = step) + step / 2
  R <- length(lats); C <- length(lons)
  S <- array(rnorm(R * C * k), dim = c(R, C, k))
  occ <- if (occupancy < 1) {
    matrix(runif(R * C) < occupancy, R, C)
  } else matrix(TRUE, R, C)
  for (kk in seq_len(k)) S[, , kk][!occ] <- 0
  list(lats = lats, lons = lons, R = R, C = C, k = k, S = S, occ = occ,
       step = step)
}

# Exact pairwise reference on occupied cell centers.
meat_pairwise <- function(g, cutoff, kernel, ncores) {
  idx <- which(g$occ)
  rr <- ((idx - 1) %% g$R) + 1
  cc <- ((idx - 1) %/% g$R) + 1
  scores <- matrix(0, length(idx), g$k)
  for (kk in seq_len(g$k)) scores[, kk] <- g$S[, , kk][idx]
  fastconley:::FastSpatialMeat(g$lats[rr], g$lons[cc], rep(1, length(idx)),
                               scores = scores, cutoff = cutoff,
                               kernel = kernel, dist_fn = "spherical",
                               ncores = ncores)
}

# ---------------------------------------------------------------------------
# Variant 1: ring_fft (exact, any kernel)
# ---------------------------------------------------------------------------

meat_ring_fft <- function(g, cutoff, kernel) {
  dlat_rad <- g$step * pi / 180
  dlam_rad <- g$step * pi / 180
  ang <- cutoff / ERAD
  rho <- min(g$R, ceiling(ang / dlat_rad) + 1L)
  lat_rad <- g$lats * pi / 180
  sphi <- sin(lat_rad); cphi <- cos(lat_rad)
  dmax_g <- min(g$C, ceiling(ang / (dlam_rad * min(cphi))) + 1L)
  npad <- stats::nextn(g$C + dmax_g)

  # Rolling buffer of per-ring score FFTs (only rings within the window are
  # alive at any time).
  Shat <- vector("list", g$R)
  get_shat <- function(r) {
    if (is.null(Shat[[r]])) {
      A <- matrix(0, npad, g$k)
      A[1:g$C, ] <- g$S[r, , ]
      Shat[[r]] <<- stats::mvfft(A)
    }
    Shat[[r]]
  }

  M <- matrix(0, g$k, g$k)
  for (r1 in seq_len(g$R)) {
    lo <- max(1L, r1 - rho); hi <- min(g$R, r1 + rho)
    acc <- matrix(0 + 0i, npad, g$k)
    any_pair <- FALSE
    for (r2 in lo:hi) {
      # dot(delta) = sin1 sin2 + cos1 cos2 cos(delta * dlam)
      base <- sphi[r1] * sphi[r2]
      cc12 <- cphi[r1] * cphi[r2]
      # max lon offset that can still accept
      coscut <- cos(ang)
      rhs <- (coscut - base) / cc12
      if (rhs > 1) next                       # no offset accepts
      dmax <- if (rhs <= -1) g$C - 1L else
        min(g$C - 1L, floor(acos(rhs) / dlam_rad) + 1L)
      delta <- 0:dmax
      dots <- base + cc12 * cos(delta * dlam_rad)
      w <- kernel_weight(dots, kernel, cutoff)
      if (all(w == 0)) next
      any_pair <- TRUE
      wpad <- numeric(npad)
      wpad[1 + delta] <- w
      if (dmax >= 1) wpad[npad + 1 - (1:dmax)] <- w[-1]
      acc <- acc + stats::fft(wpad) * get_shat(r2)
    }
    if (any_pair) {
      Cf <- Re(stats::mvfft(acc, inverse = TRUE)) / npad
      M <- M + crossprod(matrix(g$S[r1, , ], g$C, g$k),
                         Cf[1:g$C, , drop = FALSE])
    }
    # Evict ring FFTs that have left the sliding window.
    if (r1 - rho >= 1) Shat[r1 - rho] <- list(NULL)
  }
  M
}

# ---------------------------------------------------------------------------
# Variant 2: boxcar (exact, uniform kernel only)
# ---------------------------------------------------------------------------

meat_boxcar <- function(g, cutoff) {
  dlat_rad <- g$step * pi / 180
  dlam_rad <- g$step * pi / 180
  ang <- cutoff / ERAD
  coscut <- cos(ang)
  rho <- min(g$R, ceiling(ang / dlat_rad) + 1L)
  lat_rad <- g$lats * pi / 180
  sphi <- sin(lat_rad); cphi <- cos(lat_rad)

  # Prefix sums per ring: P[[r]] is (C+1) x k.
  P <- vector("list", g$R)
  for (r in seq_len(g$R)) {
    P[[r]] <- rbind(0, apply(matrix(g$S[r, , ], g$C, g$k), 2, cumsum))
  }

  cidx <- seq_len(g$C)
  M <- matrix(0, g$k, g$k)
  for (r1 in seq_len(g$R)) {
    lo <- max(1L, r1 - rho); hi <- min(g$R, r1 + rho)
    Cf <- matrix(0, g$C, g$k)
    for (r2 in lo:hi) {
      base <- sphi[r1] * sphi[r2]
      cc12 <- cphi[r1] * cphi[r2]
      rhs <- (coscut - base) / cc12
      if (rhs > 1) next
      dstar <- if (rhs <= -1) g$C - 1L else floor(acos(rhs) / dlam_rad) + 1L
      # fix up boundary using the same dot expression
      while (dstar >= 0 &&
             base + cc12 * cos(dstar * dlam_rad) < coscut) dstar <- dstar - 1L
      while (dstar + 1L <= g$C - 1L &&
             base + cc12 * cos((dstar + 1L) * dlam_rad) >= coscut) dstar <- dstar + 1L
      if (dstar < 0) next
      hi_i <- pmin(g$C, cidx + dstar)
      lo_i <- pmax(1L, cidx - dstar)
      Cf <- Cf + P[[r2]][hi_i + 1L, , drop = FALSE] -
                 P[[r2]][lo_i, , drop = FALSE]
    }
    M <- M + crossprod(matrix(g$S[r1, , ], g$C, g$k), Cf)
  }
  M
}

# ---------------------------------------------------------------------------
# Variant 3: band_fft (approximate; kernel frozen at band-center latitude)
# ---------------------------------------------------------------------------

meat_band_fft <- function(g, cutoff, kernel, B) {
  dlat_rad <- g$step * pi / 180
  dlam_rad <- g$step * pi / 180
  ang <- cutoff / ERAD
  rho <- min(g$R, ceiling(ang / dlat_rad) + 1L)
  lat_rad <- g$lats * pi / 180
  dmax <- min(g$C, ceiling(ang / (dlam_rad * min(cos(lat_rad)))) + 1L)

  nr_pad <- stats::nextn(min(g$R, B) + 2L * rho)
  nc_pad <- stats::nextn(g$C + 2L * dmax)

  M <- matrix(0, g$k, g$k)
  b0 <- 1L
  while (b0 <= g$R) {
    b1 <- min(g$R, b0 + B - 1L)
    phi_c <- mean(lat_rad[b0:b1])
    t0 <- max(1L, b0 - rho); t1 <- min(g$R, b1 + rho)
    nt <- t1 - t0 + 1L

    # 2D kernel at the band-center latitude: pair (phi_c, phi_c + dr*dlat).
    dr <- -rho:rho
    dd <- 0:dmax
    dots <- outer(sin(phi_c) * sin(phi_c + dr * dlat_rad),
                  rep(1, length(dd))) +
            outer(cos(phi_c) * cos(phi_c + dr * dlat_rad),
                  cos(dd * dlam_rad))
    W2 <- matrix(kernel_weight(dots, kernel, cutoff), length(dr), length(dd))
    Kp <- matrix(0, nr_pad, nc_pad)
    ri <- ((dr %% nr_pad) + nr_pad) %% nr_pad + 1L
    Kp[ri, 1L + dd] <- W2
    if (dmax >= 1) Kp[ri, nc_pad + 1L - (1:dmax)] <- W2[, -1, drop = FALSE]
    Khat <- stats::fft(Kp)

    out_rows <- (b0 - t0 + 1L):(b1 - t0 + 1L)
    Cf <- array(0, dim = c(b1 - b0 + 1L, g$C, g$k))
    for (kk in seq_len(g$k)) {
      A <- matrix(0, nr_pad, nc_pad)
      A[1:nt, 1:g$C] <- g$S[t0:t1, , kk]
      conv <- Re(stats::fft(stats::fft(A) * Khat, inverse = TRUE)) /
        (nr_pad * nc_pad)
      Cf[, , kk] <- conv[out_rows, 1:g$C]
    }
    Sb <- matrix(g$S[b0:b1, , ], ncol = g$k)
    Cb <- matrix(Cf, ncol = g$k)
    M <- M + crossprod(Sb, Cb)
    b0 <- b1 + 1L
  }
  M
}

# ---------------------------------------------------------------------------
# Study driver
# ---------------------------------------------------------------------------

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(v = v, t = proc.time()[["elapsed"]] - t0)
}

cat("=== C1 prototype study (", format(Sys.time()), ") ===\n\n", sep = "")

## ---- 1. Validation config: accuracy of all variants ----------------------
cat("[1] Validation configs: 300 x 400 grid at 0.05 deg (~5.6 km cells), k = 3\n")
cat("    (ring/bartlett deviations ~1e-7 are acos conditioning at short\n")
cat("    distances, inherent to spherical-bartlett, not algorithm error)\n")
for (latrange in list(c(5, 20), c(45, 60))) {
  cat(sprintf("  -- latitude %d..%d --\n", latrange[1], latrange[2]))
  g <- make_grid(latrange[1], latrange[2], 10, 30, 0.05, k = 3, seed = 11)
  for (cutoff in c(100, 300)) {
    ref_u <- meat_pairwise(g, cutoff, "uniform", ncores = 16)
    ref_b <- meat_pairwise(g, cutoff, "bartlett", ncores = 16)
    r_u <- meat_ring_fft(g, cutoff, "uniform")
    r_b <- meat_ring_fft(g, cutoff, "bartlett")
    bx <- meat_boxcar(g, cutoff)
    cat(sprintf("  cutoff %4d km: ring/unif %.2e  ring/bart %.2e  boxcar %.2e\n",
                cutoff, rel_diff(r_u, ref_u), rel_diff(r_b, ref_b),
                rel_diff(bx, ref_u)))
    for (B in c(20, 50, 100, 200)) {
      bf_u <- meat_band_fft(g, cutoff, "uniform", B)
      bf_b <- meat_band_fft(g, cutoff, "bartlett", B)
      cat(sprintf("    band B=%3d (%5.2f deg): unif %.3e  bart %.3e\n",
                  B, B * g$step, rel_diff(bf_u, ref_u), rel_diff(bf_b, ref_b)))
    }
  }
}

## ---- 2. Dense scale config: speed crossover ------------------------------
cat("\n[2] Dense scale config: 1500 x 1500 grid at 0.01 deg (~1.1 km cells),\n")
cat("    k = 4, cutoff 250 km, 100% occupancy (raster regime)\n")
gs <- make_grid(0, 15, 0, 15, 0.01, k = 4, seed = 22)
n_cells <- gs$R * gs$C
rho_c <- 250 / (111.2 * gs$step)
nnz_est <- n_cells * pi * rho_c^2 / 2
cat(sprintf("    cells = %.2fM, est. pairs = %.2e\n", n_cells / 1e6, nnz_est))

pw16 <- timed(meat_pairwise(gs, 250, "uniform", ncores = 16))
cat(sprintf("    pairwise 16c (uniform):       %7.1f s\n", pw16$t))
# 1-core projection from a calibration run on a subgrid
gcal <- make_grid(0, 4, 0, 4, 0.01, k = 4, seed = 22)
cal16 <- timed(meat_pairwise(gcal, 250, "uniform", ncores = 16))
cal1 <- timed(meat_pairwise(gcal, 250, "uniform", ncores = 1))
proj1 <- pw16$t * cal1$t / cal16$t
cat(sprintf("    pairwise 1c (projected):      %7.1f s  (calib ratio %.1fx)\n",
            proj1, cal1$t / cal16$t))

bx_s <- timed(meat_boxcar(gs, 250))
cat(sprintf("    boxcar exact (uniform, R):    %7.1f s  err %.2e\n",
            bx_s$t, rel_diff(bx_s$v, pw16$v)))
bf_s <- timed(meat_band_fft(gs, 250, "uniform", 150))
cat(sprintf("    band_fft B=150 (uniform, R):  %7.1f s  err %.2e\n",
            bf_s$t, rel_diff(bf_s$v, pw16$v)))
rf_s <- timed(meat_ring_fft(gs, 250, "uniform"))
cat(sprintf("    ring_fft exact (uniform, R):  %7.1f s  err %.2e\n",
            rf_s$t, rel_diff(rf_s$v, pw16$v)))
bfb_s <- timed(meat_band_fft(gs, 250, "bartlett", 150))
rfb_s <- timed(meat_ring_fft(gs, 250, "bartlett"))
cat(sprintf("    bartlett conv (no pairwise ref at this scale):\n"))
cat(sprintf("      band_fft B=150: %7.1f s   ring_fft: %7.1f s   band-vs-ring %.3e\n",
            bfb_s$t, rfb_s$t, rel_diff(bfb_s$v, rfb_s$v)))

## ---- 3. Sparse config: where pairwise keeps winning ----------------------
cat("\n[3] Sparse config: same grid, 5% occupancy (scattered-points regime)\n")
gsp <- make_grid(0, 15, 0, 15, 0.01, k = 4, occupancy = 0.05, seed = 33)
pw_sp <- timed(meat_pairwise(gsp, 250, "uniform", ncores = 16))
bf_sp <- timed(meat_band_fft(gsp, 250, "uniform", 150))
cat(sprintf("    pairwise 16c: %6.1f s   band_fft (R): %6.1f s   err %.2e\n",
            pw_sp$t, bf_sp$t, rel_diff(bf_sp$v, pw_sp$v)))

cat("\nDone.\n")
