# fastconley 0.5.0

Phase 1 of `notes/OPTIMIZATION_PLAN.md`: 3D cell-grid neighbor search.

## New: `neighbor = c("grid", "band")`

The spatial meat's candidate enumeration now defaults to a 3D cell grid.
Points are bucketed by their unit vectors into a cubic grid whose edge is
the unit-sphere chord equivalent of the cutoff; every supported distance is
monotone in the chord, so accepted pairs are never more than one cell apart
per axis — no pole or dateline special cases. Each row scans its own cell
plus five contiguous row ranges covering the 13 forward neighbor cells:
~3–4 candidates per accepted pair, independent of geographic extent, versus
the latitude band scan's `2*L_lon/(pi*r)`.

Both strategies call the identical per-pair accept test, so pair sets and
weights are exactly the same; results differ only by floating-point
summation order (observed <= 6e-15 relative across the validation matrix;
`neighbor = "band"` remains bitwise identical to v0.4.1). The band path is
kept for one release and will be removed in v0.6.0.

Measured single-threaded speedups (FastSpatialMeat, k = 10):

- global extent, n = 1M, 100 km: 8.6x (uniform/spherical) to 25.3x
  (bartlett/haversine)
- CONUS, n = 100k, 100 km: 4.5x; 500 km: 2.1-2.3x (that config is
  dominated by true-pair accumulation, not candidate search)
- balanced CSR build, 200k units x 4 periods, 100 km: 3.6x
- 16-thread scaling on CONUS 100k/500 km improves from ~5.3x to ~7.9x
  (less memory-bandwidth pressure from candidate streaming)

# fastconley 0.4.1

Phase 0 quick wins from `notes/OPTIMIZATION_PLAN.md` (Q1–Q6). Numerical
results are unchanged (verified bitwise against v0.4.0 on the balanced /
general / unbalanced × kernel × distance matrix at `ncores = 1`).

## Memory

- The balanced-path CSR no longer allocates its weight array for
  `kernel = "uniform"` (it was filled but never read), and column indices
  are now 32-bit. Per stored pair: 16 bytes -> 4 (uniform) / 12 (bartlett).
  A guard errors if a single period exceeds 2^32 - 1 units.

## Speed

- The haversine longitude screen uses a multiply instead of a per-pair
  divide (same accept set; the exact `a`-test is unchanged).
- The `balanced_pnl = TRUE` unit-set validation is now O(n) vectorized
  instead of two grouping passes plus per-period string materialization.

## Cleanup

- Removed the legacy `TimeDist` internal export (unused since the serial
  HAC moved to `FastSerialHacPanel`).
- README gains an optional `~/.R/Makevars` performance-build recipe.

# fastconley 0.3.0.9000

Development build for testing a faster spatial HAC path.

## Spatial meat changes

- Adds a `FastSpatialMeat()` internal Rcpp routine.
- Replaces dense `n x n` spatial distance matrices in `vcovSpHAC()` with a screened spatial edge list.
- Computes the meat via cumulative scores: `S = e * X`, `C_i = 0.5 S_i + sum_j w_ij S_j`, and `S'C + C'S`.
- Supports both `kernel = "bartlett"` and `kernel = "uniform"`.
- Supports `dist_fn` choices: haversine, spherical, chord. The upstream `flatearth` option was dropped — its formula was not symmetric in `(i, j)` and the equirectangular approximation is no longer worth a separate code path now that spherical+uniform runs as a 3D dot-product threshold with no trig.
- In balanced panels, builds one spatial edge list from the first period and reuses it for all periods.

## Performance work (most recent)

- **Templated `(dist_fn, kernel)` dispatch.** The screen + distance + kernel-weight stack is specialised at compile time for each of the six `(distance, kernel)` combinations, removing per-pair runtime branches.
- **3D unit-vector threshold for spherical and chord.** `(spherical, uniform)` pair inclusion now reduces to a 3D dot product compared against `cos(cutoff / R)` — no trig in the inner loop. `(chord, uniform)` uses squared-Euclidean against `(cutoff / R)²`. The `bartlett` variants only pay `acos`/`sqrt` on accepted pairs.
- **Row-major score buffer.** Scores `s_i = e_i · X_i` are stored row-major in a single flat `std::vector<double>` instead of being recomputed against column-major Armadillo memory on every pair.
- **Sorted-order layout.** Both the coordinate cache and the scores buffer are permuted into latitude-sorted order before the meat workers run, so every read in the hot loop is sequential — independent of how the caller laid out the input.
- **`pixel` argument** for score pre-aggregation (`R/vcovSpHAC.R`). At `pixel = 0` rows that share `(lat, lon)` within a time period are collapsed exactly; at `pixel > 0` rows are first snapped to a uniform `pixel`-km grid (speed/accuracy trade-off).
- Cross-section vs `fixest::vcov_conley` (`kernel = "uniform"`, spherical, 500 km cutoff): `fastconley` is now 3-4× faster at `pixel = 0` and 27-62× faster at `pixel = 25` across `n ∈ {5k … 100k}`, and scales better than fixest past 4 threads.

## Notes

- This is meant as a test branch/source drop, not a CRAN-ready release.
- The implementation uses RcppParallel for CSR construction and meat accumulation.
- Rcpp exports were updated manually. Running `Rcpp::compileAttributes()` is recommended after further edits.
