# fastconley 0.7.0

Workstream C2: exact grid-native meat for raster data.

## New: `method = c("auto", "pairwise", "grid")`

For the uniform kernel on a regular lat/lon lattice (raster cell centers,
gridded covariates), the within-cutoff accept set between two latitude
rings is a longitude-index interval, so the spatial meat reduces to
sliding-window sums over per-ring prefix sums — `FastGridMeat`. Cost is
O(n_ring * window * n_col * k), **independent of the pair count**, and the
accept threshold is the same dot-product constant the pairwise engine
uses, so the result is exact (agrees to FP summation order; no
approximation anywhere for natively gridded data).

- C1-study config (2.25M cells at ~1.1 km, 250 km cutoff, ~1.8e11 pairs):
  0.81 s vs 244 s for the 16-core pairwise engine — **302x**, error 5e-15.
- 9M-cell continental raster (3000 x 3000, 250 km, ~7e11 pairs): 3.5 s.
- Supports all three `dist_fn`s (monotone in the chord), multiple time
  blocks (panels apply it per period), sparse occupancy, duplicate cells,
  and is deterministic across `ncores`.

`method = "auto"` (the new default) switches to the grid engine only when
a lattice is detected, the kernel is uniform, and a flop-balance estimate
says it wins; otherwise the pairwise engine runs as before. Scattered
(non-lattice) data is unaffected. `method = "grid"` errors informatively
when its requirements are not met.

The bartlett ring-FFT variant (exact per the C1 study) is planned as a
follow-up; bartlett rasters currently stay on the pairwise engine.

# fastconley 0.6.1

Minor-backlog items M1-M4 from `notes/OPTIMIZATION_PLAN.md`.

- **Serial HAC is now O(T * k) per unit** instead of O(T^2 * k): the
  Bartlett-in-time weight decomposes over a sliding two-pointer window
  (with a per-block time shift for conditioning). A 100k-row panel with
  T = 2000 and lag = 50 computes in ~6 ms. Rows must be time-sorted within
  unit blocks (the R layer always sorts; direct callers get a clear error).
- **`vcovSpHAC.felm` gains `data =`**: when the passed frame has the same
  row count as the fit (no NAs dropped, no subset), coordinates are taken
  by direct column access with no model-frame re-evaluation; otherwise the
  frame overrides the call-recovered data in the aligned model-frame path.
- **Prep path parallelized** (grid sort via TBB parallel_sort, coordinate
  cache trig, permutation gathers) — all order-preserving, so results stay
  bit-identical across `ncores`. 1M-row global cross-section at 16 threads:
  0.46 -> 0.21 s (scaling 2.2x -> 4.8x).
- Spatial results are bitwise identical to v0.6.0; the serial meat differs
  by ~1e-15 relative from the new summation algebra.
- README benchmarks refreshed against fixest 0.14.1 at scale: 10.7-45.8x on
  cross-sections up to n = 1M, 1.4-2.9x on panel SHAC up to 1.5M obs (where
  fastconley computes ~T x more spatial work by construction). A brute-force
  arbitration shows fastconley matches the exact great-circle uniform meat
  to ~1e-15 while fixest 0.14.1's "spherical" distance deviates ~2e-2.

# fastconley 0.6.0

Phase 2 of `notes/OPTIMIZATION_PLAN.md`: memory diet, deterministic
reduction, and screen work. Results remain exact (same pairs, same
weights); summation order changed, so values differ from v0.5.0 by
<= ~5e-15 relative.

## Results are now invariant to `ncores`

All meat accumulations use a deterministic chunked reduction (fixed-size
row/block chunks, partials summed in chunk order). `ncores = 1` and
`ncores = 16` produce bit-identical matrices — the old multicore
tolerance caveat is gone.

## Memory

- The C++ entry points take a pre-computed score matrix and alias R
  memory directly (no Rcpp input copies, no intermediate unsorted score
  buffer; one gather pass builds the sorted layout). `vcovSpHAC` passes
  the aggregated scores straight through; `X`/`e` arguments remain on the
  internal wrappers for convenience. Measured peak RSS: 4M-row
  cross-section 2321 -> 1558 MB (-33%); 40k x 40 balanced panel
  1028 -> 748 MB (-27%). Same speed, identical checksums.
- New `csr_weight = c("double", "float")`: opt-in float storage for the
  balanced-path bartlett weights (8 -> 4 bytes per pair, <= ~6e-8
  relative error per weight). Default stays exact double.

## Speed

- Haversine gains a conservative 3D dot-product pre-screen (exact
  identity `a = (1 - dot)/2`); the exact a-test remains the arbiter, so
  results are bit-identical. CONUS 100k / 500 km bartlett/haversine:
  15.1 -> 10.5 s single-threaded (1.4x).
- Balanced 16-thread uniform meat: 2.30 -> 1.56 s on a 40k x 40 panel
  (better task granularity from the chunked reduction).

## Tried and reverted (documented for the record)

A unit-major T*k stacked score layout that streams the balanced CSR once
instead of once per period was implemented and benchmarked. It lost to
the period-major layout: spatial sorting makes neighbor gathers a
sliding window of k-wide rows that stays L2-resident per period; any
wider stacking pushes the window past L2 and thrashes the shared L3 at
high thread counts (2.6 s vs 2.0 s at 16 cores). The period-major
traversal stays, now deterministic and float-capable.

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
