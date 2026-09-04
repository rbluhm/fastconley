# fastconley 0.11.0

Engine extraction, a Stata port sharing the same C++ engine, and a review
round (six independent code reviews of the engine, the R layer, the Stata
command, the plugin build, the reghdfe proposal, and the validation
coverage) whose findings are fixed below.

## Corrections that change numbers

- **Clustered `felm` fits got the wrong default `ssc`.** `lfe` stores the
  cluster-adjusted degrees of freedom (number of clusters minus one) in
  `df.residual` when a fit carries a cluster specification, and the
  small-sample scale `n / df.residual` inherited that value: on a 20-row
  fixture the covariance was five times too large. The scale now uses the
  regression's own residual degrees of freedom (`N - p`), so clustered and
  unclustered fits of the same model give identical Conley covariances.
- **`vcovSpHAC.felm` ignored the `unit` and `time` arguments** and always
  used the first two absorbed fixed effects. The names are now matched
  against `names(reg$fe)` (with the model data as a fallback), a missing
  name is an error, and the first-two-FE default applies only when nothing
  was passed. A model with `| region + unit + time` used the wrong keys
  before.
- **Non-numeric time collapsed gaps in the serial HAC.** Character or factor
  time values were recoded to consecutive integers, so "2000", "2002" became
  lag 1. Values that parse as numbers now keep their numeric spacing (both
  methods, including the absorbed-FE factor in `felm` fits); a time key that
  does not parse is accepted for spatial blocking only and is an error when
  `lag_cutoff > 0`.
- **`time` without `unit` is honoured.** The fixest method silently dropped
  the time blocking when no unit was given; it now blocks by time with each
  row as its own unit, matching the Stata command. `lag_cutoff > 0` without
  a unit is an error.
- **fixest observation selection uses `fixest::obs()`**, so models extracted
  from `split =` estimations and fits with several stacked selections align
  their coordinates correctly. `fixest_multi` objects, `lean = TRUE` `feols`
  fits, and fixed-effects-only models now fail with clear messages instead
  of misleading ones.
- **Grid engine Bartlett weights use the same chord form as the pairwise
  engine** (haversine and spherical). Four raster/grid Bartlett battery
  configurations move by at most 3.3e-14 relative; everything else is
  bit-identical to 0.10.0 apart from the Bartlett haversine/spherical
  pairwise change described under "Engine".
- Edge cases that were silently wrong: exact duplicates at a zero cutoff keep
  their cross terms, antipodal points with a cutoff above the maximum
  distance are accepted, Bartlett boundary weights can no longer come out
  as -2e-16, and a raster lattice at a tiny or zero cutoff keeps its self
  terms.

## New validation and errors

- Non-finite or missing `lat`, `lon`, `unit`, `time` values, latitudes
  outside [-90, 90], and longitudes outside [-180, 360] are rejected in R
  before anything reaches C++; the engine repeats the checks (including
  scores and cutoff) for its other front-ends. A `NA` time used to split
  co-located observations into separate blocks.
- Both `lat` and `latitude` (or any two aliases) present in the data is now
  an auto-detection error, as documented.
- The fixest method accepts `data =` either as the full original data (rows
  addressed through `fixest::obs()`) or as a frame already aligned with the
  fit, such as the fitted model frame after `subset =` or NA removal. A frame
  with as many rows as the original data is always read as original data (so
  subsets that only permute rows stay correct); any other row count is an
  error.
- Unknown arguments in `...` are errors (typos such as `psd_fxi =` were
  swallowed). `maxobsmem` warns once that it is ignored; `neighbor = "band"`
  is documented as deprecated but still works.
- `lfe` fits made through the k-class path (any `kclass =` argument, even
  `kclass = 1`) are rejected: lfe stores the raw endogenous regressors for
  them instead of the projected 2SLS design, so the sandwich would be silently
  wrong (a `kclass = 1` fit gave 0.0028 where 2SLS gives 0.0097). Only OLS and
  ordinary 2SLS are supported.
- `ncores` must be a single finite positive integer. The default remains all
  logical cores, but `getOption("fastconley.ncores")` overrides it and a
  true-ish `_R_CHECK_LIMIT_CORES_` (CRAN's check convention) caps it at 2.
- Long computations can be interrupted: the engine polls a front-end hook on
  the calling thread once per block (`Rcpp::checkUserInterrupt()` in R,
  `SF_poll` in the Stata plugin) and joins its workers before raising.
- Thread-pool construction is exception-safe (a failed `std::thread`
  constructor used to terminate the process), `ncores` is capped at 1024,
  and the sub-metre-cutoff cell cap no longer collapses all points into one
  cell.
- Raster geometry is validated (`dlat`, `dlon` finite and positive, finite
  `lat0` in range): a zero longitude step used to spin forever.

## Engine

- Haversine weights are computed from the cached 3D unit vectors
  (`sin^2(theta/2) = |u_i - u_j|^2 / 4`, distance `2R asin(sqrt a)`), and the
  spherical Bartlett weight uses the same asin form instead of `acos(dot)`,
  which is ill-conditioned at small angles. No per-pair `sin`/`atan2` calls
  remain, which speeds up the pairwise engine everywhere (about 2x on Linux,
  6x on mingw-built Windows binaries, whose libm implements those functions
  in x87 microcode). Results move at the 1e-14 relative level versus 0.10.0
  on Bartlett haversine/spherical configurations.
- The C++ engine lives in a front-end-agnostic header, `src/conley_core.h`
  (namespace `conley`, `CONLEY_CORE_VERSION` 0.11.1, bumped on every change
  that alters results or entry points), shared with the Stata plugin. It
  compiles without R against plain Armadillo under C++14
  (`tests/manual/core_standalone_check.cpp` + `.R` verify bit-identity).
- RcppParallel dependency dropped: threading is a `std::thread` pool inside
  the engine with the same deterministic chunked reduction, so `ncores`
  still never affects results. `Depends: R (>= 4.0)` declared for the
  posix-threads Rtools toolchain on Windows.
- Dead code from the pre-chord screens removed; the balanced path builds its
  coordinate cache from the first period only.

## Tests and tooling

- testthat grew from 48 to 106 expectations: independent dense references
  for `fepois`/`feglm`, felm and fixest IV, the serial HAC, and pixel
  aggregation; 1-vs-2-core `identical()` checks for every engine path; and
  regression tests for each item above.
- `tests/manual/bitwise-battery.R` gained `write-ref`/`check` modes with a
  committed reference (`bitwise-battery-reference.txt`), and
  `tests/manual/core_edge_probes.cpp` exercises the edge cases in C++.
- A GitHub Actions `R CMD check --as-cran` workflow runs on Linux, macOS,
  and Windows.
- Documentation reconciled with the code: Conley (1999) throughout,
  `method = "auto"` as the default, chord versus great-circle wording, the
  balanced-panel error (not a fallback), and the pixel-aggregation distance
  bound (`sqrt(2) * pixel` for a pair, not `pixel / 2`).

## Stata port

- `stata/` holds `fastconley`, a reghdfe-family command with a pure-Mata
  engine and a compiled plugin built from `src/conley_core.h` for Linux,
  Windows, and macOS (see `stata/README.md` and `stata/CHANGELOG.md`), and
  `stata/upstream/` two alternative reghdfe proposals: a generic
  `vce(external PROVIDER, ...)` hook with fastconley as its first provider
  (primary), and a self-contained native `vce(conley ...)` patch.

CRAN preparation.

- `DESCRIPTION`: package names quoted, Conley (1999) DOI added, copyright
  holder role declared; `LICENSE` updated to match.
- `.Rbuildignore` keeps development-only files (the agent instructions,
  `notes/`, `tests/manual/`, `stata/`, `.github/`, and `cran-comments.md`)
  out of the source tarball.
- Tests: the tiny `felm` fixture in the balanced-panel validation test no
  longer emits "NaNs produced" warnings.

# fastconley 0.10.0

GLM support and documented IV support.

## `fixest::feglm()` / `fixest::fepois()` fits supported

`vcovSpHAC.fixest` now accepts GLM fits (any `feglm` family, including
`fepois`). The variance is the M-estimation sandwich `H^{-1} B H^{-1}`,
built from the maximum-likelihood score matrix and inverse Hessian that
`fixest` stores on every (non-`lean`) fit — the same construction
`fixest`'s own `vcov_conley()` uses for GLMs, verified against it at the
distance-formulation tolerance and against an exact same-distance
yardstick at ~1e-15. No estimation flag is needed (`demeaned = TRUE` is
only required for `feols`); weights, offsets, and the fixed-effect
profiling are already folded into the stored scores. Because the scores
ride through the existing engines unchanged, everything composes: pairwise
and grid/FFT engines, pixel aggregation, `ssc`, `psd_fix`, and — beyond
what `fixest` offers — the panel spatial **+ serial** HAC via
`lag_cutoff`, now available for Poisson/GLM panels. `lean = TRUE` fits
(no stored scores) and `femlm()`/`feNmlm()` fits are rejected with clear
errors.

## IV/2SLS support documented and tested

IV fits have in fact always produced the correct 2SLS Conley sandwich
through both methods — `lfe` and `fixest` store the projected
(second-stage) design in `cX` / `X_demeaned` and the structural residuals
in `residuals`, which is exactly what the score construction needs. This
is now documented and covered by the validation suite: `felm` IV and
`feols` IV agree with each other and with the exact yardstick at ~1e-15,
including weighted IV, multiple endogenous regressors, and IV panels with
serial HAC. See `tests/manual/test-glm-iv-parity.R`.

# fastconley 0.9.0

fixest feature parity: weighted fits, small-sample correction, PSD repair,
lat/lon auto-detection. **Breaking:** `ssc` and `psd_fix` default to `TRUE`
to match `fixest`'s defaults out of the box; pass
`ssc = FALSE, psd_fix = FALSE` to reproduce earlier fastconley versions and
`rbluhm/conley` bit-for-bit.

## Weighted (WLS) fits supported

`vcovSpHAC` now accepts weighted `felm()` and `feols()` fits. The meat
scores become `s_i = w_i * e_i * x_i` and the bread `(X'WX)^{-1}` — the
formula `fixest`'s own weighted Conley vcov uses (verified exactly against
it with a self-pairs-only cutoff, rel. err ~1e-15). Weights enter only the
scores and the bread, so every engine — pairwise, grid/FFT, serial HAC,
pixel aggregation, balanced CSR reuse — works unchanged. Note `lfe` stores
`sqrt(w)` on the fit; `vcovSpHAC` squares it back.

## `ssc`: small-sample correction (default `TRUE`)

Scales the variance matrix by `n / (n - K)`, with `K` counting all
estimated parameters including absorbed fixed-effect levels (taken from
the fit's residual degrees of freedom). This is exactly `fixest`'s default
Conley correction — its cluster adjustment (`G.adj` / `cluster.adj`) is a
no-op for Conley vcovs, so this one factor reproduces `fixest` defaults.
`ssc = FALSE` applies no correction, matching `rbluhm/conley` and previous
fastconley versions.

## `psd_fix`: positive semi-definite repair (default `TRUE`)

Conley spatial kernels do not guarantee a PSD variance matrix. With
`psd_fix = TRUE` (default, as in `fixest`) negative eigenvalues are
clamped to 1e-16 — the same semantics as `fixest`'s `vcov_fix` — with a
warning when the fix noticeably changed the matrix (> 1e-8). With
`psd_fix = FALSE` the matrix is returned as computed and a warning is
emitted when it is noticeably non-PSD.

## lat/lon auto-detection

`lat` and `lon` now default to `NULL` and are auto-detected from the
data's column names (`lat`/`latitude` and `lon`/`long`/`longitude`/`lng`,
case-insensitive exact matches). A message reports the pick; ambiguous or
missing matches error with instructions.

## Validation

- 24-config battery (balanced/unbalanced/cross-section x kernels x
  distances x lags x engines) bitwise-identical to v0.8.0 with
  `ssc = FALSE, psd_fix = FALSE` (the C++ engines and the default-off code
  paths are untouched).
- Weighted fits: exact match to `fixest` at self-pairs-only cutoff;
  ~1e-10 against a dense brute-force reference at 200 km for both
  kernels; felm and fixest paths agree to ~1e-10.
- New testthat coverage: `tests/testthat/test-weights-ssc.R` (14 tests).

# fastconley 0.8.0

Grid engine, part two: bartlett support (ring-FFT) and dateline wrap.

## `method = "grid"` / `"auto"` now covers `kernel = "bartlett"`

On a lattice the bartlett weight varies with the longitude offset, so the
per-ring-pair inner sum is a true 1D convolution rather than a boxcar.
`FastGridMeat` computes it via FFT (`arma::fft`): per ring pair, the
even-symmetric weight vector's (real) spectrum multiplies cached per-ring
score spectra, with one inverse FFT per target ring. Score spectra are
cached for a sliding latitude band plus a cutoff halo, and the reduction
is deterministically chunked — results remain bit-identical across
`ncores`.

Weights use the same per-distance arithmetic as the pairwise engine
(`atan2` haversine, `acos` spherical, `sqrt` chord), and the same-cell
distance is hard-set to 0, so agreement with the pairwise engine is
~1e-15 (haversine) to ~1e-12 (spherical/chord — inherent conditioning of
acos/sqrt near zero distance, not algorithm error). The `"auto"` rule
uses an FFT-aware cost model for bartlett.

- C1-study config (2.25M cells at ~1.1 km, 250 km cutoff, ~1.8e11 pairs,
  bartlett/spherical): 2.9 s vs 710 s for the 16-core pairwise engine
  (**242x**); single-threaded (20.5 s) it still beats 16-core pairwise
  35x. The kernel's evenness halves the work (T(r2,r1) = T(r1,r2)', so
  only upper-triangle ring pairs are computed, with self-ring weights
  halved and the half-meat symmetrized).

## Dateline wrap (bug fix for global rasters)

v0.7.0's grid engine clamped longitude windows at the lattice edges, so
on a raster spanning the full 360° circle it silently missed pairs that
are close "the short way" across the dateline (observed ~5e-3 relative
error on a global test raster). The engine now detects when the accept
window reaches across the dateline gap and switches both kernels to
circular windows (modular prefix-sum arcs for uniform; circular
convolution with period `n_col_full` for bartlett) — exact, validated
against the pairwise engine. When wrap would be needed but the lon step
does not tile 360° evenly (no consistent circular lattice exists),
`method = "grid"` stops with an informative error and `method = "auto"`
falls back to the pairwise engine. Non-wrapping rasters are unaffected:
results are bitwise identical to v0.7.0 (verified on the 30-config
battery).

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
