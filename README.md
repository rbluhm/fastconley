# fastconley

Fast Conley (1997) spatial HAC standard errors for `lfe::felm()` (OLS and IV/2SLS) and `fixest` models — `feols()` (OLS and IV), `feglm()`, and `fepois()` — in R.

Documentation site: **<https://rbluhm.github.io/fastconley/>** — function reference, changelog, and the [performance vignette](https://rbluhm.github.io/fastconley/articles/fastconley-performance.html).

`fastconley` is a drop-in replacement for the spatial path of [`rbluhm/conley`](https://github.com/rbluhm/conley) that scales to large cross-sections and high-dimensional regressions. `vcovSpHAC()` still accepts `felm` fits with the same call signature it always did; numerical equivalence with upstream `conley` holds at machine precision for the three supported distance functions (the upstream-only `dist_fn = "flatearth"` option was dropped, and `pixel` was added as a new optional argument — see the [Compatibility](#compatibility-with-rbluhmconley) section for the full diff). `vcovSpHAC()` is also now an S3 generic with a `fixest` method.

The original `conley` package was written by Richard Bluhm with contributions from [Darin Christensen](https://github.com/darinchristensen/conley-se). Some key speed up ideas (sorting on distances and pruning) are from [Laurent Berge](https://github.com/lrberge/fixest) and [Christian Düben](https://github.com/cdueben/conleyreg).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rbluhm/fastconley")
```

`fastconley` installs alongside the original `conley` package (different `Package:` name), so you can keep both libraries loaded in different R sessions to compare results.

### Performance build (optional)

The package ships with portable compiler flags (R's default `-O2`). The templated
inner loops benefit measurably from wider vectorization; for a machine-local
performance build, add to `~/.R/Makevars` before installing:

```make
CXXFLAGS += -O3 -march=native -funroll-loops
```

Results are numerically identical either way; this only affects speed.

## Usage

`vcovSpHAC()` is an S3 generic with methods for `lfe::felm` and `fixest` fits (`feols`, `feglm`, `fepois`).

### With `lfe::felm`

```r
library(lfe)
library(fastconley)

est <- felm(y ~ x1 + x2 | unit + time, data = panel, keepCX = TRUE)

V <- vcovSpHAC(
  est,
  unit         = "unit",
  time         = "time",
  lat          = "lat",
  lon          = "lon",
  kernel       = "bartlett",     # or "uniform"
  dist_fn      = "haversine",    # or "spherical", "chord"
  dist_cutoff  = 500,            # km
  lag_cutoff   = 0,              # serial-HAC lag (set > 0 for both)
  balanced_pnl = TRUE,           # set TRUE only if every period has the same units in the same order, with time-invariant coordinates
  ncores       = NA,             # NA = use all cores (std::thread pool in C++)
  pixel        = 0               # score-pre-aggregation cell size in km; 0 = exact dedupe only
)
```

`lfe::felm()` must be called with `keepCX = TRUE` so `cY` and `cX` are populated; without it, both `conley` and `fastconley` will fail to extract the centered design matrix.

IV/2SLS fits — the `felm` multi-part formula with `(endog ~ instruments)` — work out of the box and produce the correct 2SLS Conley sandwich: `lfe` stores the projected (second-stage) design in `cX` and the structural residuals in `residuals`, which is exactly what the score construction needs. Verified against the `fixest` IV path at machine precision (see `tests/manual/test-glm-iv-parity.R`).

### With `fixest::feols`

```r
library(fixest)
library(fastconley)

est <- feols(y ~ x1 + x2 | unit + time, data = panel, demeaned = TRUE)
```

`fixest::feols()` must be called with `demeaned = TRUE` so the centered design matrix is stored on the fit object — analogous to `keepCX = TRUE` for `felm`. The fixest method returns exactly the same VCOV as the felm method on identical data; the dispatch only changes how the centered design and residuals are extracted. Weighted fits (`feols(..., weights = ~w)` / `felm(..., weights = w)`) are supported since v0.9.0: the meat scores carry the weights and the bread uses `X'WX`, the same formula `fixest`'s own weighted Conley vcov uses. If you fit the model in one frame and the data is no longer reachable from the call site, pass `data = ` explicitly.

Since v0.9.0 three conveniences mirror `fixest`, including its defaults: `ssc = TRUE` (the default) applies the `n / (n − K)` small-sample correction — exactly `fixest`'s default Conley correction; its cluster adjustment is a no-op for Conley vcovs — `psd_fix = TRUE` (the default) clamps negative eigenvalues with `fixest::vcov_fix` semantics, and `lat` / `lon` are auto-detected from common column names when omitted. **Note the behavior change:** earlier fastconley versions and upstream `rbluhm/conley` apply no correction; pass `ssc = FALSE, psd_fix = FALSE` to reproduce their output bit-for-bit.

`fixest` accepts external covariance functions in its `vcov` argument. The most natural way to use `fastconley` in a `fixest` workflow is to define a one-argument wrapper that records the spatial settings and returns `vcovSpHAC()`:

```r
vcov_fastconley <- function(x) {
  vcovSpHAC(
    x,
    unit         = "unit",
    time         = "time",
    lat          = "lat",
    lon          = "lon",
    kernel       = "bartlett",
    dist_fn      = "haversine",
    dist_cutoff  = 500,
    lag_cutoff   = 1,
    balanced_pnl = TRUE,
    ncores       = NA,
    pixel        = 0,
    data         = panel
  )
}

summary(est, vcov = vcov_fastconley)
etable(est, vcov = vcov_fastconley)
```

The same wrapper can be passed directly at estimation time:

```r
est <- feols(
  y ~ x1 + x2 | unit + time,
  data = panel,
  demeaned = TRUE,
  vcov = vcov_fastconley
)
```

If you will print or export the same model many times, compute the covariance matrix once and pass the matrix to `fixest`:

```r
V <- vcov_fastconley(est)
summary(est, vcov = V)
etable(est, vcov = V)
```

This wires `fastconley` in as a drop-in spatial HAC for the `fixest` workflow while still computing the textbook within-period spatial + within-unit serial-HAC object (see the panel-comparison section below); it is structurally different from `fixest::vcov_conley + vcov_NW − vcov_hetero`.

IV fits (`feols(y ~ x1 | fe | x2 ~ z)`) work out of the box: with `demeaned = TRUE`, `fixest` stores the projected (second-stage) design and the structural residuals, which is exactly the 2SLS sandwich. The felm and feols IV paths agree with each other at machine precision.

### With `fixest::feglm` / `fixest::fepois` (since v0.10.0)

GLM fits need no estimation flag at all — `demeaned = TRUE` is only a `feols` requirement:

```r
est <- fepois(conflict_events ~ temp + precip | cell + year, data = panel)

V <- vcovSpHAC(
  est,
  unit = "cell", time = "year", lat = "lat", lon = "lon",
  kernel = "bartlett", dist_fn = "haversine",
  dist_cutoff = 500, lag_cutoff = 2, data = panel
)
```

The variance is the M-estimation sandwich `H⁻¹ B H⁻¹`, built from the maximum-likelihood score matrix and inverse Hessian that `fixest` stores on every fit (`lean = TRUE` fits are rejected — they carry no scores). Weights, offsets, and the fixed-effect profiling are already folded into the stored scores, and the scores ride through every engine unchanged: pairwise, grid/FFT for rasters, `pixel` aggregation, and — beyond what `fixest::vcov_conley()` offers — the panel spatial **+ serial** HAC via `lag_cutoff`, now available for Poisson/GLM panels (count outcomes à la conflict or mortality, PPML gravity). Any `feglm` family works; `femlm()`/`feNmlm()` fits are not supported.

### Matching settings to your data

Most users can leave `method = "auto"`. The explicit settings below are useful when you want the call to document the intended data layout, or when you want an early error if the data do not have the structure you expected.

| Data layout | Recommended `vcovSpHAC()` settings | Meaning |
|---|---|---|
| Scattered cross-section | Omit `unit` and `time`; use `lag_cutoff = 0`; leave `method = "auto"` or set `method = "pairwise"`. | Exact spatial Conley for ordinary point locations. |
| Repeated locations | Same as scattered data; start with `pixel = 0`. | Rows with identical coordinates are combined exactly before the covariance is computed. |
| Approximate location snapping | Set a positive `pixel`, such as `pixel = 5`, only after checking sensitivity. | Nearby locations are snapped to a kilometer grid. This can be faster, but it is an approximation. |
| Regular raster or grid | Use `method = "auto"`, or `method = "grid"` if you want to require a regular latitude-longitude lattice. | Exact Conley for complete or partly occupied regular rasters. If `method = "grid"` and the data are not a lattice, the call errors instead of falling back. |
| Balanced panel with fixed locations | Provide `unit` and `time`; set `lag_cutoff` as needed; set `balanced_pnl = TRUE`. | Spatial correlation is computed within period; serial HAC is computed within unit. This is faster when every period has the same units and locations do not move. |
| Unbalanced panel or moving locations | Provide `unit` and `time`; set `lag_cutoff` as needed; leave `balanced_pnl = FALSE`. | Uses the general panel route. This is safer when units enter, exit, duplicate within a period, or change coordinates. |

### The `pixel` argument

`pixel` controls score pre-aggregation, inspired by the same argument in `fixest::vcov_conley`. Before the pair loop runs, observations that share a (possibly pixelated) location within a time period are collapsed into a single aggregated score, so the pair loop runs on `n_cases` rows instead of all `n`.

- `pixel = 0` (default): exact-coordinate dedupe only. Two rows are collapsed iff their `(lat, lon)` match exactly. The answer is unchanged (machine precision); free win when units share coordinates (panels with time-invariant coords, point-aggregated data).
- `pixel > 0`: snap `(lat, lon)` to a uniform `pixel`-km grid before the dedupe. Nearby points collapse to a common representative. The pair distance is approximate up to roughly `pixel / 2` km — this is a speed/accuracy trade-off; pick `pixel` an order of magnitude smaller than `dist_cutoff`.

### Raster / gridded data: the grid engine

If your unit of observation is a raster cell — PRIO-GRID, gridded
population or nightlights, climate cells, anything you aggregated to a
regular lat/lon grid yourself — the spatial meat is computed by a
dedicated engine whose cost is independent of the pair count. On dense
rasters with realistic cutoffs this is two to three orders of magnitude
faster than enumerating pairs (see the benchmark sections below), with
no approximation: the same accept threshold and weights as the pairwise
engine, exact for natively gridded data.

You don't have to do anything to use it: `method = "auto"` (the default)
detects the lattice and engages the engine when a flop-balance estimate
says it wins. A typical workflow from a raster:

```r
library(terra); library(lfe); library(fastconley)

r <- rast("nightlights.tif")                      # any regular lat/lon raster
d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)    # x = lon, y = lat, holes dropped
d$gdp <- ...                                      # your cell-level variables

fit <- felm(log(lights) ~ log(gdp), data = d, keepCX = TRUE)
V <- vcovSpHAC(fit, lat = "y", lon = "x",
               kernel = "bartlett", dist_fn = "spherical",
               dist_cutoff = 250, data = d)
```

**Requirements and limits:**

- Coordinates must lie on a regular lat/lon lattice: all unique
  latitudes share a common step, all unique longitudes share a common
  step (the two steps may differ, e.g. 0.5° × 0.25°). Raster cell
  centers from `terra`/`raster` satisfy this by construction.
- **Holes are fine.** The lattice only has to describe the coordinate
  system; occupancy can be arbitrary — country masks, coastlines,
  land-only cells, missing data. Empty cells cost a little speed, never
  correctness; at extreme sparsity `auto` switches back to the pairwise
  engine on its own.
- Panels work (the engine runs per period) and need not be balanced;
  duplicates within a cell are handled.
- Both kernels (`uniform`, `bartlett`) and all three `dist_fn`s are
  supported. Rasters spanning the full 360° longitude circle wrap
  correctly across the dateline.
- **Not supported:** scattered point data (households, firms, survey
  clusters) and projected coordinates (UTM meters) — these always use
  the pairwise engine. A `pixel`-style snapping mode to let scattered
  points ride the grid engine is on the roadmap.
- Coordinates must match the lattice to within ~1e-6 of the step. If
  you computed cell centroids yourself with floating-point noise beyond
  that, detection fails silently and `auto` falls back to the (correct
  but slower) pairwise engine.

To see which engine actually ran, pass `verbose = TRUE` — it prints
`pairwise` or `grid`. If you believe your data is gridded but the grid
engine doesn't engage, run once with `method = "grid"`: instead of
falling back it errors with the specific reason (no lattice detected,
degenerate lattice, or a dateline-wrap incompatibility).

## What changed vs. `rbluhm/conley`

The spatial meat used to be assembled by repeatedly building dense `n × n` distance matrices and doing a row-by-row sandwich update. `fastconley` replaces that with:

1. **Score accumulation.** The meat is computed via `S_i = e_i X_i`, `C_i = 0.5 S_i + Σ_{j>i, d_ij ≤ cutoff} w_ij S_j`, `meat = Σ_i S_i' C_i + transpose`. No `k × n` per-row temporaries.
2. **Lat/lon candidate pruning.** Within each time block, observations are sorted by latitude. The pair loop breaks once latitude separation exceeds the cutoff and applies a conservative longitude screen before computing the exact distance.
3. **CSR neighbor lists.** Pair weights are stored as `row_ptr` / `col_idx` / `weight`, not a dense matrix. With `balanced_pnl = TRUE` the CSR graph is built once from the first period and reused across all periods.
4. **Serial HAC in C++.** The per-unit serial-correlation loop now runs entirely in C++ as a single `FastSerialHacPanel` call, instead of an R `lapply` that called into C++ once per unit.
5. **Templated `(dist_fn, kernel)` dispatch.** The per-pair screen + distance + kernel-weight stack is specialised at compile time for each of the six `(distance, kernel)` combinations, so the inner loop has no runtime branches.
6. **3D dot-product threshold for spherical/chord.** With the unit-vector cache precomputed once, `(spherical, uniform)` pair inclusion is three multiplies + two adds + one compare, no trig at all. The `bartlett` variants only pay `acos` (spherical) or `sqrt` (chord) on accepted pairs.
7. **Row-major scores + sorted-order layout.** Scores `s_i = e_i · X_i` are stored row-major in a flat buffer, and both the coord cache and the scores are permuted into latitude-sorted order before the worker runs, so every read in the hot loop is sequential.

A `CoordCache` precomputes `cos(lat)`, `sin(lat)`, and (for `dist_fn ∈ {chord, spherical}`) the 3D unit-vector coordinates once per call. CSR construction and meat accumulation are parallelised with a small `std::thread` pool inside the engine header (`src/conley_core.h`, shared with the Stata port); the chunked reduction is deterministic, so results are bit-identical for any `ncores`.

## Benchmarks

Numerical results match upstream `rbluhm/conley` to machine precision on `dist_fn ∈ {haversine, spherical, chord}` × both kernels × cross-section / balanced / unbalanced regimes; against upstream `conley` the spatial path is roughly two orders of magnitude faster (the original CSR refactor was already 150-170× on n_obs = 8,000–32,000 panels, and the recent templating + dot-product + score-permutation work adds another ~5×). The live benchmark below is against `fixest`.

### v0.5.0: cell-grid neighbor search

Since v0.5.0 candidate enumeration uses a 3D cell grid by default (`neighbor = "grid"`); the tables below were measured with the v0.4 latitude-band scan and remain its reference. On top of those numbers the grid gives (single-threaded, `FastSpatialMeat`, k = 10): **8.6–25×** on a global-extent 1M-point cross-section at 100 km, **4.5×** on CONUS 100k at 100 km, **2.1–2.3×** on CONUS 100k at 500 km (true-pair accumulation dominates there), **3.6×** on a 200k-unit balanced CSR build, and 16-thread scaling improves from ~5.3× to ~7.9×. The gain grows with geographic extent and shrinking cutoff because grid candidate work is proportional to *accepted* pairs (~3–4 candidates per accept) rather than to the latitude band (`~2·L_lon/(π·r)` per accept). Pair sets and weights are identical to the band scan; results differ only by floating-point summation order.

### v0.6.0: deterministic results + memory diet

Since v0.6.0 results are **bit-identical across `ncores` values** (deterministic chunked reduction), the C++ entry points alias R memory with zero input copies (peak RSS −27 to −33% on large runs), haversine gains an exact dot-product pre-screen (bartlett/haversine ~1.4× single-threaded), and `csr_weight = "float"` optionally halves balanced-path bartlett weight storage.

### v0.7.0: exact grid engine for rasters (`method = "auto"/"grid"`)

For the uniform kernel on a regular lat/lon lattice (raster data), the meat
reduces to exact sliding-window sums whose cost is independent of the pair
count. On 2.25M cells at ~1.1 km with a 250 km cutoff (~1.8×10¹¹ pairs):
**0.81 s vs 244 s** for the 16-core pairwise engine (302×, error 5e-15); a
9M-cell continental raster takes 3.5 s. `method = "auto"` (default) engages
it only when a lattice is detected and the kernel is uniform — scattered
data and bartlett stay on the pairwise engine. No approximation: natively
gridded data gets the exact great-circle answer.

### v0.8.0: bartlett on the grid engine (ring-FFT) + dateline wrap

The grid engine now covers `kernel = "bartlett"`: on a lattice the
per-ring-pair inner sum becomes a 1D convolution (the weight varies with
the longitude offset), computed via FFT with cached per-ring score
spectra. Same C1-study raster (2.25M cells, 250 km, ~1.8×10¹¹ pairs),
bartlett/spherical: **2.9 s vs 710 s** for the 16-core pairwise engine
(**242×**, error ~1e-13); even single-threaded (20.5 s) it beats the
16-core pairwise engine 35×. Agreement with the pairwise engine is
~1e-15 (haversine) to ~1e-12 (spherical/chord — inherent acos/sqrt
conditioning at near-zero distances, not algorithm error), and results
remain bit-identical across `ncores`.

v0.8.0 also fixes a v0.7.0 bug: on rasters spanning the full 360°
longitude circle the grid engine clamped windows at the lattice edge and
silently missed pairs across the dateline (~5e-3 relative error on a
global test raster). Windows now wrap correctly for both kernels; if wrap
would be needed but the lon step doesn't tile 360° evenly, `method =
"grid"` errors informatively and `method = "auto"` falls back to the
pairwise engine.

### vs `fixest::vcov_conley` (v0.6.1, large data)

Post-estimation only, `kernel = "uniform"`, great-circle distance, k = 10,
`pixel = 0`. Both sides fit the model once outside the timer; only the VCOV
call is wallclock-measured (min over reps). `fixest` 0.14.1 with
`distance = "spherical"`, `ssc(adj = FALSE, cluster.adj = FALSE)`,
`vcov_fix = FALSE`. Machine: i7-1360P (16 threads), 32 GB. Earlier
(v0.4-era, n ≤ 100k) tables live in the git history; at those sizes the gap
was 3–4×.

**Accuracy note.** `fastconley` matches a brute-force great-circle uniform
meat to ~1e-15. `fixest` 0.14.1's `distance = "spherical"` deviates from the
brute-force reference by ~2e-2 in the VCOV on these configs (its distance
formula flips pairs near the cutoff boundary), so the two tools' results
agree only to ~1e-2–1e-3. The workload is the same; `fastconley` is the
numerically exact one.

**Cross-section** (CONUS-sized box, uniform random coordinates):

| n | cutoff | threads | fastconley | fixest | speedup | fc / fx peak RSS |
|---:|---:|---:|---:|---:|---:|---:|
| 100,000   | 500 km | 1  | 2.53 s  | 27.2 s  | **10.7×** | 205 / 140 MB |
| 100,000   | 500 km | 16 | 0.40 s  | 4.30 s  | **10.8×** | |
| 250,000   | 500 km | 1  | 15.8 s  | 170.8 s | **10.8×** | 388 / 223 MB |
| 250,000   | 500 km | 16 | 2.34 s  | 45.8 s  | **19.6×** | |
| 500,000   | 100 km | 1  | 3.14 s  | 77.4 s  | **24.7×** | 629 / 393 MB |
| 500,000   | 100 km | 16 | 0.66 s  | 24.8 s  | **37.8×** | |
| 1,000,000 | 100 km | 1  | 12.1 s  | 300.7 s | **24.9×** | 1.20 / 0.64 GB |
| 1,000,000 | 100 km | 16 | 2.16 s  | 98.8 s  | **45.8×** | |

The gap widens with n, with shrinking cutoff, and with threads: grid
candidate enumeration is proportional to *accepted* pairs (~3–4 candidates
per accept regardless of extent), and since v0.6.1 the prep path (sort,
coordinate cache, gathers) is parallel too. `fixest` holds a smaller peak
RSS (fastconley carries the score matrix plus one sorted copy); both stay
far below dense-matrix approaches.

### Panel comparison

This is a structural comparison, not a numerical one. `fastconley` computes the textbook panel Spatial HAC = within-period spatial Conley + within-unit Newey-West serial (Hsiang 2010; Christensen & Fetzer 2015; the contract upstream `rbluhm/conley` always honored). `fixest::vcov_conley` is a cross-section object; its author's recommended panel SHAC composition is `vcov_conley + vcov_NW − vcov_hetero`, which sums over cross-period spatial pairs that fastconley excludes. The two VCOVs differ by `O(s_i·s_jᵀ)` on every same-unit cross-time pair, so the numerical answers disagree by design. The fair comparison is wallclock for "the panel SHAC each tool delivers":

Balanced panel, time-invariant coordinates, 500 km cutoff, `lag = 1`,
`kernel = "uniform"`, `dist_fn = "spherical"`, `pixel = 0` (fixest timed as
the full sum-of-3):

| n_unit | T | n_obs | threads | fastconley | fixest sum-of-3 | fc vs fx |
|---:|---:|---:|---:|---:|---:|---:|
| 50,000  | 4  | 200,000   | 1  | 2.35 s | 6.82 s | **2.9×** |
| 50,000  | 4  | 200,000   | 16 | 1.13 s | 1.05 s | 0.9× (fx faster) |
| 150,000 | 10 | 1,500,000 | 1  | 30.9 s | 60.9 s | **2.0×** |
| 150,000 | 10 | 1,500,000 | 16 | 10.9 s | 15.7 s | **1.4×** |

`fastconley` honors the within-period rule, so it runs the `n_unit × n_unit`
spatial sum once *per period* plus a serial term, while `fixest`'s
exact-coordinate dedupe collapses the panel into a single spatial pass —
fastconley does ~T× more spatial work by construction and still wins on
single-threaded wallclock at scale (the v0.4-era crossover at n_obs ≈ 80k is
long gone). The 150k × 10 run holds the balanced CSR in RAM (~3.9 GB at
500 km / 6.4×10⁸ pairs); use `csr_weight = "float"` (bartlett) or a tighter
cutoff if memory-bound. See [`inst/FIXEST_COMPARISON.md`](inst/FIXEST_COMPARISON.md) for the full structural comparison and the literature pointers.

## Compatibility with `rbluhm/conley`

Numerical results match upstream to machine precision for `dist_fn ∈ {haversine, spherical, chord}` and both kernels.

**`dist_fn = "flatearth"` was removed.** Upstream's equirectangular formula is not symmetric in `(i, j)` and produces order-dependent results; `fastconley` no longer accepts it. Use `dist_fn = "spherical"` instead — with the 3D dot-product threshold it runs without trig in the inner loop and is now no slower than the old `flatearth` path was.

**`balanced_pnl = TRUE` is now stricter.** Upstream silently produced wrong results when units had time-varying coordinates; `fastconley` errors with a clear message instead. Unequal block sizes still trigger a warning and fall back to the general path.

## Requirements

R ≥ 4.0 with a C++11 compiler (the engine uses `std::thread`, which on Windows needs the posix-threads Rtools that ship with R 4.0+). Imports `data.table` and `Rcpp` (and links against `RcppArmadillo`); `lfe` and `fixest` are suggested — you need whichever one fits your models.

## Citation

The estimator is from Conley, T. G. (1999). "GMM estimation with cross sectional dependence." *Journal of Econometrics* 92(1):1–45.

## License

MIT.
