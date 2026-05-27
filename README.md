# fastconley

Fast Conley (1997) spatial HAC standard errors for `lfe::felm()` and `fixest::feols()` models in R.

`fastconley` is a drop-in replacement for the spatial path of [`rbluhm/conley`](https://github.com/rbluhm/conley) that scales to large cross-sections and high-dimensional regressions. `vcovSpHAC()` still accepts `felm` fits with the same call signature it always did; numerical equivalence with upstream `conley` holds at machine precision for the three supported distance functions (the upstream-only `dist_fn = "flatearth"` option was dropped, and `pixel` was added as a new optional argument — see the [Compatibility](#compatibility-with-rbluhmconley) section for the full diff). `vcovSpHAC()` is also now an S3 generic with a `fixest` method.

The original `conley` package was written by Richard Bluhm with contributions from [Darin Christensen](https://github.com/darinchristensen/conley-se). Some key speed up ideas (sorting on distances and pruning) are from [Laurent Berge](https://github.com/lrberge/fixest) and [Christian Düben](https://github.com/cdueben/conleyreg).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rbluhm/fastconley")
```

`fastconley` installs alongside the original `conley` package (different `Package:` name), so you can keep both libraries loaded in different R sessions to compare results.

## Usage

`vcovSpHAC()` is an S3 generic with methods for `lfe::felm` and `fixest::feols` fits.

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
  ncores       = NA,             # NA = use all cores via RcppParallel
  pixel        = 0               # score-pre-aggregation cell size in km; 0 = exact dedupe only
)
```

`lfe::felm()` must be called with `keepCX = TRUE` so `cY` and `cX` are populated; without it, both `conley` and `fastconley` will fail to extract the centered design matrix.

### With `fixest::feols`

```r
library(fixest)
library(fastconley)

est <- feols(y ~ x1 + x2 | unit + time, data = panel, demeaned = TRUE)

V <- vcovSpHAC(
  est,
  unit         = "unit",
  time         = "time",
  lat          = "lat",
  lon          = "lon",
  kernel       = "bartlett",
  dist_fn      = "haversine",
  dist_cutoff  = 500,
  lag_cutoff   = 0,
  balanced_pnl = TRUE,
  ncores       = NA,
  pixel        = 0
)
```

`fixest::feols()` must be called with `demeaned = TRUE` so the centered design matrix is stored on the fit object — analogous to `keepCX = TRUE` for `felm`. The fixest method returns exactly the same VCOV as the felm method on identical data; the dispatch only changes how the centered design and residuals are extracted. Weighted fits (`feols(..., weights = w)`) are not currently supported. If you fit the model in one frame and the data is no longer reachable from the call site, pass `data = ` explicitly.

This wires `fastconley` in as a drop-in spatial HAC for the `fixest` workflow while still computing the textbook within-period spatial + within-unit serial-HAC object (see the panel-comparison section below); it is structurally different from `fixest::vcov_conley + vcov_NW − vcov_hetero`.

### The `pixel` argument

`pixel` controls score pre-aggregation, inspired by the same argument in `fixest::vcov_conley`. Before the pair loop runs, observations that share a (possibly pixelated) location within a time period are collapsed into a single aggregated score, so the pair loop runs on `n_cases` rows instead of all `n`.

- `pixel = 0` (default): exact-coordinate dedupe only. Two rows are collapsed iff their `(lat, lon)` match exactly. The answer is unchanged (machine precision); free win when units share coordinates (panels with time-invariant coords, point-aggregated data).
- `pixel > 0`: snap `(lat, lon)` to a uniform `pixel`-km grid before the dedupe. Nearby points collapse to a common representative. The pair distance is approximate up to roughly `pixel / 2` km — this is a speed/accuracy trade-off; pick `pixel` an order of magnitude smaller than `dist_cutoff`.

## What changed vs. `rbluhm/conley`

The spatial meat used to be assembled by repeatedly building dense `n × n` distance matrices and doing a row-by-row sandwich update. `fastconley` replaces that with:

1. **Score accumulation.** The meat is computed via `S_i = e_i X_i`, `C_i = 0.5 S_i + Σ_{j>i, d_ij ≤ cutoff} w_ij S_j`, `meat = Σ_i S_i' C_i + transpose`. No `k × n` per-row temporaries.
2. **Lat/lon candidate pruning.** Within each time block, observations are sorted by latitude. The pair loop breaks once latitude separation exceeds the cutoff and applies a conservative longitude screen before computing the exact distance.
3. **CSR neighbor lists.** Pair weights are stored as `row_ptr` / `col_idx` / `weight`, not a dense matrix. With `balanced_pnl = TRUE` the CSR graph is built once from the first period and reused across all periods.
4. **Serial HAC in C++.** The per-unit serial-correlation loop now runs entirely in C++ as a single `FastSerialHacPanel` call, instead of an R `lapply` that called into C++ once per unit.
5. **Templated `(dist_fn, kernel)` dispatch.** The per-pair screen + distance + kernel-weight stack is specialised at compile time for each of the six `(distance, kernel)` combinations, so the inner loop has no runtime branches.
6. **3D dot-product threshold for spherical/chord.** With the unit-vector cache precomputed once, `(spherical, uniform)` pair inclusion is three multiplies + two adds + one compare, no trig at all. The `bartlett` variants only pay `acos` (spherical) or `sqrt` (chord) on accepted pairs.
7. **Row-major scores + sorted-order layout.** Scores `s_i = e_i · X_i` are stored row-major in a flat buffer, and both the coord cache and the scores are permuted into latitude-sorted order before the worker runs, so every read in the hot loop is sequential.

A `CoordCache` precomputes `cos(lat)`, `sin(lat)`, and (for `dist_fn ∈ {chord, spherical}`) the 3D unit-vector coordinates once per call. CSR construction and meat accumulation are parallelised with `RcppParallel`.

## Benchmarks

Numerical results match upstream `rbluhm/conley` to machine precision on `dist_fn ∈ {haversine, spherical, chord}` × both kernels × cross-section / balanced / unbalanced regimes; against upstream `conley` the spatial path is roughly two orders of magnitude faster (the original CSR refactor was already 150-170× on n_obs = 8,000–32,000 panels, and the recent templating + dot-product + score-permutation work adds another ~5×). The live benchmark below is against `fixest`.

### vs `fixest::vcov_conley`

Post-estimation only, `kernel = "uniform"`, great-circle distance, 500 km cutoff, k=10. Both sides fit the model once outside the timer; only the VCOV call is wallclock-measured (`proc.time()[["elapsed"]]`, min of 3 reps). `fixest` flags disabled for a clean core comparison: `ssc(adj = FALSE, cluster.adj = FALSE)`, `vcov_fix = FALSE`.

**Cross-section, single-threaded** (CONUS-sized box, uniform random coordinates). `fastconley` and `fixest` compute identical VCOVs at `pixel = 0` (parity to machine precision on a brute-force reference):

| n | fc `pixel=0` | fc `pixel=25` | fixest `pixel=0` | fixest `pixel=25` |
|---:|---:|---:|---:|---:|
| 5,000   | 0.024 s | 0.015 s | 0.071 s | 0.050 s |
| 20,000  | 0.281 s | 0.083 s | 1.077 s | 0.650 s |
| 50,000  | 1.68 s  | 0.180 s | 6.68 s  | 3.91 s  |
| 100,000 | 6.92 s  | 0.250 s | 27.8 s  | 15.6 s  |

At `pixel = 0`, `fastconley` is 3-4× faster than `fixest`. At `pixel = 25`, the gap is 27-62× because per-period dedupe + the dot-product threshold + the sequential score layout compound.

**Cross-section, 4 threads** (`RcppParallel::setThreadOptions(numThreads = 4)`, `setFixest_nthreads(4)`):

| n | fc `pixel=0` | fc `pixel=25` | fixest `pixel=0` | fixest `pixel=25` |
|---:|---:|---:|---:|---:|
| 5,000   | 0.012 s | 0.010 s | 0.026 s | 0.018 s |
| 20,000  | 0.093 s | 0.039 s | 0.345 s | 0.207 s |
| 50,000  | 0.562 s | 0.124 s | 3.12 s  | 1.84 s  |
| 100,000 | 3.23 s  | 0.169 s | 12.3 s  | 7.25 s  |

**Multicore scaling**, n = 50,000, `pixel = 0`, uniform/spherical:

| ncores | fc (s) | fx (s) | fc speedup vs 1-thread | fx speedup vs 1-thread |
|---:|---:|---:|---:|---:|
| 1  | 1.75 | 6.78 | 1.00× | 1.00× |
| 2  | 0.88 | 3.73 | 1.98× | 1.82× |
| 4  | 0.72 | 3.07 | 2.44× | 2.21× |
| 8  | 0.46 | 2.37 | 3.83× | 2.86× |
| 16 | 0.33 | 1.79 | 5.27× | 3.79× |

`fastconley` scales **better than `fixest`** past two threads — the per-pair work is heavy enough that RcppParallel/TBB scheduling pays off. The fastconley/fixest ratio widens from 3.9× at 1 thread to 5.4× at 16 threads.

### Panel comparison

This is a structural comparison, not a numerical one. `fastconley` computes the textbook panel Spatial HAC = within-period spatial Conley + within-unit Newey-West serial (Hsiang 2010; Christensen & Fetzer 2015; the contract upstream `rbluhm/conley` always honored). `fixest::vcov_conley` is a cross-section object; its author's recommended panel SHAC composition is `vcov_conley + vcov_NW − vcov_hetero`, which sums over cross-period spatial pairs that fastconley excludes. The two VCOVs differ by `O(s_i·s_jᵀ)` on every same-unit cross-time pair, so the numerical answers disagree by design. The fair comparison is wallclock for "the panel SHAC each tool delivers":

`ncores = 1`, balanced panel, time-invariant coordinates, `n_time = 4`, `lag = 1`, `kernel = "uniform"`, `dist_fn = "spherical"`, `pixel = 0`:

| n_unit | n_obs | fc (s) | fixest sum-of-3 (s) | fc vs fx |
|---:|---:|---:|---:|---:|
| 500    | 2,000   | 0.014 | 0.007 | 0.5× (fx faster) |
| 2,000  | 8,000   | 0.047 | 0.021 | 0.45× |
| 8,000  | 32,000  | 0.322 | 0.199 | 0.62× |
| 20,000 | 80,000  | 1.03  | 1.13  | 1.1× (fc faster) |

`fixest` wins on small-to-medium panels because `pixel = 0` exact-coordinate dedupe collapses each unit's `n_time` rows into a single case (coordinates are time-invariant), so `fixest` runs a `n_unit × n_unit` spatial sum once; `fastconley` honors the within-period rule and runs a `n_unit × n_unit` spatial sum `n_time` times plus a serial term. The lines cross around n_obs ≈ 80,000. See [`inst/FIXEST_COMPARISON.md`](inst/FIXEST_COMPARISON.md) for the full structural comparison and the literature pointers.

## Compatibility with `rbluhm/conley`

Numerical results match upstream to machine precision for `dist_fn ∈ {haversine, spherical, chord}` and both kernels.

**`dist_fn = "flatearth"` was removed.** Upstream's equirectangular formula is not symmetric in `(i, j)` and produces order-dependent results; `fastconley` no longer accepts it. Use `dist_fn = "spherical"` instead — with the 3D dot-product threshold it runs without trig in the inner loop and is now no slower than the old `flatearth` path was.

**`balanced_pnl = TRUE` is now stricter.** Upstream silently produced wrong results when units had time-varying coordinates; `fastconley` errors with a clear message instead. Unequal block sizes still trigger a warning and fall back to the general path.

## Requirements

R ≥ 3.5 with a C++11 compiler and GNU make. Imports `data.table`, `lfe`, `Rcpp`, `RcppArmadillo`, `RcppParallel`.

## Citation

The estimator is from Conley, T. G. (1999). "GMM estimation with cross sectional dependence." *Journal of Econometrics* 92(1):1–45.

## License

MIT.
