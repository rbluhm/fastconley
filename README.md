# fastconley

Fast Conley (1997) spatial HAC standard errors for `lfe::felm()` models in R.

`fastconley` is a drop-in replacement for the spatial path of [`rbluhm/conley`](https://github.com/rbluhm/conley) that scales to large cross-sections and high-dimensional regressions. The public API (`vcovSpHAC()`) is unchanged; only the internals were rewritten.

The original `conley` package was written by Richard Bluhm with contributions from [Darin Christensen](https://github.com/darinchristensen/conley-se). Some key speed up ideas (sorting on distances and pruning) are from [Laurent Berge](https://github.com/lrberge/fixest) and [Christian Düben](https://github.com/cdueben/conleyreg).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rbluhm/fastconley")
```

`fastconley` installs alongside the original `conley` package (different `Package:` name), so you can keep both libraries loaded in different R sessions to compare results.

## Usage

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
  dist_fn      = "haversine",    # or "spherical", "chord", "flatearth"
  dist_cutoff  = 500,            # km
  lag_cutoff   = 0,              # serial-HAC lag (set > 0 for both)
  balanced_pnl = TRUE,           # set TRUE only if every period has the same units in the same order, with time-invariant coordinates
  ncores       = NA              # NA = use all cores via RcppParallel
)
```

`lfe::felm()` must be called with `keepCX = TRUE` so `cY` and `cX` are populated; without it, both `conley` and `fastconley` will fail to extract the centered design matrix.

## What changed vs. `rbluhm/conley`

The spatial meat used to be assembled by repeatedly building dense `n × n` distance matrices and doing a row-by-row sandwich update. `fastconley` replaces that with four pieces:

1. **Score accumulation.** The meat is computed via `S_i = e_i X_i`, `C_i = 0.5 S_i + Σ_{j>i, d_ij ≤ cutoff} w_ij S_j`, `meat = Σ_i S_i' C_i + transpose`. No `k × n` per-row temporaries.
2. **Lat/lon candidate pruning.** Within each time block, observations are sorted by latitude. The pair loop breaks once latitude separation exceeds the cutoff and applies a conservative longitude screen before computing the exact distance.
3. **CSR neighbor lists.** Pair weights are stored as `row_ptr` / `col_idx` / `weight`, not a dense matrix. With `balanced_pnl = TRUE` the CSR graph is built once from the first period and reused across all periods.
4. **Serial HAC in C++.** The per-unit serial-correlation loop now runs entirely in C++ as a single `FastSerialHacPanel` call, instead of an R `lapply` that called into C++ once per unit.

A `CoordCache` precomputes `cos(lat)`, `sin(lat)`, and (for `dist_fn = "chord"`) the 3D Cartesian coordinates once per call, so per-pair distance evaluation is cache-only — `chord` is now actually the cheapest of the curved-earth distances. CSR construction and meat accumulation are parallelised with `RcppParallel`.

## Benchmarks

Single-threaded, `kernel = "bartlett"`, `dist_fn = "haversine"`, `dist_cutoff = 500 km`, internal `vcovSpHAC()` time only (Rscript startup excluded).

**Spatial path** (`lag_cutoff = 0`):

| `n_unit` | `n_time` | `k` | `n_obs` | conley (s) | fastconley (s) | speedup | max VCOV diff |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 4 | 20 | 400 | 0.13 | 0.007 | 19× | 2.0e-15 |
| 800 | 4 | 20 | 3,200 | 1.49 | 0.029 | 51× | 3.1e-16 |
| 1,500 | 4 | 50 | 6,000 | 13.7 | 0.094 | 146× | 1.8e-16 |
| 3,000 | 4 | 50 | 12,000 | 52.7 | 0.305 | 173× | 7.9e-17 |
| 5,000 | 4 | 50 | 20,000 | 161 | 0.93 | 173× | 5.1e-17 |
| 8,000 | 4 | 50 | 32,000 | 373 | 2.33 | 160× | 2.8e-17 |

**Serial-HAC path** (`lag_cutoff = 1`, `balanced_pnl = TRUE`):

| `n_unit` | `n_time` | `k` | conley (s) | fastconley (s) | speedup |
|---:|---:|---:|---:|---:|---:|
| 1,500 | 4 | 20 | 3.87 | 0.068 | 57× |
| 3,000 | 4 | 20 | 12.55 | 0.195 | 64× |
| 8,000 | 4 | 20 | 78.32 | 1.179 | 66× |

The serial-only added cost (lag=1 minus lag=0) at `n_unit = 8,000, k = 20` is **1.4s** in `conley`, **0.015s** in `fastconley` — the panel-native C++ routine eliminates the per-unit R-to-C++ overhead.

**Distance functions**, `n_obs = 8,000`, `k = 5`, single-threaded internal seconds per call:

| `dist_fn` | seconds | note |
|---|---:|---|
| `flatearth` | 0.29 | 2D, cheapest |
| `chord` | 0.67 | now the fastest curved-earth distance (was the slowest before the 3D coord cache) |
| `haversine` | 0.71 | baseline |
| `spherical` | 0.84 | `acos` is more expensive than `atan2`; cached `sin(lat)` saves the per-pair sin |

Speedups grow with both `n` and `k`; max VCOV difference across the full 36-cell spatial grid is 2.06e-14 (machine precision). Memory peak at the largest case is ~0.4 GB (vs ~1.3 GB for the dense path).

## Compatibility with `rbluhm/conley`

Numerical results match upstream to machine precision for `dist_fn ∈ {haversine, spherical, chord}` and both kernels.

**Known divergence.** `dist_fn = "flatearth"` differs by `O(1e-4)` to `O(3e-3)` in absolute VCOV entries on typical panels. Upstream's `flatearth` distance is order-dependent — the cosine factor uses the latitude of whichever point appears first in input order. `fastconley` sorts by latitude before screening, so it uses the latitude of the lower-latitude point. Other distance functions are symmetric in their two arguments, so they match upstream regardless of iteration order. See [`inst/FAST_SPATIAL_NOTES.md`](inst/FAST_SPATIAL_NOTES.md).

**`balanced_pnl = TRUE` is now stricter.** Upstream silently produced wrong results when units had time-varying coordinates; `fastconley` errors with a clear message instead. Unequal block sizes still trigger a warning and fall back to the general path.

## Requirements

R ≥ 3.5 with a C++11 compiler and GNU make. Imports `data.table`, `lfe`, `Rcpp`, `RcppArmadillo`, `RcppParallel`.

## Citation

The estimator is from Conley, T. G. (1999). "GMM estimation with cross sectional dependence." *Journal of Econometrics* 92(1):1–45.

## License

MIT.
