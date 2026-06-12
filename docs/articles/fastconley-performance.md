# Conley standard errors with fastconley

When regression errors are spatially correlated — nearby observations
hit by the same local shocks — conventional standard errors understate
the uncertainty of the estimates, often badly. Conley (1999) standard
errors fix this the same way Newey-West fixes serial correlation:
covariance between observations is admitted up to a chosen distance
cutoff and zeroed beyond it. If your regressors have spatial structure
(almost everything geographic does: climate, infrastructure,
institutions, treatments rolled out by region), this is the difference
between a significant and an insignificant result far more often than is
comfortable.

`fastconley` computes Conley spatial HAC standard errors for models
estimated with
[`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
or [`lfe::felm()`](https://rdrr.io/pkg/lfe/man/felm.html) — for
cross-sections and for panels, where it implements the full Conley
(1999) spatial-*temporal* HAC. Four properties distinguish it:

- **Exact.** Distances are true great-circle distances and every
  within-cutoff pair enters with its exact kernel weight — no
  dense-matrix shortcuts, no approximate distance formulas. Against a
  brute-force reference, results agree to near machine precision (the
  dense-baseline benchmark below shows this directly).
- **Complete for panels.** With `lag_cutoff > 0` the spatial adjustment
  within each period is combined with a Newey-West serial HAC term
  within each unit — the spatial-temporal correction of Conley (1999) as
  popularized by Hsiang (2010). `fixest`’s
  [`vcov_conley()`](https://lrberge.github.io/fixest/reference/vcov_conley.html)
  is spatial-only; it makes you choose between the Conley and the
  Driscoll-Kraay/Newey-West corrections, while here they enter one
  covariance matrix together.
- **Reproducible.** Results are bit-identical across thread counts: the
  same data and settings give the same covariance matrix on a laptop
  with `ncores = 2` and a server with `ncores = 64`.
- **Fast at any size.** Small problems are instant; a
  1,000,000-observation cross-section takes about a second, and regular
  rasters with hundreds of billions of within-cutoff pairs take seconds
  rather than hours (the benchmark sections quantify this).

The usual workflow is:

1.  Estimate the model.
2.  Pick a spatial cutoff that is meaningful for the application.
3.  Pass the fitted model and coordinate columns to
    [`vcovSpHAC()`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.md)
    (when the columns have standard names like `lat` / `lon`, they are
    auto-detected and can be omitted).
4.  Use the returned covariance matrix in
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`etable()`](https://lrberge.github.io/fixest/reference/etable.html),
    or your preferred table-making code.

The examples below use `fixest`, because many users already rely on its
formula syntax and reporting tools. The same
[`vcovSpHAC()`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.md)
call also works with `felm` models when the fit was created with
`keepCX = TRUE`.

## Why It Matters: A Small Example

The quickest way to see the point is to simulate data where both the
regressor and the error have spatial structure — the textbook case where
conventional standard errors mislead — and compare the inference. We
draw 1,500 locations in a box roughly 650 km across and give both `x`
and the error a smooth spatial component with an effective range of
about 150 km.

``` r
library(fastconley)
set.seed(42)

n <- 1500
d <- data.frame(lat = runif(n, 40, 46), lon = runif(n, 5, 11))

# Smooth spatial fields via an exponential-covariance Gaussian process.
rad <- pi / 180
D <- 6371 * acos(pmin(1,
  sin(d$lat * rad) %o% sin(d$lat * rad) +
  (cos(d$lat * rad) %o% cos(d$lat * rad)) * cos(outer(d$lon, d$lon, "-") * rad)))
L <- chol(exp(-D / 150) + diag(1e-8, n))
d$x <- drop(crossprod(L, rnorm(n))) + 0.5 * rnorm(n)
e   <- drop(crossprod(L, rnorm(n))) + 0.5 * rnorm(n)
d$y <- 1 + 0.5 * d$x + e
```

``` r
library(fixest)
est <- feols(y ~ x, data = d, demeaned = TRUE)

vcov_conley <- function(x) {
  vcovSpHAC(x, lat = "lat", lon = "lon", dist_cutoff = 150,
            kernel = "bartlett", dist_fn = "spherical",
            ncores = 2, data = d)
}

se_tab <- data.frame(
  `standard errors` = c("IID", "Heteroskedasticity-robust", "Conley (150 km)"),
  `se(x)` = sapply(list(vcov(est, vcov = "iid"),
                        vcov(est, vcov = "hetero"),
                        vcov_conley(est)),
                   function(V) sqrt(diag(V))["x"]),
  check.names = FALSE
)
knitr::kable(se_tab, digits = 4, row.names = FALSE)
```

| standard errors           |  se(x) |
|:--------------------------|-------:|
| IID                       | 0.0305 |
| Heteroskedasticity-robust | 0.0308 |
| Conley (150 km)           | 0.1281 |

The coefficient is the same in all three rows; only the uncertainty
changes. Ignoring the spatial correlation here understates the standard
error by a factor of about 4.2 — enough to turn an insignificant effect
“significant”. This is the problem Conley standard errors exist to
solve, and everything below is about computing them exactly and quickly.

## Quick Start

For `fixest`, fit the model with `demeaned = TRUE`. This keeps the
centered design matrix that `fastconley` needs after estimation.

`fixest` accepts an external covariance function in its `vcov` argument.
The most convenient pattern is to define a small function that records
the coordinate columns, cutoff, kernel, and thread count, then pass that
function to [`summary()`](https://rdrr.io/r/base/summary.html),
[`etable()`](https://lrberge.github.io/fixest/reference/etable.html), or
even to
[`feols()`](https://lrberge.github.io/fixest/reference/feols.html)
itself.

``` r
library(fixest)
library(fastconley)

est <- feols(
  outcome ~ treatment + x1 + x2 | county + year,
  data = d,
  demeaned = TRUE
)

vcov_fastconley <- function(x) {
  vcovSpHAC(
    x,
    lat = "lat",
    lon = "lon",
    dist_cutoff = 100,     # kilometers
    kernel = "uniform",
    dist_fn = "spherical",
    ncores = 8,
    data = d
  )
}

summary(est, vcov = vcov_fastconley)
etable(est, vcov = vcov_fastconley)
```

The same wrapper can be used at estimation time. This is useful when you
want the fitted object to carry the requested standard errors
immediately.

``` r
est <- feols(
  outcome ~ treatment + x1 + x2 | county + year,
  data = d,
  demeaned = TRUE,
  vcov = vcov_fastconley
)
```

You can also compute the covariance matrix once and pass the matrix to
`fixest`. That is preferable when you will print or export the same
model many times and do not want the covariance function to be
recomputed.

``` r
V_conley <- vcov_fastconley(est)
summary(est, vcov = V_conley)
etable(est, vcov = V_conley)
```

For `lfe`, use `keepCX = TRUE` when fitting the model.

``` r
library(lfe)
library(fastconley)

est <- felm(
  outcome ~ treatment + x1 + x2 | county + year,
  data = d,
  keepCX = TRUE
)

V_conley <- vcovSpHAC(
  est,
  lat = "lat",
  lon = "lon",
  dist_cutoff = 100,
  kernel = "uniform",
  dist_fn = "spherical",
  ncores = 8,
  data = d
)

sqrt(diag(V_conley))
```

For panel data, provide the unit and time variables. Set
`lag_cutoff > 0` when you also want serial HAC adjustment within units.

``` r
est <- feols(
  y ~ treatment + x1 + x2 | unit + year,
  data = panel_data,
  demeaned = TRUE
)

vcov_shac <- function(x) {
  vcovSpHAC(
    x,
    unit = "unit",
    time = "year",
    lat = "lat",
    lon = "lon",
    dist_cutoff = 500,
    lag_cutoff = 1,
    kernel = "bartlett",
    dist_fn = "spherical",
    balanced_pnl = TRUE,
    ncores = 8,
    data = panel_data
  )
}

summary(est, vcov = vcov_shac)
```

Use `balanced_pnl = TRUE` only when every period contains the same units
and each unit has the same coordinates in every period. Otherwise leave
it at the default, `FALSE`.

## Choosing Settings

The most important choice is `dist_cutoff`. It should come from the
research design: the distance over which residual correlation is
plausible enough to affect inference.

| Argument              | Practical meaning                                                                                                                                                                                                                                  |
|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dist_cutoff`         | Spatial range, in kilometers. Observations farther apart do not contribute to each other’s Conley adjustment.                                                                                                                                      |
| `kernel = "uniform"`  | All observations inside the cutoff receive equal weight.                                                                                                                                                                                           |
| `kernel = "bartlett"` | Observations inside the cutoff are downweighted with distance.                                                                                                                                                                                     |
| `lag_cutoff`          | Number of time lags for serial HAC in panel data. Use `0` for spatial-only inference.                                                                                                                                                              |
| `ncores`              | Number of threads used for the covariance calculation.                                                                                                                                                                                             |
| `pixel = 0`           | Exact deduplication of identical coordinates. Positive values trade a little spatial precision for more aggregation.                                                                                                                               |
| `ssc = TRUE`          | Small-sample correction. The default scales the matrix by $n/(n - K)$, counting absorbed fixed-effect levels in $K$ — exactly `fixest`’s default Conley correction. `FALSE` applies none, reproducing `rbluhm/conley` and fastconley before 0.9.0. |
| `psd_fix = TRUE`      | Conley kernels do not formally guarantee a positive semi-definite matrix. The default clamps negative eigenvalues (same semantics as `fixest`’s `vcov_fix`); either way a warning reports a noticeably non-PSD result.                             |

Two conveniences need no argument at all. Weighted regressions
(`feols(..., weights = ~w)` or `felm(..., weights = w)`) are handled
automatically — the scores and the bread pick up the weights from the
fit, matching `fixest`’s own weighted Conley vcov. And if `lat` / `lon`
are omitted, they are auto-detected from common column names
(`lat`/`latitude`, `lon`/`long`/`longitude`/`lng`, case-insensitive),
with a message reporting the pick.

Start with `pixel = 0`. If many observations share the same coordinates,
this is exact and often faster. If the data are very large and the
application can tolerate approximation, try a positive `pixel` value
that is much smaller than `dist_cutoff`, then compare the resulting
standard errors on a representative subset.

For `dist_fn`, the three options are all exact great-circle geometries
that agree to within a fraction of a percent at sub-1000 km cutoffs:
`haversine` is the numerically best-conditioned at short distances,
`spherical` (arc-cosine) is the cheapest, and `chord` is the
straight-line distance through the sphere. Any of them is fine; just
report which one you used. `lag_cutoff` is in units of the time variable
(periods, not years unless your time variable is years).

### Cutoff Sensitivity Is Cheap — Report It

The honest answer to “what cutoff should I use?” is usually a range, and
referees increasingly ask for it. Because the covariance call is fast,
looping over cutoffs costs almost nothing — there is no reason to report
a single number when you can report the curve:

``` r
cutoffs <- c(50, 100, 150, 300, 600)
se_x <- sapply(cutoffs, function(cut) {
  V <- vcovSpHAC(est, lat = "lat", lon = "lon", dist_cutoff = cut,
                 kernel = "bartlett", dist_fn = "spherical",
                 ncores = 2, data = d)
  sqrt(diag(V))["x"]
})
knitr::kable(data.frame(`cutoff (km)` = cutoffs, `se(x)` = round(se_x, 4),
                        check.names = FALSE), row.names = FALSE)
```

| cutoff (km) |  se(x) |
|------------:|-------:|
|          50 | 0.0723 |
|         100 | 0.1079 |
|         150 | 0.1281 |
|         300 | 0.1531 |
|         600 | 0.1455 |

In this simulation the errors were built with a ~150 km correlation
range, and the standard error climbs steeply until the cutoff covers
that range, then levels off — exactly the pattern to look for in real
data. If the standard errors keep growing with the cutoff, the spatial
correlation is longer-ranged than the design assumed (or there is a
spatial trend the model should absorb).

## Matching The Call To The Data

Most users can leave `method = "auto"`. The explicit settings below are
useful when you want the call to document the intended data layout, or
when you want an early error if the data do not have the structure you
expected.

| Data layout                          | Recommended [`vcovSpHAC()`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.md) settings         | What it means                                                                                                                                                                          |
|--------------------------------------|----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Scattered cross-section              | Omit `unit` and `time`; use `lag_cutoff = 0`; leave `method = "auto"` or set `method = "pairwise"`.      | Exact spatial Conley for point locations. This is the default for ordinary latitude-longitude data.                                                                                    |
| Repeated locations                   | Same as scattered data; start with `pixel = 0`.                                                          | Rows with identical coordinates are combined exactly before the covariance is computed.                                                                                                |
| Approximate location snapping        | Set a positive `pixel`, such as `pixel = 5`, only after checking sensitivity.                            | Nearby locations are snapped to a kilometer grid. This can be faster, but it is an approximation.                                                                                      |
| Regular raster or grid               | Use `method = "auto"`, or `method = "grid"` if you want to require a regular latitude-longitude lattice. | Exact Conley for complete or partly occupied regular rasters. If `method = "grid"` and the data are not a lattice, the call errors instead of falling back.                            |
| Balanced panel with fixed locations  | Provide `unit` and `time`; set `lag_cutoff` as needed; set `balanced_pnl = TRUE`.                        | Spatial correlation is computed within each period; serial HAC is computed within unit. The balanced setting is faster when every period has the same units and locations do not move. |
| Unbalanced panel or moving locations | Provide `unit` and `time`; set `lag_cutoff` as needed; leave `balanced_pnl = FALSE`.                     | Uses the general panel route. This is safer when units enter, exit, duplicate within a period, or change coordinates.                                                                  |

A few concrete templates are often enough.

For ordinary point data:

``` r
vcov_points <- function(x) {
  vcovSpHAC(
    x,
    lat = "lat",
    lon = "lon",
    dist_cutoff = 100,
    kernel = "uniform",
    dist_fn = "spherical",
    method = "pairwise",
    ncores = 8,
    data = d
  )
}
```

For a regular raster:

``` r
vcov_raster <- function(x) {
  vcovSpHAC(
    x,
    lat = "cell_lat",
    lon = "cell_lon",
    dist_cutoff = 250,
    kernel = "bartlett",
    dist_fn = "spherical",
    method = "grid",
    ncores = 8,
    data = raster_data
  )
}
```

For a balanced panel with one serial lag:

``` r
vcov_panel <- function(x) {
  vcovSpHAC(
    x,
    unit = "unit",
    time = "year",
    lat = "lat",
    lon = "lon",
    dist_cutoff = 500,
    lag_cutoff = 1,
    kernel = "bartlett",
    dist_fn = "spherical",
    balanced_pnl = TRUE,
    ncores = 8,
    data = panel_data
  )
}
```

If you are unsure whether a raster is being recognized, run once with
`verbose = TRUE`. If you want to force the regular-raster route, use
`method = "grid"`; a failure then tells you that the coordinates are not
on a regular latitude-longitude lattice.

## What Is Estimated

Let $x_{i}$ be the centered row of regressors for observation $i$, and
let ${\widehat{u}}_{i}$ be the residual. Define the score

$$s_{i} = {\widehat{u}}_{i}x_{i}.$$

For a cross-section, the Conley meat is

$$\widehat{\Omega} = \sum\limits_{i = 1}^{n}\sum\limits_{j = 1}^{n}K\left( \frac{d_{ij}}{c} \right)1\{ d_{ij} \leq c\} s_{i}\prime s_{j},$$

where $d_{ij}$ is the great-circle distance between observations $i$ and
$j$, $c$ is `dist_cutoff`, and $K( \cdot )$ is the chosen kernel:

$$K_{\text{uniform}}(u) = 1,\qquad K_{\text{bartlett}}(u) = 1 - u,\qquad u = d_{ij}/c \leq 1.$$

The uniform kernel weights every within-cutoff pair equally; the
bartlett kernel downweights linearly with distance, which makes a
non-positive semi-definite estimate much less likely (neither spatial
kernel formally guarantees one — by default negative eigenvalues are
clamped, with the same semantics as `fixest`’s `vcov_fix`; pass
`psd_fix = FALSE` to return the matrix as computed). The covariance
matrix is the usual sandwich:

$$\widehat{V} = (X\prime X)^{- 1}\widehat{\Omega}(X\prime X)^{- 1}.$$

For weighted (WLS) fits the scores become
$s_{i} = w_{i}{\widehat{u}}_{i}x_{i}$ and the bread uses $X\prime WX$,
matching `fixest`’s weighted Conley vcov.

For panel data, the spatial part is applied within time period. If
`lag_cutoff` is positive, `fastconley` adds a serial HAC term within
unit:

$$\widehat{\Omega} = \sum\limits_{t}\sum\limits_{i:t_{i} = t}\sum\limits_{j:t_{j} = t}K\left( \frac{d_{ij}}{c} \right)1\{ d_{ij} \leq c\} s_{i}\prime s_{j} + \sum\limits_{g}\sum\limits_{t}\sum\limits_{\ell = 1}^{L}\left( 1 - \frac{\ell}{L + 1} \right)\left( s_{g,t}\prime s_{g,t - \ell} + s_{g,t - \ell}\prime s_{g,t} \right),$$

where $L$ is `lag_cutoff`. In words: spatial correlation is handled
among units observed in the same period, and serial correlation is
handled along the same unit over time. Note the serial sum starts at
$\ell = 1$: the contemporaneous ($\ell = 0$) covariance already lives in
the spatial term, which includes each observation’s own score product.
This is why the panel object is *not* the same as applying a
cross-sectional Conley formula to the stacked rows.

Two computational guarantees are worth stating precisely. First, the
sums above are evaluated **exactly**: every pair with $d_{ij} \leq c$
enters with its exact kernel weight, using true great-circle distances —
there is no sampling, no tapering beyond the kernel itself, and no
approximate distance formula. Second, the result is **deterministic**:
summation is organized in fixed-size chunks, so the same data and
settings produce a bit-identical matrix regardless of `ncores`. Both
properties hold for every engine (`pairwise` and `grid`), which is why
`method` is purely a speed choice.

By default $\widehat{V}$ is scaled by the $n/(n - K)$ small-sample
correction — the same one `fixest` applies to Conley vcovs — so out of
the box the two packages’ defaults line up. Pass `ssc = FALSE` to report
$\widehat{V}$ exactly as defined above (the convention of
`rbluhm/conley` and fastconley versions before 0.9.0).

## Why It Scales

A direct Conley calculation is often described as if it required an
$n \times n$ spatial weight matrix. That becomes impossible quickly.
With 1,000,000 observations, a dense double-precision matrix would
require about 7.3 TiB of memory before any covariance algebra.

`fastconley` avoids that object. It works with the score rows that enter
the sandwich estimator and only spends time on observations that can
affect the answer. It also uses three common data features directly:

- If rows share coordinates, their scores can be combined exactly.
- If panel locations are fixed over time, the same spatial relationships
  can be reused across periods.
- If the data are a regular latitude-longitude raster, the calculation
  is done as a grid problem instead of a point-pair problem.

You normally do not have to choose among these paths. The default
settings are intended to do the right thing for ordinary point data,
panel data, and regular rasters.

## Benchmarks

The tables below report post-estimation covariance time. The model is
fitted once outside the timer, so the numbers answer a specific
question: after the regression is estimated, how long does it take to
compute the Conley or spatial HAC covariance matrix?

The benchmark files are shipped with the package. They are not rerun
when this vignette is built.

### Machine

| item             | value                                  |
|:-----------------|:---------------------------------------|
| run date         | 2026-06-10 17:01:07                    |
| CPU              | 13th Gen Intel(R) Core(TM) i7-1360P    |
| logical CPUs     | 16                                     |
| cores per socket | 12                                     |
| threads per core | 2                                      |
| memory           | Mem: 30Gi 23Gi 6.5Gi 585Mi 2.4Gi 7.8Gi |
| operating system | Ubuntu 24.04.4 LTS                     |
| R                | R version 4.6.0 (2026-04-24)           |
| fastconley       | 0.8.0                                  |
| fixest           | 0.14.1                                 |
| lfe              | 3.1.1                                  |

### One Million Observations

This cross-section has 1,000,000 observations, 10 regressors, a 100 km
cutoff, the uniform kernel, and spherical distance. A dense spatial
weight matrix would be about 7.3 TiB.

| observations | regressors | cutoff_km | threads | seconds | speedup_vs_2_threads | dense_matrix_size |
|:-------------|-----------:|----------:|--------:|:--------|:---------------------|:------------------|
| 1,000,000    |         10 |       100 |       2 | 2.180   | 1.0x                 | 7.3 TiB           |
| 1,000,000    |         10 |       100 |       4 | 1.553   | 1.4x                 | 7.3 TiB           |
| 1,000,000    |         10 |       100 |       8 | 1.405   | 1.6x                 | 7.3 TiB           |
| 1,000,000    |         10 |       100 |      16 | 1.131   | 1.9x                 | 7.3 TiB           |

On this machine, the 16-thread run computes the covariance matrix in
about 1.131 seconds.

### Comparison With fixest

`fixest` has its own Conley covariance implementation through
[`vcov_conley()`](https://lrberge.github.io/fixest/reference/vcov_conley.html).
The comparison below uses smaller cross-sections where both methods are
comfortable to run locally. Both models are fitted once outside the
timer; only the post-estimation covariance call is timed.

| observations | cutoff_km | threads | fastconley_seconds | fixest_seconds | speedup | max_abs_difference |
|:-------------|----------:|--------:|:-------------------|:---------------|:--------|:-------------------|
| 50,000       |       100 |       8 | 0.066              | 0.369          | 5.6x    | 1.03e-07           |
| 50,000       |       500 |       1 | 1.093              | 11.516         | 10.5x   | 3.79e-07           |
| 50,000       |       500 |       8 | 0.293              | 2.434          | 8.3x    | 4.71e-07           |
| 100,000      |       100 |       8 | 0.123              | 1.140          | 9.3x    | 4.27e-08           |
| 100,000      |       500 |       8 | 0.812              | 9.595          | 11.8x   | 2.37e-07           |

The last column needs an explanation, because the two packages do not
return the same matrix: on these configurations they disagree by roughly
one percent in relative terms. The disagreement is **not** a fastconley
approximation. A brute-force reference — forming the full great-circle
distance matrix and the exact double sum, feasible at a few thousand
observations — matches `fastconley` to near machine precision (about
1e-12; see the dense-baseline table below, which performs exactly this
check). The gap comes from `fixest` 0.14.x’s internal
`distance = "spherical"` formula, which flips a small share of pairs
near the cutoff boundary relative to the exact great-circle distance.
The workload being timed is the same; `fastconley` is the numerically
exact one.

#### Matching fixest’s numbers

If you are migrating from
[`vcov_conley()`](https://lrberge.github.io/fixest/reference/vcov_conley.html)
and want to reconcile results, two things matter:

1.  **Small-sample corrections.** Since v0.9.0 the defaults already
    agree: both packages scale by $n/(n - K)$, where $K$ counts all
    estimated parameters including absorbed fixed-effect levels
    (`fixest`’s cluster adjustment is a no-op for Conley vcovs). To
    compare *uncorrected* matrices instead, pass `ssc = FALSE` to
    [`vcovSpHAC()`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.md)
    and `ssc(adj = FALSE, cluster.adj = FALSE)` to `fixest`.
2.  **Distances.** After aligning the corrections, the remaining ~1e-2
    relative difference is fixest’s approximate spherical distance, as
    described above. It shrinks with the cutoff’s distance from the mass
    of pairs and is not something a setting can remove.

On panels the two packages compute structurally different objects (see
the Balanced Panel section below), so coefficients’ standard errors
should not be expected to match even after these alignments.

### Balanced Panel

The panel scenario has 50,000 units observed for 10 periods, for 500,000
rows and 10 regressors. The covariance includes within-period spatial
Conley adjustment and one serial HAC lag within unit.

| units  | periods | observations | regressors | cutoff_km | lag_cutoff | threads | seconds | dense_one_period_size |
|:-------|--------:|:-------------|-----------:|----------:|-----------:|--------:|:--------|:----------------------|
| 50,000 |      10 | 500,000      |         10 |       500 |          1 |      16 | 2.446   | 18.6 GiB              |

This is not the same calculation as applying a cross-sectional Conley
covariance to all stacked panel rows. Spatial correlation is restricted
to observations in the same period; serial correlation is handled
separately within unit.

### Regular Raster

Regular rasters are common in climate, remote-sensing, and gridded
economic data. The benchmark below uses a 1500 by 1500
latitude-longitude lattice, or 2.25 million cells, with a 250 km cutoff.

| kernel   | cells     | cutoff_km | threads | seconds | dense_matrix_size | implied_within_cutoff_pairs |
|:---------|:----------|----------:|--------:|:--------|:------------------|:----------------------------|
| bartlett | 2,250,000 |       250 |      16 | 4.241   | 36.8 TiB          | 1.80e+11                    |
| uniform  | 2,250,000 |       250 |      16 | 0.567   | 36.8 TiB          | 1.80e+11                    |

The important point for users is that this remains an exact Conley
calculation for the raster. It is just using the regular grid structure
instead of treating the data as millions of unrelated points.

For scale: treating the same lattice as scattered points and running the
pairwise engine on all ~1.8e11 within-cutoff pairs takes minutes, not
seconds — measured during v0.8.0 development on this machine at 16
threads, about 244 s for the uniform kernel and 710 s for bartlett,
against the sub-five-second grid-engine times in the table. That
two-orders-of-magnitude gap is why `method = "auto"` exists: gridded
data should never pay the point-pair price.

### Repeated Locations

Repeated coordinates appear in panels, administrative data, firm or
household microdata linked to towns, and gridded data with multiple
observations per cell. At `pixel = 0`, `fastconley` combines only
exactly identical coordinates, so the result is unchanged apart from
floating-point summation order.

| setting  | raw_rows | coordinate_cases | seconds | speed_vs_pixel0 | max_abs_difference_vs_pixel0 |
|:---------|:---------|:-----------------|:--------|:----------------|:-----------------------------|
| pixel=0  | 100,000  | 20,000           | 0.054   | 1.0x            | 0.00e+00                     |
| pixel=10 | 100,000  | 18,580           | 0.058   | 0.9x            | 4.54e-08                     |
| pixel=25 | 100,000  | 13,042           | 0.046   | 1.2x            | 1.20e-07                     |

Positive `pixel` values are optional. They can be useful for very large
data sets with dense nearby coordinates, but they should be treated as
an approximation and checked against `pixel = 0` on a smaller sample.

### Dense Baseline

For intuition, the small benchmark below actually forms the dense
comparison object. This is feasible at a few thousand rows and already
shows why the dense approach does not scale.

| observations | regressors | cutoff_km | dense_seconds | fastconley_seconds | speedup | max_abs_difference |
|:-------------|-----------:|----------:|:--------------|:-------------------|:--------|:-------------------|
| 1,000        |          5 |       500 | 0.019         | 0.002              | 9.5x    | 2.96e-12           |
| 2,000        |          5 |       500 | 0.058         | 0.002              | 29.0x   | 5.91e-12           |
| 4,000        |          5 |       500 | 0.539         | 0.003              | 179.7x  | 8.19e-12           |

## Reproducing The Benchmarks

The quick local benchmark suite is:

``` bash
FASTCONLEY_BENCH_THREADS=8 FASTCONLEY_BENCH_REPS=2 \
  Rscript inst/benchmarks/run-fastconley-benchmarks.R
```

The large-data benchmark suite is:

``` bash
FASTCONLEY_LARGE_THREADS=2,4,8,16 FASTCONLEY_LARGE_REPS=1 \
  Rscript inst/benchmarks/run-large-data-benchmarks.R
```

The large optional `fixest` comparison is disabled by default because it
can take several minutes. To run it, add `FASTCONLEY_LARGE_FIXEST=1`.

## References

Conley, T. G. (1999). GMM estimation with cross sectional dependence.
*Journal of Econometrics*, 92(1), 1–45.

Conley, T. G. (2010). Spatial econometrics. In S. N. Durlauf & L. E.
Blume (Eds.), *Microeconometrics* (pp. 303–313). Palgrave Macmillan.

Hsiang, S. M. (2010). Temperatures and cyclones strongly associated with
economic production in the Caribbean and Central America. *Proceedings
of the National Academy of Sciences*, 107(35), 15367–15372. (Source of
the within-period spatial + within-unit serial HAC panel estimator this
package implements.)

Newey, W. K., & West, K. D. (1987). A simple, positive semi-definite,
heteroskedasticity and autocorrelation consistent covariance matrix.
*Econometrica*, 55(3), 703–708.

Bergé, L. (2018). Efficient estimation of maximum likelihood models with
multiple fixed-effects: the R package FENmlm. *CREA Discussion Papers*,
13. (The `fixest` package, whose
[`vcov_conley()`](https://lrberge.github.io/fixest/reference/vcov_conley.html)
is benchmarked above and whose “On standard-errors” vignette is the
model for the small-sample-correction discussion.)
