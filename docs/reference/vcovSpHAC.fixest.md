# Spatial HAC variance-covariance matrix for fixest models

Supports
[`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
(including IV/2SLS) and
[`fixest::feglm()`](https://lrberge.github.io/fixest/reference/feglm.html)
/
[`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html)
fits.

## Usage

``` r
# S3 method for class 'fixest'
vcovSpHAC(
  reg,
  unit = NULL,
  time = NULL,
  lat = NULL,
  lon = NULL,
  kernel = c("bartlett", "uniform"),
  dist_fn = c("haversine", "spherical", "chord"),
  dist_cutoff = NULL,
  lag_cutoff = 0,
  verbose = FALSE,
  balanced_pnl = FALSE,
  ncores = NA,
  pixel = 0,
  neighbor = c("grid", "band"),
  csr_weight = c("double", "float"),
  method = c("auto", "pairwise", "grid"),
  ssc = TRUE,
  psd_fix = TRUE,
  data = NULL,
  ...
)
```

## Arguments

- reg:

  A fitted object of class "fixest": a
  [`feols()`](https://lrberge.github.io/fixest/reference/feols.html) fit
  (including IV) with `demeaned = TRUE`, or a
  [`feglm()`](https://lrberge.github.io/fixest/reference/feglm.html) /
  [`fepois()`](https://rdrr.io/pkg/lfe/man/fepois.html) fit (any family;
  `lean = TRUE` fits are rejected because they carry no score matrix).

- unit:

  Optional name of the panel unit variable. If `NULL` the call is
  treated as a cross-section (each row is its own unit, all rows share a
  single period — no serial HAC).

- time:

  Optional name of the time variable. Ignored when `unit` is NULL.

- lat:

  Name of the latitude variable. If `NULL` (default), auto-detected from
  the data's column names. See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- lon:

  Name of the longitude variable. If `NULL` (default), auto-detected.
  See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- kernel:

  Spatial kernel, either "bartlett" or "uniform".

- dist_fn:

  Distance function, one of "haversine", "spherical", "chord".

- dist_cutoff:

  Spatial cutoff in km.

- lag_cutoff:

  Serial HAC lag cutoff.

- verbose:

  Print progress messages.

- balanced_pnl:

  Whether the panel is balanced and unit locations are time-invariant.

- ncores:

  Number of cores for the C++/RcppParallel spatial and serial routines.

- pixel:

  Score-pre-aggregation cell size, in kilometres.

- neighbor:

  Neighbor-search strategy: "grid" (default) or "band". See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- csr_weight:

  Balanced-path bartlett weight storage: "double" (default) or "float".
  See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- method:

  Spatial meat engine: "auto" (default), "pairwise", or "grid". See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- ssc:

  Small-sample correction (`n / (n - K)` when `TRUE`, the default). See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- psd_fix:

  Clamp negative eigenvalues when `TRUE` (the default). See
  [`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md).

- data:

  Optional. The model frame to draw `lat`/`lon`/ `unit`/`time` from. If
  `NULL` (default), the data is recovered from the fit's call. Pass it
  explicitly if the original data has gone out of scope, or if you want
  to override.

- ...:

  Currently unused.

## Value

A variance-covariance matrix.

## Details

For `feols`, the fit must have been called with
`feols(..., demeaned = TRUE)` so that the centered design matrix
`X_demeaned` is stored on the fit object. Weighted fits are supported
(the scores carry the weights and the bread uses \\X'WX\\, matching
`fixest`'s own weighted Conley vcov). IV fits work out of the box:
`X_demeaned` holds the projected (second-stage) design and `residuals`
the structural residuals, which is exactly the 2SLS sandwich.

For `feglm` / `fepois`, no estimation flag is needed: the variance is
the M-estimation sandwich \\H^{-1} B H^{-1}\\, built from the
maximum-likelihood score matrix and inverse Hessian that `fixest` stores
on every (non-`lean`) fit. Weights, offsets, and the fixed-effect
profiling are already folded into the stored scores. This is the same
construction `fixest`'s own
[`vcov_conley()`](https://lrberge.github.io/fixest/reference/vcov_conley.html)
uses for GLMs — but with exact great-circle distances, and with the
serial-HAC panel extension available via `lag_cutoff` (which `fixest`
does not offer for Conley vcovs).

The returned matrix can be passed to `fixest`'s `vcov` argument. For the
usual `fixest` workflow, define a one-argument wrapper such as
`function(x) vcovSpHAC(x, ...)` and pass that function to
[`summary()`](https://rdrr.io/r/base/summary.html),
[`etable()`](https://lrberge.github.io/fixest/reference/etable.html), or
`feols(vcov = )`. The wrapper keeps the coordinate names, cutoffs, panel
variables, and optional `data = ` argument together.

## Examples

``` r
if (requireNamespace("fixest", quietly = TRUE)) {
  ## feols must be fit with demeaned = TRUE (the keepCX analogue).
  set.seed(1)
  cells <- expand.grid(lat = seq(40, 49.5, by = 0.5),
                       lon = seq(-10, 9.5, by = 0.5))
  cells$x <- rnorm(nrow(cells))
  cells$y <- 0.5 * cells$x + rnorm(nrow(cells))

  fit <- fixest::feols(y ~ x, data = cells, demeaned = TRUE)
  vcov_fc <- function(x) {
    vcovSpHAC(x, lat = "lat", lon = "lon",
              kernel = "uniform", dist_fn = "spherical",
              dist_cutoff = 200, ncores = 2, data = cells)
  }
  V <- vcov_fc(fit)
  sqrt(diag(V))

  ## The same wrapper can be used directly in fixest's vcov argument.
  fit_sum <- summary(fit, vcov = vcov_fc)
  sqrt(diag(fit_sum$cov.scaled))

  ## Poisson (fepois / feglm): no demeaned = TRUE needed — the stored
  ## ML scores and inverse Hessian are used directly.
  cells$cnt <- rpois(nrow(cells), exp(0.4 * cells$x))
  fit_pois <- fixest::fepois(cnt ~ x, data = cells)
  V_pois <- vcovSpHAC(fit_pois, lat = "lat", lon = "lon",
                      kernel = "uniform", dist_fn = "spherical",
                      dist_cutoff = 200, ncores = 2, data = cells)
  sqrt(diag(V_pois))
}
#> (Intercept)           x 
#>  0.03802984  0.03316774 
```
