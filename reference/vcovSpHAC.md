# Spatial HAC variance-covariance matrix

Computes Conley (1999) spatial HAC variance-covariance matrices for
models estimated with
[`lfe::felm()`](https://rdrr.io/pkg/lfe/man/felm.html) (OLS and IV/2SLS)
or with `fixest`'s
[`feols()`](https://lrberge.github.io/fixest/reference/feols.html) (OLS
and IV),
[`feglm()`](https://lrberge.github.io/fixest/reference/feglm.html), and
[`fepois()`](https://rdrr.io/pkg/lfe/man/fepois.html). For GLM fits the
variance is the M-estimation sandwich built from the stored score matrix
and inverse Hessian. The spatial meat uses a fast CSR/cumulative-score
implementation in C++; see
[`vcovSpHAC.felm`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.felm.md)
and
[`vcovSpHAC.fixest`](https://rbluhm.github.io/fastconley/reference/vcovSpHAC.fixest.md)
for the per-method argument lists.

## Usage

``` r
vcovSpHAC(reg, ...)
```

## Arguments

- reg:

  A fitted model object.

- ...:

  Method-specific arguments.

## Value

A numeric matrix (base `"matrix"`) of dimension `k x k`, where `k` is
the number of estimated coefficients (absorbed fixed effects excluded):
the Conley spatial HAC estimate of the variance-covariance matrix of the
coefficient estimates, i.e. the sandwich `bread %*% meat %*% bread` with
the kernel-weighted spatial (and, with `lag_cutoff > 0`, serial) cross
products in the meat. Row and column names are the coefficient names of
the fit, so the matrix can be passed wherever a `vcov` is expected, e.g.
`lmtest::coeftest(reg, vcov = V)` or `summary(reg, vcov = V)` for fixest
fits, and `sqrt(diag(V))` gives the standard errors. The matrix is
symmetric; with the default `ssc = TRUE` it is scaled by `n / (n - K)`,
and with the default `psd_fix = TRUE` it is positive semi-definite.
