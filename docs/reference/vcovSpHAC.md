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
