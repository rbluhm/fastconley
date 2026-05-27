# Comparison with `fixest::vcov_conley`

This note records what was learned from a head-to-head benchmark and source
review of `fastconley::vcovSpHAC()` against `fixest::vcov_conley()` (Bergé,
v0.13.2), and documents the conceptual difference that surfaces when both are
applied to panel data.

## Numerical agreement

On a pure cross-section with `kernel = "uniform"`, `dist_fn = "haversine"`
(fastconley) vs `distance = "spherical"` (fixest), no SSC, no PSD-fix:

| n | max abs VCOV diff |
|---:|---:|
| 500 | 1.3e-18 |
| 1,500 | 6.2e-6 |
| 5,000 | 5.0e-6 |
| 60,000 | 4.0e-7 |

The residual disagreement at larger n appears to be a minor cutoff-boundary or
distance-formula difference; small enough to call the same answer in practice.

## Conceptual difference on panels

When both functions are pointed at a balanced panel with time-invariant
coordinates per unit, the VCOVs do **not** agree. The reason is structural,
not numerical.

- **`fastconley`** computes the spatial meat **within each time period only**:
  the spatial sum runs over pairs `(i, j)` with `time_i == time_j` and
  `d_geo_ij <= dist_cutoff`. Same-unit cross-time observations sit at
  geographic distance 0 but never enter the spatial sum because they live in
  different time blocks. Within-unit correlation across time is the job of a
  separate `lag_cutoff` argument (Newey-West within each unit's time series),
  which adds a serial-HAC meat to the spatial one.

- **`fixest::vcov_conley`** is defined on a cross-section. If you hand it a
  panel, every same-unit cross-time pair sits at distance 0, passes the
  cutoff, and contributes with full weight 1. There is no within-period
  restriction.

This is not a `fixest` bug. `fixest::vcov_conley()` documents Conley (1999)
as a cross-sectional spatial VCOV. For panel Spatial HAC, lrberge/fixest#177
records the intended decomposition:

```text
Spatial HAC = spatial Conley + panel Newey-West - heteroskedastic robust
```

That subtraction matters because both the spatial and serial components include
the diagonal/robust contribution.

## Which decomposition is the textbook one?

`fastconley`'s split — spatial Conley meat within each time period, plus a
Newey-West-style serial term within each unit, summed — is the standard reading
of Conley (1999) extended to panels. Supporting citations:

- **Conley (1999), *J. Econometrics* 92:1-45.** The estimator is defined for
  a cross-section indexed by location `s ∈ Z^d`. There is no time index. Any
  "panel Conley" SE is by definition an extension; the 1999 paper itself does
  not specify how repeated observations at the same point should be handled.

- **Hsiang (2010, *PNAS*), `ols_spatial_HAC.ado`.** Canonical implementation.
  The spatial loop is keyed on time period (`for ti = 1; ti <= Ntime; ti++`),
  distances computed only among same-period rows; the serial piece is a
  separate outer loop over panel units with Bartlett weights on the time gap,
  explicitly excluding same-time observations to prevent double-counting.

- **Christensen & Fetzer (2015).** R port of Hsiang's code; same structure;
  exposes `lag_cutoff` and `dist_cutoff` as independent dimensions in
  separate loops.

- **Düben, `conleyreg`** (CRAN). Manual: *"`lag_cutoff` defaults to 0, meaning
  that standard errors are only adjusted cross-sectionally."* Parallelises the
  spatial code in C++ and the serial part in R as separate code paths.

- **Bergé, `fixest`.** The user-facing `vcov_conley` documents Conley (1999)
  for cross-sections. For panel Spatial HAC, issue #177 gives the additive
  decomposition `spatial + HAC - robust`.

The fastconley behavior is also what upstream `rbluhm/conley` has always done:
the R wrapper does `lapply(timeUnique, ...)` and `Reduce("+", ...)` over
per-period meats; the C++ `XeeXhC` only ever receives one period's rows. The
fork preserved that contract — it didn't invent it.

**Verdict:** `fastconley` is on the textbook side. The cross-section result
matches `fixest::vcov_conley` to machine precision; the panel disagreement is
because plain `vcov_conley` is the wrong object for stacked panel data. For a
fair panel comparison, construct the `fixest` VCOV as
`vcov_conley + vcov_NW - vcov_hetero`.

## Benchmark recipe for panel Spatial HAC

For future benchmarks, use the same estimation sample, same centered design
matrix, and the same `unit`, `time`, `lat`, and `lon` variables on both sides.
On the `fixest` side, compute:

```r
ssc_conley <- fixest::ssc(K.adj = FALSE)
ssc_nw <- fixest::ssc(K.adj = FALSE, G.adj = FALSE)

V_spatial <- fixest::vcov_conley(
  est_fx,
  lat = ~lat,
  lon = ~lon,
  cutoff = dist_cutoff,
  pixel = pixel,
  distance = "spherical",
  ssc = ssc_conley,
  vcov_fix = FALSE
)

V_serial <- fixest::vcov_NW(
  est_fx,
  unit = ~unit,
  time = ~time,
  lag = lag_cutoff,
  ssc = ssc_nw,
  vcov_fix = FALSE
)

V_robust <- fixest::vcov_hetero(
  est_fx,
  ssc = ssc_conley,
  vcov_fix = FALSE
)

V_panel_spatial_hac <- V_spatial + V_serial - V_robust
```

Benchmark settings to keep fixed:

- `fastconley` must use `kernel = "uniform"` for apples-to-apples timing,
  because `fixest::vcov_conley()` uses a uniform spatial kernel.
- Use `distance = "spherical"` in `fixest` when comparing to
  `dist_fn = "haversine"` or `dist_fn = "spherical"` in `fastconley`;
  `fixest` defaults to `"triangular"`, which is cheaper but different.
- Use `pixel = 0` for exact/parity benchmarks. Use matched positive `pixel`
  values only for approximate speed/accuracy benchmarks.
- Pass `lag = lag_cutoff`; otherwise `fixest::vcov_NW()` chooses its own
  default lag from the number of periods.
- Disable small-sample corrections and PSD repair when comparing raw matrices:
  `ssc(K.adj = FALSE)` for Conley/robust, `ssc(K.adj = FALSE, G.adj = FALSE)`
  for panel Newey-West, and `vcov_fix = FALSE` throughout.
- Plain `fixest::vcov_conley(est_fx, ...)` on a stacked panel is still useful
  as a diagnostic, but it is not a fair panel Spatial-HAC benchmark because
  same-unit cross-time pairs enter the spatial sum at distance zero.

## Why `fixest` is faster on cross-sections

Single-thread cross-section benchmark (CONUS-sized box, k=10, 500 km cutoff,
uniform kernel, haversine/spherical) had `fixest` ~1.3-1.7× faster than
`fastconley`, gap growing with n; the gap narrows to ~1.18× at 8 threads.

The dominant reason from a source review is **score pre-aggregation**: before
the pair loop, `fixest` deduplicates rows with identical `(lat, lon)` and sums
their scores. The pair loop then runs over `n_cases` unique coordinates rather
than `N` raw rows. With the `pixel > 0` argument the same mechanism quantises
nearby points into a single representative, dropping `n_cases` further at the
cost of approximating the distance. (`fixest/src/misc_funs.cpp:124` for the
tapply-sum; `fixest/R/VCOV.R:2385-2392` for the dedupe call.)

Secondary reasons:

- `fixest` uses a plain `#pragma omp parallel for` over the outer row index.
  `fastconley` uses RcppParallel/TBB, which has higher per-task overhead.
  This is consistent with the observation that the gap narrows under many
  threads.
- `fixest` lays scores out row-major (`vcov_related.cpp:391-433`), making the
  inner K-loop tight and cache-friendly. `fastconley` reads X via Armadillo,
  which is column-major and strides across rows in the same loop.
- `fixest`'s default `distance = "triangular"` is a flat-earth
  equirectangular distance with no `sqrt`/`sin`/`asin` in the inner test —
  but this does not figure in our benchmark because we forced
  `distance = "spherical"`.

No SIMD intrinsics, no single-precision, no spatial index (KD-tree, grid).
The win is algorithm + memory layout.

## Action items

1. **Implement score pre-aggregation + a `pixel` argument** in `fastconley`.
   This is the single biggest gap. Dedup is a free win at `pixel = 0` when
   any duplicate coordinates exist; `pixel > 0` is a user-opt-in
   speed/accuracy trade-off.
2. Re-benchmark cross-section after (1).
3. Run the **panel** benchmark against the `fixest` composition
   `vcov_conley + vcov_NW - vcov_hetero`, using the parity settings above.
