# Spatial HAC variance-covariance matrix for felm() models

Spatial HAC variance-covariance matrix for felm() models

## Usage

``` r
# S3 method for class 'felm'
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
  maxobsmem = 50000L,
  data = NULL,
  ...
)
```

## Arguments

- reg:

  A fitted object of class "felm".

- unit:

  Optional name of the panel unit variable.

- time:

  Optional name of the time variable.

- lat:

  Name of the latitude variable. If `NULL` (default), it is
  auto-detected from the data's column names (`"lat"`, `"latitude"`,
  case-insensitive); a message reports the pick.

- lon:

  Name of the longitude variable. If `NULL` (default), it is
  auto-detected (`"lon"`, `"long"`, `"longitude"`, `"lng"`,
  case-insensitive).

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

  Score-pre-aggregation cell size, in kilometres. Default 0
  (exact-coordinate dedupe only). If \`pixel \> 0\`, points are snapped
  to a uniform \`pixel\`-km grid before the dedupe — a speed/accuracy
  trade-off that approximates the distance up to roughly \`pixel / 2\`.

- neighbor:

  Neighbor-search strategy for the spatial meat: "grid" (default; 3D
  cell grid, output-sensitive candidate enumeration) or "band" (latitude
  band scan, the pre-0.5.0 behavior). Both are exact and use identical
  per-pair accept tests; results agree to floating-point summation
  order.

- csr_weight:

  Storage precision for the balanced-path bartlett kernel weights:
  "double" (default, exact) or "float" (halves the per-pair weight
  memory; introduces at most ~6e-8 relative error per weight). Ignored
  for `kernel = "uniform"`, which stores no weights.

- method:

  Spatial meat engine. "pairwise" enumerates neighbor pairs (the default
  engine; works for any data). "grid" uses the exact grid-native meat —
  requires observations on a regular lat/lon lattice (e.g. raster data);
  cost is independent of the pair count, so it is dramatically faster on
  dense grids with large cutoffs. The uniform kernel uses sliding-window
  prefix sums; the bartlett kernel uses per-ring-pair FFT convolutions.
  Lattices spanning the full longitude circle wrap correctly across the
  dateline. "auto" (default) picks "grid" when it detects a lattice and
  a flop-balance estimate says it wins; both engines are exact, so the
  choice only affects speed (results agree to FP summation order, plus
  ~1e-12 acos conditioning for the bartlett spherical/chord weights).
  Pass `verbose = TRUE` to see which engine ran; if you expect gridded
  data to use the grid engine and it does not, run once with
  `method = "grid"` — it errors with the specific reason instead of
  falling back.

- ssc:

  Small-sample correction. If `TRUE` (default), the variance matrix is
  scaled by `n / (n - K)` where `K` counts all estimated parameters
  including absorbed fixed-effect levels (taken from the fit's residual
  degrees of freedom). This matches `fixest`'s default Conley correction
  (its cluster adjustment is a no-op for Conley vcovs). Pass `FALSE` for
  no correction — that reproduces `rbluhm/conley`, fastconley versions
  before 0.9.0, and `fixest` with
  `ssc(adj = FALSE, cluster.adj = FALSE)`.

- psd_fix:

  The spatial kernels do not guarantee a positive semi-definite variance
  matrix. If `TRUE` (default), negative eigenvalues are clamped (to
  1e-16, as `fixest`'s `vcov_fix` does) and a warning reports when the
  fix noticeably changed the matrix. If `FALSE`, the matrix is returned
  as computed, with a warning when it is not positive semi-definite.

- maxobsmem:

  Ignored by the fast spatial path. Kept for backward compatibility.

- data:

  Optional. The data frame to draw `lat`/`lon` from. If `NULL`
  (default), the data is recovered from the fit's call and aligned via a
  model-frame re-evaluation. Pass it explicitly if the original data has
  gone out of scope — and when no rows were dropped at fit time (no NAs,
  no `subset`), the coordinates are then taken by direct column access
  with no model-frame rebuild at all.

- ...:

  Currently unused.

## Value

A variance-covariance matrix.

## Examples

``` r
if (requireNamespace("lfe", quietly = TRUE)) {
  ## Cross-section on a regular 0.5-degree raster with holes. method =
  ## "grid" forces the exact grid engine; the default method = "auto"
  ## picks it automatically when the raster is large enough to win
  ## (on a toy example this small, pairwise is just as fast).
  set.seed(1)
  cells <- expand.grid(lat = seq(40, 49.5, by = 0.5),
                       lon = seq(-10, 9.5, by = 0.5))
  cells <- cells[sample(nrow(cells), 600), ]   # irregular occupancy
  cells$x <- rnorm(nrow(cells))
  cells$y <- 0.5 * cells$x + rnorm(nrow(cells))

  fit <- lfe::felm(y ~ x, data = cells, keepCX = TRUE)
  V <- vcovSpHAC(fit, lat = "lat", lon = "lon",
                 kernel = "bartlett", dist_fn = "spherical",
                 dist_cutoff = 200, ncores = 2, method = "grid",
                 data = cells)
  sqrt(diag(V))

  ## Panel with spatial + serial HAC (scattered points: pairwise engine)
  pnl <- data.frame(unit = rep(1:200, each = 5),
                    time = rep(1:5, times = 200),
                    lat = rep(runif(200, 40, 50), each = 5),
                    lon = rep(runif(200, -10, 10), each = 5))
  pnl$x <- rnorm(nrow(pnl))
  pnl$y <- 0.5 * pnl$x + rnorm(nrow(pnl))
  fit2 <- lfe::felm(y ~ x | unit + time, data = pnl, keepCX = TRUE)
  V2 <- vcovSpHAC(fit2, unit = "unit", time = "time",
                  lat = "lat", lon = "lon", kernel = "bartlett",
                  dist_fn = "haversine", dist_cutoff = 300,
                  lag_cutoff = 2, balanced_pnl = TRUE, ncores = 2,
                  data = pnl)
  sqrt(diag(V2))
}
#>          x 
#> 0.03549569 
```
