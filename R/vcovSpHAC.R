#' Spatial HAC variance-covariance matrix
#'
#' Computes Conley (1999) spatial HAC variance-covariance matrices for models
#' estimated with \code{lfe::felm()} (OLS and IV/2SLS) or with
#' \code{fixest}'s \code{feols()} (OLS and IV), \code{feglm()}, and
#' \code{fepois()}. For GLM fits the variance is the M-estimation sandwich
#' built from the stored score matrix and inverse Hessian. The spatial
#' meat uses a fast CSR/cumulative-score implementation in C++; see
#' \code{\link{vcovSpHAC.felm}} and \code{\link{vcovSpHAC.fixest}} for the
#' per-method argument lists.
#'
#' @param reg A fitted model object.
#' @param ... Method-specific arguments.
#' @export
vcovSpHAC <- function(reg, ...) UseMethod("vcovSpHAC")

#' @export
vcovSpHAC.default <- function(reg, ...) {
  stop("vcovSpHAC: unsupported model class '", paste(class(reg), collapse = "/"),
       "'. Pass a felm or fixest fit.")
}

#' Spatial HAC variance-covariance matrix for felm() models
#'
#' The fit must have been called with \code{felm(..., keepCX = TRUE)} so the
#' centered design matrix is stored on the object. IV/2SLS fits (the
#' \code{felm} multi-part formula with \code{(endog ~ instruments)}) work
#' out of the box: \code{lfe} stores the projected (second-stage) design in
#' \code{cX} and the structural residuals in \code{residuals}, which is
#' exactly the 2SLS sandwich. Weighted fits are supported (the scores carry
#' the weights and the bread uses \eqn{X'WX}).
#'
#' @param reg A fitted object of class "felm", including IV fits.
#' @param unit Optional name of the panel unit variable.
#' @param time Optional name of the time variable.
#' @param lat Name of the latitude variable. If \code{NULL} (default), it is
#'   auto-detected from the data's column names (\code{"lat"},
#'   \code{"latitude"}, case-insensitive); a message reports the pick.
#' @param lon Name of the longitude variable. If \code{NULL} (default), it is
#'   auto-detected (\code{"lon"}, \code{"long"}, \code{"longitude"},
#'   \code{"lng"}, case-insensitive).
#' @param kernel Spatial kernel, either "bartlett" or "uniform".
#' @param dist_fn Distance function, one of "haversine", "spherical", "chord".
#' @param dist_cutoff Spatial cutoff in km.
#' @param lag_cutoff Serial HAC lag cutoff.
#' @param verbose Print progress messages.
#' @param balanced_pnl Whether the panel is balanced and unit locations are time-invariant.
#' @param ncores Number of cores for the C++/RcppParallel spatial and serial routines.
#' @param pixel Score-pre-aggregation cell size, in kilometres. Default 0
#'   (exact-coordinate dedupe only). If `pixel > 0`, points are snapped to a
#'   uniform `pixel`-km grid before the dedupe — a speed/accuracy trade-off
#'   that approximates the distance up to roughly `pixel / 2`.
#' @param neighbor Neighbor-search strategy for the spatial meat: "grid"
#'   (default; 3D cell grid, output-sensitive candidate enumeration) or
#'   "band" (latitude band scan, the pre-0.5.0 behavior). Both are exact and
#'   use identical per-pair accept tests; results agree to floating-point
#'   summation order.
#' @param csr_weight Storage precision for the balanced-path bartlett kernel
#'   weights: "double" (default, exact) or "float" (halves the per-pair
#'   weight memory; introduces at most ~6e-8 relative error per weight).
#'   Ignored for \code{kernel = "uniform"}, which stores no weights.
#' @param method Spatial meat engine. "pairwise" enumerates neighbor pairs
#'   (the default engine; works for any data). "grid" uses the exact
#'   grid-native meat — requires observations on a regular lat/lon lattice
#'   (e.g. raster data); cost is independent of the pair count, so it is
#'   dramatically faster on dense grids with large cutoffs. The uniform
#'   kernel uses sliding-window prefix sums; the bartlett kernel uses
#'   per-ring-pair FFT convolutions. Lattices spanning the full longitude
#'   circle wrap correctly across the dateline. "auto" (default) picks
#'   "grid" when it detects a lattice and a flop-balance estimate says it
#'   wins; both engines are exact, so the choice only affects speed
#'   (results agree to FP summation order, plus ~1e-12 acos conditioning
#'   for the bartlett spherical/chord weights). Pass \code{verbose = TRUE}
#'   to see which engine ran; if you expect gridded data to use the grid
#'   engine and it does not, run once with \code{method = "grid"} — it
#'   errors with the specific reason instead of falling back.
#' @param ssc Small-sample correction. If \code{TRUE} (default), the variance
#'   matrix is scaled by \code{n / (n - K)} where \code{K} counts all
#'   estimated parameters including absorbed fixed-effect levels (taken from
#'   the fit's residual degrees of freedom). This matches \code{fixest}'s
#'   default Conley correction (its cluster adjustment is a no-op for Conley
#'   vcovs). Pass \code{FALSE} for no correction — that reproduces
#'   \code{rbluhm/conley}, fastconley versions before 0.9.0, and
#'   \code{fixest} with \code{ssc(adj = FALSE, cluster.adj = FALSE)}.
#' @param psd_fix The spatial kernels do not guarantee a positive
#'   semi-definite variance matrix. If \code{TRUE} (default), negative
#'   eigenvalues are clamped (to 1e-16, as \code{fixest}'s \code{vcov_fix}
#'   does) and a warning reports when the fix noticeably changed the matrix.
#'   If \code{FALSE}, the matrix is returned as computed, with a warning
#'   when it is not positive semi-definite.
#' @param maxobsmem Ignored by the fast spatial path. Kept for backward compatibility.
#' @param data Optional. The data frame to draw \code{lat}/\code{lon} from.
#'   If \code{NULL} (default), the data is recovered from the fit's call and
#'   aligned via a model-frame re-evaluation. Pass it explicitly if the
#'   original data has gone out of scope — and when no rows were dropped at
#'   fit time (no NAs, no \code{subset}), the coordinates are then taken by
#'   direct column access with no model-frame rebuild at all.
#' @param ... Currently unused.
#' @return A variance-covariance matrix.
#' @examples
#' if (requireNamespace("lfe", quietly = TRUE)) {
#'   ## Cross-section on a regular 0.5-degree raster with holes. method =
#'   ## "grid" forces the exact grid engine; the default method = "auto"
#'   ## picks it automatically when the raster is large enough to win
#'   ## (on a toy example this small, pairwise is just as fast).
#'   set.seed(1)
#'   cells <- expand.grid(lat = seq(40, 49.5, by = 0.5),
#'                        lon = seq(-10, 9.5, by = 0.5))
#'   cells <- cells[sample(nrow(cells), 600), ]   # irregular occupancy
#'   cells$x <- rnorm(nrow(cells))
#'   cells$y <- 0.5 * cells$x + rnorm(nrow(cells))
#'
#'   fit <- lfe::felm(y ~ x, data = cells, keepCX = TRUE)
#'   V <- vcovSpHAC(fit, lat = "lat", lon = "lon",
#'                  kernel = "bartlett", dist_fn = "spherical",
#'                  dist_cutoff = 200, ncores = 2, method = "grid",
#'                  data = cells)
#'   sqrt(diag(V))
#'
#'   ## Panel with spatial + serial HAC (scattered points: pairwise engine)
#'   pnl <- data.frame(unit = rep(1:200, each = 5),
#'                     time = rep(1:5, times = 200),
#'                     lat = rep(runif(200, 40, 50), each = 5),
#'                     lon = rep(runif(200, -10, 10), each = 5))
#'   pnl$x <- rnorm(nrow(pnl))
#'   pnl$y <- 0.5 * pnl$x + rnorm(nrow(pnl))
#'   fit2 <- lfe::felm(y ~ x | unit + time, data = pnl, keepCX = TRUE)
#'   V2 <- vcovSpHAC(fit2, unit = "unit", time = "time",
#'                   lat = "lat", lon = "lon", kernel = "bartlett",
#'                   dist_fn = "haversine", dist_cutoff = 300,
#'                   lag_cutoff = 2, balanced_pnl = TRUE, ncores = 2,
#'                   data = pnl)
#'   sqrt(diag(V2))
#' }
#' @export
vcovSpHAC.felm <- function(reg,
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
                           ...) {

  args <- validate_args(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                        neighbor, csr_weight, method, ssc, psd_fix)

  if (is.null(lat) || is.null(lon)) {
    nm_src <- if (!is.null(data)) names(data) else {
      names(tryCatch(eval(reg$call$data, environment(formula(reg))),
                     error = function(e) NULL))
    }
    cn  <- detect_coord_names(nm_src, lat, lon)
    lat <- cn$lat
    lon <- cn$lon
  }

  noFEs <- length(unit) == 0L
  if (noFEs) {
    unit <- "fe1"
    time <- "fe2"
  }

  Xvars <- names(reg$coefficients)
  if (is.null(Xvars)) Xvars <- rownames(reg$coefficients)
  if (is.null(Xvars) || length(Xvars) == 0L) {
    stop("Could not infer coefficient names from the felm object.")
  }

  # `lfe` stores FEs as factors. `as.integer(factor)` returns level codes,
  # which is what the C++ grouping needs — equality matters, not the original
  # label. Using as.numeric(as.character(...)) (the old `Fac2Num`) silently
  # produced NAs for any factor whose levels weren't numeric strings.
  N <- length(reg$cY)
  if (noFEs) {
    fe1_vec <- seq_len(N)
    fe2_vec <- rep(1L, N)
  } else {
    fe1_vec <- as.integer(reg$fe[[1L]])
    fe2_vec <- as.integer(reg$fe[[2L]])
  }
  # Coordinate recovery. With an explicitly passed `data` whose row count
  # matches the fit (nothing dropped: no NAs, no subset), take the columns
  # directly -- no model-frame re-evaluation, no second copy of the data.
  # Otherwise rebuild the model frame (aligned by rownames) as before.
  if (!is.null(data) && nrow(data) == N && is.null(reg$call$subset)) {
    for (nm in c(lat, lon)) {
      if (!nm %in% names(data)) {
        stop("Column '", nm, "' not found in `data`.")
      }
    }
    coords <- data.frame(data[[lat]], data[[lon]])
    names(coords) <- c(lat, lon)
  } else {
    coords <- expand.model.felm(model = reg, extras = c(lat, lon),
                                na.expand = TRUE, data = data)
  }

  # reg$cY is deliberately not carried along -- only the centered design
  # columns, FEs, coordinates, and residuals are used downstream.
  dt <- data.table::data.table(
    reg$cX,
    fe1 = fe1_vec,
    fe2 = fe2_vec,
    coords
  )

  data.table::setnames(dt, c("fe1", "fe2"),
    c(ifelse(!noFEs, names(reg$fe)[1L], unit),
      ifelse(!noFEs, names(reg$fe)[2L], time))
  )
  # Weighted (WLS) fits: the meat scores are s_i = w_i * e_i * x_i and the
  # bread is (X'WX)^{-1}. lfe stores sqrt(w) in reg$weights, so square it.
  # Folding w into the residual column covers both the spatial and the
  # serial-HAC meat (both build scores as e * X downstream).
  w <- if (!is.null(reg$weights)) as.numeric(reg$weights)^2 else NULL
  res <- as.numeric(reg$residuals)
  if (!is.null(w)) res <- res * w
  dt[, e := res]

  data.table::setnames(dt, c(unit, time, lat, lon),
                            c("unit", "time", "lat", "lon"))

  X <- as.matrix(dt[, Xvars, with = FALSE])
  n <- nrow(dt)
  invXX <- solve(if (is.null(w)) crossprod(X) else crossprod(X, X * w)) * n
  rm(X)

  dof_scale <- ssc_scale(ssc, n, reg$df.residual)

  vcovSpHAC_core(dt = dt, Xvars = Xvars, n = n, invXX = invXX,
                 kernel = args$kernel, dist_fn = args$dist_fn,
                 dist_cutoff = dist_cutoff, lag_cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl, ncores = args$ncores,
                 pixel = pixel, neighbor = args$neighbor,
                 csr_weight = args$csr_weight, method = args$method,
                 dof_scale = dof_scale, psd_fix = psd_fix,
                 verbose = verbose)
}

#' Spatial HAC variance-covariance matrix for fixest models
#'
#' Supports \code{fixest::feols()} (including IV/2SLS) and
#' \code{fixest::feglm()} / \code{fixest::fepois()} fits.
#'
#' For \code{feols}, the fit must have been called with
#' \code{feols(..., demeaned = TRUE)} so that the centered design matrix
#' \code{X_demeaned} is stored on the fit object. Weighted fits are supported
#' (the scores carry the weights and the bread uses \eqn{X'WX}, matching
#' \code{fixest}'s own weighted Conley vcov). IV fits work out of the box:
#' \code{X_demeaned} holds the projected (second-stage) design and
#' \code{residuals} the structural residuals, which is exactly the 2SLS
#' sandwich.
#'
#' For \code{feglm} / \code{fepois}, no estimation flag is needed: the
#' variance is the M-estimation sandwich \eqn{H^{-1} B H^{-1}}, built from
#' the maximum-likelihood score matrix and inverse Hessian that
#' \code{fixest} stores on every (non-\code{lean}) fit. Weights, offsets,
#' and the fixed-effect profiling are already folded into the stored
#' scores. This is the same construction \code{fixest}'s own
#' \code{vcov_conley()} uses for GLMs — but with exact great-circle
#' distances, and with the serial-HAC panel extension available via
#' \code{lag_cutoff} (which \code{fixest} does not offer for Conley vcovs).
#'
#' The returned matrix can be passed to \code{fixest}'s \code{vcov} argument.
#' For the usual \code{fixest} workflow, define a one-argument wrapper such as
#' \code{function(x) vcovSpHAC(x, ...)} and pass that function to
#' \code{summary()}, \code{etable()}, or \code{feols(vcov = )}. The wrapper
#' keeps the coordinate names, cutoffs, panel variables, and optional
#' \code{data = } argument together.
#'
#' @param reg A fitted object of class "fixest": a \code{feols()} fit
#'   (including IV) with \code{demeaned = TRUE}, or a \code{feglm()} /
#'   \code{fepois()} fit (any family; \code{lean = TRUE} fits are rejected
#'   because they carry no score matrix).
#' @param unit Optional name of the panel unit variable. If \code{NULL} the
#'   call is treated as a cross-section (each row is its own unit, all rows
#'   share a single period — no serial HAC).
#' @param time Optional name of the time variable. Ignored when \code{unit} is NULL.
#' @param lat Name of the latitude variable. If \code{NULL} (default),
#'   auto-detected from the data's column names. See \code{\link{vcovSpHAC.felm}}.
#' @param lon Name of the longitude variable. If \code{NULL} (default),
#'   auto-detected. See \code{\link{vcovSpHAC.felm}}.
#' @param kernel Spatial kernel, either "bartlett" or "uniform".
#' @param dist_fn Distance function, one of "haversine", "spherical", "chord".
#' @param dist_cutoff Spatial cutoff in km.
#' @param lag_cutoff Serial HAC lag cutoff.
#' @param verbose Print progress messages.
#' @param balanced_pnl Whether the panel is balanced and unit locations are time-invariant.
#' @param ncores Number of cores for the C++/RcppParallel spatial and serial routines.
#' @param pixel Score-pre-aggregation cell size, in kilometres.
#' @param neighbor Neighbor-search strategy: "grid" (default) or "band".
#'   See \code{\link{vcovSpHAC.felm}}.
#' @param csr_weight Balanced-path bartlett weight storage: "double"
#'   (default) or "float". See \code{\link{vcovSpHAC.felm}}.
#' @param method Spatial meat engine: "auto" (default), "pairwise", or
#'   "grid". See \code{\link{vcovSpHAC.felm}}.
#' @param ssc Small-sample correction (\code{n / (n - K)} when \code{TRUE},
#'   the default). See \code{\link{vcovSpHAC.felm}}.
#' @param psd_fix Clamp negative eigenvalues when \code{TRUE} (the
#'   default). See \code{\link{vcovSpHAC.felm}}.
#' @param data Optional. The model frame to draw \code{lat}/\code{lon}/
#'   \code{unit}/\code{time} from. If \code{NULL} (default), the data is
#'   recovered from the fit's call. Pass it explicitly if the original data
#'   has gone out of scope, or if you want to override.
#' @param ... Currently unused.
#' @return A variance-covariance matrix.
#' @examples
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   ## feols must be fit with demeaned = TRUE (the keepCX analogue).
#'   set.seed(1)
#'   cells <- expand.grid(lat = seq(40, 49.5, by = 0.5),
#'                        lon = seq(-10, 9.5, by = 0.5))
#'   cells$x <- rnorm(nrow(cells))
#'   cells$y <- 0.5 * cells$x + rnorm(nrow(cells))
#'
#'   fit <- fixest::feols(y ~ x, data = cells, demeaned = TRUE)
#'   vcov_fc <- function(x) {
#'     vcovSpHAC(x, lat = "lat", lon = "lon",
#'               kernel = "uniform", dist_fn = "spherical",
#'               dist_cutoff = 200, ncores = 2, data = cells)
#'   }
#'   V <- vcov_fc(fit)
#'   sqrt(diag(V))
#'
#'   ## The same wrapper can be used directly in fixest's vcov argument.
#'   fit_sum <- summary(fit, vcov = vcov_fc)
#'   sqrt(diag(fit_sum$cov.scaled))
#'
#'   ## Poisson (fepois / feglm): no demeaned = TRUE needed — the stored
#'   ## ML scores and inverse Hessian are used directly.
#'   cells$cnt <- rpois(nrow(cells), exp(0.4 * cells$x))
#'   fit_pois <- fixest::fepois(cnt ~ x, data = cells)
#'   V_pois <- vcovSpHAC(fit_pois, lat = "lat", lon = "lon",
#'                       kernel = "uniform", dist_fn = "spherical",
#'                       dist_cutoff = 200, ncores = 2, data = cells)
#'   sqrt(diag(V_pois))
#' }
#' @export
vcovSpHAC.fixest <- function(reg,
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
                             ...) {

  args <- validate_args(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                        neighbor, csr_weight, method, ssc, psd_fix)

  method_type <- reg$method_type
  if (is.null(method_type)) method_type <- "feols"
  if (!method_type %in% c("feols", "feglm")) {
    stop("vcovSpHAC.fixest supports feols() and feglm()/fepois() fits, but ",
         "this fit has method_type '", method_type, "'. femlm()/feNmlm() ",
         "fits are not supported.")
  }
  is_glm <- method_type == "feglm"

  # Use direct field access rather than coef(reg) / residuals(reg) — the
  # generic accessors in `fixest` can be 100x+ slower than field access
  # under do.call() evaluation contexts (they recompute or look up the call
  # environment), which makes vcovSpHAC inadvertently slow when wrapped.
  Xvars <- names(reg$coefficients)
  if (is_glm) {
    # M-estimation sandwich H^{-1} B H^{-1}: fixest stores the ML score
    # rows (FE profiling, weights, and offsets already folded in) and
    # cov.iid = H^{-1}. The score matrix rides through the core as the
    # "X" columns with e = 1, so every engine (pairwise, grid/FFT, serial
    # HAC, pixel aggregation) sees the same score rows it would for OLS.
    if (!is.null(reg$isBounded) && any(reg$isBounded)) {
      stop("vcovSpHAC.fixest: fits with parameters estimated at a bound ",
           "are not supported.")
    }
    cX <- reg$scores
    if (is.null(cX)) {
      stop("vcovSpHAC.fixest: the fit carries no score matrix (was it fit ",
           "with lean = TRUE?). Refit with lean = FALSE.")
    }
    cX <- as.matrix(cX)
    if (is.null(Xvars)) Xvars <- colnames(cX)
    if (ncol(cX) != length(Xvars)) {
      stop("Internal: score columns (", ncol(cX),
           ") do not match the coefficient count (", length(Xvars), ").")
    }
    colnames(cX) <- Xvars
    if (is.null(reg$cov.iid) || anyNA(reg$cov.iid)) {
      stop("vcovSpHAC.fixest: the fit's iid vcov (inverse Hessian) is ",
           "missing or contains NAs (collinearity?); cannot build the ",
           "sandwich bread.")
    }
    n <- nrow(cX)
    e <- rep(1, n)
  } else {
    if (is.null(reg$X_demeaned)) {
      stop("vcovSpHAC.fixest requires the fit to have been called with ",
           "feols(..., demeaned = TRUE) so the centered design matrix is available.")
    }
    if (is.null(Xvars)) Xvars <- colnames(reg$X_demeaned)
    cX <- reg$X_demeaned[, Xvars, drop = FALSE]
    e  <- as.numeric(reg$residuals)
    n  <- length(e)
    if (nrow(cX) != n) {
      stop("Internal: X_demeaned rows (", nrow(cX),
           ") do not match residual length (", n, ").")
    }
  }

  if (is.null(data)) {
    # `feols`'s formula environment is sometimes a captured closure that does
    # not see the caller's locals, in which case eval() walks up the search
    # path and may resolve a one-letter name like `d` to a base function.
    # Try the caller's frame first, then the formula environment.
    envs <- list(parent.frame(), environment(stats::formula(reg)))
    for (env in envs) {
      cand <- tryCatch(eval(reg$call$data, env), error = function(e) NULL)
      if (is.data.frame(cand) || data.table::is.data.table(cand)) {
        data <- cand
        break
      }
    }
    if (is.null(data)) {
      stop("Could not locate the model data on the fit's call. ",
           "Pass `data = ` to vcovSpHAC explicitly.")
    }
  }
  if (is.null(lat) || is.null(lon)) {
    cn  <- detect_coord_names(names(data), lat, lon)
    lat <- cn$lat
    lon <- cn$lon
  }
  # We index `data` by column name with `[[`, which works identically on
  # data.frame and data.table — avoid converting (a data.table-to-data.frame
  # copy at n = 8e3, k = 8 is hundreds of KB of pointless deep copy).
  for (nm in c(lat, lon)) {
    if (!nm %in% names(data)) {
      stop("Column '", nm, "' not found in the model data.")
    }
  }
  if (!is.null(unit) && !unit %in% names(data)) {
    stop("Column '", unit, "' not found in the model data.")
  }
  if (!is.null(time) && !time %in% names(data)) {
    stop("Column '", time, "' not found in the model data.")
  }

  # `obs_selection` has two stacked filters:
  #   - $subset:     positive 1-based indices into the original data, set when
  #                  the user passed `subset =` to feols.
  #   - $obsRemoved: NEGATIVE indices into the SUBSET-APPLIED view (rows dropped
  #                  for NAs, singletons, perfect collinearity).
  # The two must be applied in order: subset first, then obsRemoved.
  subset_idx  <- reg$obs_selection$subset
  removed_idx <- reg$obs_selection$obsRemoved
  noFEs <- is.null(unit)
  pick <- function(col) {
    v <- data[[col]]
    if (length(subset_idx))  v <- v[subset_idx]
    if (length(removed_idx)) v <- v[removed_idx]
    v
  }
  lat_v  <- pick(lat)
  lon_v  <- pick(lon)
  unit_v <- if (noFEs) seq_len(n) else pick(unit)
  time_v <- if (noFEs || is.null(time)) rep(1L, n) else pick(time)

  if (length(lat_v) != n) {
    stop("Row count after applying obs_selection (", length(lat_v),
         ") does not match the fitted model (", n, "). ",
         "Pass `data = ` explicitly or refit on a frame with no extra rows.")
  }

  # Weighted (WLS) feols fits: scores are s_i = w_i * e_i * x_i, bread is
  # (X'WX)^{-1}. fixest stores the weights on the original scale and
  # X_demeaned / residuals on the raw (unweighted) scale. On the GLM path
  # the stored scores already carry the weights, so nothing is folded.
  w <- if (!is_glm && !is.null(reg$weights)) as.numeric(reg$weights) else NULL
  if (!is.null(w)) e <- e * w

  dt <- data.table::data.table(
    cX,
    unit = unit_v,
    time = time_v,
    lat  = lat_v,
    lon  = lon_v,
    e    = e
  )

  # The core computes invXX %*% (XeeX / n) %*% invXX / n, i.e. it expects
  # n * bread. For GLMs the bread is the inverse Hessian fixest stores as
  # cov.iid; for (possibly weighted) feols it is (X'WX)^{-1}.
  invXX <- if (is_glm) {
    unname(reg$cov.iid) * n
  } else {
    solve(if (is.null(w)) crossprod(cX) else crossprod(cX, cX * w)) * n
  }

  dof_scale <- ssc_scale(ssc, n, reg$nobs - reg$nparams)

  vcovSpHAC_core(dt = dt, Xvars = Xvars, n = n, invXX = invXX,
                 kernel = args$kernel, dist_fn = args$dist_fn,
                 dist_cutoff = dist_cutoff, lag_cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl, ncores = args$ncores,
                 pixel = pixel, neighbor = args$neighbor,
                 csr_weight = args$csr_weight, method = args$method,
                 dof_scale = dof_scale, psd_fix = psd_fix,
                 verbose = verbose)
}

# Shared post-extraction core. `dt` must carry Xvars, unit, time, lat, lon, e
# (with any regression weights already folded into e on entry).
vcovSpHAC_core <- function(dt, Xvars, n, invXX,
                           kernel, dist_fn, dist_cutoff, lag_cutoff,
                           balanced_pnl, ncores, pixel, neighbor, csr_weight,
                           method, dof_scale = 1, psd_fix = FALSE, verbose) {

  # The FastSpatialMeat / FastSerialHacPanel C++ entry points require numeric
  # vectors for time and unit (arma::vec). Equality is the only thing they use
  # — the integer codes preserve the user's grouping while letting character
  # and factor inputs through. Do this BEFORE setorderv so the integer
  # codes match the user-supplied sort order (factor level order for factors;
  # alphabetical for raw strings).
  dt[, unit := to_group_id(dt[["unit"]])]
  dt[, time := to_group_id(dt[["time"]])]

  if (balanced_pnl) {
    data.table::setorderv(dt, c("time", "unit"))
    if (data.table::uniqueN(dt[["time"]]) > 1L) {
      # The C++ balanced path builds the CSR neighbor graph from the first
      # period and reuses it for every subsequent period. That is only valid
      # when every period contains exactly the same units (no repeats, no
      # missing units, no extras). Check that:
      #   (a) every period has the same row count,
      #   (b) no unit appears twice within a period,
      #   (c) every period sees the same unit set,
      # before checking time-invariant coordinates.
      period_n <- dt[, .N, by = "time"][["N"]]
      if (length(unique(period_n)) != 1L) {
        stop("balanced_pnl = TRUE requires each period to have the same ",
             "number of observations, but period sizes are: ",
             paste(unique(period_n), collapse = ", "), ".")
      }
      # dt is sorted by (time, unit), so each period's unit column is sorted:
      # period 1 must be duplicate-free, and every later period must repeat
      # period 1's unit sequence exactly. Together with the equal-size check
      # above this implies the same unit set in every period with no repeats.
      # O(n) with no per-period grouping or string materialization.
      n_per <- period_n[1L]
      u <- dt[["unit"]]
      u1 <- u[seq_len(n_per)]
      if (anyDuplicated(u1) > 0L) {
        stop("balanced_pnl = TRUE requires each unit to appear at most ",
             "once per period, but some (unit, time) combinations are ",
             "duplicated. Use balanced_pnl = FALSE.")
      }
      if (!isTRUE(all(u == rep(u1, times = length(period_n))))) {
        stop("balanced_pnl = TRUE requires every period to contain the ",
             "same set of units, but unit membership varies across periods ",
             "(or a unit is duplicated within a period). ",
             "Use balanced_pnl = FALSE.")
      }
      coord_var <- dt[, list(n_lat = data.table::uniqueN(lat),
                             n_lon = data.table::uniqueN(lon)),
                      by = unit]
      if (any(coord_var$n_lat > 1L) || any(coord_var$n_lon > 1L)) {
        stop("balanced_pnl = TRUE requires time-invariant coordinates per unit, ",
             "but some units have varying lat or lon across periods.")
      }
    }
  } else {
    data.table::setorderv(dt, "time")
  }

  agg <- aggregate_scores(dt, Xvars, pixel = pixel,
                          balanced_pnl = balanced_pnl, verbose = verbose)

  gi <- choose_grid_method(method, kernel, agg, dist_cutoff)
  if (verbose) {
    message("Starting fast spatial HAC meat in C++ (",
            if (is.null(gi)) "pairwise" else "grid", " engine)")
  }
  run_pairwise <- function() {
    FastSpatialMeat(
      lat = agg$lat,
      lon = agg$lon,
      time = agg$time,
      scores = agg$scores,
      cutoff = dist_cutoff,
      kernel = kernel,
      dist_fn = dist_fn,
      balanced_pnl = balanced_pnl,
      ncores = ncores,
      neighbor = neighbor,
      csr_weight = csr_weight
    )
  }
  XeeX <- if (!is.null(gi)) {
    run_grid <- function() {
      FastGridMeat(
        ring = gi$ring, col = gi$col, time = agg$time, scores = agg$scores,
        lat0 = gi$lat0, dlat = gi$dlat, dlon = gi$dlon,
        n_ring = gi$n_ring, n_col = gi$n_col, n_col_full = gi$n_col_full,
        cutoff = dist_cutoff, dist_fn = dist_fn, kernel = kernel,
        ncores = ncores
      )
    }
    if (method == "auto") {
      # The C++ wrap-feasibility check refuses lattices whose accept window
      # crosses the dateline gap when dlon does not tile 360 evenly. Under
      # "auto" that is a speed choice, not an error -- fall back.
      tryCatch(run_grid(), error = function(e) {
        if (grepl("dateline", conditionMessage(e), fixed = TRUE)) {
          run_pairwise()
        } else {
          stop(e)
        }
      })
    } else {
      run_grid()
    }
  } else {
    run_pairwise()
  }

  if (lag_cutoff > 0 && length(unique(dt[["time"]])) > 1L) {
    data.table::setorderv(dt, c("unit", "time"))

    if (verbose) message("Starting serial HAC meat")
    scores_serial <- as.matrix(dt[, Xvars, with = FALSE]) * dt[["e"]]
    XeeX_serial <- FastSerialHacPanel(
      unit = dt[["unit"]],
      time = dt[["time"]],
      cutoff = lag_cutoff,
      scores = scores_serial,
      ncores = ncores
    )

    XeeX <- XeeX + XeeX_serial
  }

  V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
  V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2
  if (dof_scale != 1) V_spatial_HAC <- V_spatial_HAC * dof_scale
  rownames(V_spatial_HAC) <- colnames(V_spatial_HAC) <- Xvars

  # The spatial kernels do not guarantee a PSD meat. Mirror fixest's
  # vcov_fix: clamp eigenvalues to 1e-16 and warn only when the fix
  # noticeably changed the matrix.
  ev <- eigen(V_spatial_HAC, symmetric = TRUE)
  if (any(ev$values <= 0)) {
    V_fixed <- tcrossprod(ev$vectors %*% diag(pmax(ev$values, 1e-16),
                                              length(ev$values)),
                          ev$vectors)
    noticeable <- max(abs(V_spatial_HAC - V_fixed)) > 1e-8
    if (psd_fix) {
      dimnames(V_fixed) <- dimnames(V_spatial_HAC)
      V_spatial_HAC <- V_fixed
      if (noticeable) {
        warning("The vcov matrix was not positive semi-definite and was ",
                "fixed by clamping negative eigenvalues.", call. = FALSE)
      }
    } else if (noticeable) {
      warning("The vcov matrix is not positive semi-definite. Pass ",
              "psd_fix = TRUE to clamp negative eigenvalues.", call. = FALSE)
    }
  }
  V_spatial_HAC
}

# Argument validation shared by both methods.
validate_args <- function(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                          neighbor = c("grid", "band"),
                          csr_weight = c("double", "float"),
                          method = c("auto", "pairwise", "grid"),
                          ssc = FALSE, psd_fix = FALSE) {
  if (!(isTRUE(ssc) || isFALSE(ssc))) {
    stop("ssc must be TRUE or FALSE.")
  }
  if (!(isTRUE(psd_fix) || isFALSE(psd_fix))) {
    stop("psd_fix must be TRUE or FALSE.")
  }
  kernel     <- match.arg(kernel,  c("bartlett", "uniform"))
  dist_fn    <- match.arg(dist_fn, c("haversine", "spherical", "chord"))
  neighbor   <- match.arg(neighbor, c("grid", "band"))
  csr_weight <- match.arg(csr_weight, c("double", "float"))
  method     <- match.arg(method, c("auto", "pairwise", "grid"))
  if (is.null(dist_cutoff) || length(dist_cutoff) != 1L ||
      !is.finite(dist_cutoff) || dist_cutoff <= 0) {
    stop("dist_cutoff must be a single positive finite number.")
  }
  if (length(lag_cutoff) != 1L || !is.finite(lag_cutoff) || lag_cutoff < 0) {
    stop("lag_cutoff must be a single non-negative finite number.")
  }
  if (length(pixel) != 1L || !is.finite(pixel) || pixel < 0) {
    stop("pixel must be a single non-negative finite number.")
  }
  if (is.na(ncores)) {
    ncores <- max(1L, parallel::detectCores(logical = TRUE))
  }
  ncores <- as.integer(max(1L, ncores))
  list(kernel = kernel, dist_fn = dist_fn, ncores = ncores, neighbor = neighbor,
       csr_weight = csr_weight, method = method)
}

# Small-sample scale factor n / max(n - K, 1), with n - K taken from the
# fit's residual degrees of freedom (the max(., 1) floor is fixest's
# convention for saturated fits). Returns 1 when ssc is FALSE or the fit
# does not expose its df (with a warning, so a silent no-op cannot be
# mistaken for a correction).
ssc_scale <- function(ssc, n, df_resid) {
  if (!isTRUE(ssc)) return(1.0)
  if (!length(df_resid) || !is.finite(df_resid)) {
    warning("ssc = TRUE, but the fit does not expose its residual degrees ",
            "of freedom; no small-sample correction applied.", call. = FALSE)
    return(1.0)
  }
  n / max(as.numeric(df_resid), 1)
}

# Guess lat/lon column names from the data when the user did not pass them.
# Case-insensitive exact matches only — a substring match could grab
# "dilation" or "elongation". Errors when no (or an ambiguous) match exists.
detect_coord_names <- function(nms, lat, lon) {
  if (is.null(nms) || !length(nms)) {
    stop("lat/lon were not supplied and the model data could not be located ",
         "to auto-detect them. Pass lat = and lon = (or data = ) explicitly.")
  }
  low <- tolower(nms)
  pick1 <- function(cands, what) {
    for (cand in cands) {
      i <- which(low == cand)
      if (length(i) == 1L) return(nms[i])
      if (length(i) > 1L) {
        stop("Auto-detection of the ", what, " column is ambiguous (several ",
             "columns named '", cand, "' up to case). Pass it explicitly.")
      }
    }
    stop("Could not auto-detect the ", what, " column (looked for: ",
         paste(cands, collapse = ", "), "). Pass it explicitly.")
  }
  if (is.null(lat)) lat <- pick1(c("lat", "latitude"), "latitude")
  if (is.null(lon)) lon <- pick1(c("lon", "long", "longitude", "lng"), "longitude")
  message("vcovSpHAC: using lat = \"", lat, "\", lon = \"", lon,
          "\" (auto-detected).")
  list(lat = lat, lon = lon)
}

# Detect whether the coordinates lie on a regular lat/lon lattice (gaps —
# missing cells — are fine; rings/columns just need a common step). Returns
# NULL or list(lat0, dlat, lon0, dlon, ring, col, n_ring, n_col) with
# 0-based integer cell indices aligned to the input rows.
detect_lonlat_grid <- function(lat, lon, tol = 1e-6) {
  ul <- sort(unique(lat))
  uo <- sort(unique(lon))
  if (length(ul) < 2L || length(uo) < 2L) return(NULL)
  dl <- diff(ul); dol <- diff(uo)
  step_l <- min(dl); step_o <- min(dol)
  if (step_l <= 0 || step_o <= 0) return(NULL)
  if (any(abs(dl / step_l - round(dl / step_l)) > tol)) return(NULL)
  if (any(abs(dol / step_o - round(dol / step_o)) > tol)) return(NULL)
  ring_d <- round((lat - ul[1L]) / step_l)
  col_d  <- round((lon - uo[1L]) / step_o)
  if (!all(is.finite(ring_d)) || !all(is.finite(col_d))) return(NULL)
  if (max(ring_d) > .Machine$integer.max - 1 ||
      max(col_d) > .Machine$integer.max - 1) return(NULL)
  if (max(abs(lat - (ul[1L] + ring_d * step_l))) > tol * step_l) return(NULL)
  if (max(abs(lon - (uo[1L] + col_d * step_o))) > tol * step_o) return(NULL)
  ring <- as.integer(ring_d)
  col  <- as.integer(col_d)
  n_col <- max(col) + 1L
  # Full-circle column count when the lon step tiles 360 degrees evenly
  # (0 otherwise). FastGridMeat uses it to wrap windows across the
  # dateline when the lattice spans the whole circle.
  cf <- 360 / step_o
  n_col_full <- if (abs(cf - round(cf)) < 1e-6 && round(cf) >= n_col) {
    as.integer(round(cf))
  } else 0L
  list(lat0 = ul[1L], dlat = step_l, lon0 = uo[1L], dlon = step_o,
       ring = ring, col = col,
       n_ring = max(ring) + 1L, n_col = n_col, n_col_full = n_col_full)
}

# Decide whether the grid-native meat applies. Returns the lattice info
# (to be passed to FastGridMeat) or NULL for the pairwise engine. Both
# engines are exact, so a "wrong" auto choice only costs speed.
choose_grid_method <- function(method, kernel, agg, dist_cutoff) {
  if (method == "pairwise") return(NULL)
  gi <- detect_lonlat_grid(agg$lat, agg$lon)
  if (is.null(gi)) {
    if (method == "grid") {
      stop("method = 'grid' requires observations on a regular lat/lon ",
           "lattice (e.g. raster cell centers); no lattice detected.")
    }
    return(NULL)
  }
  cells <- as.numeric(gi$n_ring) * as.numeric(gi$n_col)
  n <- length(agg$lat)
  if (cells > 100 * n || cells > 4e9) {
    # Degenerate lattice (e.g. two nearly-coincident coordinate values
    # implying a huge sparse grid): the dense prefix tensor would explode.
    if (method == "grid") {
      stop("method = 'grid': the implied lattice has ", format(cells),
           " cells for ", n, " observations; too sparse/degenerate. ",
           "Use method = 'pairwise'.")
    }
    return(NULL)
  }
  if (method == "grid") return(gi)
  # auto: flop-balance rule. Pairwise work ~ expected within-cutoff pairs;
  # grid work per ring pair is the window-sum width (uniform boxcar) or the
  # FFT cost in per-k units (bartlett). Prefer grid only when pairs clearly
  # dominate.
  ang <- dist_cutoff / 6371
  rho_r <- max(1, ang / (gi$dlat * pi / 180))
  cosbar <- max(0.05, mean(cos(agg$lat * pi / 180)))
  dbar <- max(1, ang / (gi$dlon * pi / 180) / cosbar)
  nb <- as.numeric(table(agg$time))
  pairs_est <- sum(nb^2) * pi * rho_r * dbar / cells / 2
  per_pair <- if (kernel == "uniform") {
    as.numeric(gi$n_col)
  } else {
    k <- ncol(agg$scores)
    npad_est <- 2^ceiling(log2(as.numeric(gi$n_col) + min(dbar, gi$n_col)))
    npad_est * (5 * log2(npad_est) + 2 * k) / k
  }
  grid_work <- length(nb) * as.numeric(gi$n_ring) * (2 * rho_r + 1) * per_pair
  if (pairs_est > 2 * grid_work) gi else NULL
}

# Map any unit/time vector to integer group codes. The FastSpatialMeat /
# FastSerialHacPanel entry points require numeric (arma::vec) input but only
# use equality, so the actual values are immaterial; we just need same labels
# to share a code. `as.integer(factor)` returns the factor's level code, which
# is what we want both for input factors and for raw character vectors.
to_group_id <- function(x) {
  if (is.integer(x) || is.numeric(x)) return(x)
  if (is.factor(x)) return(as.integer(x))
  as.integer(factor(x))
}

# Build per-time-block aggregated scores. Returns a list with element-wise
# vectors plus the score matrix ready to hand to FastSpatialMeat: lat, lon,
# time, scores.
#
# At pixel = 0 we collapse rows whose (lat, lon) match exactly. At pixel > 0
# we first snap (lat, lon) to a uniform grid whose latitude step is
# pixel / 111 km and whose longitude step is pixel / (111 * cos(lat_rep)) km,
# so cells stay roughly square at all latitudes. The representative coordinate
# for a cell is the cell centre.
aggregate_scores <- function(dt, Xvars, pixel, balanced_pnl, verbose) {
  n <- nrow(dt)
  scores <- as.matrix(dt[, Xvars, with = FALSE]) * dt[["e"]]

  if (pixel > 0) {
    lat_step <- pixel / 111.0
    lat_cell <- round(dt[["lat"]] / lat_step)
    lat_rep  <- lat_cell * lat_step
    cos_lat  <- cos(lat_rep * pi / 180)
    cos_lat[abs(cos_lat) < 1e-6] <- 1e-6  # avoid blow-up near poles
    lon_step <- pixel / (111.0 * cos_lat)
    lon_cell <- round(dt[["lon"]] / lon_step)
    lon_rep  <- lon_cell * lon_step
    key_lat <- lat_rep
    key_lon <- lon_rep
  } else {
    key_lat <- dt[["lat"]]
    key_lon <- dt[["lon"]]
  }

  agg_dt <- data.table::data.table(
    time = dt[["time"]],
    lat  = key_lat,
    lon  = key_lon,
    scores
  )
  data.table::setnames(agg_dt, colnames(scores), Xvars)

  agg <- agg_dt[, lapply(.SD, sum), by = c("time", "lat", "lon"), .SDcols = Xvars]
  # Preserve the spatial path's expectations: balanced wants (time, <stable
  # within-period order>) — we use (time, lat, lon); general just wants
  # contiguous time blocks.
  data.table::setorderv(agg, c("time", "lat", "lon"))

  n_agg <- nrow(agg)
  if (verbose && n_agg < n) {
    message(sprintf("Score pre-aggregation: %d rows -> %d cases (pixel = %g km)",
                    n, n_agg, pixel))
  }

  # Balanced path additionally requires identical row count per period and a
  # stable within-period ordering. If the dedupe broke that invariant (some
  # periods lost points the others kept), fall back to the general path by
  # warning and returning unaggregated rows.
  if (balanced_pnl) {
    counts <- agg[, .N, by = "time"][["N"]]
    if (length(unique(counts)) != 1L) {
      warning("Pre-aggregation produced different case counts per period; ",
              "skipping aggregation. Pass balanced_pnl = FALSE to silence.")
      return(list(
        lat = dt[["lat"]], lon = dt[["lon"]], time = dt[["time"]],
        scores = scores
      ))
    }
  }

  list(
    lat = agg[["lat"]], lon = agg[["lon"]], time = agg[["time"]],
    scores = as.matrix(agg[, Xvars, with = FALSE])
  )
}

expand.model.felm <- function(model, extras, envir = environment(formula(model)),
                              na.expand = FALSE, data = NULL) {

  topaste <- c(names(model$fe), names(model$clustervar), extras)
  fescluext <- parse(text = paste("~", paste(topaste, collapse = "+")))[[1L]]

  if (is.null(data)) data <- eval(model$call$data, envir)
  ff <- foo ~ bar + baz

  ff[[2L]] <- parse(text = paste("~", model$lhs))[[1L]][[2L]]
  ff[[3L]][[2L]] <- formula(model)[[2L]]
  ff[[3L]][[3L]] <- fescluext[[2L]]

  if (!na.expand) {
    naa <- model$call$na.action
    subset <- model$call$subset
    rval <- eval(call("model.frame", ff, data = data, subset = subset,
                      na.action = naa), envir)[, extras]
  } else {
    subset <- model$call$subset
    rval <- eval(call("model.frame", ff, data = data, subset = subset,
                      na.action = I), envir)
    oldmf <- model.frame(model)
    keep <- match(rownames(oldmf), rownames(rval))
    rval <- rval[keep, extras]
    class(rval) <- "data.frame"
  }
  rval
}
