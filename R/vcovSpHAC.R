#' Spatial HAC variance-covariance matrix
#'
#' Computes Conley (1999) spatial HAC variance-covariance matrices for models
#' estimated with \code{lfe::felm()} or \code{fixest::feols()}. The spatial
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
#' @param reg A fitted object of class "felm".
#' @param unit Optional name of the panel unit variable.
#' @param time Optional name of the time variable.
#' @param lat Name of the latitude variable.
#' @param lon Name of the longitude variable.
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
#' @param maxobsmem Ignored by the fast spatial path. Kept for backward compatibility.
#' @param data Optional. The data frame to draw \code{lat}/\code{lon} from.
#'   If \code{NULL} (default), the data is recovered from the fit's call and
#'   aligned via a model-frame re-evaluation. Pass it explicitly if the
#'   original data has gone out of scope — and when no rows were dropped at
#'   fit time (no NAs, no \code{subset}), the coordinates are then taken by
#'   direct column access with no model-frame rebuild at all.
#' @param ... Currently unused.
#' @return A variance-covariance matrix.
#' @export
vcovSpHAC.felm <- function(reg,
                           unit = NULL,
                           time = NULL,
                           lat,
                           lon,
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
                           maxobsmem = 50000L,
                           data = NULL,
                           ...) {

  args <- validate_args(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                        neighbor, csr_weight)

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
  dt[, e := as.numeric(reg$residuals)]

  data.table::setnames(dt, c(unit, time, lat, lon),
                            c("unit", "time", "lat", "lon"))

  X <- as.matrix(dt[, Xvars, with = FALSE])
  n <- nrow(dt)
  invXX <- solve(crossprod(X)) * n
  rm(X)

  vcovSpHAC_core(dt = dt, Xvars = Xvars, n = n, invXX = invXX,
                 kernel = args$kernel, dist_fn = args$dist_fn,
                 dist_cutoff = dist_cutoff, lag_cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl, ncores = args$ncores,
                 pixel = pixel, neighbor = args$neighbor,
                 csr_weight = args$csr_weight, verbose = verbose)
}

#' Spatial HAC variance-covariance matrix for fixest::feols models
#'
#' Requires the fit to have been called with \code{feols(..., demeaned = TRUE)}
#' so that the centered design matrix \code{X_demeaned} is stored on the fit
#' object. Weighted fits are not supported.
#'
#' @param reg A fitted object of class "fixest".
#' @param unit Optional name of the panel unit variable. If \code{NULL} the
#'   call is treated as a cross-section (each row is its own unit, all rows
#'   share a single period — no serial HAC).
#' @param time Optional name of the time variable. Ignored when \code{unit} is NULL.
#' @param lat Name of the latitude variable.
#' @param lon Name of the longitude variable.
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
#' @param data Optional. The model frame to draw \code{lat}/\code{lon}/
#'   \code{unit}/\code{time} from. If \code{NULL} (default), the data is
#'   recovered from the fit's call. Pass it explicitly if the original data
#'   has gone out of scope, or if you want to override.
#' @param ... Currently unused.
#' @return A variance-covariance matrix.
#' @export
vcovSpHAC.fixest <- function(reg,
                             unit = NULL,
                             time = NULL,
                             lat,
                             lon,
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
                             data = NULL,
                             ...) {

  args <- validate_args(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                        neighbor, csr_weight)

  if (is.null(reg$X_demeaned)) {
    stop("vcovSpHAC.fixest requires the fit to have been called with ",
         "feols(..., demeaned = TRUE) so the centered design matrix is available.")
  }
  if (!is.null(reg$weights)) {
    stop("vcovSpHAC.fixest does not currently support weighted fits.")
  }

  # Use direct field access rather than coef(reg) / residuals(reg) — the
  # generic accessors in `fixest` can be 100x+ slower than field access
  # under do.call() evaluation contexts (they recompute or look up the call
  # environment), which makes vcovSpHAC inadvertently slow when wrapped.
  Xvars <- names(reg$coefficients)
  if (is.null(Xvars)) Xvars <- colnames(reg$X_demeaned)
  cX    <- reg$X_demeaned[, Xvars, drop = FALSE]
  e     <- as.numeric(reg$residuals)
  n     <- length(e)
  if (nrow(cX) != n) {
    stop("Internal: X_demeaned rows (", nrow(cX),
         ") do not match residual length (", n, ").")
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

  dt <- data.table::data.table(
    cX,
    unit = unit_v,
    time = time_v,
    lat  = lat_v,
    lon  = lon_v,
    e    = e
  )

  invXX <- solve(crossprod(cX)) * n

  vcovSpHAC_core(dt = dt, Xvars = Xvars, n = n, invXX = invXX,
                 kernel = args$kernel, dist_fn = args$dist_fn,
                 dist_cutoff = dist_cutoff, lag_cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl, ncores = args$ncores,
                 pixel = pixel, neighbor = args$neighbor,
                 csr_weight = args$csr_weight, verbose = verbose)
}

# Shared post-extraction core. `dt` must carry Xvars, unit, time, lat, lon, e.
vcovSpHAC_core <- function(dt, Xvars, n, invXX,
                           kernel, dist_fn, dist_cutoff, lag_cutoff,
                           balanced_pnl, ncores, pixel, neighbor, csr_weight,
                           verbose) {

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

  if (verbose) message("Starting fast spatial HAC meat in C++")
  XeeX <- FastSpatialMeat(
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
  rownames(V_spatial_HAC) <- colnames(V_spatial_HAC) <- Xvars
  V_spatial_HAC
}

# Argument validation shared by both methods.
validate_args <- function(kernel, dist_fn, dist_cutoff, lag_cutoff, pixel, ncores,
                          neighbor = c("grid", "band"),
                          csr_weight = c("double", "float")) {
  kernel     <- match.arg(kernel,  c("bartlett", "uniform"))
  dist_fn    <- match.arg(dist_fn, c("haversine", "spherical", "chord"))
  neighbor   <- match.arg(neighbor, c("grid", "band"))
  csr_weight <- match.arg(csr_weight, c("double", "float"))
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
       csr_weight = csr_weight)
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
