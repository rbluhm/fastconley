#' Spatial HAC variance-covariance matrix for felm() models
#'
#' Computes Conley-style spatial HAC variance-covariance matrices for models
#' estimated with lfe::felm(). This development version uses a fast spatial
#' CSR/cumulative-score implementation instead of dense distance matrices.
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
#'   (exact-coordinate dedupe only — observations that share the same
#'   `(lat, lon)` within a time period are collapsed into a single
#'   aggregated score before the pair loop runs). If `pixel > 0`, points are
#'   snapped to a uniform `pixel`-km grid before the dedupe, so nearby points
#'   are collapsed too — this is a speed/accuracy trade-off that approximates
#'   the distance up to roughly `pixel / 2`. Pre-aggregation only changes
#'   speed; with `pixel = 0` the answer is identical to the un-aggregated
#'   computation up to floating-point order. With `pixel > 0` the answer is
#'   approximate. Inspired by `fixest`'s `pixel` argument.
#' @param maxobsmem Ignored by the fast spatial path. Kept for backward compatibility.
#' @return A variance-covariance matrix.
#' @export
vcovSpHAC <- function(reg,
                      unit = NULL,
                      time = NULL,
                      lat,
                      lon,
                      kernel = c("bartlett", "uniform"),
                      dist_fn = c("haversine", "spherical", "chord"),
                      dist_cutoff,
                      lag_cutoff = 0,
                      verbose = FALSE,
                      balanced_pnl = FALSE,
                      ncores = NA,
                      pixel = 0,
                      maxobsmem = 50000L) {

  kernel <- match.arg(kernel)
  dist_fn <- match.arg(dist_fn)

  if (missing(dist_cutoff) || length(dist_cutoff) != 1L || !is.finite(dist_cutoff) || dist_cutoff <= 0) {
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

  Fac2Num <- function(x) {
    as.numeric(as.character(x))
  }

  noFEs <- length(unit) == 0L

  if (inherits(reg, "felm")) {
    if (noFEs) {
      unit <- "fe1"
      time <- "fe2"
    }

    Xvars <- names(reg$coefficients)
    if (is.null(Xvars)) Xvars <- rownames(reg$coefficients)
    if (is.null(Xvars) || length(Xvars) == 0L) {
      stop("Could not infer coefficient names from the felm object.")
    }

    dt <- data.table::data.table(
      reg$cY,
      reg$cX,
      fe1 = ifelse(rep(!noFEs, length(reg$cY)), Fac2Num(reg$fe[[1L]]), seq_along(reg$cY)),
      fe2 = ifelse(rep(!noFEs, length(reg$cY)), Fac2Num(reg$fe[[2L]]), rep(1, length(reg$cY))),
      expand.model.felm(model = reg, extras = c(lat, lon), na.expand = TRUE)
    )

    data.table::setnames(dt, c("fe1", "fe2"),
      c(ifelse(!noFEs, names(reg$fe)[1L], unit),
        ifelse(!noFEs, names(reg$fe)[2L], time))
    )
    dt[, e := as.numeric(reg$residuals)]
  } else {
    stop("Model class not recognized. This function currently expects a felm object.")
  }

  n <- nrow(dt)

  orig_names <- c(unit, time, lat, lon)
  new_names <- c("unit", "time", "lat", "lon")
  data.table::setnames(dt, orig_names, new_names)

  # Spatial C++ routine requires time blocks to be contiguous. In balanced panels
  # the reusable CSR graph also assumes the within-period row order is stable, so
  # we order by unit within time.
  if (balanced_pnl) {
    data.table::setorderv(dt, c("time", "unit"))
    # The fast balanced path builds the CSR neighbor graph from the first
    # period and reuses it for every period. That is only valid if each
    # unit's coordinates are time-invariant.
    if (data.table::uniqueN(dt[["time"]]) > 1L) {
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

  X <- as.matrix(dt[, Xvars, with = FALSE])
  e <- dt[["e"]]
  invXX <- solve(crossprod(X)) * n

  # Score pre-aggregation per time block: collapse rows that share a
  # (possibly pixelated) location into a single row whose "score" is the sum
  # of the constituents' scores. At pixel = 0 this only collapses exact
  # duplicates; at pixel > 0 nearby points are snapped to a common grid
  # representative. The balanced-panel path's CSR-reuse assumption (identical
  # row count and order per period) is preserved by the dedupe key, which
  # always includes time. Both kernels weigh distance 0 at exactly 1, so
  # collapsing same-location rows is exact at pixel = 0.
  agg <- aggregate_scores(dt, Xvars, pixel = pixel,
                          balanced_pnl = balanced_pnl, verbose = verbose)

  if (verbose) message("Starting fast spatial HAC meat in C++")
  XeeX <- FastSpatialMeat(
    lat = agg$lat,
    lon = agg$lon,
    time = agg$time,
    X = agg$X,
    e = agg$e,
    cutoff = dist_cutoff,
    kernel = kernel,
    dist_fn = dist_fn,
    balanced_pnl = balanced_pnl,
    ncores = ncores
  )

  if (lag_cutoff > 0 && length(unique(dt[["time"]])) > 1L) {
    data.table::setorderv(dt, c("unit", "time"))

    if (verbose) message("Starting serial HAC meat")
    X_serial <- as.matrix(dt[, Xvars, with = FALSE])
    XeeX_serial <- FastSerialHacPanel(
      unit = dt[["unit"]],
      time = dt[["time"]],
      cutoff = lag_cutoff,
      X = X_serial,
      e = dt[["e"]],
      ncores = ncores
    )

    XeeX <- XeeX + XeeX_serial
  }

  V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
  V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2
  rownames(V_spatial_HAC) <- colnames(V_spatial_HAC) <- Xvars
  V_spatial_HAC
}

# Build per-time-block aggregated scores. Returns a list with element-wise
# vectors/matrices ready to hand to FastSpatialMeat: lat, lon, time, X, e.
#
# At pixel = 0 we collapse rows whose (lat, lon) match exactly. At pixel > 0
# we first snap (lat, lon) to a uniform grid whose latitude step is
# pixel / 111 km and whose longitude step is pixel / (111 * cos(lat_rep)) km,
# so cells stay roughly square at all latitudes. The representative coordinate
# for a cell is the cell centre.
aggregate_scores <- function(dt, Xvars, pixel, balanced_pnl, verbose) {
  n <- nrow(dt)
  X <- as.matrix(dt[, Xvars, with = FALSE])
  e <- dt[["e"]]
  scores <- X * e

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
        X = X, e = e
      ))
    }
  }

  X_agg <- as.matrix(agg[, Xvars, with = FALSE])
  list(
    lat = agg[["lat"]], lon = agg[["lon"]], time = agg[["time"]],
    X = X_agg,
    e = rep(1.0, n_agg)
  )
}

expand.model.felm <- function(model, extras, envir = environment(formula(model)),
                              na.expand = FALSE) {

  topaste <- c(names(model$fe), names(model$clustervar), extras)
  fescluext <- parse(text = paste("~", paste(topaste, collapse = "+")))[[1L]]

  data <- eval(model$call$data, envir)
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
