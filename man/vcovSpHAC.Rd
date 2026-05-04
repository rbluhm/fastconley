\name{vcovSpHAC}
\alias{vcovSpHAC}
\title{Spatial HAC variance-covariance matrix for felm models}
\usage{
vcovSpHAC(reg, unit = NULL, time = NULL, lat, lon,
          kernel = c("bartlett", "uniform"),
          dist_fn = c("haversine", "spherical", "chord", "flatearth"),
          dist_cutoff, lag_cutoff = 0, verbose = FALSE,
          balanced_pnl = FALSE, ncores = NA, maxobsmem = 50000L)
}
\arguments{
\item{reg}{A model estimated by \code{lfe::felm()}.}
\item{unit}{Unit identifier for panel data. If omitted, a cross-section is assumed.}
\item{time}{Time identifier for panel data. If omitted, a cross-section is assumed.}
\item{lat, lon}{Names of latitude and longitude variables in the model data.}
\item{kernel}{Spatial kernel: \code{"bartlett"} or \code{"uniform"}.}
\item{dist_fn}{Distance function.}
\item{dist_cutoff}{Spatial cutoff in kilometers.}
\item{lag_cutoff}{Serial lag cutoff.}
\item{verbose}{Whether to print progress messages.}
\item{balanced_pnl}{Whether the panel is balanced and unit ordering is identical across periods. If true, the spatial neighbor graph is computed once and reused.}
\item{ncores}{Number of threads for RcppParallel.}
\item{maxobsmem}{Kept for backward compatibility; the fast spatial path no longer uses dense distance matrices.}
}
\description{
Computes Conley-style spatial HAC variance-covariance matrices for \code{felm()} models. The fast spatial path uses score accumulation, latitude/longitude screening, and CSR neighbor lists rather than dense distance matrices.
}
