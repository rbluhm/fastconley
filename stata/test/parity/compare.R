# Cross-language parity harness, step 3.
#   FASTCONLEY_PARITY_TOL=1e-8 Rscript stata/test/parity/compare.R OUT_DIR
# Reports max relative differences between the R package's vcov/coefficients
# and Stata's e(V)/e(b); exits 1 if either exceeds the environment tolerance.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("OUT_DIR is required")
out <- args[1]
tol <- as.numeric(Sys.getenv("FASTCONLEY_PARITY_TOL", "1e-8"))
if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("FASTCONLEY_PARITY_TOL must be a positive finite number")
cfg <- read.csv(file.path(out, "configs.csv"), stringsAsFactors = FALSE)
rd <- function(f) unname(as.matrix(read.csv(f, header = FALSE)))
worst <- 0; bad <- 0L
for (nm in cfg$name) {
  VR <- rd(file.path(out, paste0(nm, "_R.csv"))); VS <- rd(file.path(out, paste0(nm, "_S.csv")))
  bR <- rd(file.path(out, paste0(nm, "_Rb.csv"))); bS <- rd(file.path(out, paste0(nm, "_Sb.csv")))
  rel <- max(abs(VR - VS)) / max(abs(VR), .Machine$double.eps)
  relb <- max(abs(c(bR) - c(bS))) / max(abs(bR), .Machine$double.eps)
  ok <- rel <= tol && relb <= tol
  if (!ok) bad <- bad + 1L
  worst <- max(worst, rel)
  cat(sprintf("%-42s max V rel diff %.2e   max b rel diff %.2e   %s\n", nm, rel, relb, if (ok) "ok" else "FAIL"))
}
cat(sprintf("%d configs, worst max V rel diff %.2e, tolerance %.0e, %d failures\n", nrow(cfg), worst, tol, bad))
quit(status = if (bad) 1L else 0L)
