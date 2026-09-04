# Cross-language parity harness, step 3.
#   Rscript stata/test/parity/compare.R OUT_DIR [tolerance]
# Reports max relative difference between the R package's vcov and the Stata
# command's e(V) (slopes only) per configuration; exits 1 if any exceeds tol.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; tol <- if (length(args) > 1) as.numeric(args[2]) else 1e-8
cfg <- read.csv(file.path(out, "configs.csv"), stringsAsFactors = FALSE)
rd <- function(f) unname(as.matrix(read.csv(f, header = FALSE)))
worst <- 0; bad <- 0L
for (nm in cfg$name) {
  VR <- rd(file.path(out, paste0(nm, "_R.csv"))); VS <- rd(file.path(out, paste0(nm, "_S.csv")))
  bR <- rd(file.path(out, paste0(nm, "_Rb.csv"))); bS <- rd(file.path(out, paste0(nm, "_Sb.csv")))
  rel <- max(abs(VR - VS)) / max(abs(VR))
  relb <- max(abs(c(bR) - c(bS))) / max(abs(bR))
  ok <- rel <= tol
  if (!ok) bad <- bad + 1L
  worst <- max(worst, rel)
  cat(sprintf("%-42s V rel diff %.2e   b rel diff %.2e   %s\n", nm, rel, relb, if (ok) "ok" else "FAIL"))
}
cat(sprintf("%d configs, worst V rel diff %.2e, tolerance %.0e, %d failures\n", nrow(cfg), worst, tol, bad))
quit(status = if (bad) 1L else 0L)
