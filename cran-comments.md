## Resubmission

This 0.11.1 replaces fastconley 0.11.0, submitted on 2026-09-07 and still in
the incoming queue. It fixes two regressions of 0.11.0 in the felm method
(a silently changed default for panel blocking with absorbed fixed effects,
and a fragile recovery of coordinate columns from the model call) reported by
a downstream replication pipeline; see NEWS.md. No other changes.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The spell checker flags HAC, which is the standard abbreviation for
  heteroscedasticity and autocorrelation consistent.

## Test environments

* local Ubuntu 24.04, R 4.6.1, `R CMD check --as-cran` (0.11.1 tarball):
  new-submission note plus two machine-local notes (the distribution's
  `-mno-omit-leaf-frame-pointer` compiler flag; no remote clock check)
* GitHub Actions `R CMD check --as-cran --no-manual`, R release, with
  `_R_CHECK_LIMIT_CORES_=true`: ubuntu-latest, macos-latest, windows-latest,
  all OK (run for commit c548b97 / tag v0.11.0)
* win-builder R-devel (R Under development 2026-09-04 r90492 ucrt,
  x86_64-w64-mingw32, 0.11.0 tarball): Status: 1 NOTE (new submission;
  spell check on HAC) (https://win-builder.r-project.org/Q64nQgiG8KyN)
* win-builder R-release (R 4.6.1 ucrt, x86_64-w64-mingw32, 0.11.0 tarball):
  Status: 1 NOTE (new submission; spell check on HAC)
  (https://win-builder.r-project.org/K4ENW148dbpl)
* macOS builder (Apple M1, macOS 26.6, R 4.6.1 patched, clang 1700, 0.11.0
  tarball, 2026-09-05): Status: OK, 0 notes
  (https://mac.R-project.org/macbuilder/results/1788638975-a75c4a580b33fc44/)
* Sanitizers: the C++ engine (all code in `src/conley_core.h`) built with
  gcc `-fsanitize=address,undefined` and exercised through its edge-case
  probes, the deterministic golden cases, and a 100,000-point benchmark at
  1 and 4 threads: no reports.

## Notes

* The package compiles C++ via Rcpp and RcppArmadillo; threading uses
  `std::thread` (linked with `-pthread`), so there is no RcppParallel or
  OpenMP dependency and no `SystemRequirements`.
* All examples, tests, and vignette chunks that run at check time use at
  most 2 threads; the default thread count honours `_R_CHECK_LIMIT_CORES_`.
* The test suite runs in about 15 seconds.
* No reverse dependencies.
