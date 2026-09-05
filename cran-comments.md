## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The spell checker flags HAC, which is the standard abbreviation for
  heteroscedasticity and autocorrelation consistent.

## Test environments

* local Ubuntu 24.04, R 4.6.1, `R CMD check --as-cran` (0.11.0 tarball):
  new-submission note plus two machine-local notes (the distribution's
  `-mno-omit-leaf-frame-pointer` compiler flag; no remote clock check)
* GitHub Actions `R CMD check --as-cran --no-manual`, R release, with
  `_R_CHECK_LIMIT_CORES_=true`: ubuntu-latest, macos-latest, windows-latest,
  all OK (run for commit c548b97 / tag v0.11.0)
* win-builder R-devel and R-release (0.11.0 tarball uploaded 2026-09-05):
  RESULTS_PENDING
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
