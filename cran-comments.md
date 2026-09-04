## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The spell checker flags HAC, which is the standard abbreviation for
  heteroscedasticity and autocorrelation consistent.

## Test environments

* local Ubuntu 24.04, R 4.6.1, `R CMD check --as-cran`
* win-builder, R-devel (Windows, gcc 14.3): OK, 1 note (new submission)
* macOS builder (Apple M1, macOS 14 SDK, clang 17, R 4.6.1): OK, 0 notes

## Notes

* The package compiles C++ via Rcpp and RcppArmadillo; threading uses
  `std::thread` (linked with `-pthread`), so there is no RcppParallel or
  OpenMP dependency and no `SystemRequirements`.
* All examples, tests, and vignette chunks that run at check time use at
  most 2 threads.
* No reverse dependencies.
