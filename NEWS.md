# fastconley 0.3.0.9000

Development build for testing a faster spatial HAC path.

## Spatial meat changes

- Adds a `FastSpatialMeat()` internal Rcpp routine.
- Replaces dense `n x n` spatial distance matrices in `vcovSpHAC()` with a screened spatial edge list.
- Computes the meat via cumulative scores: `S = e * X`, `C_i = 0.5 S_i + sum_j w_ij S_j`, and `S'C + C'S`.
- Supports both `kernel = "bartlett"` and `kernel = "uniform"`.
- Supports the existing distance function choices: haversine, spherical, chord, and flatearth.
- In balanced panels, builds one spatial edge list from the first period and reuses it for all periods.

## Notes

- This is meant as a test branch/source drop, not a CRAN-ready release.
- The implementation uses RcppParallel for CSR construction and meat accumulation.
- Rcpp exports were updated manually. Running `Rcpp::compileAttributes()` is recommended after further edits.
