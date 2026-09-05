# fastconley for Stata changelog

## Unreleased

- `fastconley, version` returns its fields in `r()`; `stata/bench/` holds the benchmark suite behind `stata/PERFORMANCE.md` (the R vignette's six sections timed with the plugin, the Mata fallback, and acreg).

- `e(vce_seconds)` reports the wall-clock time of the covariance step alone (Mata timer slot 100), the quantity the benchmarks in `stata/PERFORMANCE.md` compare with the R package.

- Shipped plugins rebuilt by CI from engine 0.11.1 (build fd27197): stripped, exporting only `stata_call`/`pginit`, GLIBC floor 2.25 on Linux, ad-hoc signed universal macOS binary, installed under their platform names by the .pkg.

- Added `fastconley_reghdfe_vce`, the first provider for the proposed generic reghdfe `vce(external PROVIDER, ...)` hook, with standalone-identical Mata/plugin numerics and provider metadata.
- Made the generic external-VCE hook the primary upstream proposal; retained the native pure-Mata `vce(conley ...)` patch as an alternative.
- Optimized the pure-Mata uniform kernel with an all-accepted tile path, conditional verbose pair counting, and reuse of per-left-tile work vectors (about 10-20% faster on 100,000 points).
- Mata uniform and Bartlett/chord tiles clamp the squared chord to 4 so exact antipodes are accepted when the cutoff covers the whole sphere (matches the compiled engine).

## 0.2.0 - 2026-09-04

- Added pweight preparation, strict numeric-option validation, grouped-expression-safe IV parsing, IV regressor collinearity handling, no-absorb IV constants, IV prediction, and fuller IV stored results.
- Added rank-safe Wald reporting, bounded Mata tile workspace, stable near-cutoff Mata acceptance, tiny-cutoff cell-key support, plugin row guards, loader diagnostics, and `fastconley, version`.
- Added platform-specific plugin installation, macOS console platform codes, license and notice files, regression tests, and a one-command R-Stata parity driver.

## 0.1.0 - 2026-09-03

- Initial Stata port with reghdfe integration, Mata and compiled engines, OLS/2SLS, panel serial HAC, pixel aggregation, and raster-grid acceleration.
