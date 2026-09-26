# fastconley for Stata changelog

## Unreleased

- Compiled-engine path: one plugin call (`vce`, engine 0.11.3) now does the whole row preparation on the raw sample rows, in C++: balanced-panel validation, aggregation of identical (time, lat, lon) keys, lattice detection and the `method(auto)` rule, the spatial meat, and the (unit, time) sort for the serial meat. Mata no longer sorts the sample three times, runs `uniqrows()` on the latitudes to reject scattered data, or copies the rows into temporary variables twice for `lag()`. Covariance step on one million scattered points at 8 threads: 1.75 to 1.09 s with `method(pairwise)` and 2.15 to 1.12 s with the default `method(auto)`; 200,000-row Bartlett panel with `lag(2)`: 0.70 to 0.22 s (`balanced`) and 0.59 to 0.32 s (unbalanced); pair-dominated runs are unchanged. e(V) is bit-identical to the Mata-prepared path whenever no two sample rows share a (time, location) key and no two share a (unit, time) key (21 configurations checked against the 0.11.1 plugin, including the grid engine with dateline wrap, balanced panels, IV, and the reghdfe hook). Where rows share a location, their scores are summed in input order rather than in the order of Mata's unstable `order()`, so e(V) can differ in the last bit (1e-17 relative in the `pixel()` checks). `pixel()` snapping stays in Mata (`fastconley_pixel_keys()`), so the snapped coordinates are unchanged. Specification errors found by the plugin (balanced-panel validation, `method(grid)` without a usable lattice) still exit with r(3498). The Mata engine keeps its own preparation. `engine(plugin)` now requires at most 2,147,483,647 sample observations (previously prepared rows after aggregation).
- Shared engine 0.11.3: the per-row k x k meat update walks the accumulator column by column (bit-identical; about 1.2x faster spatial meat at k = 100).
- `fastconley.mata`: `fastconley_prepare()` is split into `fastconley_prepare_scores()` (bread and scores) and the existing `fastconley_prepare_rows()`; the upstream full patch's generated `Conley.mata` is regenerated accordingly.
- Shared engine 0.11.2 (`src/conley_core.h`, bit-identical results): a branch-free pair loop with register-held accumulation (about 1.5-2x faster pair work), the same accumulation for the balanced-panel CSR stream, a CSR build that no longer evaluates Bartlett weights to count neighbours, and a parallel radix sort for the cell order. The shipped plugins are rebuilt by CI from engine 0.11.2 (build bf83654, stata-plugin run 36054420796) and both ados expect 0.11.2. Tested in Stata 18 on Linux against the 0.11.1 plugins: e(V) bit-identical on 21 configurations (scattered, global, balanced and unbalanced panels with lag(), raster grid with dateline wrap, pixel(), IV, weights, both kernels, the reghdfe hook provider), test_basic.do with the plugin required, the R-Stata parity driver, and both upstream tests pass; covariance step 1.4-2.0x faster single-threaded on 100,000 scattered points, about 1.2x at one million, where preparation dominates. The version moved because `spatial_meat()` gained an optional score row map, which the plugin does not use.
- Mata engine: the cell dictionary (a live Mata associative array) is no longer in scope during the pair traversal; cells are resolved through a direct-address numeric lookup (bounded at 64 MiB) or a precomputed neighbour table. A live associative array makes every Mata function call in scope roughly a hundred times slower, which is why the fallback could not finish one million points; the traversal order and results are unchanged.
- Mata engine: default `tile(512)` (was 1024); 512 is fastest for both kernels at 50k-100k points because the tile-sized work matrices stay cache resident. Only the accumulation order changes.
- Preparation: the coordinate merge returns immediately when every row is its own location, the period count is computed once and shared, the unused `Ua*Ub'` product is skipped below 200 km, and raster detection rejects on latitude before sorting longitude.

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
