# fastconley for Stata: performance

This is the Stata counterpart of the R package's performance vignette. The same six benchmark sections are run through the `fastconley` command with the compiled plugin at several thread counts and with the pure-Mata fallback, and, where its O(n²) loop is feasible, through acreg (Colella, Lalive, Sakalli and Thoenig, 2019, "Inference with Arbitrary Clustering", IZA Discussion Paper 12584; on SSC), the established Stata implementation of Conley standard errors. The tables were produced by `stata/bench/run_bench.sh` on the machine below and rendered by `stata/bench/render_performance.py`; they are not regenerated automatically.

## What is timed

- **fastconley (plugin, Mata)**: `e(vce_seconds)`, the covariance step alone, after reghdfe has partialled out the fixed effects and solved the regression. This is the same quantity the R vignette reports for `vcovSpHAC()` (post-estimation only). The whole-command time including reghdfe is in the results file (`total_seconds`) and is typically 0.03 to 0.05 s longer on the small cases and a few seconds longer at one million rows.
- **acreg**: the whole command, because acreg estimates and corrects in one pass and exposes no separate timer. On these configurations its regression is a negligible part of the total.
- **R (same machine)**: the numbers shipped with the R package (`inst/benchmarks/`), run on the same CPU on 2026-06-10 with fastconley 0.8.0; the engine has changed since (chord-form weights, RcppParallel replaced by a std::thread pool with a serial sort), so the R column is a historical reference on the same hardware, not a current measurement.
- Every fastconley call uses `nossc nopsdfix` so that its covariance is comparable to acreg, which applies no small-sample correction. The relative-difference columns compare the slope block of `e(V)` with the plugin result on the same data using Stata's `mreldif`, which divides each element's difference by one plus the reference value; for covariance entries far below one that is close to an absolute difference, so treat those columns as evidence of agreement, not as a norm-relative error comparable to the R package's checks. Each fastconley time is the repetition with the smaller covariance time (of two), with its own total.

## Machine

| item | value |
|---|---|
| run date | 2026-09-05T17:27:50Z; plugin rows re-timed 2026-09-26 (engine 0.11.3, build 5f7cfdd; plugin rows only, same machine; Mata and acreg rows are from run_date) |
| CPU | 13th Gen Intel(R) Core(TM) i7-1360P |
| logical CPUs | 16 |
| memory | 30 GB |
| OS | Ubuntu 24.04.4 LTS |
| Stata | 18 IC (4 licensed / 16 cores) |
| fastconley ado / engine build | 0.2.0 / 5f7cfdd |
| acreg | December 2020 (1.1.0) |
| plugin thread counts | 1, 4, 8, 16 |

## Dense baseline

Small cross-sections (5 regressors, 500 km, uniform kernel, spherical distance) where the R vignette checks the engine against an explicit dense weight matrix. acreg runs comfortably here.

| observations | plugin (8 thr.) | Mata | acreg | R fastconley (8 thr.) | R dense | Mata vs plugin | acreg vs plugin |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.001 | 0.012 | 0.185 | 0.002 | 0.019 | 1.4e-18 | 6.5e-05 |
| 2,000 | 0.003 | 0.023 | 0.406 | 0.002 | 0.058 | 7.6e-19 | 3.1e-05 |
| 4,000 | 0.003 | 0.050 | 1.14 | 0.003 | 0.539 | 2.2e-19 | 1.0e-05 |

The acreg column differs from fastconley at the 1e-5 level because acreg measures distance on an equirectangular plane (111 km per degree of latitude, scaled by the cosine of the latitude for longitude) rather than on the sphere; a few pairs near the cutoff boundary change status.

## Scattered cross-sections

10 regressors, uniform kernel, spherical distance, points drawn uniformly over the contiguous United States. The R columns are the vignette's `vcovSpHAC()` and `fixest::vcov_conley()` times at the same thread count where the vignette ran one (1 and 8 threads).

| observations | cutoff km | threads | plugin | Mata (1 thr.) | acreg | R fastconley | R fixest | plugin vs Mata |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50,000 | 100 | 1 | 0.041 | 0.673 | 179 | n/a | n/a | 16x |
| 50,000 | 100 | 4 | 0.031 |  |  | n/a | n/a | 22x |
| 50,000 | 100 | 8 | 0.027 |  |  | 0.066 | 0.369 | 25x |
| 50,000 | 100 | 16 | 0.026 |  |  | n/a | n/a | 26x |
| 50,000 | 500 | 1 | 0.329 | 3.16 | 181 | 1.09 | 11.5 | 9.6x |
| 50,000 | 500 | 4 | 0.141 |  |  | n/a | n/a | 22x |
| 50,000 | 500 | 8 | 0.095 |  |  | 0.293 | 2.43 | 33x |
| 50,000 | 500 | 16 | 0.071 |  |  | n/a | n/a | 44x |
| 100,000 | 100 | 1 | 0.131 | 1.45 | 777 | n/a | n/a | 11x |
| 100,000 | 100 | 4 | 0.079 |  |  | n/a | n/a | 18x |
| 100,000 | 100 | 8 | 0.069 |  |  | 0.123 | 1.14 | 21x |
| 100,000 | 100 | 16 | 0.057 |  |  | n/a | n/a | 25x |
| 100,000 | 500 | 1 | 1.27 | 17.3 | 782 | n/a | n/a | 14x |
| 100,000 | 500 | 4 | 0.428 |  |  | n/a | n/a | 40x |
| 100,000 | 500 | 8 | 0.313 |  |  | 0.812 | 9.60 | 55x |
| 100,000 | 500 | 16 | 0.243 |  |  | n/a | n/a | 71x |

acreg's Mata loop touches every pair of observations once per observation, so its time grows with n² regardless of the cutoff, while fastconley's cell grid only visits candidate pairs within the cutoff. The Mata fallback of fastconley uses the same cell grid, so it stays proportional to the number of pairs but runs single-threaded in interpreted Mata.

## Balanced panel with serial HAC

10,000 units observed in 4 periods (40,000 rows), unit and time fixed effects absorbed, 5 regressors, 500 km, uniform kernel, and a one-lag serial Bartlett term (`lag(1) balanced`). fastconley follows the Hsiang (2010) convention used by the R package: contemporaneous spatial correlation within each period plus own-unit serial correlation. acreg's `hac` option instead applies the product of the temporal Bartlett weight and the spatial indicator to every cross-unit pair, a different estimator, so the two are timed but not compared numerically.

| method | threads | seconds | R (8 thr.) |
|---|---:|---:|---:|
| plugin | 1 | 0.045 |  |
| plugin | 4 | 0.029 |  |
| plugin | 8 | 0.026 | 0.328 |
| plugin | 16 | 0.025 |  |
| Mata | 1 | 0.312 | |
| acreg (`hac lag(1) pfe1 pfe2`) | 1 | 75.3 | |
| R fixest composition (conley + NW - hetero) | 8 | | 0.123 |

## Regular raster

A 180 × 180 latitude/longitude lattice (32,400 cells at 0.05°), 3 regressors, 250 km. The plugin has a dedicated grid engine (prefix sums for the uniform kernel, per-ring FFT convolutions for Bartlett) that the Mata fallback does not have; both are exact, so the grid and pairwise results agree to roundoff.

| kernel | plugin grid (8 thr.) | plugin pairwise (8 thr.) | Mata pairwise | acreg | R grid | R pairwise | grid vs pairwise |
|---|---:|---:|---:|---:|---:|---:|---:|
| uniform | 0.011 | 0.097 | 7.12 | 65.1 | 0.030 | 0.219 | 8.1e-20 |
| bartlett | 0.020 | 0.242 | 9.85 | 69.3 | 0.055 | 0.548 | 6.4e-20 |

## Repeated locations and pixel aggregation

100,000 rows on 20,000 distinct locations (5 rows each), 5 regressors, 250 km, Bartlett kernel. `pixel(0)` merges rows with identical coordinates exactly; `pixel(10)` and `pixel(25)` snap coordinates to a 10 km or 25 km lattice first, which is approximate. acreg has no aggregation and is not attempted at this size.

| pixel km | plugin (8 thr.) | Mata | R fastconley (8 thr.) | vs exact (pixel 0) |
|---:|---:|---:|---:|---:|
| 0 | 0.038 | 0.516 | 0.054 | 0 |
| 10 | 0.047 | 0.526 | 0.058 | 5.7e-08 |
| 25 | 0.046 | 0.373 | 0.046 | 1.3e-07 |

## One million observations

A global cross-section of 1,000,000 points (latitude -55 to 70, all longitudes), 10 regressors, 100 km, uniform kernel. A dense weight matrix would be about 7.3 TiB.

| threads | plugin | Mata | R fastconley |
|---:|---:|---:|---:|
| 1 | 1.16 | 15.4 | n/a |
| 4 | 0.949 |  | 1.55 |
| 8 | 0.757 |  | 1.41 |
| 16 | 0.874 |  | 1.13 |

The Mata fallback agrees with the plugin to 3.5e-20 relative on this case.
acreg is not attempted here: its cost grows with n², and its whole-command time at 100,000 points is already about 700 s.

## Fixed preparation cost

`cutoff(-1)` keeps only the diagonal of the meat, so `e(vce_seconds)` then measures everything except the pair enumeration and accumulation. With the plugin (engine 0.11.3 and later) that is: computing the scores and the bread in Mata, copying the raw sample rows into temporary variables, the plugin reading them through Stata's plugin interface, its own row preparation in C++ (sorting, merging identical coordinates, the lattice check, the unit-time sort for `lag()`), the engine's coordinate cache, sort, and score gather (a negative cutoff runs the engine's band path, whose sort is more expensive than the cell grid's), and assembling the sandwich. The Mata column is the fallback, which prepares the rows in Mata. The runs use data generated with a different seed than the timed runs.

| observations | plugin (8 thr.) | Mata |
|---:|---:|---:|
| 100,000 | 0.063 | 0.081 |
| 1,000,000 | 0.961 | 1.15 |

Until engine 0.11.2 the ado prepared the rows in Mata: three sorts of the sample, a `uniqrows()` of the latitudes to reject scattered data under `method(auto)`, and a second copy of the rows for `lag()`, which made this fixed work most of the covariance time at one million rows. Engine 0.11.3 moved that preparation into the plugin, which halved the one-million-row time. What remains fixed on the Stata side is the data movement: in a timer split of the one-million-row case at 8 threads, computing scores and bread took about 0.17 s and copying the rows into temporary variables about 0.25 s, and the plugin reads those 13 million values one call at a time through the plugin interface. Because the no-pairs run takes the band path, it can exceed the full cell-grid run; read it as an upper bound.

## Reading the numbers

- The plugin is the same C++ engine as the R package, so plugin and R times differ only by the front-end (Stata tempvars versus R memory aliasing), by build flags (the plugin uses -O3, R its default -O2), and by the engine changes since the R numbers were recorded.
- Thread scaling flattens beyond 8 threads on this 12-core, 16-thread laptop (hybrid P/E cores and memory bandwidth), and at one million rows moving the data between Stata and the plugin is a large share of the time (see the previous section). Stata's own licence (MP with 4 cores here) does not limit the plugin's `threads()`.
- On the scattered cross-sections the Mata fallback is 10 to 16 times slower than the single-threaded plugin (the two columns come from runs three weeks apart on the same machine, see the run date above) but has the same complexity, so it remains usable up to a million observations. It is what `engine(auto)` uses when no plugin is available for the platform.
- acreg and fastconley agree to about 1e-5 on the uniform kernel and 1e-6 on Bartlett at these cutoffs; the residual is acreg's planar distance approximation, not a difference in the estimator.

## Reproducing

```bash
# from the repository root; needs reghdfe, ftools, require, and acreg on the adopath
bash stata/bench/run_bench.sh stata/bench/results     # ~1-3 hours, mostly acreg
python3 stata/bench/render_performance.py stata/bench/results stata/PERFORMANCE.md
```

`run_bench.sh` accepts `THREADS`, `ACREG_MAX_N`, `ACREG_TIMEOUT`, `LARGE_N`, and `REPS` in the environment; the driver `stata/bench/bench_vignette.do` can also be run directly for one section with the `BENCH_*` globals it documents. Results are one CSV row per timed call with the machine and software versions attached.

