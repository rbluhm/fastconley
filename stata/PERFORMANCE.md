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
| run date | 2026-09-05T06:10:44Z |
| CPU | 13th Gen Intel(R) Core(TM) i7-1360P |
| logical CPUs | 16 |
| memory | 30 GB |
| OS | Ubuntu 24.04.4 LTS |
| Stata | 18 IC (4 licensed / 16 cores) |
| fastconley ado / engine build | 0.2.0 / fd27197 |
| acreg | December 2020 (1.1.0) |
| plugin thread counts | 1, 4, 8, 16 |

## Dense baseline

Small cross-sections (5 regressors, 500 km, uniform kernel, spherical distance) where the R vignette checks the engine against an explicit dense weight matrix. acreg runs comfortably here.

| observations | plugin (8 thr.) | Mata | acreg | R fastconley (8 thr.) | R dense | Mata vs plugin | acreg vs plugin |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.003 | 0.012 | 0.126 | 0.002 | 0.019 | 1.5e-18 | 6.5e-05 |
| 2,000 | 0.004 | 0.024 | 0.347 | 0.002 | 0.058 | 7.6e-19 | 3.1e-05 |
| 4,000 | 0.008 | 0.057 | 1.09 | 0.003 | 0.539 | 2.7e-19 | 1.0e-05 |

The acreg column differs from fastconley at the 1e-5 level because acreg measures distance on an equirectangular plane (111 km per degree of latitude, scaled by the cosine of the latitude for longitude) rather than on the sphere; a few pairs near the cutoff boundary change status.

## Scattered cross-sections

10 regressors, uniform kernel, spherical distance, points drawn uniformly over the contiguous United States. The R columns are the vignette's `vcovSpHAC()` and `fixest::vcov_conley()` times at the same thread count where the vignette ran one (1 and 8 threads).

| observations | cutoff km | threads | plugin | Mata (1 thr.) | acreg | R fastconley | R fixest | plugin vs Mata |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50,000 | 100 | 1 | 0.106 | 1.63 | 151 | n/a | n/a | 15x |
| 50,000 | 100 | 4 | 0.081 |  |  | n/a | n/a | 20x |
| 50,000 | 100 | 8 | 0.075 |  |  | 0.066 | 0.369 | 22x |
| 50,000 | 100 | 16 | 0.071 |  |  | n/a | n/a | 23x |
| 50,000 | 500 | 1 | 0.778 | 4.01 | 153 | 1.09 | 11.5 | 5.2x |
| 50,000 | 500 | 4 | 0.308 |  |  | n/a | n/a | 13x |
| 50,000 | 500 | 8 | 0.217 |  |  | 0.293 | 2.43 | 18x |
| 50,000 | 500 | 16 | 0.159 |  |  | n/a | n/a | 25x |
| 100,000 | 100 | 1 | 0.295 | 2.62 | 692 | n/a | n/a | 8.9x |
| 100,000 | 100 | 4 | 0.199 |  |  | n/a | n/a | 13x |
| 100,000 | 100 | 8 | 0.182 |  |  | 0.123 | 1.14 | 14x |
| 100,000 | 100 | 16 | 0.154 |  |  | n/a | n/a | 17x |
| 100,000 | 500 | 1 | 3.17 | 18.8 | 760 | n/a | n/a | 5.9x |
| 100,000 | 500 | 4 | 1.03 |  |  | n/a | n/a | 18x |
| 100,000 | 500 | 8 | 1.03 |  |  | 0.812 | 9.60 | 18x |
| 100,000 | 500 | 16 | 0.877 |  |  | n/a | n/a | 21x |

acreg's Mata loop touches every pair of observations once per observation, so its time grows with n² regardless of the cutoff, while fastconley's cell grid only visits candidate pairs within the cutoff. The Mata fallback of fastconley uses the same cell grid, so it stays proportional to the number of pairs but runs single-threaded in interpreted Mata.

## Balanced panel with serial HAC

10,000 units observed in 4 periods (40,000 rows), unit and time fixed effects absorbed, 5 regressors, 500 km, uniform kernel, and a one-lag serial Bartlett term (`lag(1) balanced`). fastconley follows the Hsiang (2010) convention used by the R package: contemporaneous spatial correlation within each period plus own-unit serial correlation. acreg's `hac` option instead applies the product of the temporal Bartlett weight and the spatial indicator to every cross-unit pair, a different estimator, so the two are timed but not compared numerically.

| method | threads | seconds | R (8 thr.) |
|---|---:|---:|---:|
| plugin | 1 | 0.174 |  |
| plugin | 4 | 0.157 |  |
| plugin | 8 | 0.144 | 0.328 |
| plugin | 16 | 0.147 |  |
| Mata | 1 | 0.313 | |
| acreg (`hac lag(1) pfe1 pfe2`) | 1 | 62.4 | |
| R fixest composition (conley + NW - hetero) | 8 | | 0.123 |

## Regular raster

A 180 × 180 latitude/longitude lattice (32,400 cells at 0.05°), 3 regressors, 250 km. The plugin has a dedicated grid engine (prefix sums for the uniform kernel, per-ring FFT convolutions for Bartlett) that the Mata fallback does not have; both are exact, so the grid and pairwise results agree to roundoff.

| kernel | plugin grid (8 thr.) | plugin pairwise (8 thr.) | Mata pairwise | acreg | R grid | R pairwise | grid vs pairwise |
|---|---:|---:|---:|---:|---:|---:|---:|
| uniform | 0.047 | 0.209 | 7.94 | 72.8 | 0.030 | 0.219 | 8.1e-20 |
| bartlett | 0.071 | 0.480 | 8.96 | 76.5 | 0.055 | 0.548 | 6.4e-20 |

## Repeated locations and pixel aggregation

100,000 rows on 20,000 distinct locations (5 rows each), 5 regressors, 250 km, Bartlett kernel. `pixel(0)` merges rows with identical coordinates exactly; `pixel(10)` and `pixel(25)` snap coordinates to a 10 km or 25 km lattice first, which is approximate. acreg has no aggregation and is not attempted at this size.

| pixel km | plugin (8 thr.) | Mata | R fastconley (8 thr.) | vs exact (pixel 0) |
|---:|---:|---:|---:|---:|
| 0 | 0.108 | 0.512 | 0.054 | 0 |
| 10 | 0.116 | 0.473 | 0.058 | 5.7e-08 |
| 25 | 0.103 | 0.389 | 0.046 | 1.3e-07 |

## One million observations

A global cross-section of 1,000,000 points (latitude -55 to 70, all longitudes), 10 regressors, 100 km, uniform kernel. A dense weight matrix would be about 7.3 TiB.

| threads | plugin | Mata | R fastconley |
|---:|---:|---:|---:|
| 1 | 2.80 | n/a | n/a |
| 4 | 2.18 |  | 1.55 |
| 8 | 2.12 |  | 1.41 |
| 16 | 2.02 |  | 1.13 |

The Mata fallback did not finish this case: aborted: exceeded the 5400 s cap (plugin: 2.0-2.8 s). Its cell grid keeps the pair count proportional to the data, but interpreted Mata pays a fixed cost per cell pair, and at one million points over the whole globe there are hundreds of thousands of occupied cells.
acreg is not attempted here: its cost grows with n², and its whole-command time at 100,000 points is already about 700 s.

## Fixed preparation cost

`cutoff(-1)` keeps only the diagonal of the meat, so `e(vce_seconds)` then measures everything except the pair enumeration and accumulation: reading the sample into Mata, sorting by period, merging identical coordinates, marshalling scores into temporary variables for the plugin, the engine's own coordinate cache, sort, and score gather (a negative cutoff still runs the engine's band path), and assembling the sandwich. It is therefore an upper bound on the Stata-side preparation, measured on data generated with a different seed than the timed runs.

| observations | plugin (8 thr.) | Mata |
|---:|---:|---:|
| 100,000 | 0.184 | 0.177 |
| 1,000,000 | 2.29 | 1.98 |

At one million rows this fixed work is most of the covariance time reported above, which is why the plugin's thread scaling looks flat there: the pair work itself is a fraction of a second at 16 threads. Splitting the fixed part between the Mata preparation and the engine's own sort and gather, and trimming the Mata side (skipping the coordinate merge when locations are unique, streaming scores to the plugin without temporary variables), is the obvious next optimisation for very large samples.

## Reading the numbers

- The plugin is the same C++ engine as the R package, so plugin and R times differ only by the front-end (Stata tempvars versus R memory aliasing), by build flags (the plugin uses -O3, R its default -O2), and by the engine changes since the R numbers were recorded.
- Thread scaling flattens beyond 8 threads on this 12-core, 16-thread laptop (hybrid P/E cores and memory bandwidth), and at one million rows the fixed preparation cost dominates (see the previous section). Stata's own licence (MP with 4 cores here) does not limit the plugin's `threads()`.
- The Mata fallback is 5 to 15 times slower than the single-threaded plugin but has the same complexity, so it remains usable at a few hundred thousand observations. It is what `engine(auto)` uses when no plugin is available for the platform.
- acreg and fastconley agree to about 1e-5 on the uniform kernel and 1e-6 on Bartlett at these cutoffs; the residual is acreg's planar distance approximation, not a difference in the estimator.

## Reproducing

```bash
# from the repository root; needs reghdfe, ftools, require, and acreg on the adopath
bash stata/bench/run_bench.sh stata/bench/results     # ~1-3 hours, mostly acreg
python3 stata/bench/render_performance.py stata/bench/results stata/PERFORMANCE.md
```

`run_bench.sh` accepts `THREADS`, `ACREG_MAX_N`, `ACREG_TIMEOUT`, `LARGE_N`, and `REPS` in the environment; the driver `stata/bench/bench_vignette.do` can also be run directly for one section with the `BENCH_*` globals it documents. Results are one CSV row per timed call with the machine and software versions attached.

