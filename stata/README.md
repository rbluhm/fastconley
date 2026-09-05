# fastconley for Stata

Conley (1999) spatial HAC standard errors for linear models with many absorbed
fixed effects, built on [reghdfe](https://github.com/sergiocorreia/reghdfe).
This is the Stata port of the [fastconley](https://github.com/rbluhm/fastconley)
R package and reproduces `vcovSpHAC()` on the same fit (validated by
`test/parity/`).

```stata
ssc install ftools
ssc install reghdfe     // >= 6.12.5; reghdfe 6.13+ also needs: ssc install require
ssc install require
net install fastconley, from("https://raw.githubusercontent.com/rbluhm/fastconley/main/stata/src/")

fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300)
fastconley y x1 x2, absorb(id year) lat(lat) lon(lon) cutoff(500) ///
    kernel(uniform) unit(id) time(year) lag(2) balanced
fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)   // 2SLS
```

OLS retains reghdfe's output, `e()` surface, footnote, and prediction. IV
models use the `(endog = instruments)` syntax and are estimated by 2SLS on
the partialled-out data; they post the core fit/VCE results listed in
`help fastconley`, plus `predict` support for `xb` and saved residuals. With
`noabsorb`, IV includes `_cons` unless `noconstant` is specified.

`time()` alone blocks the spatial covariance by time. Serial `lag()` requires
both `unit()` and `time()`. Numeric-looking string times retain their numeric
scale (`"2000"` to `"2002"` is a gap of two); other strings can define spatial
blocks but are rejected for positive serial lags.

The Stata command covers linear OLS and 2SLS models. The R package's GLM
support, including `feglm`/`fepois`, has no Stata analogue.

## Option mapping

| Stata | R `vcovSpHAC()` |
|-------|-----------------|
| `cutoff(#)` | `dist_cutoff` |
| `kernel(bartlett|uniform)` | `kernel` |
| `distance(haversine|spherical|chord)` | `dist_fn` |
| `unit()`, `time()`, `lag(#)` | `unit`, `time`, `lag_cutoff` |
| `balanced` | `balanced_pnl = TRUE` |
| `pixel(#)` | `pixel` |
| `nossc` | `ssc = FALSE` |
| `nopsdfix` | `psd_fix = FALSE` |

Coming from `acreg`: `latitude()`/`longitude()` are the same names,
`distcutoff()` is `cutoff()`, `lagcutoff()` is `lag()`, `id()`/`time()` are
`unit()`/`time()`, `bartlett` is the default kernel (use `kernel(uniform)` for
acreg's default), and fixed effects go in `absorb()` with no limit on their
number.

## Engines

Two interchangeable engines compute the meat; `e(engine)` records which ran.

- **plugin** (default when it loads): the R package's C++ engine
  (`src/conley_core.h`) compiled as a Stata plugin. Multithreaded
  (`threads()`, default `c(processors_mach)`, the machine's processor count;
  the plugin's threads are independent of the Stata licence, so Stata/SE and
  /BE get the same speed-up as /MP), streaming (no pair list in memory),
  and the only engine with the exact raster **grid engine**
  (`method(auto|pairwise|grid)`), which is pair-count independent on regular
  lat/lon lattices and wraps across the dateline. `fastconley.pkg` installs
  the binary under its platform name: `fastconley_linux64.plugin`,
  `fastconley_win64.plugin`, or universal `fastconley_macosx.plugin` (including
  console macOS). A generic `fastconley.plugin` remains a local-install
  fallback; the engine version is checked on load.
- **mata**: pure Mata (cell-grid neighbour search at cell-pair granularity
  with tiled dense blocks, per-period stacking for balanced panels, vectorised
  serial HAC). Always available; `engine(auto)` falls back to it when no
  binary loads.

The 21-configuration parity harness (including three IV models) observed a
worst relative V difference of 2.0e-11, in the unbalanced serial lag-2 case;
the default pass tolerance is 1e-8.

Indicative timings measured 2026-09-04 on a 4-core Stata/MP Linux laptop
(16 logical CPUs), 100,000 scattered points with one absorbed FE, `k = 4`
(including `_cons`), and a 500 km cutoff: plugin Bartlett 4.8 s / 1.5 s /
0.8 s at 1 / 4 / 16
threads, uniform 2.5 s / 0.8 s / 0.5 s; Mata (best of two runs in one Stata
instance) 21.8 s Bartlett/haversine and 23.0 s uniform/spherical. A balanced
40,000 x 3 panel at 500 km with one serial lag: plugin
about 1 s, Mata 5 s with `balanced`. A 216,000-cell 0.1-degree raster at
300 km: grid engine 0.7 s, pairwise 1.0 s. For comparison, `acreg` needs
1.6 s for 5,000 points where `fastconley` needs 0.1 s, and its cost grows
with the square of the sample.

Bartlett distances in both engines come from the chord between the points'
unit vectors (`2R asin(|u_i - u_j| / 2)`), which needs no per-pair `sin` or
`atan2` and stays accurate at small angles. The Mata engine builds the chord
from coordinate differences below 200 km and from the dot product above. At
the 200-km switch, the measured bound is 5e-12 relative in squared chord and
5e-10 km absolute in distance.

The Mata `tile()` default is 1024. It has a hard cap of 8192 and is also
rejected when the conservative `5 * tile^2 * 8` workspace estimate exceeds
1 GiB. Run `fastconley, version` for the ado version, expected engine version,
loaded file/build, status, and loader-attempt return codes.

The model Wald F uses the covariance before PSD clamping. It is stored as
missing when the slope block is rank deficient, so the eigenvalue floor cannot
manufacture an enormous statistic.

For `cutoff(-1)`, unweighted, aweighted, and pweighted VCEs reproduce
`reghdfe, vce(robust)`. With fweights the Conley score outer product is
`(w e x)(w e x)'`: the `w` copies are spatially identical points. That is the
appropriate Conley interpretation but differs from robust's
`w (e x)(e x)'`.

Coordinates stored as `float` (Stata's default) carry about 6e-8 relative
rounding noise. The command loosens its lattice-detection tolerance for float
variables and the grid engine then treats the lattice as exact, which agrees
with the pairwise engine to about 1e-12; storing coordinates as `double`
makes the two identical.

### Building the plugin

`stata/plugin/` holds the SPI 3.0 front-end and a Makefile with `linux`,
`macos` (universal Mach-O), and `windows` (mingw-w64, static) targets; only the
Armadillo headers are needed (`make linux ARMA_INC=/path/to/include`). The
GitHub Actions workflow `.github/workflows/stata-plugin.yml` builds all three
and uploads them as artifacts; copy them into `stata/src/` before a release.
The engine version stamped into the binary (`CONLEY_CORE_VERSION`) must match
the ado, otherwise the ado ignores the binary and uses Mata.

## Tests

From the repository root, with Stata in batch mode:

```bash
stata-mp -b do stata/test/test_basic.do
stata/test/parity/run_parity.sh
```

The parity driver runs `gen_reference.R`, Stata with both Mata and plugin
engines, and `compare.R`. It checks the Stata log because Stata batch mode can
return shell status 0 after a do-file error. Override the comparison tolerance
with `FASTCONLEY_PARITY_TOL`; the default is `1e-8`.

Note that `reghdfe` drops singleton groups by default and fixest does not;
the parity harness passes `keepsingletons`.

On Windows, run the same do-files with `StataSE-64.exe /e do stata\test\test_basic.do`
(or `StataMP-64.exe`) from the repository root in PowerShell or cmd. From Git
Bash the `/e` flag is rewritten as a path; use `//e` or set
`MSYS_NO_PATHCONV=1`. The log lands in the working directory.

Reference timings from a Windows 11 laptop (Stata/SE, i7-1265U) before the
chord-form change: plugin 106 s single-threaded / 15 s with `threads(8)`,
Mata 169 s. That machine traced the single-thread gap to mingw-w64's x87
`sin`/`atan2`; the chord form removes those calls (46.6 s to 7.1 s on its
standalone benchmark) and CI keeps building with mingw.

## Release checklist

1. Push (or run the `stata-plugin` workflow manually) and download the three
   artifacts into `stata/src/` as `fastconley_linux64.plugin`,
   `fastconley_win64.plugin`, `fastconley_macosx.plugin`. The Linux one can
   also be built locally with `make -C stata/plugin linux`.
2. Bump `*! version` in `fastconley.ado`, `{* *! version` in the help file,
   and `Distribution-Date` in `fastconley.pkg`; keep `CONLEY_CORE_VERSION`
   in `src/conley_core.h` and the single expected-version constant at the top of `LoadPlugin` in sync
   with the R package version when the engine changes.
3. From the repository root: `stata-mp -b do stata/test/test_basic.do`, the
   parity harness on both engines, and `stata/upstream/test_upstream.do` if
   the upstream patch is being refreshed.
4. SSC: zip the contents of `stata/src/` (ado, mata, sthlp, pkg, toc, the
   three plugins) and send them to the SSC maintainer as described at
   https://ideas.repec.org/c/boc/bocode/ (package name `fastconley`).
   GitHub installs keep working from `stata/src/` at any time.

## Performance

`stata/PERFORMANCE.md` is the Stata counterpart of the R package's performance vignette: the same six benchmark sections (dense baseline, 50k/100k cross-sections, a 10,000-unit panel with serial HAC, a 180 x 180 raster, repeated locations with pixel aggregation, one million points) timed through `fastconley` with the compiled plugin at 1/4/8/16 threads and with the Mata fallback, next to acreg where its O(n^2) loop is feasible and next to the R numbers recorded on the same machine. `e(vce_seconds)` is the covariance-step time each run reports. Regenerate with `bash stata/bench/run_bench.sh` and `python3 stata/bench/render_performance.py`.

## Upstream proposal

The primary proposal in `stata/upstream/` is a small generic hook for
reghdfe 6.15.0. With that patch installed, the same OLS fit can be requested
without making reghdfe own the Conley algorithm:

```stata
reghdfe y x1 x2, absorb(region) ///
    vce(external fastconley, lat(lat) lon(lon) cutoff(300))
```

`fastconley_reghdfe_vce.ado` is the first provider. It uses the same Mata or
compiled engine, option validation, small-sample factor, and PSD repair as the
standalone command; reghdfe retains its native posting, replay, `predict`,
`test`, and `margins` behavior. The hook is OLS-only. The older full patch,
which adds a pure-Mata `vce(conley ...)` implementation directly to reghdfe,
remains in `stata/upstream/` as an alternative. See that directory's README
for both trade-offs and reproducible tests.
