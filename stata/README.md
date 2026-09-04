# fastconley for Stata

Conley (1999) spatial HAC standard errors for linear models with many absorbed
fixed effects, built on [reghdfe](https://github.com/sergiocorreia/reghdfe).
This is the Stata port of the [fastconley](https://github.com/rbluhm/fastconley)
R package and reproduces `vcovSpHAC()` on the same fit (validated by
`test/parity/`).

```stata
net install fastconley, from("https://raw.githubusercontent.com/rbluhm/fastconley/main/stata/src/")
ssc install reghdfe     // >= 6.12.5, plus ftools

fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300)
fastconley y x1 x2, absorb(id year) lat(lat) lon(lon) cutoff(500) ///
    kernel(uniform) unit(id) time(year) lag(2) balanced
fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)   // 2SLS
```

Everything reghdfe reports is reported (`e()`, footnote, `predict`); only the
variance, `e(F)`, and `e(vcetype)` change. IV models use the `(endog =
instruments)` syntax and are estimated by 2SLS on the partialled-out data;
point estimates match ivreghdfe/ivreg2 and, with a negative cutoff, so does
their robust variance. See `help fastconley`.

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
  (`threads()`, default `c(processors)`), streaming (no pair list in memory),
  and the only engine with the exact raster **grid engine**
  (`method(auto|pairwise|grid)`), which is pair-count independent on regular
  lat/lon lattices and wraps across the dateline. `fastconley.pkg` installs
  the binary matching the platform (`g` lines: Linux x86-64, Windows x86-64,
  macOS universal); the version is checked on load.
- **mata**: pure Mata (cell-grid neighbour search at cell-pair granularity
  with tiled dense blocks, per-period stacking for balanced panels, vectorised
  serial HAC). Always available; `engine(auto)` falls back to it when no
  binary loads.

Both agree with each other and with the R package to about 1e-14
(`test/parity/`, 21 configurations on both engines, including three IV models).

Indicative timings on a 4-core Stata/MP laptop (16 logical CPUs), 100,000
scattered points with a 500 km cutoff: plugin Bartlett 10 s / 3.0 s / 1.6 s at
1 / 4 / 16 threads, uniform 2.5 s / 0.8 s / 0.5 s; Mata 25 s Bartlett, 8 s
uniform. A balanced 40,000 x 3 panel at 500 km with one serial lag: plugin
1.0 s, Mata 5 s with `balanced`. A 216,000-cell 0.1-degree raster at 300 km:
grid engine 0.7 s, pairwise 1.0 s. For comparison, `acreg` needs 1.6 s for
5,000 points where `fastconley` needs 0.1 s, and its cost grows with the
square of the sample.

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
stata-mp -b do stata/test/test_basic.do          # smoke, robust-equivalence, engines, raster, panel, errors
Rscript stata/test/parity/gen_reference.R OUT     # R reference (fixest + fastconley)
stata-mp -b do stata/test/parity/run_stata.do    # set globals PARITY_DIR and PARITY_ENGINE (mata|plugin)
Rscript stata/test/parity/compare.R OUT 1e-8
```

Note that `reghdfe` drops singleton groups by default and fixest does not;
the parity harness passes `keepsingletons`.

## Release checklist

1. Push (or run the `stata-plugin` workflow manually) and download the three
   artifacts into `stata/src/` as `fastconley_linux64.plugin`,
   `fastconley_win64.plugin`, `fastconley_macosx.plugin`. The Linux one can
   also be built locally with `make -C stata/plugin linux`.
2. Bump `*! version` in `fastconley.ado`, `{* *! version` in the help file,
   and `Distribution-Date` in `fastconley.pkg`; keep `CONLEY_CORE_VERSION`
   in `src/conley_core.h` and the expected version in `LoadPlugin` in sync
   with the R package version when the engine changes.
3. From the repository root: `stata-mp -b do stata/test/test_basic.do`, the
   parity harness on both engines, and `stata/upstream/test_upstream.do` if
   the upstream patch is being refreshed.
4. SSC: zip the contents of `stata/src/` (ado, mata, sthlp, pkg, toc, the
   three plugins) and send them to the SSC maintainer as described at
   https://ideas.repec.org/c/boc/bocode/ (package name `fastconley`).
   GitHub installs keep working from `stata/src/` at any time.

## Upstream proposal

`stata/upstream/` contains a ready-to-review contribution adding
`vce(conley latvar lonvar, cutoff(#) ...)` to reghdfe itself (pure Mata,
same engine); see its README.
