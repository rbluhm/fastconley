# Proposal: `vce(conley ...)` for reghdfe

This directory holds a ready-to-review contribution that adds Conley (1999)
spatial HAC standard errors to reghdfe as a native `vce()` type, using the
pure-Mata engine of the `fastconley` Stata command. It has not been submitted
yet; the files here are the patch and its test.

```stata
reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300))
reghdfe y x1 x2, absorb(id year) vce(conley lat lon, cutoff(500) kernel(uniform) lag(2) unit(id) time(year))
```

Syntax: `vce(conley latvar lonvar, cutoff(#) [kernel(bartlett|uniform)
dist(haversine|spherical|chord) lag(#) unit(varname) time(varname) pixel(#)
nossc nopsdfix])`. Defaults: Bartlett kernel, haversine distance, no serial
lags, the same `N/(N-K-df_a)` small-sample factor as `vce(robust)`, and the
positive-semi-definite fix reghdfe already applies to cluster variances.

## Files

- `Conley.mata`: new file for `current-code/`. `reghdfe_vce_conley(S, sol, D,
  X, w, vce_mode)` follows `reghdfe_vce_dkraay()` line by line (scores as
  residual times weights, sandwich `D M D q`, `reghdfe_fix_psd`), plus the
  engine functions (cell-grid neighbour search with tiled dense blocks,
  vectorised serial HAC, aggregation of coincident points), derived from
  `stata/src/fastconley.mata` with a `reghdfe_conley_` prefix.
- `reghdfe-vce-conley.patch` (against reghdfe 6.14.1 `current-code/`):
  `reghdfe.ado` (a `ParseConley` subprogram, markout of the coordinate and
  panel variables, filling the HDFE fields, `compact` keep list),
  `FE.mata` (fields and defaults), `Regression.mata` (whitelist and dispatch),
  `Solution.mata` (fields, `e(vcetype) = "Conley"`, `e(title3)`,
  `e(conley_*)`), `reghdfe_header.ado` (cutoff line), `reghdfe.mata`
  (include), `reghdfe.pkg`. `build.py` regenerates `src/` as usual.
- `ftools-vce-conley.patch` (against ftools 2.50.0): `ms_parse_vce.ado`
  accepts `conley` with two coordinate arguments and returns the sub-options
  through `s(vceextra)` instead of rejecting them.
- `test_upstream.do`: runs the patched reghdfe against `fastconley` (cross
  section, uniform/chord, aweights with pixel aggregation, `noabsorb`, panel
  with two serial lags, error paths); all agree to better than 1e-12.

## How it was validated

The engine is the Mata engine of the `fastconley` Stata command, which
matches the fastconley R package to about 1e-14 on 21 configurations
(kernels, distances, weights, dateline and pole geometries, raster lattices,
balanced and unbalanced panels with serial lags, IV) and reproduces
`reghdfe, vce(robust)` exactly when the cutoff excludes every cross term.

## Notes for the maintainers

- No compiled code, no new dependencies. Cost on 100,000 scattered points
  with a 500 km cutoff is about 25 s (Bartlett) or 8 s (uniform) on a 4-core
  Stata/MP laptop; `fastconley` additionally ships a compiled engine that
  does the same in 1.6 s, which is why the standalone command exists.
- With fweights an observation with weight w counts as w identical,
  perfectly correlated observations at the same location (the Driscoll-Kraay
  convention in reghdfe), so `vce(robust)` is reproduced only under aweights,
  pweights, and unweighted fits.
- `unit()`/`time()` are explicit rather than taken from `xtset`, because
  the serial lag is measured in the units of the time variable and many
  Conley applications are cross-sections or unbalanced panels; taking
  `_dta[_TStvar]` as a default when set would be a one-line change.
