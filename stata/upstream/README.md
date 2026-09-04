# Proposal: native `vce(conley ...)` for reghdfe

This directory contains a reviewable proposal for Conley (1999) spatial HAC
standard errors in reghdfe's native OLS estimator. It has not been submitted.

```stata
reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300))
reghdfe y x1 x2, absorb(id year) ///
    vce(conley lat lon, cutoff(500) kernel(uniform) lag(2) unit(id) time(year))
```

The grammar is `vce(conley latvar lonvar, cutoff(#)
[kernel(bartlett|uniform) distance(haversine|spherical|chord) lag(#)
unit(varname) time(varname) pixel(#) balanced nossc nopsdfix])`. `con` and
`conl` are accepted abbreviations. String suboption values are normalized to
lower case. All Conley syntax-validation failures return `r(198)`.

Defaults are Bartlett weights, haversine distance, no serial lags, no pixel
aggregation, robust-style `N/(N-df_m-df_a)` small-sample scaling, and PSD
repair. The repair follows fastconley/fixest semantics: eigenvalues at or below
zero are replaced by `1e-16`, with a note only when the matrix changes by more
than `1e-8`; `nopsdfix` leaves it unchanged and reports a noticeable failure.

## Scope and limitations

Native `vce(conley)` is OLS-only. ivreghdfe rejects VCE types other than
unadjusted, robust, and cluster before calling reghdfe, so this patch does not
make Conley reachable for IV. The standalone fastconley command is separately
validated for IV and raster/FFT workloads; neither claim applies to this native
integration, whose engine is pairwise pure Mata.

With `group()` or `individual()`, latitude, longitude, unit, and time variables
are included in `ValidateGroups`; they must be constant within group. The
`balanced` suboption validates equal period sizes, identical unit membership,
and time-invariant coordinates, while using the same exact general covariance
path. Numeric-looking string times keep their numeric spacing; other string
keys are grouped in sorted-value order.

For fweights, Conley treats an observation of weight `w` as `w` colocated,
perfectly correlated observations. Consequently the negative-cutoff diagnostic
used by standalone development tests is not equal to reghdfe's robust VCE for
fweights. Unweighted, aweight, and pweight score construction follows
reghdfe's normalized-weight conventions.

## Files and design

- `make_conley_mata.py` embeds the reghdfe front end and discovers every
  top-level `fastconley_*` function in `stata/src/fastconley.mata`. It applies
  the `reghdfe_conley_` namespace, rewrites standalone-only diagnostics, and
  deterministically generates `Conley.mata` without reading the old output.
- `reghdfe-vce-conley.patch` adds `Conley.mata`, local pre-parsing before stock
  `ms_parse_vce`, a generic `vce_extra_keepvars` field used by compact loading,
  full `e(conley_*)` posting and replay display, help, version/date/changelog
  updates, and `test/part1/conley.do` in the upstream certification layout.
- `test_upstream.do` compares the patched command with standalone fastconley
  and exercises parsing, replay/posting, compact and `poolsize(1)`, aweights,
  fweights, native pweights, panel lags, stale `s()` state, PSD repair, and
  group-coordinate validation.

No ftools patch is needed. reghdfe recognizes `con`/`conl` before calling
`ms_parse_vce`, parses the complete Conley clause locally, hands stock ftools a
valid robust VCE, and then restores `vcetype="conley"`. The declared dependency
therefore remains stock ftools 2.50.0.

## Rebuild and validation

The patches are pinned to:

- reghdfe `29b2b203f534413369d675d0c0f674d61291f536` (tag 6.14.1 and the
  verified current `master` at preparation time);
- ftools `7b3663e49ea5c5b81638c55be29edf416e68e8b7` (`*! version 2.50.0
  09jan2026`), unmodified.

From the fastconley repository root, regenerate after every engine change:

```bash
python3 stata/upstream/make_conley_mata.py
git diff --exit-code -- stata/upstream/Conley.mata
```

To validate the upstream patch, copy the pinned repositories to scratch, then:

```bash
git -C reghdfe apply --check /path/to/reghdfe-vce-conley.patch
git -C reghdfe apply /path/to/reghdfe-vce-conley.patch
cd reghdfe/current-code && python3 build.py
```

Install the resulting reghdfe `src/` and the unmodified pinned ftools `src/`
into a scratch PLUS, add fastconley's `stata/src` to the adopath, set global
`UPSTREAM_PLUS`, and run `stata-mp -b do stata/upstream/test_upstream.do` from
the fastconley root. Stata batch can return shell status 0 after a do-file
error, so the log must contain `test_upstream.do: all checks passed`.

`test_upstream.do` primarily checks reghdfe integration because both commands
use the same Mata engine. It does not independently validate the distance
formula, compiled plugin, IV, raster/FFT dispatch, cross-platform numerics, or
large-data performance. Those belong to fastconley's R/Stata parity and engine
test suites. The upstream-native `test/part1/conley.do` has no fastconley
dependency and independently reconstructs the PSD floor from the raw VCE.

## Likely maintainer concerns

The largest review question is ownership: this proposal adds a sizeable
estimator-specific algorithm to reghdfe, including generated source and a new
public option/result surface. Other concerns are the pure-Mata runtime,
long-term synchronization with fastconley, the OLS-only boundary at ivreghdfe,
the statistical assumptions behind cutoff/kernel choices, and whether a
provider interface would be easier to maintain than native Conley code.

An alternative is a generic external-VCE hook, for example
`vce(external:<provider> ...)`. reghdfe would call the provider with `S`, `sol`,
`D`, `X`, `w`, and `vce_mode`; the provider would update `sol.V` and declare the
raw variables it needs through `vce_extra_keepvars` before partial-out. That
keeps compact/team-FE preservation generic while leaving Conley parsing,
documentation, algorithm ownership, IV/raster claims, and optional compiled
acceleration in fastconley. A production hook would also need a stable callback
contract, namespace/error rules, capability/version negotiation, posting and
replay conventions, and tests showing that external providers cannot corrupt
reghdfe state.
