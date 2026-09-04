# Primary proposal: a generic external-VCE hook for reghdfe

This directory contains two alternative, unsubmitted proposals against
reghdfe 6.14.1. The primary proposal is the small generic hook in
`reghdfe-external-vce-hook.patch`; fastconley is its first provider.

```stata
reghdfe y x1 x2, absorb(region) ///
    vce(external fastconley, lat(lat) lon(lon) cutoff(300))
```

The hook keeps estimator-specific algorithms out of reghdfe. reghdfe owns
sample construction, partialling out, OLS, compact-mode preservation, result
posting, replay, and postestimation. fastconley continues to own Conley syntax,
validation, small-sample and PSD rules, the pure-Mata engine, compiled plugin,
and its statistical documentation. A second provider can use the same hook
without another reghdfe parser or solver change.

## External-provider contract

The grammar is:

```stata
vce(external PROVIDER [, provider_options])
vce(ext PROVIDER [, provider_options])
```

Before ftools' `ms_parse_vce` runs, reghdfe recognizes the external clause,
verifies that `PROVIDER_reghdfe_vce` is on the adopath, and gives
`ms_parse_vce` a valid robust VCE. Missing or unavailable providers return
`r(198)` with the expected program name.

The provider is an s-class Stata program with two calls:

1. Before partial-out, reghdfe calls
   `PROVIDER_reghdfe_vce, keepvars provider_options`. The provider returns a
   simple varlist in `s(keepvars)` and may return the standard-error column
   title in `s(vcetype)`; the default title is `External`. reghdfe adds those
   variables to the estimation sample, compact/pool loading, and
   `ValidateGroups`. Thus provider variables must be constant within
   `group()` when team/individual FEs are used.
2. After `reghdfe_solve_ols`, reghdfe calls
   `PROVIDER_reghdfe_vce, compute provider_options` with `HDFE` in Mata scope.
   The provider must replace `HDFE.solution.V` with the full covariance in
   coefficient units, including `_cons` when reported and including its own
   small-sample adjustment. It may return whitespace-delimited `name=value`
   lists in `s(post_scalars)` and `s(post_macros)`; reghdfe posts those names
   to `e()`.

Both calls use `capture noisily`, so the provider's diagnostic remains visible.
A provider error is normalized to `r(498)` followed by
`vce provider PROVIDER failed`. After a successful compute, reghdfe recomputes
the model Wald F and posts `e(vce)=external`, the provider title in
`e(vcetype)`, `e(vce_provider)`, and the raw `e(vce_options)`. Replay,
`predict`, `test`, and `margins` remain reghdfe-native.

The hook is OLS-only. ivreghdfe rejects VCE types it does not know before this
path can run. A provider owns all numerical and statistical validation,
small-sample scaling, and covariance repair; reghdfe deliberately does not
interpret those choices.

## fastconley provider

`../src/fastconley_reghdfe_vce.ado` implements both calls. It accepts:

```text
lat() lon() cutoff() kernel() dist() lag() unit() time() balanced pixel()
engine() method() threads() tile() nossc nopsdfix verbose
```

The provider includes the standalone command's shared Mata routines and loads
the same installed plugin binary under a provider-private Stata program name.
On compute it reconstructs the same
weighted design, bread, scores, constant extension, DoF factor, spatial and
serial meat, engine dispatch, and PSD repair as standalone fastconley. The
provider therefore retains plugin acceleration and a pure-Mata fallback;
reghdfe contains no Conley-specific code.

## Files

- `reghdfe-external-vce-hook.patch`: primary patch, including help, version
  6.15.0 metadata, and `test/part1/external-vce.do`, whose tiny dummy provider
  makes the upstream certification test independent of fastconley.
- `test_hook.do`: end-to-end hook/provider tests run from this repository.
- `reghdfe-vce-conley.patch`: alternative full native patch. It adds
  `vce(conley lat lon, ...)`, the generated `Conley.mata`, help, posting, and
  native certification tests directly to reghdfe.
- `Conley.mata` and `make_conley_mata.py`: generated engine and generator used
  only by the full-patch alternative.
- `test_upstream.do`: integration suite for the full-patch alternative.

No ftools patch is needed by either proposal. Both pre-parse their new grammar
inside reghdfe and use unmodified ftools 2.50.0.

## Pinned sources and clean application

Both reghdfe patches target
`29b2b203f534413369d675d0c0f674d61291f536` (tag 6.14.1, also current master
when prepared). ftools is pinned at
`7b3663e49ea5c5b81638c55be29edf416e68e8b7` (2.50.0). The patches are
alternatives and must be applied to separate fresh copies, never stacked.

From the fastconley repository root:

```bash
mkdir -p .scratch
cp -a /path/to/pinned/reghdfe .scratch/reghdfe-hook
cp -a /path/to/pinned/reghdfe .scratch/reghdfe-full
cp -a /path/to/pinned/ftools .scratch/ftools

git -C .scratch/reghdfe-hook apply --check "$PWD/stata/upstream/reghdfe-external-vce-hook.patch"
git -C .scratch/reghdfe-hook apply "$PWD/stata/upstream/reghdfe-external-vce-hook.patch"
(cd .scratch/reghdfe-hook/current-code && python3 build.py)

git -C .scratch/reghdfe-full apply --check "$PWD/stata/upstream/reghdfe-vce-conley.patch"
git -C .scratch/reghdfe-full apply "$PWD/stata/upstream/reghdfe-vce-conley.patch"
(cd .scratch/reghdfe-full/current-code && python3 build.py)
```

Create two scratch PLUS directories. Copy the corresponding built reghdfe
`src/`, stock ftools `src/`, and `require.ado`/`require.sthlp` into each. Add
this repository's `stata/src` to the adopath. A scratch runner can then set:

```stata
global UPSTREAM_PLUS "/absolute/path/to/.scratch/plus-hook"
global FULL_PATCH_PLUS "/absolute/path/to/.scratch/plus-full"
do stata/upstream/test_hook.do
```

Run it from the repository root with `stata-mp -b do .scratch/run_hook.do`.
Stata batch mode can return shell status 0 after a do-file error, so inspect
the log and require the final line `test_hook.do: all checks passed`.

`test_hook.do` checks external fastconley against the standalone command at a
`1e-15` relative threshold for a cross section, aweights plus pixel
aggregation, a balanced panel with lag 2 and unit/time keys, and
uniform/chord/nossc/nopsdfix. It also checks every printed difference,
provider/raw-option posting, replay, `predict`, `test`, `margins`, compact,
`poolsize(1)`, compact plus `poolsize(1)`, missing-provider `r(198)`, visible
provider failures normalized to `r(498)`, dummy-provider dispatch, group
validation, and equality to the alternative native `vce(conley)` patch.

The repository's ordinary Stata gate remains:

```bash
stata-mp -b do stata/test/test_basic.do
```

The shipped plugins currently identify as engine 0.11.0 while the ado expects
0.11.1, so `engine(auto)` falling back to Mata in that test is expected.

## Alternative: native `vce(conley ...)`

The full patch offers a single-command, dependency-free Conley experience once
reghdfe is installed. Its option surface and stored results are visible in
reghdfe help, and its pure-Mata engine can be reviewed and tested in one
repository. It does not provide the compiled plugin or raster/FFT engine and
does not make Conley available to ivreghdfe.

Its trade-off is ownership: roughly nine hundred generated Mata lines plus
Conley parsing, validation, posting, documentation, and tests become reghdfe
maintenance. The generated copy must remain synchronized with fastconley, and
future algorithm changes require coordinated upstream updates. Validate it
separately with `test_upstream.do`; regenerate `Conley.mata` only for this
alternative via `python3 stata/upstream/make_conley_mata.py`.

## Maintainer review questions

For the generic hook, the central questions are whether the two-call ABI is
stable enough, whether retaining the trimmed standardized design for an
external run is an acceptable temporary memory cost, how strictly provider
`e()` names should be namespaced, and whether executing third-party ado code
inside reghdfe needs capability/version negotiation beyond `which` and clear
failure isolation. Its benefit is a small generic surface with no statistical
algorithm to maintain.

For the full patch, the central questions are the size and generated nature of
the added engine, pure-Mata performance, long-term synchronization, and
reghdfe assuming responsibility for Conley defaults and validation. Its
benefit is a self-contained native option with no provider contract.
