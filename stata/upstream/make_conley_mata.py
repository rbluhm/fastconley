#!/usr/bin/env python3
"""Regenerate the native reghdfe Conley source from fastconley's Mata engine.

The reghdfe-specific front end is deliberately embedded here: Conley.mata is
an output, never an input. Every top-level fastconley_* function found in
stata/src/fastconley.mata is copied in source order and namespaced with the
reghdfe_conley_ prefix.
"""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "stata" / "src" / "fastconley.mata"
OUTPUT = ROOT / "stata" / "upstream" / "Conley.mata"


FRONT_END = r'''// --------------------------------------------------------------------------
// Conley (1999) spatial HAC standard errors: vce(conley latvar lonvar, ...)
// --------------------------------------------------------------------------
// Arbitrary correlation between observations closer than a cutoff (km),
// with Bartlett (1 - d/cutoff) or uniform weights, and optionally
// within-unit serial correlation up to lag() periods (Bartlett weights in
// time). The meat is M + M' with M = sum_i s_i' c_i,
// c_i = 0.5 s_i + sum_{j > i, d_ij <= cutoff} w_ij s_j, on scores
// s_i = (w_i e_i) x_i. The pure-Mata engine below is generated from the
// fastconley Stata package and namespaced for reghdfe.
//
// References:
// - Conley (1999), "GMM estimation with cross sectional dependence",
//   Journal of Econometrics
// - Hsiang (2010), PNAS (the panel spatial + serial convention)

mata:

// Time keys: numeric variables are used as-is; string time values that all
// parse as numbers keep their numeric scale (so "2000" and "2002" stay two
// units apart for the serial HAC); other strings are equality-coded and only
// define spatial blocks, so they are rejected when lag() is positive. This
// mirrors fastconley and its R original.
real colvector reghdfe_conley_read_time(`FixedEffects' S, string scalar varname)
{
	string colvector raw

	if (!st_isstrvar(varname)) return(st_data(S.sample, varname))
	raw = st_sdata(S.sample, varname)
	if (reghdfe_conley_numeric_string(raw)) return(strtoreal(raw))
	if (S.conley_lag > 0) {
		_error(198, "lag() requires numeric time; nonnumeric string time values only define spatial blocks")
	}
	return(reghdfe_conley_group_codes(raw))
}

// Unit keys always keep string identity ("01" and "1" are distinct units).
real colvector reghdfe_conley_read_unit(`FixedEffects' S, string scalar varname)
{
	if (!st_isstrvar(varname)) return(st_data(S.sample, varname))
	return(reghdfe_conley_group_codes(st_sdata(S.sample, varname)))
}


void reghdfe_conley_validate_balanced(real colvector lat,
	                                  real colvector lon,
	                                  real colvector time,
	                                  real colvector unit)
{
	real colvector p, sorted_time, sorted_unit, sorted_lat, sorted_lon, first_unit
	real matrix info
	real scalar n, T, n_per, t

	n = rows(lat)
	p = order((time, unit), (1, 2))
	sorted_time = time[p]
	sorted_unit = unit[p]
	sorted_lat = lat[p]
	sorted_lon = lon[p]
	T = rows(uniqrows(sorted_time))
	if (T <= 1) return
	info = panelsetup(sorted_time, 1)
	n_per = n / T
	if (n_per != floor(n_per) | any(info[., 2] - info[., 1] :+ 1 :!= n_per))
		_error(3498, "vce(conley): balanced requires equally sized periods")
	first_unit = sorted_unit[|1 \ n_per|]
	if (rows(uniqrows(first_unit)) != n_per)
		_error(3498, "vce(conley): balanced requires unique units within period")
	for (t = 2; t <= T; t++) {
		if (any(sorted_unit[|(t-1)*n_per+1 \ t*n_per|] :!= first_unit))
			_error(3498, "vce(conley): balanced requires the same units in every period")
		if (any(sorted_lat[|(t-1)*n_per+1 \ t*n_per|] :!= sorted_lat[|1 \ n_per|]) |
		    any(sorted_lon[|(t-1)*n_per+1 \ t*n_per|] :!= sorted_lon[|1 \ n_per|]))
			_error(3498, "vce(conley): balanced requires time-invariant coordinates")
	}
}


`Void' reghdfe_vce_conley(`FixedEffects' S,
                          `Solution' sol,
                          `Matrix' D,
                          `Matrix' X,
                          `Variable' w,
                          `String' vce_mode)
{
	`Matrix'                Scores, M, Sagg, V_user, scale
	`Vector'                resid, lat, lon, time, unit, alat, alon, atime
	`RowVector'             stdev_x
	`Real'                  dof_adj, stdev_y
	`Boolean'               noticeable
	`Integer'               n

	assert_msg(S.conley_lat != "" & S.conley_lon != "", "vce(conley) requires latitude and longitude variables")
	assert_msg(!missing(S.conley_cutoff) & S.conley_cutoff > 0, "vce(conley) requires cutoff(#) > 0")

	// Scores exactly as in reghdfe_vce_dkraay: residuals times normalized weights.
	resid = S.weight_type != "" ? sol.resid :* w : sol.resid
	Scores = (sol.report_constant ? (X, J(rows(X), 1, 1)) : X) :* resid
	n = rows(Scores)

	lat = st_data(S.sample, S.conley_lat)
	lon = st_data(S.sample, S.conley_lon)
	time = S.conley_time == "" ? J(n, 1, 1) : reghdfe_conley_read_time(S, S.conley_time)
	unit = S.conley_unit == "" ? (1::n) : reghdfe_conley_read_unit(S, S.conley_unit)
	if (S.conley_balanced) reghdfe_conley_validate_balanced(lat, lon, time, unit)

	if (S.verbose > 0) {
		printf("{txt}# Estimating Conley spatial HAC Variance-Covariance Matrix\n\n")
		printf("{txt}   - Kernel: {res}%s{txt}; distance: {res}%s{txt}; cutoff: {res}%g{txt} km; lag: {res}%g{txt}\n",
		       S.conley_kernel, S.conley_dist, S.conley_cutoff, S.conley_lag)
	}

	// Aggregate identical (time, lat, lon) points, optionally after pixel snap.
	alat = lat; alon = lon; atime = time; Sagg = Scores
	reghdfe_conley_aggregate(alat, alon, atime, Sagg, S.conley_pixel)
	M = reghdfe_conley_spatial_meat(alat, alon, atime, Sagg, 1, S.conley_cutoff,
	                                S.conley_kernel, S.conley_dist, 512, S.verbose > 0)
	if (S.conley_lag > 0 & rows(uniqrows(time)) > 1) {
		M = M + reghdfe_conley_serial_meat(unit, time, Scores, S.conley_lag)
	}

	// Same small-sample factor as vce(robust); nossc disables it.
	dof_adj = sol.N / (sol.N - S.df_a - sol.df_m)
	if (vce_mode == "vce_asymptotic") dof_adj = sol.N / (sol.N - 1)
	if (!S.conley_ssc) dof_adj = 1
	if (S.verbose > 0 & vce_mode != "vce_asymptotic") {
		printf("{txt}   - Small-sample-adjustment: q = %g\n", dof_adj)
	}

	sol.V = D * M * D * dof_adj
	_makesymmetric(sol.V)
	// reghdfe calls VCE providers before undoing standardization. Clamp in
	// user units, as standalone fastconley and the R package do, then convert
	// back so reghdfe's generic post-VCE scaling produces V_user exactly.
	if (sol.is_standardized) {
		stdev_y = sol.stdevs[1]
		stdev_x = sol.K ? sol.stdevs[|2 \ cols(sol.stdevs)|] : J(1, 0, .)
		if (sol.report_constant) stdev_x = stdev_x, 1
		stdev_x = stdev_x :/ stdev_y
		scale = stdev_x' * stdev_x
		V_user = sol.V :/ scale
		noticeable = reghdfe_conley_psd_fix(V_user, S.conley_psdfix)
		sol.V = V_user :* scale
	}
	else noticeable = reghdfe_conley_psd_fix(sol.V, S.conley_psdfix)
	if (noticeable) {
		if (S.conley_psdfix) printf("{txt}note: the Conley vcov was not positive semi-definite; negative eigenvalues were clamped (nopsdfix to disable)\n")
		else printf("{txt}warning: the Conley vcov is not positive semi-definite (PSD fix disabled)\n")
	}

	sol.conley_lat = S.conley_lat
	sol.conley_lon = S.conley_lon
	sol.conley_unit = S.conley_unit
	sol.conley_time = S.conley_time
	sol.conley_cutoff = S.conley_cutoff
	sol.conley_lag = S.conley_lag
	sol.conley_pixel = S.conley_pixel
	sol.conley_kernel = S.conley_kernel
	sol.conley_dist = S.conley_dist
	sol.conley_ssc = S.conley_ssc
	sol.conley_psdfix = S.conley_psdfix
	sol.conley_balanced = S.conley_balanced
	sol.conley_engine = "mata"
	if (S.verbose > 0) printf("\n")
}

'''


FUNCTION_RE = re.compile(
    r"(?m)^[A-Za-z][A-Za-z0-9_ ]*\b(fastconley_[A-Za-z0-9_]+)\s*\("
)


def closing_brace(text: str, opening: int) -> int:
    """Return one past the brace matching text[opening], ignoring comments/strings."""
    depth = 0
    i = opening
    in_string = False
    in_line_comment = False
    in_block_comment = False
    while i < len(text):
        pair = text[i : i + 2]
        char = text[i]
        if in_line_comment:
            if char == "\n":
                in_line_comment = False
        elif in_block_comment:
            if pair == "*/":
                in_block_comment = False
                i += 1
        elif in_string:
            if char == '"':
                in_string = False
        elif pair == "//":
            in_line_comment = True
            i += 1
        elif pair == "/*":
            in_block_comment = True
            i += 1
        elif char == '"':
            in_string = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise ValueError("unterminated Mata function")


def extract_engine(source: str) -> tuple[str, list[str]]:
    matches = list(FUNCTION_RE.finditer(source))
    if not matches:
        raise ValueError(f"no fastconley_* functions found in {SOURCE}")
    names = [match.group(1) for match in matches]
    if len(names) != len(set(names)):
        raise ValueError("duplicate fastconley_* function definition")

    # Copy one contiguous region so comments and helper declarations between
    # functions survive. The endpoints themselves are discovered, not hard-coded.
    first = matches[0].start()
    separator = source.rfind("// ---------------------------------------------------------------------------", 0, first)
    if separator >= 0:
        first = separator
    last_open = source.find("{", matches[-1].end())
    if last_open < 0:
        raise ValueError(f"missing body for {names[-1]}")
    last = closing_brace(source, last_open)
    engine = source[first:last].strip() + "\n"

    # Every function-like fastconley identifier in the extracted engine must be
    # defined there; this catches new helpers whose declaration shape we missed.
    referenced = set(re.findall(r"\b(fastconley_[A-Za-z0-9_]+)\s*\(", engine))
    missing = sorted(referenced.difference(names))
    if missing:
        raise ValueError("unrecognized fastconley_* function definitions: " + ", ".join(missing))

    engine = re.sub(r"\bfastconley_", "reghdfe_conley_", engine)
    engine = engine.replace('"fastconley:', '"reghdfe vce(conley):')
    engine = engine.replace("use engine(plugin)", "use a larger cutoff")
    engine = engine.replace('"method(grid): ', '"reghdfe vce(conley): ')
    engine = engine.replace("; use method(pairwise)", "")
    engine = engine.replace("Mata engine", "native Mata code")
    return engine, names


def main() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    engine, names = extract_engine(source)
    output = (
        FRONT_END
        + "// ---------------------------------------------------------------------------\n"
        + "// Engine derived from stata/src/fastconley.mata by make_conley_mata.py.\n"
        + "// Do not edit this generated file; edit the source engine or generator.\n"
        + "// ---------------------------------------------------------------------------\n\n"
        + engine
        + "\nend\n"
    )
    OUTPUT.write_text(output, encoding="utf-8")
    print(f"Conley.mata regenerated: {output.count(chr(10))} lines; {len(names)} functions")


if __name__ == "__main__":
    main()
