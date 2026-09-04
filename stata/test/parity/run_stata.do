* Cross-language parity harness, step 2 (Stata side).
*   stata-mp -b do stata/test/parity/run_stata.do   (from the repository root;
*   optional args: OUT_DIR ENGINE)
* For every configuration in configs.csv: load the .dta, run fastconley with
* nossc nopsdfix keepsingletons tolerance(1e-14), write e(V) (slopes only)
* and e(b) as CSV for compare.R.
clear all
adopath ++ "stata/src"
args parity_dir parity_engine
if ("`parity_dir'" != "") global PARITY_DIR "`parity_dir'"
if ("`parity_engine'" != "") global PARITY_ENGINE "`parity_engine'"
local env_dir : environment FASTCONLEY_PARITY_DIR
local env_engine : environment FASTCONLEY_PARITY_ENGINE
if ("`env_dir'" != "") global PARITY_DIR "`env_dir'"
if ("`env_engine'" != "") global PARITY_ENGINE "`env_engine'"
if ("$PARITY_DIR" == "") global PARITY_DIR "stata/test/parity/out"
if ("$PARITY_ENGINE" == "") global PARITY_ENGINE "mata"

capture program drop _fc_write
program define _fc_write
	args matname fname
	mata: _fc_write_mat("`matname'", "`fname'")
end
mata:
void _fc_write_mat(string scalar matname, string scalar fname)
{
	real matrix M
	real scalar fh, i, j
	string scalar line
	M = st_matrix(matname)
	unlink(fname)
	fh = fopen(fname, "w")
	for (i = 1; i <= rows(M); i++) {
		line = ""
		for (j = 1; j <= cols(M); j++) line = line + sprintf("%s%21.17g", j > 1 ? "," : "", M[i, j])
		fput(fh, line)
	}
	fclose(fh)
}
end

import delimited using "$PARITY_DIR/configs.csv", clear varnames(1) stringcols(_all)
local ncfg = _N
forvalues i = 1/`ncfg' {
	local name = name[`i']
	local data = data[`i']
	local sopts = sopts[`i']
	local cfg`i' `"`name'|`data'|`sopts'"'
}
forvalues i = 1/`ncfg' {
	gettoken name rest : cfg`i', parse("|")
	gettoken bar rest : rest, parse("|")
	gettoken data rest : rest, parse("|")
	gettoken bar sopts : rest, parse("|")
	di as text _n "==== `name'"
	use "$PARITY_DIR/`data'.dta", clear
	* weights, if any, are written as a leading [aw=w] token in sopts
	local wexp
	if (strpos(`"`sopts'"', "[") == 1) {
		gettoken wexp sopts : sopts
	}
	* method(grid) is plugin-only; with the Mata engine the grid configs run pairwise
	if ("$PARITY_ENGINE" == "mata") local sopts : subinstr local sopts "method(grid)" "method(pairwise)"
	* IV configurations: x2 (cs) or xe (bal) instrumented
	local model y x1 x2
	if (strpos("`name'", "csiv") == 1) local model y x1 (x2 = z1 z2)
	if (strpos("`name'", "baliv") == 1) local model y x1 (xe = z)
	fastconley `model' `wexp', `sopts' nossc nopsdfix keepsingletons tolerance(1e-14) engine($PARITY_ENGINE)
	assert e(engine) == "$PARITY_ENGINE"
	matrix V = e(V)
	matrix V = V[1..2, 1..2]
	matrix b = e(b)
	matrix b = b[1, 1..2]
	_fc_write V "$PARITY_DIR/`name'_S.csv"
	_fc_write b "$PARITY_DIR/`name'_Sb.csv"
}
di as result _n "run_stata.do: `ncfg' configurations written (engine $PARITY_ENGINE)"
