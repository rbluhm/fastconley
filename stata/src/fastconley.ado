*! version 0.1.0 03sep2026  fastconley: Conley (1999) spatial HAC standard errors for reghdfe models
*! Richard Bluhm. Port of the fastconley R package (https://github.com/rbluhm/fastconley).

program define fastconley, eclass
	version 14.1

	* Replay
	if replay() {
		Replay `0'
		exit
	}

	* Estimate; clean up Mata objects whether or not it succeeds
	cap noi Estimate `0'
	loc rc = c(rc)
	cap mata: mata drop HDFE
	cap mata: mata drop fc_*
	cap mata: mata drop hdfe_*
	if (`rc') exit `rc'
end


program Estimate, eclass
	syntax anything(equalok) [if] [in] [fw aw pw/], ///
		LATitude(varname numeric) LONgitude(varname numeric) CUToff(real) ///
		[Absorb(string) NOAbsorb ///
		 KERnel(string) DISTance(string) LAG(real 0) UNIT(varname) TIME(varname) ///
		 BALanced PIXel(real 0) ENGine(string) METHod(string) THReads(integer 0) TILE(integer 1024) ///
		 NEIGHbor(string) CSRweight(string) ///
		 NOSSC NOPSDfix Verbose ///
		 RESiduals(name) noCONStant ///
		 noHEADer noTABLE noFOOTnote * ]

	* ---- model: "depvar indepvars" or IV "depvar exog (endog = instruments)" --
	loc iv = strpos(`"`anything'"', "(") > 0
	if (`iv') {
		gettoken depvar rest : anything
		loc p1 = strpos(`"`rest'"', "(")
		loc p2 = strpos(`"`rest'"', ")")
		if (`p1' == 0 | `p2' < `p1' | strpos(substr(`"`rest'"', `p2' + 1, .), "(")) {
			di as error "IV syntax is: depvar [exog] (endog = instruments)"
			exit 198
		}
		loc exog = strtrim(stritrim(substr(`"`rest'"', 1, `p1' - 1) + " " + substr(`"`rest'"', `p2' + 1, .)))
		loc inside = substr(`"`rest'"', `p1' + 1, `p2' - `p1' - 1)
		gettoken endog inst : inside, parse("=")
		gettoken eq inst : inst, parse("=")
		loc endog = strtrim("`endog'")
		loc inst = strtrim("`inst'")
		if ("`eq'" != "=" | "`endog'" == "" | "`inst'" == "") {
			di as error "IV syntax is: depvar [exog] (endog = instruments)"
			exit 198
		}
		foreach l in depvar exog endog inst {
			if (regexm(" ``l'' ", " (i|c|b|o|ib[0-9]+)\.") | strpos("``l''", "#")) {
				di as error "factor-variable terms are not supported in IV models; absorb them or use indicator variables"
				exit 198
			}
			if ("``l''" != "") {
				tsunab `l' : ``l''
			}
		}
		loc allvars `depvar' `exog' `endog' `inst'
	}
	else loc allvars `anything'
	loc 0 `"`allvars' `if' `in'"'
	if ("`weight'" != "") loc 0 `"`0' [`weight'=`exp']"'
	syntax varlist(fv ts numeric) [if] [in] [fw aw pw/]
	if (`iv') {
		gettoken depvar : varlist
	}

	* ---- dependencies ------------------------------------------------------
	cap which reghdfe
	if (c(rc)) {
		di as error "fastconley requires reghdfe (ssc install reghdfe)"
		exit 111
	}
	ms_get_version reghdfe, min_version("6.12.5")

	* ---- options -----------------------------------------------------------
	if ("`kernel'" == "") loc kernel bartlett
	if (!inlist("`kernel'", "bartlett", "uniform")) {
		di as error "kernel() must be bartlett or uniform"
		exit 198
	}
	if ("`distance'" == "") loc distance haversine
	if (!inlist("`distance'", "haversine", "spherical", "chord")) {
		di as error "distance() must be haversine, spherical, or chord"
		exit 198
	}
	if (`cutoff' == 0) {
		di as error "cutoff() must be positive (km); a negative cutoff gives the heteroskedasticity-only meat"
		exit 198
	}
	if (`lag' < 0) {
		di as error "lag() must be >= 0"
		exit 198
	}
	if (`lag' > 0 & ("`time'" == "" | "`unit'" == "")) {
		di as error "lag() requires unit() and time()"
		exit 198
	}
	if ("`balanced'" != "" & ("`time'" == "" | "`unit'" == "")) {
		di as error "balanced requires unit() and time()"
		exit 198
	}
	if (`pixel' < 0) {
		di as error "pixel() must be >= 0"
		exit 198
	}
	if ("`engine'" == "") loc engine auto
	if (!inlist("`engine'", "auto", "mata", "plugin")) {
		di as error "engine() must be auto, mata, or plugin"
		exit 198
	}
	if ("`method'" == "") loc method auto
	if (!inlist("`method'", "auto", "pairwise", "grid")) {
		di as error "method() must be auto, pairwise, or grid"
		exit 198
	}
	if ("`neighbor'" == "") loc neighbor grid
	if (!inlist("`neighbor'", "grid", "band")) {
		di as error "neighbor() must be grid or band"
		exit 198
	}
	if ("`csrweight'" == "") loc csrweight double
	if (!inlist("`csrweight'", "double", "float")) {
		di as error "csrweight() must be double or float"
		exit 198
	}
	if (`threads' <= 0) loc threads = c(processors)
	if (`tile' < 16) {
		di as error "tile() must be >= 16"
		exit 198
	}
	if ("`absorb'" == "" & "`noabsorb'" == "") {
		di as error "absorb() or noabsorb is required"
		exit 198
	}
	if ("`absorb'" != "" & "`noabsorb'" != "") {
		di as error "absorb() and noabsorb are mutually exclusive"
		exit 198
	}
	loc absopt = cond("`absorb'" != "", "absorb(`absorb')", "noabsorb")
	* Split display options from the reghdfe pass-through options
	_get_diopts diopts options, `options'
	loc ssc = ("`nossc'" == "")
	loc psdfix = ("`nopsdfix'" == "")
	loc verbose = ("`verbose'" != "")
	loc report_constant = ("`constant'" != "noconstant")
	loc balanced_flag = ("`balanced'" != "")
	if ("`weight'" != "") loc wexp "[`weight'=`exp']"

	* ---- engine ------------------------------------------------------------
	if ("`engine'" != "mata") {
		LoadPlugin
		if (r(ok)) {
			loc engine plugin
			loc plugin_build `r(build)'
		}
		else if ("`engine'" == "plugin") {
			di as error "engine(plugin): `r(why)'"
			exit 198
		}
		else {
			loc engine mata
			if (`verbose') di as text "note: compiled engine not available (`r(why)'); using the Mata engine"
		}
	}
	if ("`engine'" == "mata") {
		if ("`method'" == "grid") {
			di as error "method(grid) requires engine(plugin); the Mata engine has no raster engine"
			exit 198
		}
		if ("`neighbor'" != "grid" | "`csrweight'" != "double") {
			di as text "note: neighbor() and csrweight() apply to the plugin engine only"
		}
	}

	* ---- sample ------------------------------------------------------------
	marksample touse
	markout `touse' `latitude' `longitude'
	if ("`unit'`time'" != "") markout `touse' `unit' `time', strok
	qui count if `touse'
	if (r(N) == 0) error 2000

	* ---- 1. reghdfe: build the HDFE object and partial out, stop before regressing
	qui reghdfe `varlist' `wexp' if `touse', `absopt' noregress vce(unadjusted) `constant' `options'

	* The regression sample may have shrunk (singletons); rebuild touse from HDFE.sample
	qui replace `touse' = 0
	mata: st_store(HDFE.sample, "`touse'", J(rows(HDFE.sample), 1, 1))

	* ---- 2. keep the partialled-out data before the solver trims it
	mata: fc_Xstd = cols(HDFE.solution.data) > 1 ? HDFE.solution.data[|1, 2 \ rows(HDFE.solution.data), cols(HDFE.solution.data)|] : J(rows(HDFE.solution.data), 0, .)
	mata: fc_status0 = HDFE.solution.indepvar_status
	mata: fc_stdevs = HDFE.solution.is_standardized ? HDFE.solution.stdevs : J(1, cols(HDFE.solution.means), 1)
	mata: st_local("wt", HDFE.weight_type)
	mata: st_local("fc_N", strofreal("`wt'" == "fweight" ? quadsum(HDFE.weights) : rows(HDFE.sample), "%21.17g"))
	* Regression weights exactly as reghdfe_solve_ols builds them
	if ("`wt'" == "") mata: fc_w = J(0, 1, .)
	else if ("`wt'" == "fweight") mata: fc_w = HDFE.weights
	else mata: fc_w = HDFE.weights * (`fc_N' / quadsum(HDFE.weights))

	if (`iv') {
		* ---- 3'. 2SLS on the partialled-out data (no reghdfe solver involved)
		mata: fc_names = select(HDFE.solution.varlist, HDFE.solution.indepvar_status :== 0)
		mata: fc_data = HDFE.solution.data :* fc_stdevs
		mata: fastconley_iv_prepare(fc_data, fc_names, tokens("`exog'"), tokens("`endog'"), tokens("`inst'"), fc_w, `verbose')
		mata: st_local("kk", strofreal(fc_kk))
		mata: st_local("df_a", strofreal(HDFE.df_a))
		loc df_m = `kk'
		loc df_r = `fc_N' - `df_m' - `df_a'
		if (`df_r' <= 0) {
			di as error "no residual degrees of freedom"
			exit 2001
		}
		loc dof_adj = cond(`ssc', `fc_N' / `df_r', 1)
	}
	else {
		* ---- 3. reghdfe's own solver: b, residuals, df, collinearity handling, R2s
		mata: HDFE.solution.report_constant = HDFE.has_intercept & `report_constant'
		mata: reghdfe_solve_ols(HDFE, HDFE.solution, "vce_small")
		* the solver trims stdevs/means to the kept (non-collinear) regressors
		mata: fc_stdevs = HDFE.solution.is_standardized ? HDFE.solution.stdevs : J(1, cols(HDFE.solution.means), 1)
		mata: fc_tmpN = (HDFE.weight_type == "aweight" | HDFE.weight_type == "pweight") ? HDFE.solution.N : HDFE.solution.sumweights
	}

	* ---- 4. Conley VCE ----------------------------------------------------

	* Coordinates and panel identifiers on the regression sample
	mata: fc_lat = st_data(HDFE.sample, "`latitude'")
	mata: fc_lon = st_data(HDFE.sample, "`longitude'")
	if ("`time'" == "") mata: fc_time = J(rows(HDFE.sample), 1, 1)
	else {
		cap confirm string variable `time'
		if (!c(rc)) mata: fc_time = fastconley_group(st_sdata(HDFE.sample, "`time'"))
		else mata: fc_time = st_data(HDFE.sample, "`time'")
	}
	if ("`unit'" == "") mata: fc_unit = (1::rows(HDFE.sample))
	else {
		cap confirm string variable `unit'
		if (!c(rc)) mata: fc_unit = fastconley_group(st_sdata(HDFE.sample, "`unit'"))
		else mata: fc_unit = st_data(HDFE.sample, "`unit'")
	}

	if (`iv') {
		mata: fc_time = fc_time
		mata: fc_unit = fc_unit
		mata: fastconley_prepare_rows(fc_lat, fc_lon, fc_time, fc_unit, `balanced_flag', `pixel', `verbose')
	}
	else {
		mata: fastconley_prepare(fc_Xstd, fc_status0, HDFE.solution.indepvar_status, ///
			fc_stdevs, HDFE.solution.means, HDFE.solution.resid, fc_w, HDFE.solution.report_constant, ///
			fc_tmpN, fc_lat, fc_lon, fc_time, fc_unit, `balanced_flag', `pixel', `verbose')
		mata: st_local("kk", strofreal(fc_kk))
		mata: st_local("dof_adj", strofreal(`ssc' ? HDFE.solution.N / (HDFE.solution.N - HDFE.solution.df_m - HDFE.df_a) : 1, "%21.17g"))
	}
	mata: st_local("n_sp", strofreal(rows(fc_sp_S)))
	mata: st_local("sp_balanced", strofreal(fc_sp_balanced))
	mata: st_local("n_periods", strofreal(rows(uniqrows(fc_time))))
	loc method_used pairwise

	if ("`engine'" == "plugin") {
		* Hand the prepared rows to the compiled engine through temporary
		* variables (rows 1..n_sp; n_sp <= _N since aggregation only merges).
		loc use_grid 0
		if ("`method'" != "pairwise") {
			* float coordinates carry ~6e-8 relative rounding noise: loosen the lattice test
			loc grid_tol 1e-6
			loc lat_type : type `latitude'
			loc lon_type : type `longitude'
			if ("`lat_type'" == "float" | "`lon_type'" == "float") loc grid_tol 1e-3
			mata: fc_use_grid = fastconley_choose_grid("`method'", "`kernel'", ///
				fc_sp_lat, fc_sp_lon, fc_sp_time, fc_kk, `cutoff', `grid_tol')
			mata: st_local("use_grid", strofreal(fc_use_grid))
		}
		tempvar v_lat v_lon v_time
		loc svars
		forvalues j = 1/`kk' {
			tempvar s`j'
			loc svars `svars' `s`j''
		}
		foreach v in `v_lat' `v_lon' `v_time' `svars' {
			qui gen double `v' = .
		}
		mata: st_store((1::rows(fc_sp_S)), tokens("`v_lat' `v_lon' `v_time' `svars'"), ///
			(fc_sp_lat, fc_sp_lon, fc_sp_time, fc_sp_S))
		tempname M sc_cutoff sc_lag sc_lat0 sc_dlat sc_dlon
		scalar `sc_cutoff' = `cutoff'
		scalar `sc_lag' = `lag'
		matrix `M' = J(`kk', `kk', 0)
		loc done 0
		if (`use_grid') {
			tempvar v_ring v_col
			qui gen double `v_ring' = .
			qui gen double `v_col' = .
			mata: st_store((1::rows(fc_sp_S)), tokens("`v_ring' `v_col'"), (fc_ring, fc_col))
			mata: st_numscalar("`sc_lat0'", fc_grid[1])
			mata: st_numscalar("`sc_dlat'", fc_grid[2])
			mata: st_numscalar("`sc_dlon'", fc_grid[4])
			mata: st_local("gargs", "`sc_lat0' `sc_dlat' `sc_dlon' " + invtokens(strofreal(fc_grid[(5, 6, 7)], "%21.17g")))
			if (`verbose') di as text "# Conley spatial meat (plugin, grid engine, `kernel' kernel, `distance' distance, cutoff `cutoff' km, `threads' threads)"
			cap plugin call fastconley_plugin `v_ring' `v_col' `v_time' `svars' in 1/`n_sp', ///
				grid `gargs' `sc_cutoff' `distance' `kernel' `threads' `M'
			if (c(rc) == 0) {
				loc done 1
				loc method_used grid
			}
			else if ("`method'" == "auto" & strpos(`"`fc_plugin_error'"', "dateline")) {
				if (`verbose') di as text "   - lattice does not tile the dateline; falling back to the pairwise engine"
			}
			else {
				exit 198
			}
		}
		if (!`done') {
			if (`verbose') di as text "# Conley spatial meat (plugin, pairwise engine, `kernel' kernel, `distance' distance, cutoff `cutoff' km, `threads' threads)"
			plugin call fastconley_plugin `v_lat' `v_lon' `v_time' `svars' in 1/`n_sp', ///
				spatial `sc_cutoff' `kernel' `distance' `sp_balanced' `threads' `neighbor' `csrweight' `M'
			if ("`fc_unbalanced_fallback'" == "1") {
				di as text "note: balanced requested but period sizes differ after aggregation; general path used"
			}
		}
		mata: fc_meat = st_matrix("`M'")
		if (`lag' > 0 & `n_periods' > 1) {
			if (`verbose') di as text "# Serial HAC meat (plugin, lag cutoff `lag')"
			tempvar v_unit
			qui gen double `v_unit' = .
			mata: fc_p = order((fc_unit, fc_time), (1, 2))
			mata: st_store((1::rows(fc_S)), tokens("`v_unit' `v_time' `svars'"), ///
				(fc_unit[fc_p], fc_time[fc_p], fc_S[fc_p, .]))
			mata: st_local("n_full", strofreal(rows(fc_S)))
			matrix `M' = J(`kk', `kk', 0)
			plugin call fastconley_plugin `v_unit' `v_time' `svars' in 1/`n_full', serial `sc_lag' `threads' `M'
			mata: fc_meat = fc_meat + st_matrix("`M'")
		}
	}
	else {
		mata: fc_meat = fastconley_meat_mata(fc_sp_lat, fc_sp_lon, fc_sp_time, fc_sp_S, fc_sp_balanced, ///
			`cutoff', "`kernel'", "`distance'", `tile', `verbose')
		if (`lag' > 0 & `n_periods' > 1) {
			if (`verbose') di as text "# Serial HAC meat (Mata engine, lag cutoff `lag')"
			mata: fc_meat = fc_meat + fastconley_serial_meat(fc_unit, fc_time, fc_S, `lag')
		}
	}
	mata: fc_V = fastconley_assemble(fc_D, fc_meat, `dof_adj', `psdfix')
	if ("`fc_psd_noticeable'" == "1") {
		if (`psdfix') di as text "note: the Conley vcov was not positive semi-definite; negative eigenvalues were clamped (nopsdfix to disable)"
		else di as text "warning: the Conley vcov is not positive semi-definite (psd fix disabled)"
	}

	* ---- 5. post ------------------------------------------------------------
	tempname b V
	ereturn clear
	if (`iv') {
		* 2SLS: our own b/V and the statistics reghdfe_header displays; the
		* absorbed-FE footnote still comes from the HDFE object
		mata: st_local("xnames", invtokens(fc_xnames))
		mata: st_matrix("`b'", fc_b')
		mata: st_matrix("`V'", fc_V)
		matrix colnames `b' = `xnames'
		matrix colnames `V' = `xnames'
		matrix rownames `V' = `xnames'
		ereturn post `b' `V', esample(`touse') depname(`depvar')
		if ("`residuals'" != "") {
			mata: HDFE.save_variable("`residuals'", fc_resid, "Residuals")
			ereturn local resid "`residuals'"
		}
		ereturn scalar N = `fc_N'
		ereturn scalar df_m = `df_m'
		ereturn scalar df_r = `df_r'
		ereturn scalar rank = `df_m'
		mata: st_numscalar("e(rss)", fc_rss)
		mata: st_numscalar("e(tss)", HDFE.solution.tss[1])
		mata: st_numscalar("e(tss_within)", fc_tss_within)
		ereturn scalar mss = e(tss) - e(rss)
		ereturn scalar r2 = 1 - e(rss) / e(tss)
		ereturn scalar r2_within = 1 - e(rss) / e(tss_within)
		ereturn scalar rmse = sqrt(e(rss) / e(df_r))
		mata: st_numscalar("e(F)", fastconley_wald_F(fc_b, fc_V, fc_kk, `df_m'))
		ereturn local depvar "`depvar'"
		ereturn local indepvars "`xnames'"
		ereturn local instd "`endog'"
		ereturn local insts "`exog' `inst'"
		ereturn local exexog "`inst'"
		ereturn local inexog "`exog'"
		ereturn local title "HDFE 2SLS regression"
		ereturn local model "iv"
		ereturn local marginsnotok "Residuals SCore"
		if ("`wt'" != "") {
			ereturn local wtype "`wt'"
			ereturn local wexp "= `exp'"
		}
		mata: HDFE.post_footnote()
	}
	else {
		mata: HDFE.solution.V = fc_V
		mata: HDFE.solution.F = fastconley_wald_F(HDFE.solution.b, HDFE.solution.V, HDFE.solution.K, HDFE.solution.df_m)
		mata: HDFE.solution.expand_results("`b'", "`V'", HDFE.verbose)
		mata: st_local("depvar", HDFE.solution.depvar)
		mata: st_local("indepvars", invtokens(HDFE.solution.fullindepvars))
		if ("`indepvars'" != "") {
			matrix colnames `b' = `indepvars'
			matrix colnames `V' = `indepvars'
			matrix rownames `V' = `indepvars'
			_ms_findomitted `b' `V'
			ereturn post `b' `V', esample(`touse') buildfvinfo depname(`depvar')
		}
		else {
			ereturn post, esample(`touse') buildfvinfo depname(`depvar')
		}
		if ("`residuals'" != "") {
			mata: HDFE.save_variable("`residuals'", HDFE.solution.resid, "Residuals")
			mata: HDFE.solution.residuals_varname = "`residuals'"
		}
		mata: HDFE.solution.vcetype = "robust"   // Solution::post() accepts only its own list; overridden below
		mata: HDFE.solution.post()
		mata: HDFE.post_footnote()
	}

	ereturn local cmd "fastconley"
	ereturn local cmdline `"fastconley `0'"'
	ereturn local vcetype "Conley"
	ereturn local vce "conley"
	ereturn local title3 "Conley spatial HAC standard errors"
	ereturn local estat_cmd ""
	ereturn scalar iv = `iv'
	ereturn local conley_kernel "`kernel'"
	ereturn local conley_dist "`distance'"
	ereturn local engine "`engine'"
	ereturn local method "`method_used'"
	if ("`engine'" == "plugin") {
		ereturn local engine_build "`plugin_build'"
		ereturn scalar threads = `threads'
	}
	ereturn local latvar "`latitude'"
	ereturn local lonvar "`longitude'"
	if ("`unit'" != "") ereturn local unitvar "`unit'"
	if ("`time'" != "") ereturn local timevar "`time'"
	ereturn scalar conley_cutoff = `cutoff'
	ereturn scalar conley_lag = `lag'
	ereturn scalar conley_pixel = `pixel'
	ereturn scalar conley_balanced = `balanced_flag'
	ereturn scalar ssc = `ssc'
	ereturn scalar psd_fix = `psdfix'
	ereturn scalar dof_adj = `dof_adj'

	Replay, `header' `table' `footnote' `diopts'
end


* The compiled engine is loaded once, at the bottom of this file, when the
* ado is loaded (a plugin program cannot be re-defined from inside a
* program). r(ok) = 1 when it is available and its engine version matches.
program LoadPlugin, rclass
	if ("$FASTCONLEY_PLUGIN_FILE" == "") {
		return scalar ok = 0
		return local why "no plugin binary for this platform on the adopath"
		exit
	}
	cap plugin call fastconley_plugin, check
	if (c(rc)) {
		return scalar ok = 0
		return local why "$FASTCONLEY_PLUGIN_FILE failed the version check"
		exit
	}
	if ("$FASTCONLEY_ENGINE_VERSION" != "0.11.0") {
		return scalar ok = 0
		return local why "$FASTCONLEY_PLUGIN_FILE has engine version $FASTCONLEY_ENGINE_VERSION, this ado expects 0.11.0"
		exit
	}
	return scalar ok = 1
	return local file "$FASTCONLEY_PLUGIN_FILE"
	return local version "$FASTCONLEY_ENGINE_VERSION"
	return local build "$FASTCONLEY_ENGINE_BUILD"
end


program Replay, rclass
	syntax [, noHEADer noTABLE noFOOTnote *]
	if (`"`e(cmd)'"' != "fastconley") error 301
	_get_diopts options, `options'
	if ("`header'" == "") reghdfe_header
	if ("`header'" == "" & "`table'" == "") di ""
	if ("`table'" == "") _coef_table, `options'
	return add
	if ("`footnote'" == "") reghdfe_footnote
	if ("`footnote'" == "") {
		loc lagtxt = cond(e(conley_lag) > 0, ", serial lag cutoff " + string(e(conley_lag)), "")
		di as text "Conley SEs: " as res "`e(conley_kernel)'" as text " kernel, " ///
			as res "`e(conley_dist)'" as text " distance, cutoff " as res string(e(conley_cutoff)) as text " km`lagtxt'" ///
			as text " (`e(engine)' engine, `e(method)')"
	}
end

* Mata: reghdfe's classes must be compiled in this ado's context for the
* mata: one-liners above to see the HDFE object (same pattern as ivreghdfe).
cap findfile "reghdfe.mata"
if (_rc) {
	di as error "fastconley requires reghdfe (ssc install reghdfe)"
	exit 9
}
include "reghdfe.mata", adopath
include "fastconley.mata", adopath

* Compiled engine (optional): try the installed name first, then the
* platform-specific file shipped in the source tree. Loaded once per session;
* engine(auto) falls back to Mata when nothing loads. (Top-level ado code
* must stay flat: no nested braces.)
cap program drop fastconley_plugin
global FASTCONLEY_PLUGIN_FILE
local fc_platform_plugin fastconley_linux64.plugin
if (strpos("`c(os)'", "Windows")) local fc_platform_plugin fastconley_win64.plugin
if (inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac")) local fc_platform_plugin fastconley_macosx.plugin
cap program fastconley_plugin, plugin using("fastconley.plugin")
if (!c(rc)) global FASTCONLEY_PLUGIN_FILE "fastconley.plugin"
if ("$FASTCONLEY_PLUGIN_FILE" == "") cap program fastconley_plugin, plugin using("`fc_platform_plugin'")
if ("$FASTCONLEY_PLUGIN_FILE" == "" & !c(rc)) global FASTCONLEY_PLUGIN_FILE "`fc_platform_plugin'"
