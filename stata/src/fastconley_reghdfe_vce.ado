*! version 0.2.0 04sep2026  external-VCE provider for reghdfe 6.15.0+

program define fastconley_reghdfe_vce, sclass
	version 14.1
	syntax, [KEEPVARS COMPUTE] ///
		LATitude(varname numeric) LONgitude(varname numeric) CUToff(real) ///
		[ KERnel(string) DISTance(string) LAG(real 0) UNIT(varname) TIME(varname) ///
		  BALanced PIXel(real 0) ENGine(string) METHod(string) ///
		  THReads(real 0) TILE(real 1024) NOSSC NOPSDfix Verbose ]

	if (("`keepvars'" != "") + ("`compute'" != "") != 1) {
		di as error "fastconley reghdfe provider requires exactly one of keepvars or compute"
		exit 198
	}
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
	if (missing(`cutoff')) {
		di as error "cutoff() may not be missing"
		exit 198
	}
	if (`cutoff' == 0) {
		di as error "cutoff() must be positive (km); a negative cutoff gives the heteroskedasticity-only meat"
		exit 198
	}
	if (missing(`lag') | `lag' < 0 | `lag' != floor(`lag')) {
		di as error "lag() must be a nonnegative integer"
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
	if (missing(`pixel') | `pixel' < 0) {
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
	if (missing(`threads') | `threads' != floor(`threads')) {
		di as error "threads() must be an integer"
		exit 198
	}
	if (`threads' <= 0) loc threads = c(processors_mach)
	if (missing(`tile') | `tile' != floor(`tile')) {
		di as error "tile() must be an integer"
		exit 198
	}
	if (`tile' < 16 | `tile' > 8192) {
		di as error "tile() must be between 16 and 8192"
		exit 198
	}
	if (5 * `tile'^2 * 8 > 1073741824) {
		di as error "tile() would require more than 1 GiB of estimated Mata workspace"
		exit 198
	}

	loc balanced_flag = ("`balanced'" != "")
	loc ssc = ("`nossc'" == "")
	loc psdfix = ("`nopsdfix'" == "")
	loc verbose_flag = ("`verbose'" != "")

	if ("`keepvars'" != "") {
		sreturn clear
		sreturn local keepvars "`latitude' `longitude' `unit' `time'"
		sreturn local vcetype "Conley"
		exit
	}

	sreturn clear
	loc engine_request `engine'
	if ("`engine'" != "mata") {
		FastconleyRHLoad
		if (r(ok)) {
			loc engine plugin
			loc plugin_build `r(build)'
		}
		else if ("`engine'" == "plugin") {
			di as error "engine(plugin): `r(why)'"
			FastconleyRHDiag
			exit 198
		}
		else {
			loc engine mata
			if (`verbose_flag') {
				di as text "note: compiled engine not available (`r(why)'); using the Mata engine"
				FastconleyRHDiag
			}
		}
	}
	if ("`engine'" == "mata" & "`method'" == "grid") {
		di as error "method(grid) requires engine(plugin); the Mata engine has no raster engine"
		exit 198
	}

	* Reconstruct exactly the standalone command's post-solve bread, scores,
	* sample keys, and prepared spatial rows from the external-hook HDFE object.
	mata: fc_Xstd = cols(HDFE.solution.data) > 1 ? HDFE.solution.data[|1, 2 \ rows(HDFE.solution.data), cols(HDFE.solution.data)|] : J(rows(HDFE.solution.data), 0, .)
	mata: fc_status = J(1, cols(HDFE.solution.data), 0)
	mata: fc_stdevs = HDFE.solution.is_standardized ? HDFE.solution.stdevs : J(1, cols(HDFE.solution.means), 1)
	mata: fc_tmpN = (HDFE.weight_type == "aweight" | HDFE.weight_type == "pweight") ? HDFE.solution.N : HDFE.solution.sumweights
	mata: fc_w = HDFE.weight_type == "" ? J(0, 1, .) : (HDFE.weight_type == "fweight" ? HDFE.weights : HDFE.weights * (HDFE.solution.N / quadsum(HDFE.weights)))
	mata: fc_lat = st_data(HDFE.sample, "`latitude'")
	mata: fc_lon = st_data(HDFE.sample, "`longitude'")
	if ("`time'" == "") mata: fc_time = J(rows(HDFE.sample), 1, 1)
	else {
		cap confirm string variable `time'
		if (!c(rc)) {
			mata: st_local("fc_numeric_time", strofreal(fastconley_numeric_string(st_sdata(HDFE.sample, "`time'"))))
			if (`lag' > 0 & !`fc_numeric_time') {
				di as error "lag() requires numeric time; nonnumeric string time values only define spatial blocks"
				exit 198
			}
			mata: fc_time = fastconley_group(st_sdata(HDFE.sample, "`time'"))
		}
		else mata: fc_time = st_data(HDFE.sample, "`time'")
	}
	if ("`unit'" == "") mata: fc_unit = (1::rows(HDFE.sample))
	else {
		cap confirm string variable `unit'
		if (!c(rc)) mata: fc_unit = fastconley_group_codes(st_sdata(HDFE.sample, "`unit'"))
		else mata: fc_unit = st_data(HDFE.sample, "`unit'")
	}
	mata: fastconley_prepare(fc_Xstd, fc_status, fc_status, fc_stdevs, HDFE.solution.means, ///
		HDFE.solution.resid, fc_w, HDFE.solution.report_constant, fc_tmpN, fc_lat, fc_lon, ///
		fc_time, fc_unit, `balanced_flag', `pixel', `verbose_flag')
	mata: st_local("kk", strofreal(fc_kk))
	mata: st_local("dof_adj", strofreal(`ssc' ? HDFE.solution.N / (HDFE.solution.N - HDFE.solution.df_m - HDFE.df_a) : 1, "%21.17g"))
	mata: st_local("n_sp", strofreal(rows(fc_sp_S)))
	mata: st_local("n_full", strofreal(rows(fc_S)))
	mata: st_local("sp_balanced", strofreal(fc_sp_balanced))
	mata: st_local("n_periods", strofreal(rows(uniqrows(fc_time))))

	loc method_used pairwise
	if ("`engine'" == "plugin" & (`n_sp' > 2147483647 | (`lag' > 0 & `n_full' > 2147483647))) {
		if ("`engine_request'" == "plugin") {
			di as error "engine(plugin) supports at most 2,147,483,647 prepared rows; use engine(mata)"
			exit 198
		}
		loc engine mata
		if ("`method'" == "grid") loc method pairwise
		if (`verbose_flag') di as text "note: prepared sample exceeds the plugin row-index limit; using the Mata pairwise engine"
	}

	if ("`engine'" == "plugin") {
		loc use_grid 0
		if ("`method'" != "pairwise") {
			loc grid_tol 1e-6
			loc lat_type : type `latitude'
			loc lon_type : type `longitude'
			if ("`lat_type'" == "float" | "`lon_type'" == "float") loc grid_tol 1e-3
			mata: fc_use_grid = fastconley_choose_grid("`method'", "`kernel'", fc_sp_lat, fc_sp_lon, fc_sp_time, fc_kk, `cutoff', `grid_tol')
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
		mata: st_store((1::rows(fc_sp_S)), tokens("`v_lat' `v_lon' `v_time' `svars'"), (fc_sp_lat, fc_sp_lon, fc_sp_time, fc_sp_S))
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
			if (`verbose_flag') di as text "# Conley spatial meat (plugin, grid engine, `kernel' kernel, `distance' distance, cutoff `cutoff' km, `threads' threads)"
			loc fc_plugin_error
			cap plugin call fastconley_rh_plugin `v_ring' `v_col' `v_time' `svars' in 1/`n_sp', ///
				grid `gargs' `sc_cutoff' `distance' `kernel' `threads' `M'
			loc grid_rc = c(rc)
			if (`grid_rc' == 0) {
				loc done 1
				loc method_used grid
			}
			else if ("`method'" == "auto" & strpos(`"`fc_plugin_error'"', "dateline")) {
				if (`verbose_flag') di as text "   - lattice does not tile the dateline; falling back to the pairwise engine"
			}
			else {
				if (`"`fc_plugin_error'"' != "") di as error "fastconley plugin: `fc_plugin_error'"
				else di as error "fastconley plugin grid call failed with return code `grid_rc'"
				exit 198
			}
		}
		if (!`done') {
			if (`verbose_flag') di as text "# Conley spatial meat (plugin, pairwise engine, `kernel' kernel, `distance' distance, cutoff `cutoff' km, `threads' threads)"
			plugin call fastconley_rh_plugin `v_lat' `v_lon' `v_time' `svars' in 1/`n_sp', ///
				spatial `sc_cutoff' `kernel' `distance' `sp_balanced' `threads' grid double `M'
			if ("`fc_unbalanced_fallback'" == "1") di as text "note: balanced requested but period sizes differ after aggregation; general path used"
		}
		mata: fc_meat = st_matrix("`M'")
		if (`lag' > 0 & `n_periods' > 1) {
			if (`verbose_flag') di as text "# Serial HAC meat (plugin, lag cutoff `lag')"
			tempvar v_unit
			qui gen double `v_unit' = .
			mata: fc_p = order((fc_unit, fc_time), (1, 2))
			mata: st_store((1::rows(fc_S)), tokens("`v_unit' `v_time' `svars'"), (fc_unit[fc_p], fc_time[fc_p], fc_S[fc_p, .]))
			matrix `M' = J(`kk', `kk', 0)
			plugin call fastconley_rh_plugin `v_unit' `v_time' `svars' in 1/`n_full', serial `sc_lag' `threads' `M'
			mata: fc_meat = fc_meat + st_matrix("`M'")
		}
	}
	else {
		mata: fc_meat = fastconley_meat_mata(fc_sp_lat, fc_sp_lon, fc_sp_time, fc_sp_S, fc_sp_balanced, ///
			`cutoff', "`kernel'", "`distance'", `tile', `verbose_flag')
		if (`lag' > 0 & `n_periods' > 1) {
			if (`verbose_flag') di as text "# Serial HAC meat (Mata engine, lag cutoff `lag')"
			mata: fc_meat = fc_meat + fastconley_serial_meat(fc_unit, fc_time, fc_S, `lag')
		}
	}
	mata: HDFE.solution.V = fastconley_assemble(fc_D, fc_meat, `dof_adj', `psdfix')
	if ("`fc_psd_noticeable'" == "1") {
		if (`psdfix') di as text "note: the Conley vcov was not positive semi-definite; negative eigenvalues were clamped (nopsdfix to disable)"
		else di as text "warning: the Conley vcov is not positive semi-definite (psd fix disabled)"
	}

	loc post_scalars "conley_cutoff=`cutoff' conley_lag=`lag' conley_pixel=`pixel' conley_balanced=`balanced_flag' ssc=`ssc' psd_fix=`psdfix' dof_adj=`dof_adj'"
	loc post_macros "conley_kernel=`kernel' conley_dist=`distance' engine=`engine' method=`method_used' latvar=`latitude' lonvar=`longitude'"
	if ("`unit'" != "") loc post_macros "`post_macros' unitvar=`unit'"
	if ("`time'" != "") loc post_macros "`post_macros' timevar=`time'"
	if ("`engine'" == "plugin") {
		loc post_scalars "`post_scalars' threads=`threads'"
		loc post_macros "`post_macros' engine_build=`plugin_build'"
	}
	sreturn local post_scalars "`post_scalars'"
	sreturn local post_macros "`post_macros'"
end


program FastconleyRHLoad, rclass
	loc expected_engine_version "0.11.1"
	return local expected "`expected_engine_version'"
	if ("$FASTCONLEY_RH_PLUGIN_FILE" == "") {
		return scalar ok = 0
		return local why "no plugin binary for this platform on the adopath"
		exit
	}
	global FASTCONLEY_ENGINE_VERSION
	global FASTCONLEY_ENGINE_BUILD
	cap plugin call fastconley_rh_plugin, check
	loc check_rc = c(rc)
	if (`check_rc') {
		return scalar ok = 0
		return local why "$FASTCONLEY_RH_PLUGIN_FILE failed the version check (return code `check_rc')"
		exit
	}
	if ("$FASTCONLEY_ENGINE_VERSION" != "`expected_engine_version'") {
		return scalar ok = 0
		return local why "$FASTCONLEY_RH_PLUGIN_FILE has engine version $FASTCONLEY_ENGINE_VERSION, this ado expects `expected_engine_version'"
		exit
	}
	return scalar ok = 1
	return local file "$FASTCONLEY_RH_PLUGIN_FILE"
	return local version "$FASTCONLEY_ENGINE_VERSION"
	return local build "$FASTCONLEY_ENGINE_BUILD"
end


program FastconleyRHDiag
	di as text "plugin loader attempts: $FASTCONLEY_RH_PLUGIN_TRIED"
	di as text "plugin loader return codes: $FASTCONLEY_RH_PLUGIN_RCS"
end


* Compile the same class and engine sources used by the standalone command.
cap findfile "reghdfe.mata"
if (_rc) {
	di as error "fastconley provider requires reghdfe 6.15.0+"
	exit 9
}
include "reghdfe.mata", adopath
include "fastconley.mata", adopath

* A plugin program is private to its ado file, so the provider needs its own
* loader name even though it uses the same installed binary as fastconley.
cap program drop fastconley_rh_plugin
global FASTCONLEY_RH_PLUGIN_FILE
global FASTCONLEY_RH_PLUGIN_TRIED
global FASTCONLEY_RH_PLUGIN_RCS
local fc_rh_platform_plugin fastconley_linux64.plugin
if ("`c(os)'" == "Windows") local fc_rh_platform_plugin fastconley_win64.plugin
if ("`c(os)'" == "MacOSX") local fc_rh_platform_plugin fastconley_macosx.plugin
cap program fastconley_rh_plugin, plugin using("`fc_rh_platform_plugin'")
local fc_rh_platform_rc = c(rc)
global FASTCONLEY_RH_PLUGIN_TRIED "`fc_rh_platform_plugin'"
global FASTCONLEY_RH_PLUGIN_RCS "`fc_rh_platform_rc'"
if (`fc_rh_platform_rc' == 0) global FASTCONLEY_RH_PLUGIN_FILE "`fc_rh_platform_plugin'"
local fc_rh_try_generic = ("$FASTCONLEY_RH_PLUGIN_FILE" == "")
if (`fc_rh_try_generic') cap program fastconley_rh_plugin, plugin using("fastconley.plugin")
local fc_rh_generic_rc = cond(`fc_rh_try_generic', c(rc), .)
if (`fc_rh_try_generic') global FASTCONLEY_RH_PLUGIN_TRIED "$FASTCONLEY_RH_PLUGIN_TRIED fastconley.plugin"
if (`fc_rh_try_generic') global FASTCONLEY_RH_PLUGIN_RCS "$FASTCONLEY_RH_PLUGIN_RCS `fc_rh_generic_rc'"
if (`fc_rh_try_generic' & `fc_rh_generic_rc' == 0) global FASTCONLEY_RH_PLUGIN_FILE "fastconley.plugin"
