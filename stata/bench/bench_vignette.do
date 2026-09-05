* bench_vignette.do -- Stata counterpart of the fastconley R vignette benchmarks
* (inst/benchmarks/run-fastconley-benchmarks.R and run-large-data-benchmarks.R).
*
* Mirrors the six sections of the R suite (dense baseline, scattered
* cross-sections, balanced panel with serial HAC, regular raster, repeated
* locations / pixel aggregation, one million observations) and times, for
* each configuration, the compiled plugin at several thread counts, the
* pure-Mata fallback, and -- where its O(n^2) loop is feasible -- acreg
* (Colella, Lalive, Sakalli, Thoenig). The quantity reported for fastconley
* is e(vce_seconds), the covariance step alone (what the R vignette times);
* total_seconds is the whole command including reghdfe's partial-out.
* acreg has no separate covariance timer, so its total time is reported.
*
* Controlled through globals (see run_bench.sh):
*   BENCH_SECTION    dense | xsection | panel | raster | pixel | large | overhead
*   BENCH_METHODS    subset of "plugin mata acreg" (default "plugin mata")
*   BENCH_THREADS    thread counts for the plugin (default "1 4 8 16")
*   BENCH_XS_CASES   xsection cases as n:cutoff pairs (default all four)
*   BENCH_LARGE_N    observations for the large section (default 1000000)
*   BENCH_OUT        CSV to append to (default stata/bench/results/bench.csv)
*   BENCH_REPS       repetitions, the minimum is kept (default 2)
* Run from the repository root: stata-mp -b do stata/bench/bench_vignette.do
version 14.1
clear all
set more off
set matsize 800
adopath ++ "stata/src"

if ("$BENCH_METHODS" == "") global BENCH_METHODS "plugin mata"
if ("$BENCH_THREADS" == "") global BENCH_THREADS "1 4 8 16"
if ("$BENCH_OUT" == "") global BENCH_OUT "stata/bench/results/bench.csv"
if ("$BENCH_REPS" == "") global BENCH_REPS 2
if ("$BENCH_LARGE_N" == "") global BENCH_LARGE_N 1000000
if ("$BENCH_XS_CASES" == "") global BENCH_XS_CASES "50000:100 50000:500 100000:100 100000:500"

* ---------------------------------------------------------------------------
* Bookkeeping: versions, CSV writer, timing helpers
* ---------------------------------------------------------------------------
fastconley, version
local ado_version = r(ado_version)
local engine_build = r(engine_build)
local engine_status = r(status)
local acreg_version ""
cap which acreg
if (!_rc) {
	cap findfile acreg.ado
	if (!_rc) {
		tempname fh
		file open `fh' using "`r(fn)'", read text
		file read `fh' line
		while (r(eof) == 0 & "`acreg_version'" == "") {
			if (strpos(`"`line'"', "*! Version") == 1) local acreg_version = trim(subinstr(`"`line'"', "*! Version", "", 1))
			file read `fh' line
		}
		file close `fh'
	}
}
local run_date = trim("`c(current_date)'") + " " + "`c(current_time)'"
local stata_version = "`c(stata_version)' `c(flavor)' (`c(processors)' licensed / `c(processors_mach)' cores)"

cap confirm file "$BENCH_OUT"
if (_rc) {
	tempname fh
	file open `fh' using "$BENCH_OUT", write text replace
	file write `fh' "section,benchmark,method,engine,threads,n_obs,n_unit,n_time,k,cutoff_km,kernel,dist,lag,pixel,vce_seconds,total_seconds,rel_diff_vs_plugin,notes,run_date,stata_version,ado_version,engine_build,acreg_version" _n
	file close `fh'
}

program define bench_row
	syntax, section(string) benchmark(string) method(string) engine(string) threads(string) ///
		n_obs(string) k(string) cutoff(string) kernel(string) dist(string) ///
		[n_unit(string) n_time(string) lag(string) pixel(string) vce(string) total(string) rel(string) notes(string) meta(string)]
	if ("`n_unit'" == "") local n_unit "NA"
	if ("`n_time'" == "") local n_time "NA"
	if ("`lag'" == "") local lag 0
	if ("`pixel'" == "") local pixel 0
	if ("`vce'" == "") local vce "NA"
	if ("`rel'" == "") local rel "NA"
	tempname fh
	file open `fh' using "$BENCH_OUT", write text append
	file write `fh' `"`section',`benchmark',`method',`engine',`threads',`n_obs',`n_unit',`n_time',`k',`cutoff',`kernel',`dist',`lag',`pixel',`vce',`total',`rel',"`notes'",`meta'"' _n
	file close `fh'
	di as text "  [`section'] `benchmark' | `method' engine=`engine' threads=`threads' | vce=`vce' s total=`total' s rel=`rel'"
end
local meta `""`run_date'","`stata_version'","`ado_version'","`engine_build'","`acreg_version'""'

* time_fastconley: runs `cmd' BENCH_REPS times, keeps the minimum e(vce_seconds)
* and total wall time, and leaves e(V) of the last run in memory.
program define time_fastconley, rclass
	syntax, cmd(string)
	local best_vce = .
	local best_tot = .
	forvalues r = 1/$BENCH_REPS {
		timer clear 90
		timer on 90
		qui `cmd'
		timer off 90
		qui timer list 90
		local tot = r(t90)
		local v = e(vce_seconds)
		if (`v' < `best_vce') local best_vce = `v'
		if (`tot' < `best_tot') local best_tot = `tot'
	}
	return scalar vce = `best_vce'
	return scalar total = `best_tot'
end

program define time_acreg, rclass
	syntax, cmd(string)
	timer clear 91
	timer on 91
	qui `cmd'
	timer off 91
	qui timer list 91
	return scalar total = r(t91)
end

* slope-block relative difference between two e(V)-style matrices
program define slope_reldiff, rclass
	args A B kk
	tempname A2 B2
	matrix `A2' = `A'[1..`kk', 1..`kk']
	matrix `B2' = `B'[1..`kk', 1..`kk']
	return scalar rel = mreldif(`A2', `B2')
end

local methods "$BENCH_METHODS"
local has_plugin = ("`engine_status'" == "ready")
if (strpos(" `methods' ", " plugin ") & !`has_plugin') {
	di as error "plugin engine not available (`engine_status'); plugin rows will be skipped"
}
local tmax = 1
foreach t of numlist $BENCH_THREADS {
	if (`t' > `tmax') local tmax = `t'
}

* ---------------------------------------------------------------------------
* Data generators (same shapes as the R scripts; Stata RNG, seeds per section)
* ---------------------------------------------------------------------------
program define make_xsection
	syntax, n(integer) k(integer) seed(integer) [lat0(real 25) lat1(real 49) lon0(real -125) lon1(real -67)]
	clear
	set seed `seed'
	set obs `n'
	gen double lat = runiform(`lat0', `lat1')
	gen double lon = runiform(`lon0', `lon1')
	gen double y = rnormal()
	forvalues j = 1/`k' {
		gen double x`j' = rnormal()
		qui replace y = y + (0.1 + 0.4 * (`j' - 1) / max(`k' - 1, 1)) * x`j'
	}
end

* ---------------------------------------------------------------------------
* 1. Dense baseline: n = 1000, 2000, 4000; k = 5; cutoff 500; uniform/spherical
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "dense") {
	local k = 5
	local cutoff = 500
	foreach n in 1000 2000 4000 {
		make_xsection, n(`n') k(`k') seed(`=1000 + `n'')
		local rhs ""
		forvalues j = 1/`k' {
			local rhs "`rhs' x`j'"
		}
		local base "y `rhs', noabsorb lat(lat) lon(lon) cutoff(`cutoff') kernel(uniform) dist(spherical) nossc nopsdfix"
		tempname Vp
		if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`tmax')")
			local vce = r(vce)
			local tot = r(total)
			matrix `Vp' = e(V)
			bench_row, section(dense) benchmark("dense uniform spherical") method(fastconley) engine(plugin) threads(`tmax') n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(0) notes("reference") meta(`meta')
		}
		if (strpos(" `methods' ", " mata ")) {
			time_fastconley, cmd("fastconley `base' engine(mata)")
			local vce = r(vce)
			local tot = r(total)
			local rel "NA"
			cap slope_reldiff e(V) `Vp' `k'
			if (!_rc) local rel = r(rel)
			bench_row, section(dense) benchmark("dense uniform spherical") method(fastconley) engine(mata) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(`rel') notes("") meta(`meta')
		}
		if (strpos(" `methods' ", " acreg ") & "`acreg_version'" != "") {
			if (!strpos(" `methods' ", " plugin ") & `has_plugin') {
				qui fastconley `base' engine(plugin) threads(`tmax')
				matrix `Vp' = e(V)
			}
			time_acreg, cmd("acreg y `rhs', spatial latitude(lat) longitude(lon) dist(`cutoff')")
			local tot = r(total)
			local rel "NA"
			cap slope_reldiff e(V) `Vp' `k'
			if (!_rc) local rel = r(rel)
			bench_row, section(dense) benchmark("dense uniform spherical") method(acreg) engine(acreg) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(planar) total(`tot') rel(`rel') notes("planar equirectangular distance; whole command timed") meta(`meta')
		}
	}
}

* ---------------------------------------------------------------------------
* 2. Scattered cross-sections: k = 10; uniform/spherical; cases n:cutoff
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "xsection") {
	local k = 10
	local ii = 0
	foreach case in $BENCH_XS_CASES {
		local ++ii
		local n = substr("`case'", 1, strpos("`case'", ":") - 1)
		local cutoff = substr("`case'", strpos("`case'", ":") + 1, .)
		make_xsection, n(`n') k(`k') seed(`=200 + `ii'')
		local rhs ""
		forvalues j = 1/`k' {
			local rhs "`rhs' x`j'"
		}
		local label "n=`n' cutoff=`cutoff'"
		local base "y `rhs', noabsorb lat(lat) lon(lon) cutoff(`cutoff') kernel(uniform) dist(spherical) nossc nopsdfix method(pairwise)"
		tempname Vp
		local have_ref = 0
		if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
			foreach t of numlist $BENCH_THREADS {
				time_fastconley, cmd("fastconley `base' engine(plugin) threads(`t')")
				local vce = r(vce)
				local tot = r(total)
				if (!`have_ref') {
					matrix `Vp' = e(V)
					local have_ref = 1
				}
				bench_row, section(cross_section) benchmark("`label'") method(fastconley) engine(plugin) threads(`t') n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(0) notes("") meta(`meta')
			}
		}
		if (strpos(" `methods' ", " mata ")) {
			time_fastconley, cmd("fastconley `base' engine(mata)")
			local vce = r(vce)
			local tot = r(total)
			local rel "NA"
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(cross_section) benchmark("`label'") method(fastconley) engine(mata) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(`rel') notes("") meta(`meta')
		}
		if (strpos(" `methods' ", " acreg ") & "`acreg_version'" != "") {
			if (!`have_ref' & `has_plugin') {
				qui fastconley `base' engine(plugin) threads(`tmax')
				matrix `Vp' = e(V)
				local have_ref = 1
			}
			time_acreg, cmd("acreg y `rhs', spatial latitude(lat) longitude(lon) dist(`cutoff')")
			local tot = r(total)
			local rel "NA"
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(cross_section) benchmark("`label'") method(acreg) engine(acreg) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(planar) total(`tot') rel(`rel') notes("planar equirectangular distance; whole command timed") meta(`meta')
		}
	}
}

* ---------------------------------------------------------------------------
* 3. Balanced panel with serial HAC: 10000 units x 4 periods, k = 5, 500 km, lag 1
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "panel") {
	local n_unit = 10000
	local n_time = 4
	local k = 5
	local cutoff = 500
	local lag = 1
	clear
	set seed 300
	set obs `n_unit'
	gen int unit = _n
	gen double lat = runiform(25, 49)
	gen double lon = runiform(-125, -67)
	gen double unit_fe = rnormal()
	expand `n_time'
	bysort unit: gen int time = _n
	sort time unit
	gen double y = unit_fe + rnormal()
	forvalues t = 1/`n_time' {
		local tfe = rnormal()
		qui replace y = y + `tfe' if time == `t'
	}
	local rhs ""
	forvalues j = 1/`k' {
		gen double x`j' = rnormal()
		qui replace y = y + (0.1 + 0.4 * (`j' - 1) / (`k' - 1)) * x`j'
		local rhs "`rhs' x`j'"
	}
	local n = _N
	local base "y `rhs', absorb(unit time) lat(lat) lon(lon) cutoff(`cutoff') kernel(uniform) dist(spherical) unit(unit) time(time) lag(`lag') balanced nossc nopsdfix"
	tempname Vp
	local have_ref = 0
	if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
		foreach t of numlist $BENCH_THREADS {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`t')")
			local vce = r(vce)
			local tot = r(total)
			if (!`have_ref') {
				matrix `Vp' = e(V)
				local have_ref = 1
			}
			bench_row, section(panel) benchmark("balanced panel SHAC") method(fastconley) engine(plugin) threads(`t') n_obs(`n') n_unit(`n_unit') n_time(`n_time') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) lag(`lag') vce(`vce') total(`tot') rel(0) notes("lag=1; spatial (contemporaneous) + own-unit serial Bartlett") meta(`meta')
		}
	}
	if (strpos(" `methods' ", " mata ")) {
		time_fastconley, cmd("fastconley `base' engine(mata)")
		local vce = r(vce)
		local tot = r(total)
		local rel "NA"
		if (`have_ref') {
			slope_reldiff e(V) `Vp' `k'
			local rel = r(rel)
		}
		bench_row, section(panel) benchmark("balanced panel SHAC") method(fastconley) engine(mata) threads(1) n_obs(`n') n_unit(`n_unit') n_time(`n_time') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) lag(`lag') vce(`vce') total(`tot') rel(`rel') notes("") meta(`meta')
	}
	if (strpos(" `methods' ", " acreg ") & "`acreg_version'" != "") {
		time_acreg, cmd("acreg y `rhs', spatial latitude(lat) longitude(lon) dist(`cutoff') id(unit) time(time) lag(`lag') hac pfe1(unit) pfe2(time)")
		local tot = r(total)
		local rel "NA"
		if (`have_ref') {
			slope_reldiff e(V) `Vp' `k'
			local rel = r(rel)
		}
		bench_row, section(panel) benchmark("balanced panel SHAC") method(acreg) engine(acreg) threads(1) n_obs(`n') n_unit(`n_unit') n_time(`n_time') k(`k') cutoff(`cutoff') kernel(uniform) dist(planar) lag(`lag') total(`tot') rel(`rel') notes("acreg hac: spatial x temporal product weights (different convention); planar distance; whole command timed") meta(`meta')
	}
}

* ---------------------------------------------------------------------------
* 4. Regular raster: 180 x 180 lattice (step 0.05 deg from 35N, 10W), k = 3, 250 km
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "raster") {
	local n_side = 180
	local k = 3
	local cutoff = 250
	clear
	set seed 400
	set obs `=`n_side' * `n_side''
	gen int ilat = mod(_n - 1, `n_side')
	gen int ilon = floor((_n - 1) / `n_side')
	gen double lat = 35 + 0.05 * ilat
	gen double lon = -10 + 0.05 * ilon
	gen double x1 = rnormal()
	gen double x2 = rnormal()
	gen double x3 = rnormal()
	gen double y = 0.3 * x1 - 0.2 * x2 + 0.1 * x3 + rnormal()
	local n = _N
	foreach kern in uniform bartlett {
		local base "y x1 x2 x3, noabsorb lat(lat) lon(lon) cutoff(`cutoff') kernel(`kern') dist(spherical) nossc nopsdfix"
		tempname Vp
		local have_ref = 0
		if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`tmax') method(grid)")
			local vce = r(vce)
			local tot = r(total)
			matrix `Vp' = e(V)
			local have_ref = 1
			bench_row, section(grid) benchmark("regular lat/lon lattice") method(grid) engine(plugin) threads(`tmax') n_obs(`n') k(`k') cutoff(`cutoff') kernel(`kern') dist(spherical) vce(`vce') total(`tot') rel(0) notes("180 x 180 lattice; grid engine") meta(`meta')
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`tmax') method(pairwise)")
			local vce = r(vce)
			local tot = r(total)
			slope_reldiff e(V) `Vp' `k'
			bench_row, section(grid) benchmark("regular lat/lon lattice") method(pairwise) engine(plugin) threads(`tmax') n_obs(`n') k(`k') cutoff(`cutoff') kernel(`kern') dist(spherical) vce(`vce') total(`tot') rel(`r(rel)') notes("3D cell-grid pair enumeration") meta(`meta')
		}
		if (strpos(" `methods' ", " mata ")) {
			time_fastconley, cmd("fastconley `base' engine(mata)")
			local vce = r(vce)
			local tot = r(total)
			local rel "NA"
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(grid) benchmark("regular lat/lon lattice") method(pairwise) engine(mata) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(`kern') dist(spherical) vce(`vce') total(`tot') rel(`rel') notes("Mata has no grid engine") meta(`meta')
		}
		if (strpos(" `methods' ", " acreg ") & "`acreg_version'" != "") {
			local bopt = cond("`kern'" == "bartlett", "bartlett", "")
			if (!`have_ref' & `has_plugin') {
				qui fastconley `base' engine(plugin) threads(`tmax')
				matrix `Vp' = e(V)
				local have_ref = 1
			}
			time_acreg, cmd("acreg y x1 x2 x3, spatial latitude(lat) longitude(lon) dist(`cutoff') `bopt'")
			local tot = r(total)
			local rel "NA"
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(grid) benchmark("regular lat/lon lattice") method(acreg) engine(acreg) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(`kern') dist(planar) total(`tot') rel(`rel') notes("planar equirectangular distance; whole command timed") meta(`meta')
		}
	}
}

* ---------------------------------------------------------------------------
* 5. Repeated locations: 20000 unique points x 5 rows, k = 5, 250 km, bartlett
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "pixel") {
	local n_unique = 20000
	local repeats = 5
	local k = 5
	local cutoff = 250
	clear
	set seed 500
	set obs `n_unique'
	gen double lat = runiform(25, 49)
	gen double lon = runiform(-125, -67)
	expand `repeats'
	gen double y = rnormal()
	local rhs ""
	forvalues j = 1/`k' {
		gen double x`j' = rnormal()
		qui replace y = y + (0.1 + 0.4 * (`j' - 1) / (`k' - 1)) * x`j'
		local rhs "`rhs' x`j'"
	}
	local n = _N
	tempname Vp
	local have_ref = 0
	foreach px in 0 10 25 {
		local base "y `rhs', noabsorb lat(lat) lon(lon) cutoff(`cutoff') kernel(bartlett) dist(spherical) pixel(`px') nossc nopsdfix"
		if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`tmax')")
			local vce = r(vce)
			local tot = r(total)
			if (`px' == 0) {
				matrix `Vp' = e(V)
				local have_ref = 1
			}
			local rel = 0
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(pixel) benchmark("repeated-coordinate cross-section") method(pixel=`px') engine(plugin) threads(`tmax') n_obs(`n') k(`k') cutoff(`cutoff') kernel(bartlett) dist(spherical) pixel(`px') vce(`vce') total(`tot') rel(`rel') notes("`=cond(`px' == 0, "exact coordinate dedupe", "approximate snapping")'; rel vs pixel=0") meta(`meta')
		}
		if (strpos(" `methods' ", " mata ")) {
			time_fastconley, cmd("fastconley `base' engine(mata)")
			local vce = r(vce)
			local tot = r(total)
			local rel "NA"
			if (`have_ref') {
				slope_reldiff e(V) `Vp' `k'
				local rel = r(rel)
			}
			bench_row, section(pixel) benchmark("repeated-coordinate cross-section") method(pixel=`px') engine(mata) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(bartlett) dist(spherical) pixel(`px') vce(`vce') total(`tot') rel(`rel') notes("rel vs plugin pixel=0") meta(`meta')
		}
	}
	if (strpos(" `methods' ", " acreg ") & "`acreg_version'" != "") {
		if (!`have_ref' & `has_plugin') {
			qui fastconley `base' engine(plugin) threads(`tmax')
			matrix `Vp' = e(V)
			local have_ref = 1
		}
		time_acreg, cmd("acreg y `rhs', spatial latitude(lat) longitude(lon) dist(`cutoff') bartlett")
		local tot = r(total)
		local rel "NA"
		if (`have_ref') {
			slope_reldiff e(V) `Vp' `k'
			local rel = r(rel)
		}
		bench_row, section(pixel) benchmark("repeated-coordinate cross-section") method(acreg) engine(acreg) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(bartlett) dist(planar) total(`tot') rel(`rel') notes("no aggregation; planar distance; whole command timed") meta(`meta')
	}
}

* ---------------------------------------------------------------------------
* 6. Large cross-section: BENCH_LARGE_N observations, k = 10, 100 km, uniform
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "large") {
	local n = $BENCH_LARGE_N
	local k = 10
	local cutoff = 100
	make_xsection, n(`n') k(`k') seed(600) lat0(-55) lat1(70) lon0(-180) lon1(180)
	local rhs ""
	forvalues j = 1/`k' {
		local rhs "`rhs' x`j'"
	}
	local base "y `rhs', noabsorb lat(lat) lon(lon) cutoff(`cutoff') kernel(uniform) dist(spherical) nossc nopsdfix method(pairwise)"
	tempname Vp
	local have_ref = 0
	if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
		foreach t of numlist $BENCH_THREADS {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`t')")
			local vce = r(vce)
			local tot = r(total)
			if (!`have_ref') {
				matrix `Vp' = e(V)
				local have_ref = 1
			}
			bench_row, section(large_xsection) benchmark("n=`n' cutoff=`cutoff'") method(fastconley) engine(plugin) threads(`t') n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(0) notes("lat in [-55,70], lon in [-180,180]") meta(`meta')
		}
	}
	if (strpos(" `methods' ", " mata ")) {
		time_fastconley, cmd("fastconley `base' engine(mata)")
		local vce = r(vce)
		local tot = r(total)
		local rel "NA"
		if (`have_ref') {
			slope_reldiff e(V) `Vp' `k'
			local rel = r(rel)
		}
		bench_row, section(large_xsection) benchmark("n=`n' cutoff=`cutoff'") method(fastconley) engine(mata) threads(1) n_obs(`n') k(`k') cutoff(`cutoff') kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(`rel') notes("") meta(`meta')
	}
}
* ---------------------------------------------------------------------------
* 7. Fixed overhead: cutoff(-1) keeps only the diagonal, so e(vce_seconds) is
*    the preparation cost (sorting, aggregation, marshalling) without any pairs
* ---------------------------------------------------------------------------
if ("$BENCH_SECTION" == "overhead") {
	local k = 10
	foreach n in 100000 $BENCH_LARGE_N {
		make_xsection, n(`n') k(`k') seed(`=700 + (`n' > 100000)') lat0(-55) lat1(70) lon0(-180) lon1(180)
		local rhs ""
		forvalues j = 1/`k' {
			local rhs "`rhs' x`j'"
		}
		local base "y `rhs', noabsorb lat(lat) lon(lon) cutoff(-1) kernel(uniform) dist(spherical) nossc nopsdfix method(pairwise)"
		if (`has_plugin' & strpos(" `methods' ", " plugin ")) {
			time_fastconley, cmd("fastconley `base' engine(plugin) threads(`tmax')")
			local vce = r(vce)
			local tot = r(total)
			bench_row, section(overhead) benchmark("n=`n' no pairs") method(fastconley) engine(plugin) threads(`tmax') n_obs(`n') k(`k') cutoff(-1) kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(0) notes("cutoff(-1): diagonal only; preparation cost") meta(`meta')
		}
		if (strpos(" `methods' ", " mata ")) {
			time_fastconley, cmd("fastconley `base' engine(mata)")
			local vce = r(vce)
			local tot = r(total)
			bench_row, section(overhead) benchmark("n=`n' no pairs") method(fastconley) engine(mata) threads(1) n_obs(`n') k(`k') cutoff(-1) kernel(uniform) dist(spherical) vce(`vce') total(`tot') rel(NA) notes("cutoff(-1): diagonal only; preparation cost") meta(`meta')
		}
	}
}
di as result "bench_vignette.do: section $BENCH_SECTION done"
exit
