* Integration checks for reghdfe's generic external-VCE hook and fastconley.
* Run from the fastconley repository root with UPSTREAM_PLUS set to the
* hook-patched install and FULL_PATCH_PLUS to the alternative native patch.
version 14.1
clear all
set more off
if ("$UPSTREAM_PLUS" == "") {
	di as error "set global UPSTREAM_PLUS to the hook-patched scratch PLUS"
	exit 198
}
if ("$FULL_PATCH_PLUS" == "") {
	di as error "set global FULL_PATCH_PLUS to the full-patch scratch PLUS"
	exit 198
}
sysdir set PLUS "$UPSTREAM_PLUS"
adopath ++ "stata/src"
which ftools
which ms_parse_vce
which reghdfe
which fastconley_reghdfe_vce

program define dummy_reghdfe_vce, sclass
	syntax, [KEEPVARS COMPUTE] KEEP(varname numeric) [SCALE(real 2) FAIL]
	sreturn clear
	assert ("`keepvars'" != "") + ("`compute'" != "") == 1
	if ("`keepvars'" != "") {
		sreturn local keepvars "`keep'"
		sreturn local vcetype "Dummy"
		exit
	}
	if ("`fail'" != "") {
		di as error "dummy provider requested failure"
		exit 459
	}
	mata: dummy_reghdfe_compute(`scale', "`keep'")
	sreturn local post_scalars "dummy_scale=`scale'"
	sreturn local post_macros "dummy_status=ok"
end

* --------------------------------------------------------------------------
* Cross section: exact standalone parity, posting/replay/postestimation,
* compact/poolsize preservation, aweights/pixels, and uniform/chord/nossc.
* --------------------------------------------------------------------------
set seed 5
set obs 1200
gen double lat = runiform(25, 49)
gen double lon = runiform(-124, -67)
gen int region = ceil(runiform() * 5)
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double aw = runiform(.5, 2)
gen double y = .5 * x1 - .3 * x2 + rnormal()
save ".scratch/test_hook_cross.dta", replace

reghdfe y x1 x2, absorb(region) vce(external fastconley, lat(lat) lon(lon) cutoff(300))
matrix Vhook = e(V)
assert "`e(vce)'" == "external" & "`e(vcetype)'" == "Conley"
assert "`e(vce_provider)'" == "fastconley"
assert "`e(vce_options)'" == "lat(lat) lon(lon) cutoff(300)"
assert e(conley_cutoff) == 300 & e(conley_lag) == 0 & e(conley_pixel) == 0
assert e(conley_balanced) == 0 & e(ssc) == 1 & e(psd_fix) == 1
assert "`e(conley_kernel)'" == "bartlett" & "`e(conley_dist)'" == "haversine"
assert "`e(engine)'" == "mata" & "`e(method)'" == "pairwise"
preserve
clear
svmat double Vhook
save ".scratch/test_hook_Vhook.dta", replace
restore

reghdfe
predict double xb_hook, xb
assert !missing(xb_hook) if e(sample)
test x1 = x2
margins, dydx(x1)

fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
matrix Vstandalone = e(V)
scalar __diff = mreldif(Vhook, Vstandalone)
display as result "external fastconley vs standalone, cross section: " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15

foreach opts in "compact" "poolsize(1)" "compact poolsize(1)" {
	reghdfe y x1 x2, absorb(region) vce(external fastconley, lat(lat) lon(lon) cutoff(300) engine(mata)) `opts'
	matrix Vcurrent = e(V)
	scalar __diff = mreldif(Vhook, Vcurrent)
	display as result "external fastconley `opts' difference: " %21.17g scalar(__diff)
	assert scalar(__diff) < 1e-15
}

reghdfe y x1 x2 [aw=aw], absorb(region) vce(external fastconley, lat(lat) lon(lon) cutoff(300) pixel(20) engine(mata))
matrix Vhook_aw = e(V)
fastconley y x1 x2 [aw=aw], absorb(region) lat(lat) lon(lon) cutoff(300) pixel(20) engine(mata)
matrix Vstandalone = e(V)
scalar __diff = mreldif(Vhook_aw, Vstandalone)
display as result "external fastconley vs standalone, aweight + pixel: " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15

reghdfe y x1 x2, absorb(region) vce(external fastconley, lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(chord) nossc nopsdfix engine(mata) method(pairwise) threads(2) tile(64))
matrix Vhook_uniform = e(V)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) distance(chord) nossc nopsdfix engine(mata) method(pairwise) threads(2) tile(64)
matrix Vstandalone = e(V)
scalar __diff = mreldif(Vhook_uniform, Vstandalone)
display as result "external fastconley vs standalone, uniform/chord/nossc: " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15

* Dummy provider: independent hook check plus provider metadata and failure.
quietly include "reghdfe.mata", adopath
mata:
void dummy_reghdfe_compute(real scalar scale, string scalar keep)
{
	external class FixedEffects scalar HDFE
	assert(rows(st_data(HDFE.sample, keep)) == rows(HDFE.sample))
	HDFE.solution.V = HDFE.solution.V * scale
}
end
reghdfe y x1 x2, absorb(region) vce(robust)
matrix Vrobust = e(V)
matrix Vtwice_robust = 2 * Vrobust
reghdfe y x1 x2, absorb(region) vce(ext dummy, keep(lat) scale(2)) compact poolsize(1)
matrix Vdummy = e(V)
scalar __diff = mreldif(Vdummy, Vtwice_robust)
display as result "dummy provider vs twice robust: " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15
assert "`e(vce)'" == "external" & "`e(vcetype)'" == "Dummy"
assert e(dummy_scale) == 2 & "`e(dummy_status)'" == "ok"

capture noisily reghdfe y x1 x2, absorb(region) vce(external)
assert _rc == 198
capture noisily reghdfe y x1 x2, absorb(region) vce(external provider_that_does_not_exist)
assert _rc == 198
capture noisily reghdfe y x1 x2, absorb(region) vce(external fastconley, lat(lat) lon(lon) cutoff(.))
assert _rc == 498
capture noisily reghdfe y x1 x2, absorb(region) vce(external dummy, keep(lat) fail)
assert _rc == 498

* --------------------------------------------------------------------------
* Balanced panel with serial lag 2.
* --------------------------------------------------------------------------
clear
set seed 17
set obs 80
gen int unit = _n
gen double lat = runiform(35, 60)
gen double lon = runiform(-10, 30)
expand 5
bys unit: gen int time = _n
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double y = .5 * x1 - .3 * x2 + rnormal()
reghdfe y x1 x2, absorb(unit time) vce(external fastconley, lat(lat) lon(lon) cutoff(300) lag(2) unit(unit) time(time) balanced engine(mata))
matrix Vhook_panel = e(V)
fastconley y x1 x2, absorb(unit time) lat(lat) lon(lon) cutoff(300) lag(2) unit(unit) time(time) balanced engine(mata)
matrix Vstandalone = e(V)
scalar __diff = mreldif(Vhook_panel, Vstandalone)
display as result "external fastconley vs standalone, panel lag(2): " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15

* Team/group outcomes: provider keepvars must be constant within group.
clear
set obs 40
gen int team = ceil(_n / 2)
gen int person = _n
gen double lat = 35 + team / 100
gen double lon = -80 + team / 100
gen double x = team
gen double y = 2 * team + mod(team, 3)
reghdfe y x, noabsorb group(team) vce(external fastconley, lat(lat) lon(lon) cutoff(100) engine(mata))
replace lat = lat + 1 if mod(_n, 2) == 0
capture noisily reghdfe y x, absorb(person) group(team) individual(person) vce(external fastconley, lat(lat) lon(lon) cutoff(100) engine(mata))
assert _rc == 498
display as result "team-FE nonconstant provider keepvars rejected with r(498)"

* --------------------------------------------------------------------------
* The alternative full patch must produce the same native Conley covariance
* on the saved cross section. clear all forces reghdfe/Mata to reload.
* --------------------------------------------------------------------------
clear all
sysdir set PLUS "$FULL_PATCH_PLUS"
adopath ++ "stata/src"
use ".scratch/test_hook_cross.dta", clear
which reghdfe
reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300))
matrix Vnative = e(V)
preserve
use ".scratch/test_hook_Vhook.dta", clear
mkmat Vhook*, matrix(Vhook_saved)
restore
scalar __diff = mreldif(Vhook_saved, Vnative)
display as result "external fastconley vs alternative native vce(conley): " %21.17g scalar(__diff)
assert scalar(__diff) < 1e-15

display as result _n "test_hook.do: all checks passed"
