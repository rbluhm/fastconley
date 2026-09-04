* Integration checks for native reghdfe vce(conley) against the standalone
* fastconley command. Run from the fastconley repository root after setting
* global UPSTREAM_PLUS to a scratch PLUS containing patched reghdfe and stock
* ftools 2.50.0.
version 14.1
clear all
set more off
if ("$UPSTREAM_PLUS" != "") sysdir set PLUS "$UPSTREAM_PLUS"
adopath ++ "stata/src"
which ftools
which ms_parse_vce
which reghdfe
set seed 5

* --------------------------------------------------------------------------
* Cross section: parity, posting, abbreviations, compact, and weights.
* --------------------------------------------------------------------------
set obs 3000
gen double lat = runiform(25, 49)
gen double lon = runiform(-124, -67)
gen int region = ceil(runiform() * 5)
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double w = runiform(0.5, 2)
gen int fw = ceil(3 * w)
gen double y = 0.5 * x1 - 0.3 * x2 + rnormal()

reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300))
assert "`e(vce)'" == "conley" & "`e(vcetype)'" == "Conley"
assert "`e(conley_lat)'" == "lat" & "`e(conley_lon)'" == "lon"
assert "`e(conley_unit)'" == "" & "`e(conley_time)'" == ""
assert e(conley_cutoff) == 300 & e(conley_lag) == 0 & e(conley_pixel) == 0
assert "`e(conley_kernel)'" == "bartlett" & "`e(conley_distance)'" == "haversine"
assert "`e(conley_dist)'" == "haversine" & "`e(conley_engine)'" == "mata"
assert e(conley_ssc) == 1 & e(conley_psdfix) == 1 & e(conley_balanced) == 0
matrix Vbase = e(V)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vbase, Vcurrent)
display as result "reghdfe vce(conley) vs fastconley: " scalar(__diff)
assert scalar(__diff) < 1e-12

reghdfe y x1 x2, absorb(region) vce(con lat lon, cutoff(300) kernel(UNIFORM) distance(CHORD) nossc nopsdfix) verbose(1)
assert "`e(conley_kernel)'" == "uniform" & "`e(conley_distance)'" == "chord"
assert e(conley_ssc) == 0 & e(conley_psdfix) == 0
matrix Vcon = e(V)
reghdfe y x1 x2, absorb(region) vce(conl lat lon, cutoff(300) kernel(uniform) distance(chord) nossc nopsdfix)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vcon, Vcurrent)
display as result "con/conl abbreviation difference: " scalar(__diff)
assert scalar(__diff) < 1e-14

reghdfe y x1 x2 [aw = w], absorb(region) vce(conley lat lon, cutoff(300) pixel(20))
matrix Vaw = e(V)
fastconley y x1 x2 [aw = w], absorb(region) lat(lat) lon(lon) cutoff(300) pixel(20) engine(mata)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vaw, Vcurrent)
display as result "aweight + pixel: " scalar(__diff)
assert scalar(__diff) < 1e-12

reghdfe y x1 x2 [fw = fw], absorb(region) vce(conley lat lon, cutoff(300))
matrix Vfw = e(V)
fastconley y x1 x2 [fw = fw], absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vfw, Vcurrent)
display as result "fweight: " scalar(__diff)
assert scalar(__diff) < 1e-12

* This specifically verifies that local pre-parsing reaches Conley with stock
* ftools; standalone fastconley's pweight preparatory call is tested elsewhere.
reghdfe y x1 x2 [pw = w], absorb(region) vce(conley lat lon, cutoff(300))
mata: assert(!missing(st_matrix("e(V)")))
display as result "pweight native Conley: finite e(V)"

reghdfe y x1 x2, noabsorb vce(conley lat lon, cutoff(300))
matrix Vnoabsorb = e(V)
fastconley y x1 x2, noabsorb lat(lat) lon(lon) cutoff(300) engine(mata)
matrix Vcurrent = e(V)
assert mreldif(Vnoabsorb, Vcurrent) < 1e-12

reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300)) compact
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vbase, Vcurrent)
display as result "compact difference: " scalar(__diff)
assert scalar(__diff) < 1e-12
reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300)) compact poolsize(1)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vbase, Vcurrent)
display as result "compact poolsize(1) difference: " scalar(__diff)
assert scalar(__diff) < 1e-12
reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300)) poolsize(1)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vbase, Vcurrent)
display as result "poolsize(1) difference: " scalar(__diff)
assert scalar(__diff) < 1e-12

* Every Conley user-syntax failure is normalized to r(198).
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat, cutoff(300))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon)" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(.))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) kernel(bad))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) lag(2))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) lag(.) unit(region) time(region))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) pixel(.))" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) balanced)" == 198

* --------------------------------------------------------------------------
* Panel lags, balanced validation, compact loading, replay, and stale s().
* --------------------------------------------------------------------------
clear
set seed 17
set obs 80
gen int unit = _n
gen double ulat = runiform(35, 60)
gen double ulon = runiform(-10, 30)
expand 5
bys unit: gen int time = _n
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double y = 0.5 * x1 - 0.3 * x2 + rnormal()

reghdfe y x1 x2, absorb(unit time) vce(conley ulat ulon, cutoff(300) lag(2) unit(unit) time(time) balanced)
matrix Vpanel = e(V)
assert e(conley_lag) == 2 & e(conley_balanced) == 1
assert "`e(conley_unit)'" == "unit" & "`e(conley_time)'" == "time"
reghdfe
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced engine(mata)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vpanel, Vcurrent)
display as result "panel lag(2): " scalar(__diff)
assert scalar(__diff) < 1e-12

reghdfe y x1 x2, absorb(unit time) vce(conley ulat ulon, cutoff(300))
assert e(conley_lag) == 0 & e(conley_pixel) == 0 & e(conley_balanced) == 0
assert "`e(conley_unit)'" == "" & "`e(conley_time)'" == ""
assert "`e(conley_kernel)'" == "bartlett" & "`e(conley_distance)'" == "haversine"
display as result "option-sparse call cleared prior panel configuration"

reghdfe y x1 x2, absorb(unit time) vce(conley ulat ulon, cutoff(300) lag(2) unit(unit) time(time)) compact
matrix Vcurrent = e(V)
assert mreldif(Vpanel, Vcurrent) < 1e-12
reghdfe y x1 x2, absorb(unit time) vce(conley ulat ulon, cutoff(300) lag(2) unit(unit) time(time)) compact poolsize(1)
matrix Vcurrent = e(V)
assert mreldif(Vpanel, Vcurrent) < 1e-12

* --------------------------------------------------------------------------
* Default PSD repair must equal standalone fastconley's 1e-16 clamp.
* --------------------------------------------------------------------------
clear
set seed 1
set obs 20
gen double lat = runiform(40, 42)
gen double lon = runiform(-2, 2)
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double x3 = rnormal()
gen double y = rnormal()
reghdfe y x1 x2 x3, noabsorb vce(conley lat lon, cutoff(150) kernel(uniform) nossc nopsdfix)
mata: symeigensystem(st_matrix("e(V)"), __psd_X=., __psd_lambda=.); st_numscalar("raw_min_eigen", min(__psd_lambda))
assert scalar(raw_min_eigen) < -1e-8
reghdfe y x1 x2 x3, noabsorb vce(conley lat lon, cutoff(150) kernel(uniform) nossc)
matrix Vpsd = e(V)
fastconley y x1 x2 x3, noabsorb lat(lat) lon(lon) cutoff(150) kernel(uniform) nossc engine(mata)
matrix Vcurrent = e(V)
scalar __diff = mreldif(Vpsd, Vcurrent)
display as result "PSD-clamped difference: " scalar(__diff)
assert scalar(__diff) < 1e-12
mata: symeigensystem(st_matrix("Vpsd"), __psd_X=., __psd_lambda=.); assert(min(__psd_lambda) >= 0)

* --------------------------------------------------------------------------
* group()/individual() reduction may not choose arbitrary coordinates.
* --------------------------------------------------------------------------
clear
set obs 40
gen int team = ceil(_n / 2)
gen int person = _n
gen double lat = 35 + team / 100
gen double lon = -80 + team / 100
gen double x = team
gen double y = 2 * team + mod(team, 3)
reghdfe y x, noabsorb group(team) vce(conley lat lon, cutoff(100))
replace lat = lat + 1 if mod(_n, 2) == 0
capture noisily reghdfe y x, noabsorb group(team) individual(person) absorb(person) vce(conley lat lon, cutoff(100))
assert _rc == 498
display as result "team-FE nonconstant coordinates rejected with r(498)"


* --------------------------------------------------------------------------
* Key handling: string units keep identity ("01" != "1"); numeric-looking
* string time keeps its scale; nonnumeric string time is rejected with lag().
* --------------------------------------------------------------------------
clear
set obs 400
set seed 17
gen int unum = ceil(_n / 20)
gen double lat = runiform(30, 45)
gen double lon = runiform(-110, -80)
bysort unum (lat): gen int t = _n
gen double x = rnormal()
gen double y = 0.5 * x + rnormal()
gen str2 ustr = string(unum, "%02.0f")
replace ustr = "1" if unum == 1 & t > 10
gen int unum2 = unum
replace unum2 = 99 if unum == 1 & t > 10
reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(unum2) time(t) lag(2))
matrix __Vnum = e(V)
reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(ustr) time(t) lag(2))
matrix __Vstr = e(V)
scalar __diff = mreldif(__Vstr, __Vnum)
display as result "string unit identity (01 vs 1 distinct): " scalar(__diff)
assert scalar(__diff) < 1e-12
gen str4 tstr = string(2000 + 2 * t)
gen int tnum = 2000 + 2 * t
reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(unum) time(tnum) lag(2))
matrix __Vt1 = e(V)
reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(unum) time(tstr) lag(2))
matrix __Vt2 = e(V)
scalar __diff = mreldif(__Vt1, __Vt2)
display as result "numeric-looking string time keeps scale: " scalar(__diff)
assert scalar(__diff) < 1e-12
gen str5 tbad = "T" + string(t)
rcof `"reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(unum) time(tbad) lag(2))"' == 198
reghdfe y x, noabsorb vce(conley lat lon, cutoff(300) unit(unum) time(tbad))
display as result "nonnumeric string time rejected with lag(), accepted without"

display as result _n "test_upstream.do: all checks passed"
