* Validates the vce(conley) proposal for reghdfe against the fastconley
* Stata command (same engine, so results must agree to summation order).
* Run from the fastconley repository root with a PLUS directory holding the
* patched reghdfe + ftools:  global UPSTREAM_PLUS "/path/to/plus"
clear all
if ("$UPSTREAM_PLUS" != "") sysdir set PLUS "$UPSTREAM_PLUS"
adopath ++ "stata/src"
which reghdfe
set seed 5

* cross-section, one absorbed FE
set obs 3000
gen lat = runiform(25, 49)
gen lon = runiform(-124, -67)
gen region = ceil(runiform() * 5)
gen x1 = rnormal()
gen x2 = rnormal()
gen w = runiform(0.5, 2)
gen y = 0.5 * x1 - 0.3 * x2 + rnormal()

reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300))
assert "`e(vce)'" == "conley" & "`e(vcetype)'" == "Conley"
assert e(conley_cutoff) == 300 & "`e(conley_kernel)'" == "bartlett"
matrix Vr = e(V)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
di "reghdfe vce(conley) vs fastconley: " mreldif(Vr, e(V))
assert mreldif(Vr, e(V)) < 1e-12

reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) kernel(uniform) dist(chord) nossc) verbose(1)
matrix Vr = e(V)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(chord) nossc engine(mata)
assert mreldif(Vr, e(V)) < 1e-12

reghdfe y x1 x2 [aw = w], absorb(region) vce(conley lat lon, cutoff(300) pixel(20))
matrix Vr = e(V)
fastconley y x1 x2 [aw = w], absorb(region) lat(lat) lon(lon) cutoff(300) pixel(20) engine(mata)
di "aweight + pixel: " mreldif(Vr, e(V))
assert mreldif(Vr, e(V)) < 1e-12

reghdfe y x1 x2, noabsorb vce(conley lat lon, cutoff(300))
matrix Vr = e(V)
fastconley y x1 x2, noabsorb lat(lat) lon(lon) cutoff(300) engine(mata)
assert mreldif(Vr, e(V)) < 1e-12

* errors
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat, cutoff(300))" == 9
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon)" == 198
rcof "reghdfe y x1 x2, absorb(region) vce(conley lat lon, cutoff(300) lag(2))" == 9

* panel with serial lags
clear
set obs 400
gen unit = _n
gen ulat = runiform(35, 60)
gen ulon = runiform(-10, 30)
expand 5
bys unit: gen time = _n
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 0.5 * x1 - 0.3 * x2 + rnormal()
reghdfe y x1 x2, absorb(unit time) vce(conley ulat ulon, cutoff(300) lag(2) unit(unit) time(time))
matrix Vr = e(V)
assert e(conley_lag) == 2
reghdfe                                   // replay
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) engine(mata)
di "panel lag(2): " mreldif(Vr, e(V))
assert mreldif(Vr, e(V)) < 1e-12
di as result _n "test_upstream.do: all checks passed"
