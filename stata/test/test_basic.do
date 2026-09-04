* Smoke + consistency tests for the Stata fastconley command (Mata engine).
* Run from the repository root:  stata-mp -b do stata/test/test_basic.do
clear all
adopath ++ "stata/src"
set seed 20260903

* ---- cross-section with one absorbed FE ------------------------------------
set obs 3000
gen lat = runiform(25, 49)
gen lon = runiform(-124, -67)
gen region = ceil(runiform() * 5)
gen x1 = rnormal()
gen x2 = rnormal()
gen w = runiform(0.5, 2)
gen y = 0.5 * x1 - 0.3 * x2 + rnormal()

fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300)
assert e(cmd) == "fastconley"
assert e(vcetype) == "Conley"
assert e(N) == 3000
assert e(conley_cutoff) == 300
matrix V1 = e(V)
local default_engine = e(engine)
di "default engine: `default_engine'"
fastconley                                  // replay

* ---- engines: plugin (if present) and Mata must agree ---------------------
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
assert e(engine) == "mata"
matrix Vm = e(V)
cap fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(plugin) verbose
local have_plugin = (_rc == 0)
if (`have_plugin') {
	assert e(engine) == "plugin"
	assert e(method) == "pairwise"
	matrix Vp = e(V)
	di "plugin vs mata: " mreldif(Vp, Vm)
	assert mreldif(Vp, Vm) < 1e-12
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(plugin) neighbor(band) threads(2)
	assert mreldif(e(V), Vm) < 1e-12
	assert e(threads) == 2
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(plugin) kernel(uniform) dist(chord)
	matrix Vpu = e(V)
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata) kernel(uniform) dist(chord)
	assert mreldif(e(V), Vpu) < 1e-12
}
else di as text "plugin not available on this platform; skipping plugin checks"
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata) method(grid)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(nope)" == 198

* uniform kernel, spherical distance, verbose
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(spherical) verbose

* ---- cutoff(-1) = heteroskedasticity-only meat: must equal reghdfe vce(robust)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
matrix Vc = e(V)
scalar Fc = e(F)
reghdfe y x1 x2, absorb(region) vce(robust)
matrix Vr = e(V)
di "mreldif(V conley cutoff(-1), V reghdfe robust) = " mreldif(Vc, Vr)
assert mreldif(Vc, Vr) < 1e-10
assert reldif(Fc, e(F)) < 1e-10

* same with noabsorb (reported _cons goes through the extended bread)
fastconley y x1 x2, noabsorb lat(lat) lon(lon) cutoff(-1) nopsdfix
matrix Vc = e(V)
reghdfe y x1 x2, noabsorb vce(robust)
matrix Vr = e(V)
di "noabsorb: mreldif = " mreldif(Vc, Vr)
assert mreldif(Vc, Vr) < 1e-10

* same with aweights
fastconley y x1 x2 [aw = w], absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
matrix Vc = e(V)
reghdfe y x1 x2 [aw = w], absorb(region) vce(robust)
matrix Vr = e(V)
di "aweight: mreldif = " mreldif(Vc, Vr)
assert mreldif(Vc, Vr) < 1e-10

* nossc: dof_adj must be 1 and V scaled by (N - df_m - df_a)/N
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) nossc
assert e(dof_adj) == 1
matrix V2 = e(V) * (e(N) / (e(N) - e(df_m) - e(df_a)))
di "ssc scale check: " mreldif(V1, V2)
assert mreldif(V1, V2) < 1e-12

* pixel aggregation runs and stays close (25 km snap on a 300 km cutoff)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) pixel(25)
matrix Vp = e(V)
di "pixel(25) vs exact: " mreldif(Vp, V1)
assert mreldif(Vp, V1) < 0.05

* factor variables, collinear regressor, residuals()
gen x3 = 2 * x1
fastconley y x1 x2 x3 i.region, absorb(region) lat(lat) lon(lon) cutoff(300) resid(res)
assert _b[x3] == 0
assert !missing(_se[x1]) & _se[x1] > 0
qui sum res
assert abs(r(mean)) < 1e-6

* string lat/lon rejected, missing option rejected
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(gaussian)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) lag(2)" == 198
rcof "fastconley y x1 x2, lat(lat) lon(lon) cutoff(300)" == 198

* ---- raster lattice: grid engine (plugin) must equal pairwise -------------
if (`have_plugin') {
	clear
	set obs 4800
	gen ring = mod(_n - 1, 60)
	gen col = floor((_n - 1) / 60)
	gen lat = 35 + 0.5 * ring
	gen lon = -20 + 0.5 * col
	drop if runiform() < 0.25
	gen region = ceil(runiform() * 3)
	gen x1 = rnormal()
	gen x2 = rnormal()
	gen y = 0.5 * x1 - 0.3 * x2 + rnormal()
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(pairwise) engine(plugin)
	matrix Vpw = e(V)
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(grid) engine(plugin) verbose
	assert e(method) == "grid"
	di "grid vs pairwise (bartlett): " mreldif(e(V), Vpw)
	assert mreldif(e(V), Vpw) < 1e-10
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(pairwise) kernel(uniform) dist(spherical) engine(plugin)
	matrix Vpw = e(V)
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(grid) kernel(uniform) dist(spherical) engine(plugin)
	assert mreldif(e(V), Vpw) < 1e-10
	* southern hemisphere, negative longitudes: lattice origin is negative
	replace lat = -60 + 0.5 * ring
	replace lon = -70 + 0.5 * col
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(pairwise) engine(plugin)
	matrix Vpw = e(V)
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(grid) engine(plugin)
	assert e(method) == "grid"
	assert mreldif(e(V), Vpw) < 1e-10
	replace lat = 35 + 0.5 * ring
	replace lon = -20 + 0.5 * col
	* auto picks grid on a dense lattice with a large cutoff
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(800) kernel(uniform) engine(plugin)
	di "method chosen by auto at 800 km: " e(method)
	* scattered points: method(grid) must error, auto must run pairwise
	replace lat = lat + runiform() * 0.01
	rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) method(grid) engine(plugin)" == 3498
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(250) engine(plugin)
	assert e(method) == "pairwise"
}

* ---- small cutoff: the Mata engine must keep full precision (chord form) ----
* points within ~6 km of each other, 1 km cutoff: 1 - dot would lose ~8 digits
clear
set obs 3000
gen double lat = 48 + runiform() * 0.05
gen double lon = 11 + runiform() * 0.08
gen region = ceil(runiform() * 3)
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 0.5 * x1 - 0.3 * x2 + rnormal()
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(1) engine(mata) verbose
matrix Vsm = e(V)
if (`have_plugin') {
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(1) engine(plugin)
	di "1 km cutoff, plugin vs mata: " mreldif(e(V), Vsm)
	assert mreldif(e(V), Vsm) < 1e-12
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(1) dist(spherical) engine(mata)
	matrix Vsm = e(V)
	fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(1) dist(spherical) engine(plugin)
	di "1 km cutoff spherical, plugin vs mata: " mreldif(e(V), Vsm)
	assert mreldif(e(V), Vsm) < 1e-12
}

* ---- IV / 2SLS -------------------------------------------------------------
clear
set obs 4000
gen lat = runiform(25, 49)
gen lon = runiform(-124, -67)
gen region = ceil(runiform() * 5)
gen z1 = rnormal()
gen z2 = rnormal()
gen u = rnormal()
gen x1 = rnormal()
gen x2 = 0.7 * z1 - 0.4 * z2 + 0.5 * u + rnormal()
gen y = 0.5 * x1 - 0.3 * x2 + u + rnormal()
gen w = runiform(0.5, 2)
* cutoff(-1) must reproduce ivreghdfe's robust 2SLS variance (same N/(N-K-df_a) factor)
fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
assert e(iv) == 1
assert e(N) == 4000 & e(df_m) == 2
assert "`e(instd)'" == "x2" & "`e(inexog)'" == "x1"
matrix Vf = e(V)
local bx1 = _b[x1]
local bx2 = _b[x2]
cap which ivreghdfe
local have_ivreghdfe = (_rc == 0)
if (`have_ivreghdfe') {
	ivreghdfe y x1 (x2 = z1 z2), absorb(region) robust
	assert reldif(_b[x1], `bx1') < 1e-12 & reldif(_b[x2], `bx2') < 1e-12
	matrix Vi = e(V)
	di "IV cutoff(-1) vs ivreghdfe robust: " reldif(Vf[1,1], Vi["x1","x1"]) " " reldif(Vf[2,2], Vi["x2","x2"]) " " reldif(Vf[1,2], Vi["x1","x2"])
	assert reldif(Vf[1,1], Vi["x1","x1"]) < 1e-10
	assert reldif(Vf[2,2], Vi["x2","x2"]) < 1e-10
	assert reldif(Vf[1,2], Vi["x1","x2"]) < 1e-10
}
else di as text "ivreghdfe not installed; skipping the 2SLS cross-check"
* weighted
fastconley y x1 (x2 = z1 z2) [aw = w], absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
matrix Vf = e(V)
local bx2 = _b[x2]
if (`have_ivreghdfe') {
	ivreghdfe y x1 (x2 = z1 z2) [aw = w], absorb(region) robust
	assert reldif(_b[x2], `bx2') < 1e-12
	assert reldif(Vf[2,2], e(V)["x2","x2"]) < 1e-10
}
* spatial, both engines, replay, residuals
fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata)
matrix Vm = e(V)
if (`have_plugin') {
	fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300) engine(plugin) resid(ivres)
	di "IV spatial plugin vs mata: " mreldif(e(V), Vm)
	assert mreldif(e(V), Vm) < 1e-12
	qui sum ivres
	assert abs(r(mean)) < 1e-6
}
fastconley
* just identified, two endogenous, no exogenous
fastconley y (x1 x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)
assert e(df_m) == 2
* errors: syntax, underidentified, factor variables in IV lists
rcof "fastconley y x1 (x2 = z1 z2" == 198
rcof "fastconley y x1 (x2 x1 = z1), absorb(region) lat(lat) lon(lon) cutoff(300)" == 3498
rcof "fastconley y i.region (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)" == 198

* ---- panel: spatial + serial HAC, balanced vs general path ----------------
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
gen str8 sunit = "u" + string(unit)

fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) engine(mata)
matrix Vg = e(V)
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced verbose engine(mata)
matrix Vb = e(V)
di "balanced vs general: " mreldif(Vg, Vb)
assert mreldif(Vg, Vb) < 1e-10
if (`have_plugin') {
	fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced engine(plugin) verbose
	di "plugin balanced vs mata: " mreldif(e(V), Vg)
	assert mreldif(e(V), Vg) < 1e-12
	fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) engine(plugin) csrweight(float)
	assert mreldif(e(V), Vg) < 1e-12
	fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced engine(plugin) csrweight(float) pixel(10)
	di "plugin balanced float pixel10 vs mata: " mreldif(e(V), Vg)
	assert mreldif(e(V), Vg) < 1e-3
}
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(sunit) time(time) lag(2) balanced
matrix Vs = e(V)
assert mreldif(Vg, Vs) < 1e-10

* IV in a panel with serial lags, balanced vs general
gen z = rnormal()
gen xe = 0.6 * z + 0.4 * rnormal() + 0.3 * (y - 0.5 * x1 + 0.3 * x2)
fastconley y x1 (xe = z), absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) engine(mata)
matrix Vivg = e(V)
fastconley y x1 (xe = z), absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced
di "IV panel balanced vs general: " mreldif(e(V), Vivg)
assert mreldif(e(V), Vivg) < 1e-10
drop z xe
* lag(0) with a panel equals the pure spatial meat over time blocks
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time)
matrix V0 = e(V)
assert mreldif(V0, Vg) > 1e-4

* unbalanced: drop rows; balanced must error, general must run
drop if runiform() < 0.2
rcof "fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) balanced" == 3498
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2)
assert e(N) < 2000

di as result _n "test_basic.do: all checks passed"
