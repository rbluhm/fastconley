* Smoke + consistency tests for the Stata fastconley command (Mata and plugin).
* Run from the repository root:  stata-mp -b do stata/test/test_basic.do
clear all
adopath ++ "stata/src"
set seed 20260903

which fastconley
fastconley, version

* ---- cross-section with one absorbed FE ------------------------------------
set obs 3000
gen lat = runiform(25, 49)
gen lon = runiform(-124, -67)
gen region = ceil(runiform() * 5)
gen x1 = rnormal()
gen x2 = rnormal()
gen w = runiform(0.5, 2)
gen fw = ceil(3 * w)
gen pw = w
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
else {
	if ("$FASTCONLEY_REQUIRE_PLUGIN" != "") {
		di as error "FASTCONLEY_REQUIRE_PLUGIN is set but the compiled engine did not load:"
		fastconley, version
		error 498
	}
	di as text "plugin not available on this platform; skipping plugin checks"
}
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(mata) method(grid)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) engine(nope)" == 198

* uniform kernel, spherical distance, verbose
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(spherical) engine(mata) verbose
matrix Vverbose = e(V)
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(spherical) engine(mata)
di "uniform Mata verbose vs non-verbose: " mreldif(e(V), Vverbose)
assert mreldif(e(V), Vverbose) <= 1e-14

* Uniform Mata meat against an independent dense pair loop. The small-cutoff
* geometry occupies one cell whose bounding box is wider than the cutoff, so
* the ordinary acceptance-matrix path is used. At the large cutoff all 300
* points and every tile bounding box fit inside the cutoff, exercising the
* all-accepted path on every eligible off-diagonal tile pair. Diagonal self
* tiles deliberately retain the original strict-upper/0.5-diagonal path.
quietly include "stata/src/fastconley.mata"
mata:
void fastconley_test_uniform_dense()
{
	real scalar n, k, cc, cutoff, edge, q, lo, threshold, i, j, a, rd, junk
	real scalar ti, ti2, tj, tj2, tile
	real colvector idx, lat, lon, latr, ycoord, zcoord, ea, eb
	real matrix S, U, got, ref, ordinary, fasttile, Ua, Ub, A, W, WS
	real rowvector si, sj

	n = 300
	k = 3
	idx = 1::n
	S = (sin(idx:/7), cos(idx:/11), (mod(idx, 13):-6):/7)
	for (cc = 1; cc <= 2; cc++) {
		if (cc == 1) {
			cutoff = 100
			edge = 2 * sin(cutoff / (2 * 6371))
			q = floor(1 / edge)
			lo = q * edge - 1
			ycoord = lo :+ edge :* (0.1 :+ 0.8 :* mod(idx, 2))
			zcoord = lo :+ edge :* (0.1 :+ 0.8 :* mod(floor((idx:-1):/2), 2))
			latr = asin(zcoord)
			lat = latr :* (180 / pi())
			lon = asin(ycoord :/ cos(latr)) :* (180 / pi())
		}
		else {
			cutoff = 7000
			lat = 25 :+ mod(idx:*37, 2400) :/ 100
			lon = -124 :+ mod(idx:*83, 5700) :/ 100
		}
		latr = lat :* (pi() / 180)
		U = (cos(latr) :* cos(lon :* (pi()/180)),
		     cos(latr) :* sin(lon :* (pi()/180)), sin(latr))
		threshold = 4 * sin(min((pi(), cutoff / 6371)) / 2)^2
		ref = J(k, k, 0)
		for (i = 1; i <= n; i++) {
			si = S[i, .]
			ref = ref + si' * si
			for (j = i + 1; j <= n; j++) {
				a = sum((U[i, .] - U[j, .]) :^ 2)
				if (a <= threshold) {
					sj = S[j, .]
					ref = ref + si' * sj + sj' * si
				}
			}
		}
		got = fastconley_spatial_meat(lat, lon, J(n, 1, 1), S, 1,
		                                  cutoff, "uniform", "spherical",
		                                  cc == 1 ? 512 : 37, 0)
		rd = mreldif(got, ref)
		if (cc == 1) st_numscalar("fc_uniform_dense_small", rd)
		else {
			st_numscalar("fc_uniform_dense_large", rd)
			// Drive fastconley_cell_pair() as a cross-cell pair and reproduce
			// its pre-optimization acceptance-matrix math at the same tile size.
			// Both 150-row halves span the region, so their box bound is below
			// this cutoff and every cross-tile pair takes the fast branch.
			fasttile = J(k, k, 0)
			junk = fastconley_cell_pair(fasttile, U, S, 1, k, 1, 150,
			                                151, 300, 0, cos(cutoff/6371),
			                                cutoff, "uniform", "spherical", 37)
			ordinary = J(k, k, 0)
			tile = 37
			for (ti = 1; ti <= 150; ti = ti + tile) {
				ti2 = min((150, ti + tile - 1))
				Ua = U[|ti, 1 \ ti2, 3|]
				ea = J(rows(Ua), 1, 1)
				for (tj = 151; tj <= 300; tj = tj + tile) {
					tj2 = min((300, tj + tile - 1))
					Ub = U[|tj, 1 \ tj2, 3|]
					eb = J(rows(Ub), 1, 1)
					A = (Ua[.,1]*eb' - ea*Ub[.,1]') :^ 2 +
					    (Ua[.,2]*eb' - ea*Ub[.,2]') :^ 2 +
					    (Ua[.,3]*eb' - ea*Ub[.,3]') :^ 2
					W = (A :<= threshold)
					WS = W * S[|tj, 1 \ tj2, k|]
					ordinary = ordinary + S[|ti, 1 \ ti2, k|]' * WS
				}
			}
			st_numscalar("fc_uniform_fast_vs_ordinary", mreldif(fasttile, ordinary))
		}
	}
}
fastconley_test_uniform_dense()
end
di "uniform dense reference, small cutoff: " scalar(fc_uniform_dense_small)
di "uniform dense reference, all-accepted cutoff: " scalar(fc_uniform_dense_large)
di "uniform all-accepted vs ordinary tile path: " scalar(fc_uniform_fast_vs_ordinary)
assert scalar(fc_uniform_dense_small) < 1e-13
assert scalar(fc_uniform_dense_large) < 1e-13
assert scalar(fc_uniform_fast_vs_ordinary) <= 1e-14

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

* fweights run (the Conley duplicate-location semantics intentionally differ
* from robust); pweights both run and reproduce robust at cutoff(-1)
quietly summarize fw, meanonly
local fN = r(sum)
fastconley y x1 x2 [fw = fw], absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
assert e(N) == `fN'
assert !missing(_se[x1])
fastconley y x1 x2 [pw = pw], absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
matrix Vc = e(V)
reghdfe y x1 x2 [pw = pw], absorb(region) vce(robust)
matrix Vr = e(V)
di "pweight: mreldif = " mreldif(Vc, Vr)
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
assert mreldif(Vp, V1) < 1e-4

* factor variables, collinear regressor, residuals()
gen x3 = 2 * x1
fastconley y x1 x2 x3 i.region, absorb(region) lat(lat) lon(lon) cutoff(300) resid(res)
assert _b[x3] == 0
assert !missing(_se[x1]) & _se[x1] > 0
qui sum res
assert abs(r(mean)) < 1e-6

* grouped time-series/factor expressions are not mistaken for IV clauses
gen t = _n
gen region2 = mod(region, 3) + 1
tsset t
fastconley y L(1/2).x1, noabsorb lat(lat) lon(lon) cutoff(300)
assert e(N) == 2998
fastconley y c.(x1 x2) i.(region region2), noabsorb lat(lat) lon(lon) cutoff(300)
assert e(N) == 3000

* absorb(..., savefe) follows reghdfe's store_alphas flow
capture drop __hdfe1__
fastconley y x1 x2, absorb(region, savefe) lat(lat) lon(lon) cutoff(300)
confirm numeric variable __hdfe1__, exact
drop __hdfe1__

* string lat/lon rejected, missing option rejected
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(gaussian)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) lag(2)" == 198
rcof "fastconley y x1 x2, lat(lat) lon(lon) cutoff(300)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(.)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) lag(.)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) pixel(.)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) tile(.)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) threads(.)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) lag(1.5)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) threads(1.5)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) tile(8193)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) tile(6000)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) vce(robust)" == 198
rcof "fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) cluster(region)" == 198

* The un-clamped slope covariance is rank zero when every point is mutually
* correlated under a uniform kernel, so PSD flooring must not manufacture F.
replace lat = 40
replace lon = -90
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) engine(mata)
assert missing(e(F))

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
	* A regular 7-degree longitude lattice cannot tile 360 degrees. The forced
	* plugin failure must print fc_plugin_error before returning 198.
	preserve
	clear
	set obs 204
	gen ring = mod(_n - 1, 4)
	gen col = floor((_n - 1) / 4)
	gen double lat = -10 + 7 * ring
	gen double lon = -175 + 7 * col
	gen region = mod(_n, 3) + 1
	gen x1 = rnormal()
	gen y = x1 + rnormal()
	cap noi fastconley y x1, absorb(region) lat(lat) lon(lon) cutoff(1500) method(grid) engine(plugin)
	local grid_fail_rc = _rc
	assert `grid_fail_rc' == 198
	restore
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
fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(.01) engine(mata)
assert e(N) == 3000 & !missing(_se[x1])
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

* ---- near-cutoff uniform acceptance: stable squared-coordinate test --------
clear
set obs 4000
gen block = ceil(_n / 2)
gen byte second = mod(_n, 2) == 0
gen double lat = 0
gen double lon = second * ((100 + cond(mod(block, 2), -1e-9, 1e-9)) / 6371) * 180 / _pi
gen x1 = rnormal()
gen y = .4 * x1 + rnormal()
fastconley y x1, noabsorb lat(lat) lon(lon) cutoff(100) kernel(uniform) time(block) engine(mata)
matrix Vboundary = e(V)
if (`have_plugin') {
	fastconley y x1, noabsorb lat(lat) lon(lon) cutoff(100) kernel(uniform) time(block) engine(plugin)
	di "near-cutoff uniform plugin vs mata: " mreldif(e(V), Vboundary)
	assert mreldif(e(V), Vboundary) < 1e-12
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
fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix residuals(ivres0)
assert e(iv) == 1
assert e(N) == 4000 & e(df_m) == 2
assert "`e(instd)'" == "x2" & "`e(inexog)'" == "x1"
assert "`e(insts)'" == "x1 z1 z2" & "`e(exogr)'" == "x1"
assert "`e(predict)'" == "fastconley_p"
assert e(df_r) < e(N) & e(rank) == 2 & e(rmse) > 0
assert strpos(`"`e(cmdline)'"', "cutoff(-1)") > 0
assert strpos(`"`e(cmdline)'"', "residuals(ivres0)") > 0
predict double ivxb0, xb
predict double ivrp0, residuals
assert reldif(ivrp0, ivres0) < 1e-15 if e(sample)
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

* Regressor collinearity is removed before instrument construction and posted
* as an omitted coefficient with a zero V row/column.
gen x3 = 2 * x1
fastconley y x1 x3 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(-1) nopsdfix
assert _b[x3] == 0 & _se[x3] == 0
local iv_colnames : colfullnames e(b)
assert strpos("`iv_colnames'", "o.x3") > 0
matrix Vfcoll = e(V)
if (`have_ivreghdfe') {
	ivreghdfe y x1 x3 (x2 = z1 z2), absorb(region) robust
	assert _b[x3] == 0
	assert reldif(Vfcoll["x1","x1"], e(V)["x1","x1"]) < 1e-10
	assert reldif(Vfcoll["x2","x2"], e(V)["x2","x2"]) < 1e-10
}
drop x3

* noabsorb IV includes _cons and matches ivregress 2sls, including robust V;
* noconstant removes it and matches the corresponding no-constant model.
fastconley y x1 (x2 = z1 z2), noabsorb lat(lat) lon(lon) cutoff(-1) nossc nopsdfix residuals(ivres_na)
matrix bfn = e(b)
matrix Vfn = e(V)
assert colsof(bfn) == 3 & !missing(_b[_cons]) & e(rank) == 3
predict double ivxb_na, xb
predict double ivrp_na, residuals
assert reldif(ivrp_na, ivres_na) < 1e-15 if e(sample)
ivregress 2sls y x1 (x2 = z1 z2), vce(robust)
matrix Viv = e(V)
assert reldif(bfn[1,"x1"], _b[x1]) < 1e-12
assert reldif(bfn[1,"x2"], _b[x2]) < 1e-12
assert reldif(bfn[1,"_cons"], _b[_cons]) < 1e-12
di "noabsorb IV vs ivregress robust V (nossc): " reldif(Vfn["x1","x1"], Viv["x1","x1"])
foreach a in x1 x2 _cons {
	foreach b in x1 x2 _cons {
		assert reldif(Vfn["`a'","`b'"], Viv["`a'","`b'"]) < 1e-10
	}
}
fastconley y x1 (x2 = z1 z2), noabsorb noconstant lat(lat) lon(lon) cutoff(-1) nossc nopsdfix
matrix bfn = e(b)
matrix Vfn = e(V)
assert colsof(bfn) == 2 & e(rank) == 2
ivregress 2sls y x1 (x2 = z1 z2), noconstant vce(robust)
matrix Viv = e(V)
assert reldif(bfn[1,"x1"], _b[x1]) < 1e-12
assert reldif(bfn[1,"x2"], _b[x2]) < 1e-12
foreach a in x1 x2 {
	foreach b in x1 x2 {
		assert reldif(Vfn["`a'","`b'"], Viv["`a'","`b'"]) < 1e-10
	}
}

* IV savefe uses the same reghdfe store_alphas flow as OLS.
capture drop __hdfe1__
fastconley y x1 (x2 = z1 z2), absorb(region, savefe) lat(lat) lon(lon) cutoff(300)
confirm numeric variable __hdfe1__, exact
drop __hdfe1__

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
rcof "fastconley y x1 (x2 u = z1), absorb(region) lat(lat) lon(lon) cutoff(300)" == 3498
rcof "fastconley y i.region (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)" == 198

* ---- panel: spatial + serial HAC, balanced vs general path ----------------
* Reset after optional plugin-only tests so their random draws cannot change
* the deterministic panel fixture when the shipped plugin is unavailable.
set seed 20260904
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
gen double timegap = 2000 + 2 * (time - 1)
gen str8 stimegap = string(timegap, "%9.0f")
gen str8 stimelabel = "t" + string(time)

* time() alone blocks the spatial covariance; it is not ignored
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) time(time)
matrix Vtimeonly = e(V)
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300)
assert mreldif(e(V), Vtimeonly) > 1e-6

* numeric-looking string time keeps gaps; arbitrary labels cannot drive lags
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(timegap) lag(2)
matrix Vgap = e(V)
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(stimegap) lag(2)
assert mreldif(e(V), Vgap) < 1e-12
rcof "fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(stimelabel) lag(2)" == 198

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
quietly count
local unbal_N = r(N)
fastconley y x1 x2, absorb(unit time) lat(ulat) lon(ulon) cutoff(300) unit(unit) time(time) lag(2) keepsingletons
matrix Vunbal = e(V)
di "final unbalanced raw N/final N/V11/V22: " `unbal_N' " " e(N) " " Vunbal[1,1] " " Vunbal[2,2]
di %21.17g Vunbal[1,1] " " %21.17g Vunbal[2,2]
assert e(N) == 1619 & `unbal_N' == 1620
assert reldif(Vunbal[1,1], .0009170766510187081) < 1e-12
assert reldif(Vunbal[2,2], .0009104991328768975) < 1e-12
assert !missing(e(F))

* ---- exact antipodes: a cutoff covering the whole sphere must accept them
* (squared chord 4 + O(1e-15) is clamped; matches the C++ engine) ---------
mata:
	ap_lat = (28.597969631664455 \ -28.597969631664455)
	ap_lon = (126.38749223202467 \ 126.38749223202467 - 180)
	ap_t = (1 \ 1)
	ap_S = (1 \ 2)
	foreach d in ("haversine", "spherical", "chord") {
		ap_M = fastconley_spatial_meat(ap_lat, ap_lon, ap_t, ap_S, 1, 25000, "uniform", d, 1024, 0)
		if (abs(ap_M[1,1] - 9) > 1e-12) _error(9, "antipodes rejected by uniform " + d)
		ap_M = fastconley_spatial_meat(ap_lat, ap_lon, ap_t, ap_S, 1, 300, "uniform", d, 1024, 0)
		if (abs(ap_M[1,1] - 5) > 1e-12) _error(9, "self terms wrong for uniform " + d)
	}
	ap_M = fastconley_spatial_meat(ap_lat, ap_lon, ap_t, ap_S, 1, 2 * 12742, "bartlett", "chord", 1024, 0)
	if (abs(ap_M[1,1] - 7) > 1e-9) _error(9, "antipode bartlett/chord weight wrong")
end
di as result "exact antipodes accepted at whole-sphere cutoffs (uniform x 3 distances, bartlett chord)"

di as result _n "test_basic.do: all checks passed"
