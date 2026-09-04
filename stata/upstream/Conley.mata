// --------------------------------------------------------------------------
// Conley (1999) spatial HAC standard errors: vce(conley latvar lonvar, ...)
// --------------------------------------------------------------------------
// Arbitrary correlation between observations closer than a cutoff (km),
// with Bartlett (1 - d/cutoff) or uniform weights, and optionally
// within-unit serial correlation up to lag() periods (Bartlett weights in
// time). The meat is M + M' with M = sum_i s_i' c_i,
// c_i = 0.5 s_i + sum_{j > i, d_ij <= cutoff} w_ij s_j, on scores
// s_i = (w_i e_i) x_i, computed with a 3D cell-grid neighbour search at
// cell-pair granularity and dense tiled blocks (pure Mata). Identical
// (to floating-point summation order) to the fastconley R and Stata
// packages, whose engine this file is derived from.
//
// References:
// - Conley (1999), "GMM estimation with cross sectional dependence",
//   Journal of Econometrics
// - Hsiang (2010), PNAS (the panel spatial + serial convention)

mata:

`Void' reghdfe_vce_conley(`FixedEffects' S,
                          `Solution' sol,
                          `Matrix' D,
                          `Matrix' X,
                          `Variable' w,
                          `String' vce_mode)
{
	`Matrix'                Scores, M, Sagg
	`Vector'                resid, lat, lon, time, unit, alat, alon, atime
	`Real'                  dof_adj
	`Integer'               n

	assert_msg(S.conley_lat != "" & S.conley_lon != "", "vce(conley) requires latitude and longitude variables")
	assert_msg(!missing(S.conley_cutoff) & S.conley_cutoff > 0, "vce(conley) requires cutoff(#) > 0")

	// Scores exactly as in reghdfe_vce_dkraay: residuals times (normalized) weights
	resid = S.weight_type != "" ? sol.resid :* w : sol.resid
	Scores = (sol.report_constant ? (X, J(rows(X), 1, 1)) : X) :* resid
	n = rows(Scores)

	lat = st_data(S.sample, S.conley_lat)
	lon = st_data(S.sample, S.conley_lon)
	if (S.conley_time == "") time = J(n, 1, 1)
	else time = st_isstrvar(S.conley_time) ? reghdfe_conley_group(st_sdata(S.sample, S.conley_time)) : st_data(S.sample, S.conley_time)
	if (S.conley_unit == "") unit = (1::n)
	else unit = st_isstrvar(S.conley_unit) ? reghdfe_conley_group(st_sdata(S.sample, S.conley_unit)) : st_data(S.sample, S.conley_unit)

	if (S.verbose > 0) {
		printf("{txt}# Estimating Conley spatial HAC Variance-Covariance Matrix\n\n")
		printf("{txt}   - Kernel: {res}%s{txt}; distance: {res}%s{txt}; cutoff: {res}%g{txt} km; lag: {res}%g{txt}\n",
		       S.conley_kernel, S.conley_dist, S.conley_cutoff, S.conley_lag)
	}

	// Aggregate identical (time, lat, lon) points [+ pixel snap], then the spatial meat
	alat = lat; alon = lon; atime = time; Sagg = Scores
	reghdfe_conley_aggregate(alat, alon, atime, Sagg, S.conley_pixel)
	M = reghdfe_conley_spatial_meat(alat, alon, atime, Sagg, 1, S.conley_cutoff,
	                                S.conley_kernel, S.conley_dist, 1024, S.verbose > 0)
	if (S.conley_lag > 0 & rows(uniqrows(time)) > 1) {
		M = M + reghdfe_conley_serial_meat(unit, time, Scores, S.conley_lag)
	}

	// Same small-sample factor as vce(robust); nossc disables it
	dof_adj = sol.N / (sol.N - S.df_a - sol.df_m)
	if (vce_mode == "vce_asymptotic") dof_adj = sol.N / (sol.N - 1)
	if (!S.conley_ssc) dof_adj = 1
	if (S.verbose > 0 & vce_mode != "vce_asymptotic") {
		printf("{txt}   - Small-sample-adjustment: q = %g\n", dof_adj)
	}

	sol.V = D * M * D * dof_adj
	_makesymmetric(sol.V)
	if (S.conley_psdfix) (void) reghdfe_fix_psd(sol.V, sol.report_constant)

	sol.conley_cutoff = S.conley_cutoff
	sol.conley_lag = S.conley_lag
	sol.conley_kernel = S.conley_kernel
	sol.conley_dist = S.conley_dist
	if (S.verbose > 0) printf("\n")
}

end

// ---------------------------------------------------------------------------
// Engine (derived from fastconley.mata of the fastconley Stata package by
// make_conley_mata.py; do not edit here). Plain Mata, no reghdfe types.
// ---------------------------------------------------------------------------
mata:

// ---------------------------------------------------------------------------
real colvector reghdfe_conley_group(transmorphic colvector x)
{
	real colvector p, codes, chg
	transmorphic colvector xs
	real scalar n
	n = rows(x)
	if (n == 0) return(J(0, 1, .))
	p = order(x, 1)
	xs = x[p]
	chg = J(n, 1, 1)
	if (n > 1) chg[|2 \ n|] = (xs[|2 \ n|] :!= xs[|1 \ n-1|])
	codes = J(n, 1, .)
	codes[p] = runningsum(chg)
	return(codes)
}

// ---------------------------------------------------------------------------
void reghdfe_conley_aggregate(real colvector lat, real colvector lon,
                          real colvector time, real matrix S, real scalar pixel)
{
	real colvector klat, klon, p, coslat, lon_step, grp
	real matrix keys, info
	real scalar lat_step, n
	n = rows(S)
	if (pixel > 0) {
		lat_step = pixel / 111
		klat = round(lat :/ lat_step) :* lat_step
		coslat = cos(klat :* pi() / 180)
		coslat = coslat :* (abs(coslat) :>= 1e-6) :+ 1e-6 :* (abs(coslat) :< 1e-6)
		lon_step = pixel :/ (111 :* coslat)
		klon = round(lon :/ lon_step) :* lon_step
	}
	else {
		klat = lat
		klon = lon
	}
	keys = (time, klat, klon)
	p = order(keys, (1, 2, 3))
	keys = keys[p, .]
	S = S[p, .]
	// group id = running count of rows whose (time, lat, lon) differs from the previous row
	grp = J(n, 1, 1)
	if (n > 1) grp[|2 \ n|] = (rowsum(keys[|2, 1 \ n, 3|] :!= keys[|1, 1 \ n-1, 3|]) :> 0)
	info = panelsetup(runningsum(grp), 1)
	if (rows(info) == n) {           // nothing to merge
		time = keys[., 1]; lat = keys[., 2]; lon = keys[., 3]
		return
	}
	S = panelsum(S, info)
	time = keys[info[., 1], 1]
	lat = keys[info[., 1], 2]
	lon = keys[info[., 1], 3]
}

// ---------------------------------------------------------------------------
real matrix reghdfe_conley_spatial_meat(real colvector lat, real colvector lon,
                                    real colvector time, real matrix S,
                                    real scalar T, real scalar cutoff,
                                    string scalar kernel, string scalar dist,
                                    real scalar tile, real scalar verbose)
{
	real scalar n, k, ang, edge, G, coscut, b, c, kk, nb, lo, hi, ncell, ca, cb, o
	real scalar cx, cy, cz, nx, ny, nz, u, npairs
	real colvector latr, lonr, key, p, ucid
	real matrix U, M, tinfo, cinfo, offs, cell
	transmorphic A

	n = rows(S)
	k = cols(S) / T
	M = J(k, k, 0)
	if (n == 0) return(M)
	if (cutoff < 0) {
		// Only the 0.5 S_i diagonal survives: M + M' = sum over periods of S_t' S_t.
		for (b = 1; b <= T; b++) {
			M = M + quadcross(S[|1, (b-1)*k+1 \ n, b*k|], S[|1, (b-1)*k+1 \ n, b*k|])
		}
		return(M)
	}
	if (dist == "chord") ang = 2 * asin(min((1, cutoff / (2 * 6371))))
	else ang = cutoff / 6371
	if (ang > pi()) ang = pi()
	coscut = cos(ang)
	edge = 2 * sin(ang / 2)
	if (edge <= 0) _error(3498, "fastconley: cutoff must be positive (or negative for the no-pairs sentinel)")
	G = ceil(2 / edge) + 1
	if (G^3 >= 2^53) _error(3498, "fastconley: cutoff too small for the Mata engine's cell grid; use engine(plugin)")

	latr = lat :* (pi() / 180)
	lonr = lon :* (pi() / 180)
	U = (cos(latr) :* cos(lonr), cos(latr) :* sin(lonr), sin(latr))
	cell = floor((U :+ 1) :/ edge)
	cell = cell :* (cell :>= 0)
	cell = cell :* (cell :<= G-1) :+ (G-1) :* (cell :> G-1)
	key = (cell[., 1] :* G :+ cell[., 2]) :* G :+ cell[., 3]

	// Sort rows by (time, cell key): cells are contiguous within a time block.
	p = order((time, key), (1, 2))
	U = U[p, .]
	S = S[p, .]
	key = key[p]
	tinfo = panelsetup(time[p], 1)
	nb = rows(tinfo)

	// The 13 lexicographically-forward neighbor offsets (dx, dy, dz).
	offs = (0,0,1 \ 0,1,-1 \ 0,1,0 \ 0,1,1 \ 1,-1,-1 \ 1,-1,0 \ 1,-1,1 \
	        1,0,-1 \ 1,0,0 \ 1,0,1 \ 1,1,-1 \ 1,1,0 \ 1,1,1)

	npairs = 0
	for (b = 1; b <= nb; b++) {
		lo = tinfo[b, 1]
		hi = tinfo[b, 2]
		cinfo = panelsetup(key[|lo \ hi|], 1) :+ (lo - 1)   // absolute rows
		ncell = rows(cinfo)
		ucid = key[cinfo[., 1]]
		A = asarray_create("real", 1)
		for (c = 1; c <= ncell; c++) asarray(A, ucid[c], c)
		for (ca = 1; ca <= ncell; ca++) {
			u = ucid[ca]
			cz = mod(u, G)
			cy = mod(floor(u / G), G)
			cx = floor(u / (G * G))
			npairs = npairs + reghdfe_conley_cell_pair(M, U, S, T, k, cinfo[ca, 1], cinfo[ca, 2],
			                                       cinfo[ca, 1], cinfo[ca, 2], 1,
			                                       coscut, cutoff, kernel, dist, tile)
			for (o = 1; o <= 13; o++) {
				nx = cx + offs[o, 1]; ny = cy + offs[o, 2]; nz = cz + offs[o, 3]
				if (nx < 0 | nx >= G | ny < 0 | ny >= G | nz < 0 | nz >= G) continue
				u = (nx * G + ny) * G + nz
				if (!asarray_contains(A, u)) continue
				cb = asarray(A, u)
				npairs = npairs + reghdfe_conley_cell_pair(M, U, S, T, k, cinfo[ca, 1], cinfo[ca, 2],
				                                       cinfo[cb, 1], cinfo[cb, 2], 0,
				                                       coscut, cutoff, kernel, dist, tile)
			}
		}
	}
	if (verbose) printf("{txt}   - spatial meat: %g rows, %g time blocks, %g accepted pairs (cell edge %g)\n",
	                    n, nb, npairs, edge)
	return(M + M')
}

// One cell pair (rows ia..ib against ja..jb), tiled. `self` = same cell:
// tiles on the diagonal use the strict upper triangle with 0.5 on the
// diagonal, exactly the C++ C_i = 0.5 S_i + sum_{j>i}. Returns the number of
// accepted (i<j) pairs for the verbose summary.
real scalar reghdfe_conley_cell_pair(real matrix M, real matrix U, real matrix S,
                                 real scalar T, real scalar k,
                                 real scalar ia, real scalar ib,
                                 real scalar ja, real scalar jb, real scalar self,
                                 real scalar coscut, real scalar cutoff,
                                 string scalar kernel, string scalar dist,
                                 real scalar tile)
{
	real scalar ti, ti2, tj, tj2, t, npairs, c0, c1
	real matrix Ua, Ub, Dot, W, WS, acc, d, a
	real colvector ea, eb

	npairs = 0
	for (ti = ia; ti <= ib; ti = ti + tile) {
		ti2 = min((ib, ti + tile - 1))
		Ua = U[|ti, 1 \ ti2, 3|]
		for (tj = (self ? ti : ja); tj <= jb; tj = tj + tile) {
			tj2 = min((jb, tj + tile - 1))
			Ub = U[|tj, 1 \ tj2, 3|]
			Dot = Ua * Ub'
			if (max(Dot) < coscut) {
				// No pair within the cutoff; only a self tile still needs its diagonal.
				if (!(self & tj == ti)) continue
			}
			if (kernel == "uniform") {
				// Every supported distance is monotone in the chord, so the
				// accept test is the dot threshold alone (the C++ pair_weight
				// uniform specializations use the same criterion). The dot
				// product's rounding only matters exactly at the boundary.
				acc = (Dot :>= coscut)
				W = acc
			}
			else {
				// Bartlett needs distances. The squared chord |u_i - u_j|^2 equals
				// 2 - 2 Dot, but that difference cancels at small angles: its
				// relative error is ~2e-16 / (theta^2 / 2), i.e. 1e-12 at 200 km
				// and 1e-8 at 1 km. Below 200 km the chord is therefore built
				// from coordinate differences (full precision, eight extra
				// passes per tile); above, the dot form is exact to 1e-12 and
				// ~40% cheaper. (Outer products with ones vectors: Mata's colon
				// operators do not broadcast a column against a row.)
				if (cutoff < 200) {
					ea = J(rows(Ua), 1, 1)
					eb = J(rows(Ub), 1, 1)
					a = (Ua[., 1] * eb' - ea * Ub[., 1]') :^ 2 + (Ua[., 2] * eb' - ea * Ub[., 2]') :^ 2 + (Ua[., 3] * eb' - ea * Ub[., 3]') :^ 2
				}
				else {
					a = 2 :- 2 :* Dot
					a = a :* (a :> 0)
				}
				if (dist == "chord") {
					d = 6371 :* sqrt(a)
				}
				else {
					// great-circle angle 2 asin(chord / 2) (haversine and
					// spherical are the same distance; this form is well
					// conditioned where acos(Dot) is not)
					a = a :/ 4
					a = a :* (a :<= 1) :+ (a :> 1)
					d = (2 * 6371) :* asin(sqrt(a))
				}
				acc = (d :<= cutoff)
				W = acc :* (1 :- d :/ cutoff)
			}
			if (self & tj == ti) {
				npairs = npairs + (sum(acc) - rows(acc)) / 2
				W = uppertriangle(W, 0.5)
			}
			else npairs = npairs + sum(acc)
			// M += Sa' W Sb, per stacked period.
			WS = W * S[|tj, 1 \ tj2, T*k|]
			for (t = 1; t <= T; t++) {
				c0 = (t-1)*k + 1; c1 = t*k
				M = M + S[|ti, c0 \ ti2, c1|]' * WS[|1, c0 \ rows(WS), c1|]
			}
		}
	}
	return(npairs)
}

// ---------------------------------------------------------------------------
real matrix reghdfe_conley_serial_meat(real colvector unit, real colvector time,
                                   real matrix S, real scalar L)
{
	real scalar n, k, s
	real colvector p, u, t, dt, same, w
	real matrix M, Ss

	n = rows(S)
	k = cols(S)
	M = J(k, k, 0)
	if (n < 2 | L <= 0) return(M)
	p = order((unit, time), (1, 2))
	u = unit[p]; t = time[p]; Ss = S[p, .]
	for (s = 1; s < n; s++) {
		dt = t[|s+1 \ n|] - t[|1 \ n-s|]
		same = (u[|s+1 \ n|] :== u[|1 \ n-s|])
		if (!any(same :& (dt :<= L))) break
		w = same :* (dt :> 0) :* (dt :<= L) :* (1 :- dt :/ (L + 1))
		M = M + quadcross(Ss[|1, 1 \ n-s, k|], w, Ss[|s+1, 1 \ n, k|])
	}
	return(M + M')
}

end
