*! fastconley.mata -- Conley (1999) spatial HAC engine in pure Mata.
*! Mirrors the C++ engine of the fastconley R package (src/conley_core.h):
*! same score-accumulation meat M + M' with M = sum_i S_i' C_i and
*! C_i = 0.5 S_i + sum_{j>i, d_ij <= cutoff} w_ij S_j, same distances, same
*! kernels, same serial (within-unit) Bartlett HAC. Neighbor search runs at
*! cell-pair granularity (3D unit-vector cell grid, edge = chord of the
*! cutoff) with tiled dense BLAS inside each cell pair. Nothing here refers to
*! reghdfe; the ado passes plain matrices.

mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// Grouping without ftools: integer codes 1..L in sorted-value order for a
// string or real column (R's as.integer(factor(x))).
// ---------------------------------------------------------------------------
real colvector fastconley_group(transmorphic colvector x)
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
// Score pre-aggregation: sum scores over identical (time, lat, lon) keys,
// optionally after snapping to a pixel-km grid (same arithmetic as the R
// package's aggregate_scores). Returns rows sorted by (time, lat, lon).
// ---------------------------------------------------------------------------
void fastconley_aggregate(real colvector lat, real colvector lon,
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
// Spatial meat. S may hold T stacked periods side by side (n x T*k) sharing
// the coordinates of one period (balanced panel); T = 1 is the general
// case with `time` marking contiguous blocks after sorting. cutoff < 0 is
// the no-pairs sentinel (heteroskedasticity-only meat, sum_i s_i s_i').
// ---------------------------------------------------------------------------
real matrix fastconley_spatial_meat(real colvector lat, real colvector lon,
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
			npairs = npairs + fastconley_cell_pair(M, U, S, T, k, cinfo[ca, 1], cinfo[ca, 2],
			                                       cinfo[ca, 1], cinfo[ca, 2], 1,
			                                       coscut, cutoff, kernel, dist, tile)
			for (o = 1; o <= 13; o++) {
				nx = cx + offs[o, 1]; ny = cy + offs[o, 2]; nz = cz + offs[o, 3]
				if (nx < 0 | nx >= G | ny < 0 | ny >= G | nz < 0 | nz >= G) continue
				u = (nx * G + ny) * G + nz
				if (!asarray_contains(A, u)) continue
				cb = asarray(A, u)
				npairs = npairs + fastconley_cell_pair(M, U, S, T, k, cinfo[ca, 1], cinfo[ca, 2],
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
real scalar fastconley_cell_pair(real matrix M, real matrix U, real matrix S,
                                 real scalar T, real scalar k,
                                 real scalar ia, real scalar ib,
                                 real scalar ja, real scalar jb, real scalar self,
                                 real scalar coscut, real scalar cutoff,
                                 string scalar kernel, string scalar dist,
                                 real scalar tile)
{
	real scalar ti, ti2, tj, tj2, t, npairs, c0, c1
	real matrix Ua, Dot, W, WS, acc, d, a

	npairs = 0
	for (ti = ia; ti <= ib; ti = ti + tile) {
		ti2 = min((ib, ti + tile - 1))
		Ua = U[|ti, 1 \ ti2, 3|]
		for (tj = (self ? ti : ja); tj <= jb; tj = tj + tile) {
			tj2 = min((jb, tj + tile - 1))
			Dot = Ua * U[|tj, 1 \ tj2, 3|]'
			if (max(Dot) < coscut) {
				// No pair within the cutoff; only a self tile still needs its diagonal.
				if (!(self & tj == ti)) continue
			}
			if (kernel == "uniform") {
				// Every supported distance is monotone in the chord, so the
				// accept test is the dot threshold alone (the C++ pair_weight
				// uniform specializations use the same criterion).
				acc = (Dot :>= coscut)
				W = acc
			}
			else {
				if (dist == "haversine") {
					a = (1 :- Dot) :/ 2
					a = a :* (a :> 0)
					a = a :* (a :< 1) :+ (a :>= 1)
					d = (2 * 6371) :* asin(sqrt(a))
				}
				else if (dist == "spherical") {
					a = Dot :* (abs(Dot) :<= 1) :+ (Dot :> 1) :- (Dot :< -1)
					d = 6371 :* acos(a)
				}
				else {
					a = 2 :- 2 :* Dot
					a = a :* (a :> 0)
					d = 6371 :* sqrt(a)
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
// Serial (within-unit) Bartlett HAC meat with bandwidth L+1: pairs in the
// same unit with 0 < t_j - t_i <= L get weight 1 - dt/(L+1). Any row order.
// ---------------------------------------------------------------------------
real matrix fastconley_serial_meat(real colvector unit, real colvector time,
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

// ---------------------------------------------------------------------------
// PSD fix: clamp eigenvalues to 1e-16 (fixest / fastconley convention).
// Returns 1 when the clamp changed V by more than 1e-8.
// ---------------------------------------------------------------------------
real scalar fastconley_psd_fix(real matrix V, real scalar apply)
{
	real matrix Vs, X, Vf
	real rowvector lambda
	real scalar noticeable
	Vs = (V + V') / 2
	symeigensystem(Vs, X = ., lambda = .)
	if (min(lambda) > 0) return(0)
	lambda = lambda :* (lambda :> 1e-16) :+ 1e-16 :* (lambda :<= 1e-16)
	Vf = X * diag(lambda) * X'
	Vf = (Vf + Vf') / 2
	noticeable = max(abs(Vs - Vf)) > 1e-8
	if (apply) V = Vf
	return(noticeable)
}

// ---------------------------------------------------------------------------
// Raster lattice detection and the pairwise-vs-grid choice (ports of the R
// package's detect_lonlat_grid() / choose_grid_method()). On success
// fastconley_detect_grid() fills the globals fc_grid = (lat0, dlat, lon0,
// dlon, n_ring, n_col, n_col_full) and fc_ring / fc_col (0-based indices).
// ---------------------------------------------------------------------------
// `tol` is relative to the lattice step: 1e-6 for double coordinates (the R
// package's value); the ado passes 1e-3 when a coordinate variable is stored
// as float, whose rounding noise (~6e-8 relative) would otherwise defeat the
// exact-multiple test. The step is refined as range / number of steps so
// that float noise in the smallest gap does not accumulate along the ring.
real scalar fastconley_detect_grid(real colvector lat, real colvector lon, real scalar tol)
{
	real colvector ul, uo, dl, dol, ring_d, col_d
	real scalar step_l, step_o, n_col, cf, n_col_full, nsteps
	external real rowvector fc_grid
	external real colvector fc_ring, fc_col
	ul = uniqrows(lat); uo = uniqrows(lon)
	if (rows(ul) < 2 | rows(uo) < 2) return(0)
	dl = ul[|2 \ rows(ul)|] - ul[|1 \ rows(ul)-1|]
	dol = uo[|2 \ rows(uo)|] - uo[|1 \ rows(uo)-1|]
	step_l = min(dl); step_o = min(dol)
	if (step_l <= 0 | step_o <= 0) return(0)
	nsteps = round((ul[rows(ul)] - ul[1]) / step_l)
	if (nsteps >= 1) step_l = (ul[rows(ul)] - ul[1]) / nsteps
	nsteps = round((uo[rows(uo)] - uo[1]) / step_o)
	if (nsteps >= 1) step_o = (uo[rows(uo)] - uo[1]) / nsteps
	if (any(abs(dl :/ step_l - round(dl :/ step_l)) :> tol)) return(0)
	if (any(abs(dol :/ step_o - round(dol :/ step_o)) :> tol)) return(0)
	ring_d = round((lat :- ul[1]) :/ step_l)
	col_d = round((lon :- uo[1]) :/ step_o)
	if (max(ring_d) > 2^31 - 2 | max(col_d) > 2^31 - 2) return(0)
	if (max(abs(lat :- (ul[1] :+ ring_d :* step_l))) > tol * step_l) return(0)
	if (max(abs(lon :- (uo[1] :+ col_d :* step_o))) > tol * step_o) return(0)
	n_col = max(col_d) + 1
	cf = 360 / step_o
	n_col_full = (abs(cf - round(cf)) < 1e-6 & round(cf) >= n_col) ? round(cf) : 0
	fc_grid = (ul[1], step_l, uo[1], step_o, max(ring_d) + 1, n_col, n_col_full)
	fc_ring = ring_d
	fc_col = col_d
	return(1)
}

real scalar fastconley_choose_grid(string scalar method, string scalar kernel,
                                   real colvector lat, real colvector lon,
                                   real colvector time, real scalar kk, real scalar cutoff,
                                   real scalar tol)
{
	real scalar cells, n, ang, rho_r, cosbar, dbar, pairs_est, per_pair, npad, grid_work
	real colvector nb
	real matrix info
	external real rowvector fc_grid
	if (method == "pairwise") return(0)
	if (!fastconley_detect_grid(lat, lon, tol)) {
		// (_error() messages must stay short)
		if (method == "grid") _error(3498, "method(grid): no regular lat/lon lattice detected")
		return(0)
	}
	cells = fc_grid[5] * fc_grid[6]
	n = rows(lat)
	if (cells > 100 * n | cells > 4e9) {
		if (method == "grid") _error(3498, "method(grid): lattice too sparse or degenerate; use method(pairwise)")
		return(0)
	}
	if (method == "grid") return(1)
	ang = cutoff / 6371
	rho_r = max((1, ang / (fc_grid[2] * pi() / 180)))
	cosbar = max((0.05, mean(cos(lat :* pi() / 180))))
	dbar = max((1, ang / (fc_grid[4] * pi() / 180) / cosbar))
	info = panelsetup(time, 1)
	nb = info[., 2] - info[., 1] :+ 1
	pairs_est = sum(nb :^ 2) * pi() * rho_r * dbar / cells / 2
	if (kernel == "uniform") per_pair = fc_grid[6]
	else {
		npad = 2 ^ ceil(log(fc_grid[6] + min((dbar, fc_grid[6]))) / log(2))
		per_pair = npad * (5 * log(npad) / log(2) + 2 * kk) / kk
	}
	grid_work = rows(nb) * fc_grid[5] * (2 * rho_r + 1) * per_pair
	return(pairs_est > 2 * grid_work)
}

// ---------------------------------------------------------------------------
// Step 1 of the VCE: bread, scores, and engine inputs from reghdfe-style
// pieces (all plain matrices; see the ado for how they are pulled out of the
// HDFE object after reghdfe_solve_ols). Fills the globals
//   fc_D        bread (K [+1] square, _cons last), fc_kk = its dimension
//   fc_S        scores (n x kk) in sample order, with fc_time / fc_unit
//   fc_sp_*     spatial-engine input rows: fc_sp_lat, fc_sp_lon, fc_sp_time,
//               fc_sp_S, aggregated over identical (time, lat, lon) [+ pixel
//               snap] and sorted by (time, lat, lon); or, when `balanced`
//               was requested and the aggregation broke the balance, the
//               unaggregated rows sorted by (time, unit)
//   fc_sp_balanced  1 when the balanced path applies to fc_sp_* rows
// Inputs:
//   Xstd     partialled-out regressors as returned by partial_out (n x K0,
//            standardized; the pre-solve columns with status 0)
//   status0  indepvar_status before the solve (1 + K_full), status1 after
//   stdevs   sol.stdevs after the solve (1 + K), means likewise (standardized)
//   resid    sol.resid after the solve (unstandardized)
//   w        regression weights as solve_ols builds them (fw raw, aw/pw
//            normalized to N, 1 if none), or J(0,1,.)
//   tmp_N    as in reghdfe_extend_b_and_xx
// ---------------------------------------------------------------------------
void fastconley_prepare(real matrix Xstd, real rowvector status0,
                        real rowvector status1, real rowvector stdevs,
                        real rowvector means, real colvector resid,
                        real colvector w, real scalar report_constant,
                        real scalar tmp_N,
                        real colvector lat, real colvector lon,
                        real colvector time, real colvector unit,
                        real scalar balanced, real scalar pixel, real scalar verbose)
{
	real rowvector ok0, keep, means_x, side
	real matrix X, inv_xx
	real colvector e
	real scalar K, n, corner
	external real matrix fc_D, fc_S
	external real colvector fc_time, fc_unit
	external real scalar fc_kk

	n = rows(Xstd)
	ok0 = selectindex(status0[|2 \ cols(status0)|] :== 0)
	keep = (cols(ok0) ? status1[ok0 :+ 1] :== 0 : J(1, 0, .))
	X = cols(ok0) ? select(Xstd, keep) : J(n, 0, .)
	K = cols(X)
	if (K != cols(stdevs) - 1) _error(3498, "fastconley: internal error (kept regressors do not match stdevs)")
	if (K) X = X :* stdevs[|2 \ K+1|]

	// Bread: reghdfe_rmcoll's inverse, extended for _cons as in
	// reghdfe_extend_b_and_xx (partitioned inverse using the means).
	if (rows(w)) inv_xx = K ? invsym(quadcross(X, w, X)) : J(0, 0, .)
	else inv_xx = K ? invsym(quadcross(X, X)) : J(0, 0, .)
	if (report_constant) {
		means_x = K ? means[|2 \ K+1|] :* stdevs[|2 \ K+1|] : J(1, 0, .)
		corner = (1 / tmp_N) + (K ? means_x * inv_xx * means_x' : 0)
		side = K ? -means_x * inv_xx : J(1, 0, .)
		fc_D = (inv_xx, side' \ side, corner)
		X = (K ? X :+ means_x : J(n, 0, .)), J(n, 1, 1)
	}
	else fc_D = inv_xx
	fc_kk = cols(X)
	if (fc_kk == 0) _error(3498, "fastconley: no coefficients to compute a variance for")

	// Scores s_i = (w_i e_i) x_i, exactly the weights reghdfe's robust meat uses.
	e = rows(w) ? resid :* w : resid
	fc_S = X :* e
	fc_time = time
	fc_unit = unit

	fastconley_prepare_rows(lat, lon, time, unit, balanced, pixel, verbose)
}

// Step 1b (shared by OLS and IV): balanced-panel validation, aggregation of
// identical (time, lat, lon) keys, and the spatial-engine input rows, from
// the globals fc_S / fc_time / fc_unit set by the caller.
void fastconley_prepare_rows(real colvector lat, real colvector lon,
                             real colvector time, real colvector unit,
                             real scalar balanced, real scalar pixel, real scalar verbose)
{
	real colvector p, alat, alon, atime, ucheck, uprev
	real scalar n, T, n_per, t
	external real matrix fc_S, fc_sp_S
	external real colvector fc_sp_lat, fc_sp_lon, fc_sp_time
	external real scalar fc_sp_balanced

	n = rows(fc_S)
	// Balanced-panel validation (same conditions and messages as the R package).
	if (balanced) {
		p = order((time, unit), (1, 2))
		T = rows(uniqrows(time))
		if (T > 1) {
			atime = time[p]
			ucheck = unit[p]
			n_per = n / T
			if (n_per != floor(n_per) | any(panelsetup(atime, 1)[., 2] - panelsetup(atime, 1)[., 1] :+ 1 :!= n_per))
				_error(3498, "balanced requires each period to have the same number of observations")
			uprev = ucheck[|1 \ n_per|]
			if (rows(uniqrows(uprev)) != n_per)
				_error(3498, "balanced requires each unit to appear at most once per period")
			for (t = 2; t <= T; t++) {
				if (any(ucheck[|(t-1)*n_per+1 \ t*n_per|] :!= uprev))
					_error(3498, "balanced requires every period to contain the same set of units")
			}
			alat = lat[p]; alon = lon[p]
			for (t = 2; t <= T; t++) {
				if (any(alat[|(t-1)*n_per+1 \ t*n_per|] :!= alat[|1 \ n_per|]) |
				    any(alon[|(t-1)*n_per+1 \ t*n_per|] :!= alon[|1 \ n_per|]))
					_error(3498, "balanced requires time-invariant coordinates per unit")
			}
		}
	}

	// Pre-aggregation of exact-duplicate (time, lat, lon) keys (+ pixel snap).
	fc_sp_lat = lat; fc_sp_lon = lon; fc_sp_time = time; fc_sp_S = fc_S
	fastconley_aggregate(fc_sp_lat, fc_sp_lon, fc_sp_time, fc_sp_S, pixel)
	if (verbose & rows(fc_sp_S) < n) printf("{txt}   - score pre-aggregation: %g rows -> %g cases (pixel = %g km)\n", n, rows(fc_sp_S), pixel)
	fc_sp_balanced = 0
	if (balanced & rows(uniqrows(fc_sp_time)) > 1) {
		T = rows(uniqrows(fc_sp_time))
		n_per = rows(fc_sp_S) / T
		if (n_per == floor(n_per) & !any(panelsetup(fc_sp_time, 1)[., 2] - panelsetup(fc_sp_time, 1)[., 1] :+ 1 :!= n_per)) {
			fc_sp_balanced = 1
		}
		else {
			// Aggregation broke the balance (R warns and uses the raw rows).
			if (verbose) printf("{txt}   - pre-aggregation produced different case counts per period; using unaggregated rows\n")
			p = order((time, unit), (1, 2))
			fc_sp_lat = lat[p]; fc_sp_lon = lon[p]; fc_sp_time = time[p]; fc_sp_S = fc_S[p, .]
			fc_sp_balanced = 1
		}
	}
}

// ---------------------------------------------------------------------------
// IV / 2SLS. `data` holds the partialled-out, unstandardized columns named
// in `names` (depvar first). X = (exog, endog), Z = (exog, instruments):
//   Xhat = Z (Z'WZ)^-1 Z'WX,  b = (Xhat'WX)^-1 Xhat'Wy,  e = y - X b,
//   bread D = (Xhat'WX)^-1,   scores s_i = (w_i e_i) xhat_i,
// which is the projected-design sandwich the R package uses for lfe/fixest
// IV fits. Regressors or instruments dropped as collinear with the fixed
// effects are simply absent from `names`; an absent endogenous regressor is
// an error. Fills fc_D, fc_S, fc_kk, fc_b, fc_resid, fc_rss, fc_tss_within
// and the string vector fc_xnames.
// ---------------------------------------------------------------------------
void fastconley_iv_prepare(real matrix data, string rowvector names,
                           string rowvector exog, string rowvector endog,
                           string rowvector inst, real colvector w, real scalar verbose)
{
	real colvector y, e
	real matrix X, Z, ZZ, iZZ, Xhat, XhX, iXhX
	real rowvector xi, zi, keepz
	real scalar j, k
	string scalar v
	external real matrix fc_D, fc_S
	external real colvector fc_b, fc_resid
	external real scalar fc_kk, fc_rss, fc_tss_within
	external string rowvector fc_xnames

	y = data[., 1]
	xi = J(1, 0, .)
	zi = J(1, 0, .)
	fc_xnames = J(1, 0, "")
	for (j = 1; j <= cols(exog); j++) {
		v = exog[j]
		k = fastconley_findname(names, v)
		if (k) {
			xi = (xi, k)
			zi = (zi, k)
			fc_xnames = (fc_xnames, v)
		}
		else if (verbose) printf("{txt}   - exogenous regressor %s dropped (collinear with the fixed effects)\n", v)
	}
	for (j = 1; j <= cols(endog); j++) {
		v = endog[j]
		k = fastconley_findname(names, v)
		if (!k) _error(3498, "endogenous regressor " + v + " is collinear with the fixed effects")
		xi = (xi, k)
		fc_xnames = (fc_xnames, v)
	}
	for (j = 1; j <= cols(inst); j++) {
		v = inst[j]
		k = fastconley_findname(names, v)
		if (k) zi = (zi, k)
		else if (verbose) printf("{txt}   - instrument %s dropped (collinear with the fixed effects)\n", v)
	}
	if (!cols(xi)) _error(3498, "no regressors left")
	if (cols(zi) < cols(xi)) _error(3498, "underidentified: fewer instruments than regressors")
	X = data[., xi]
	Z = data[., zi]

	ZZ = rows(w) ? quadcross(Z, w, Z) : quadcross(Z, Z)
	iZZ = invsym(ZZ)
	if (diag0cnt(iZZ)) {
		// drop collinear instruments (e.g. duplicates); the projection is unchanged
		keepz = selectindex(diagonal(iZZ)' :!= 0)
		if (verbose) printf("{txt}   - %g collinear instrument(s) dropped\n", cols(Z) - cols(keepz))
		Z = Z[., keepz]
		ZZ = rows(w) ? quadcross(Z, w, Z) : quadcross(Z, Z)
		iZZ = invsym(ZZ)
		if (cols(Z) < cols(X)) _error(3498, "underidentified after dropping collinear instruments")
	}
	Xhat = Z * (iZZ * (rows(w) ? quadcross(Z, w, X) : quadcross(Z, X)))
	XhX = rows(w) ? quadcross(Xhat, w, X) : quadcross(Xhat, X)
	iXhX = invsym(XhX)
	if (diag0cnt(iXhX)) _error(3498, "collinear regressors or weak identification; cannot invert Xhat'X")
	fc_b = iXhX * (rows(w) ? quadcross(Xhat, w, y) : quadcross(Xhat, y))
	fc_resid = y - X * fc_b
	fc_D = (iXhX + iXhX') / 2
	e = rows(w) ? fc_resid :* w : fc_resid
	fc_S = Xhat :* e
	fc_kk = cols(X)
	fc_rss = rows(w) ? quadcross(fc_resid, w, fc_resid) : quadcross(fc_resid, fc_resid)
	fc_tss_within = rows(w) ? quadcross(y, w, y) : quadcross(y, y)
}

real scalar fastconley_findname(string rowvector names, string scalar v)
{
	real scalar j
	for (j = 1; j <= cols(names); j++) if (names[j] == v) return(j)
	return(0)
}

// Step 2 (Mata engine): spatial meat on the prepared rows, stacking periods
// side by side when the balanced path applies.
real matrix fastconley_meat_mata(real colvector lat, real colvector lon,
                                 real colvector time, real matrix S,
                                 real scalar balanced, real scalar cutoff,
                                 string scalar kernel, string scalar dist,
                                 real scalar tile, real scalar verbose)
{
	real scalar T, n_per, kk, t
	real matrix Sstack
	real colvector plat, plon, ptime
	T = 1
	plat = lat; plon = lon; ptime = time; Sstack = S
	if (balanced & rows(uniqrows(time)) > 1) {
		T = rows(uniqrows(time))
		n_per = rows(S) / T
		kk = cols(S)
		Sstack = J(n_per, T * kk, .)
		for (t = 1; t <= T; t++) Sstack[|1, (t-1)*kk+1 \ n_per, t*kk|] = S[|(t-1)*n_per+1, 1 \ t*n_per, kk|]
		plat = lat[|1 \ n_per|]; plon = lon[|1 \ n_per|]; ptime = J(n_per, 1, 1)
	}
	if (verbose) printf("{txt}# Conley spatial meat (Mata engine, %s kernel, %s distance, cutoff %g km%s)\n",
	                    kernel, dist, cutoff, T > 1 ? sprintf(", %g stacked periods", T) : "")
	return(fastconley_spatial_meat(plat, plon, ptime, Sstack, T, cutoff, kernel, dist, tile, verbose))
}

// Step 3: sandwich, symmetrize, small-sample factor, PSD fix.
real matrix fastconley_assemble(real matrix D, real matrix meat, real scalar dof_adj, real scalar psdfix)
{
	real matrix V
	real scalar fixed
	V = D * meat * D * dof_adj
	V = (V + V') / 2
	fixed = fastconley_psd_fix(V, psdfix)
	st_local("fc_psd_noticeable", strofreal(fixed))
	return(V)
}

// Wald F for the joint significance of the slopes, from the new V (copy of
// the block after the VCE dispatch in reghdfe_solve_ols).
real scalar fastconley_wald_F(real colvector b, real matrix V, real scalar K, real scalar df_m)
{
	real rowvector idx
	real matrix inv_V
	if (!K) return(.)
	idx = 1..K
	inv_V = invsym(V[idx, idx])
	if (diag0cnt(inv_V)) return(.)
	return(b[idx]' * inv_V * b[idx] / df_m)
}

end
