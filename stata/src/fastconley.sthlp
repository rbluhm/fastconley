{smcl}
{* *! version 0.1.0 03sep2026}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ivreghdfe" "help ivreghdfe"}{...}
{viewerjumpto "Syntax" "fastconley##syntax"}{...}
{viewerjumpto "Description" "fastconley##description"}{...}
{viewerjumpto "Options" "fastconley##options"}{...}
{viewerjumpto "Remarks" "fastconley##remarks"}{...}
{viewerjumpto "Examples" "fastconley##examples"}{...}
{viewerjumpto "Stored results" "fastconley##results"}{...}
{viewerjumpto "References" "fastconley##references"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:fastconley} {hline 2}}Conley (1999) spatial HAC standard errors for linear models with many fixed effects{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:fastconley}
{depvar} [{indepvars}]
[{cmd:(}{it:endogvars} {cmd:=} {it:instruments}{cmd:)}]
{ifin}
{it:{weight}}
{cmd:,}
{opth lat:itude(varname)}
{opth lon:gitude(varname)}
{opt cut:off(#)}
{{opth a:bsorb(varlist)}|{opt noa:bsorb}}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth a:bsorb(varlist)}}categorical variables to absorb, as in {help reghdfe}{p_end}
{synopt :{opt noa:bsorb}}no fixed effects (constant only){p_end}
{synopt :{opt nocons:tant}}do not report the constant{p_end}

{syntab:Conley kernel}
{synopt :{opth lat:itude(varname)}}latitude in decimal degrees{p_end}
{synopt :{opth lon:gitude(varname)}}longitude in decimal degrees{p_end}
{synopt :{opt cut:off(#)}}spatial cutoff in kilometres{p_end}
{synopt :{opt ker:nel(string)}}{opt bartlett} (default) or {opt uniform}{p_end}
{synopt :{opt dist:ance(string)}}{opt haversine} (default), {opt spherical}, or {opt chord}{p_end}
{synopt :{opth unit(varname)}}panel unit identifier{p_end}
{synopt :{opth time(varname)}}time identifier{p_end}
{synopt :{opt lag(#)}}serial (within-unit) Bartlett lag cutoff; default 0{p_end}
{synopt :{opt bal:anced}}the panel is balanced with time-invariant coordinates (faster){p_end}
{synopt :{opt pix:el(#)}}snap points to a #-km grid before computing (approximation; default 0){p_end}

{syntab:SE adjustments}
{synopt :{opt nossc}}no small-sample correction ({it:N}/({it:N}-{it:K}-{it:df_a}) is applied by default){p_end}
{synopt :{opt nopsd:fix}}do not clamp negative eigenvalues of the vcov{p_end}

{syntab:Engine}
{synopt :{opt eng:ine(string)}}{opt auto} (default), {opt mata}, or {opt plugin}{p_end}
{synopt :{opt meth:od(string)}}{opt auto} (default), {opt pairwise}, or {opt grid} (raster lattices; plugin only){p_end}
{synopt :{opt thr:eads(#)}}threads for the plugin engine; default {cmd:c(processors_mach)}{p_end}
{synopt :{opt neigh:bor(string)}}{opt grid} (default) or {opt band} candidate search (plugin only){p_end}
{synopt :{opt csr:weight(string)}}{opt double} (default) or {opt float} neighbour-weight storage on the balanced path (plugin only){p_end}
{synopt :{opt tile(#)}}tile size for the Mata engine's dense blocks; default 1024{p_end}
{synopt :{opt v:erbose}}report engine progress{p_end}

{syntab:Other}
{synopt :{opth res:iduals(newvar)}}save residuals{p_end}
{synopt :{it:reghdfe options}}any other option is passed to {help reghdfe} (e.g. {opt tol:erance()}, {opt keepsin:gletons}){p_end}
{synopt :{it:display options}}{opt nohead:er}, {opt notable}, {opt nofoot:note}, and {help estimation options##display_options:display options}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}{opt aweight}s, {opt fweight}s, and {opt pweight}s are allowed; see {help weight}.
Weights enter the scores as in {cmd:reghdfe}'s robust variance ({it:w}{it:e}{it:x}); with
{opt fweight}s an observation therefore stands for {it:w} identical, perfectly correlated
observations at the same location, so a negative cutoff reproduces {cmd:vce(robust)}
only for {opt aweight}s, {opt pweight}s, and unweighted fits.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables and time-series operators.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:fastconley} fits a linear regression with arbitrary absorbed fixed effects
(via {help reghdfe}) and reports Conley (1999) spatial HAC standard errors:
the variance of the estimated coefficients allows arbitrary correlation
between observations closer than {opt cutoff()} kilometres (weighted by the
chosen kernel) and, in panels, arbitrary serial correlation within units up
to {opt lag()} periods (Bartlett weights). It is a port of the
{browse "https://github.com/rbluhm/fastconley":fastconley} R package and
reproduces its {cmd:vcovSpHAC()} results on the same fit: same kernels,
distances, panel semantics, aggregation of coincident points, and small-sample
convention.

{pstd}
The estimation output, {cmd:e()} results, {cmd:predict}, and the absorbed-FE
footnote are those of {cmd:reghdfe}; only the variance matrix, {cmd:e(F)},
and {cmd:e(vcetype)} differ.

{pstd}
With {cmd:(}{it:endogvars} {cmd:=} {it:instruments}{cmd:)} the model is
estimated by two-stage least squares on the partialled-out variables
(the same 2SLS sandwich the R package applies to {cmd:lfe} and {cmd:fixest} IV
fits): the bread is the inverse cross-product of the projected regressors and
the scores are the structural residuals times the projected regressors.
Point estimates equal those of {help ivreghdfe} and {cmd:ivreg2}, and with a
negative cutoff the variance equals their {cmd:robust} variance. Factor-variable
terms are not allowed in the IV lists (absorb them or use indicators), no
constant is reported, and {cmd:predict} is not available after IV fits.


{marker options}{...}
{title:Options}

{dlgtab:Conley kernel}

{phang}
{opth latitude(varname)} and {opth longitude(varname)} give the coordinates in
decimal degrees. Observations with identical (time, latitude, longitude) are
aggregated before the spatial sum, which is exact.

{phang}
{opt cutoff(#)} is the distance in kilometres beyond which the spatial kernel
is zero. As an undocumented convenience for testing, a negative cutoff drops
every cross-observation term and reproduces {cmd:reghdfe, vce(robust)}.

{phang}
{opt kernel(bartlett)} weights a pair at distance {it:d} by 1 - {it:d}/cutoff;
{opt kernel(uniform)} weights every pair within the cutoff by 1.

{phang}
{opt distance()} chooses the great-circle formula: {opt haversine} (default),
{opt spherical} (arc-cosine), or {opt chord} (straight line through the
sphere), all with Earth radius 6371 km.

{phang}
{opth unit(varname)} and {opth time(varname)} define the panel. Spatial
correlation is allowed within each period; with {opt lag(#)} > 0, serial
correlation within each unit is also allowed for time differences up to #
(in the units of {it:time}, so {it:time} should be numeric and evenly spaced;
string identifiers are accepted for {it:unit}).

{phang}
{opt balanced} asserts that every period contains the same units with
time-invariant coordinates. The neighbour structure is then computed once and
reused across periods. The command verifies the assertion and errors if it
fails.

{phang}
{opt pixel(#)} snaps coordinates to a #-km grid before aggregating, which can
reduce the number of distinct points dramatically at the cost of approximating
distances by up to about #/2 km.

{dlgtab:SE adjustments}

{phang}
{opt nossc} disables the small-sample factor {it:N}/({it:N}-{it:K}-{it:df_a}),
where {it:df_a} counts the absorbed fixed effects as {cmd:reghdfe} does. The
default matches {cmd:reghdfe, vce(robust)} and the R package's
{cmd:ssc = TRUE}; {opt nossc} matches {cmd:ssc = FALSE}.

{phang}
{opt nopsdfix} disables the positive-semi-definite fix. Spatial kernels do not
guarantee a positive-semi-definite variance; by default negative eigenvalues
are clamped to 1e-16 and a note is printed when that changes the matrix
noticeably.

{dlgtab:Engine}

{phang}
{opt engine()} selects the computational engine. {opt plugin} is the compiled
engine shared with the R package (multithreaded, streaming, and the only one
with the raster engine); it is installed automatically on Linux, Windows, and
macOS and loaded once per session. {opt mata} is the pure-Mata implementation,
always available. {opt auto} (the default) uses the plugin when it loads and
its engine version matches, and Mata otherwise. Both engines give the same
answer to floating-point summation order; {cmd:e(engine)} records which ran.

{phang}
{opt method()} chooses between the pairwise engine and the exact grid engine
for observations on a regular latitude/longitude lattice (raster cell
centres). The grid engine's cost does not depend on the number of pairs, so
it wins on dense rasters with large cutoffs. {opt auto} detects a lattice and
uses a flop-balance rule; {opt grid} forces it and errors when no lattice is
detected; {opt pairwise} never uses it. Lattices that span the full
longitude circle wrap across the dateline. {cmd:e(method)} records the choice.
Plugin only. Lattice detection tolerates the rounding noise of coordinates
stored as {cmd:float}; store them as {cmd:double} for exact agreement between
the grid and pairwise engines.

{phang}
{opt threads(#)} sets the plugin's thread count (results do not depend on
it). The default is the machine's processor count, {cmd:c(processors_mach)},
not the licensed count: the plugin's threads are independent of the Stata
licence, so Stata/SE and /BE users get the full speed-up as well. {opt neighbor()} and {opt csrweight()} are the plugin's candidate-search
strategy and neighbour-weight precision, mirroring the R package's
{cmd:neighbor} and {cmd:csr_weight}; the defaults are exact and fastest.

{phang}
{opt tile(#)} sets the dense block size the Mata engine uses inside each pair
of neighbouring cells; larger tiles are faster but use more memory
(about 5 x #^2 x 8 bytes). Both engines compute Bartlett distances from the
chord between unit vectors, accurate to better than 1e-12 relative at any
cutoff; the Mata engine builds the chord from coordinate differences below
200 km, where the dot-product form would lose precision.


{marker remarks}{...}
{title:Remarks}

{pstd}
The variance is {it:V} = {it:D M D q}, with {it:D} the inverse of the
weighted cross-product of the demeaned regressors (extended for the reported
constant exactly as {cmd:reghdfe} does), {it:M} the kernel-weighted sum of
score cross-products, and {it:q} the small-sample factor. Scores fold the
regression weights into the residual as {cmd:reghdfe}'s robust variance does.

{pstd}
Correspondence with the R package ({cmd:vcovSpHAC()} arguments):
{opt cutoff()} = {cmd:dist_cutoff}, {opt kernel()} = {cmd:kernel},
{opt distance()} = {cmd:dist_fn}, {opt lag()} = {cmd:lag_cutoff},
{opt balanced} = {cmd:balanced_pnl = TRUE}, {opt pixel()} = {cmd:pixel},
{opt nossc} = {cmd:ssc = FALSE}, {opt nopsdfix} = {cmd:psd_fix = FALSE},
{opt method()} = {cmd:method}, {opt neighbor()} = {cmd:neighbor}, {opt csrweight()} = {cmd:csr_weight},
{opt threads()} = {cmd:ncores}.


{marker examples}{...}
{title:Examples}

{pstd}Cross-section, one absorbed fixed effect, 300 km Bartlett kernel{p_end}
{phang2}{cmd:. fastconley y x1 x2, absorb(region) lat(lat) lon(lon) cutoff(300)}{p_end}

{pstd}Panel with unit and time effects, 500 km uniform kernel and 2 serial lags{p_end}
{phang2}{cmd:. fastconley y x1 x2, absorb(id year) lat(lat) lon(lon) cutoff(500) kernel(uniform) unit(id) time(year) lag(2) balanced}{p_end}

{pstd}Weighted, no fixed effects, spherical distance{p_end}
{phang2}{cmd:. fastconley y x1 x2 [aw = pop], noabsorb lat(lat) lon(lon) cutoff(200) dist(spherical)}{p_end}

{pstd}Instrumental variables: x2 instrumented by z1 and z2{p_end}
{phang2}{cmd:. fastconley y x1 (x2 = z1 z2), absorb(region) lat(lat) lon(lon) cutoff(300)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:fastconley} stores everything {help reghdfe##results:reghdfe} stores, with
{cmd:e(cmd)} = {cmd:fastconley}, {cmd:e(vcetype)} = {cmd:Conley}, {cmd:e(vce)} = {cmd:conley}, and in addition:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(conley_cutoff)}}spatial cutoff (km){p_end}
{synopt:{cmd:e(conley_lag)}}serial lag cutoff{p_end}
{synopt:{cmd:e(conley_pixel)}}pixel size (km){p_end}
{synopt:{cmd:e(conley_balanced)}}1 if {opt balanced} was used{p_end}
{synopt:{cmd:e(ssc)}}1 if the small-sample factor was applied{p_end}
{synopt:{cmd:e(psd_fix)}}1 if the PSD fix was enabled{p_end}
{synopt:{cmd:e(dof_adj)}}small-sample factor applied{p_end}
{synopt:{cmd:e(threads)}}threads used (plugin only){p_end}
{synopt:{cmd:e(iv)}}1 for a 2SLS fit{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(conley_kernel)}}kernel{p_end}
{synopt:{cmd:e(conley_dist)}}distance function{p_end}
{synopt:{cmd:e(engine)}}engine used ({cmd:mata} or {cmd:plugin}){p_end}
{synopt:{cmd:e(method)}}spatial engine used ({cmd:pairwise} or {cmd:grid}){p_end}
{synopt:{cmd:e(engine_build)}}plugin build identifier (plugin only){p_end}
{synopt:{cmd:e(latvar)}}, {cmd:e(lonvar)}, {cmd:e(unitvar)}, {cmd:e(timevar)}}variables used{p_end}
{synopt:{cmd:e(instd)}}, {cmd:e(insts)}}instrumented regressors and instruments (IV fits; {cmd:e(title)} is "HDFE 2SLS regression"){p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Conley, T. G. 1999. GMM estimation with cross sectional dependence.
{it:Journal of Econometrics} 92: 1-45.

{phang}
Correia, S. 2017. Linear models with high-dimensional fixed effects: an
efficient and feasible estimator. Working paper.

{phang}
Hsiang, S. M. 2010. Temperatures and cyclones strongly associated with
economic production in the Caribbean and Central America.
{it:PNAS} 107: 15367-15372.


{title:Author}

{pstd}
Richard Bluhm. Issues and source: {browse "https://github.com/rbluhm/fastconley"}.
