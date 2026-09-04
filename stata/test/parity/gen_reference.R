# Cross-language parity harness, step 1 (R side).
#   Rscript stata/test/parity/gen_reference.R OUT_DIR
# Writes seeded datasets as .dta plus the fastconley R package's Conley vcov
# (fixest fit, ssc = FALSE, psd_fix = FALSE, tight demeaning) for each
# configuration in configs.csv. Step 2 is run_stata.do, step 3 compare.R.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(fastconley); library(fixest); library(haven) })
set.seed(20260904)

mk_xy <- function(n) { x1 <- rnorm(n); x2 <- rnorm(n)
  data.frame(x1 = x1, x2 = x2, y = 0.5 * x1 - 0.3 * x2 + rnorm(n)) }

# cs: cross-section, one-way FE
n <- 3000
cs <- cbind(mk_xy(n), lat = runif(n, 25, 49), lon = runif(n, -124, -67),
            region = sample(1:5, n, TRUE), w = runif(n, 0.5, 2))
# world: dateline + poles, one-way FE
n <- 4000
world <- cbind(mk_xy(n), lat = c(runif(n * 0.6, -85, 85), runif(n * 0.4, 60, 89.9)),
               lon = runif(n, -180, 180), region = sample(1:4, n, TRUE))
# bal: balanced panel, two-way FE
n_u <- 400; T_ <- 5
bal <- data.frame(unit = rep(seq_len(n_u), each = T_), time = rep(seq_len(T_), n_u),
                  lat = rep(runif(n_u, 35, 60), each = T_), lon = rep(runif(n_u, -10, 30), each = T_))
bal <- cbind(bal, mk_xy(nrow(bal)))
# unbal: 20% rows dropped (no singleton units by construction: T = 5, drop <= 2 per unit)
keep <- ave(runif(nrow(bal)), bal$unit, FUN = function(u) rank(u) > 2 | runif(length(u)) > 0.5)
unbal <- bal[keep == 1, ]

# ras: regular 0.5-degree lattice with holes (raster engine), one-way FE
ras <- expand.grid(lat = 35 + 0.5 * (0:59), lon = -20 + 0.5 * (0:79))
ras <- ras[-sample(nrow(ras), round(0.25 * nrow(ras))), ]
ras <- cbind(ras, mk_xy(nrow(ras)), region = sample(1:3, nrow(ras), TRUE))
# wrap: 1-degree lattice spanning the full longitude circle (dateline wrap)
wrap <- expand.grid(lat = 60 + (0:11), lon = -180 + (0:359))
wrap <- cbind(wrap, mk_xy(nrow(wrap)), region = sample(1:3, nrow(wrap), TRUE))
# IV: x2 endogenous with instruments z1 z2 (cs), xe with instrument z (bal)
n <- nrow(cs); cs$z1 <- rnorm(n); cs$z2 <- rnorm(n); u <- rnorm(n)
cs$x2 <- 0.7 * cs$z1 - 0.4 * cs$z2 + 0.5 * u + rnorm(n); cs$y <- 0.5 * cs$x1 - 0.3 * cs$x2 + u + rnorm(n)
n <- nrow(bal); bal$z <- rnorm(n); u <- rnorm(n)
bal$xe <- 0.6 * bal$z + 0.5 * u + rnorm(n); bal$y <- 0.5 * bal$x1 - 0.3 * bal$xe + u + rnorm(n)
write_dta(cs, file.path(out, "cs.dta")); write_dta(world, file.path(out, "world.dta"))
write_dta(ras, file.path(out, "ras.dta")); write_dta(wrap, file.path(out, "wrap.dta"))
write_dta(bal, file.path(out, "bal.dta")); write_dta(unbal, file.path(out, "unbal.dta"))

f_cs    <- feols(y ~ x1 + x2 | region, data = cs, demeaned = TRUE, fixef.tol = 1e-11)
f_csw   <- feols(y ~ x1 + x2 | region, data = cs, weights = ~w, demeaned = TRUE, fixef.tol = 1e-11)
f_world <- feols(y ~ x1 + x2 | region, data = world, demeaned = TRUE, fixef.tol = 1e-11)
f_ras   <- feols(y ~ x1 + x2 | region, data = ras, demeaned = TRUE, fixef.tol = 1e-11)
f_wrap  <- feols(y ~ x1 + x2 | region, data = wrap, demeaned = TRUE, fixef.tol = 1e-11)
f_bal   <- feols(y ~ x1 + x2 | unit + time, data = bal, demeaned = TRUE, fixef.tol = 1e-11)
f_csiv  <- feols(y ~ x1 | region | x2 ~ z1 + z2, data = cs, demeaned = TRUE, fixef.tol = 1e-11)
f_csivw <- feols(y ~ x1 | region | x2 ~ z1 + z2, data = cs, weights = ~w, demeaned = TRUE, fixef.tol = 1e-11)
f_baliv <- feols(y ~ x1 | unit + time | xe ~ z, data = bal, demeaned = TRUE, fixef.tol = 1e-11)
f_unbal <- feols(y ~ x1 + x2 | unit + time, data = unbal, demeaned = TRUE, fixef.tol = 1e-11)

# name | data | fixest fit | R args | Stata options (stata-side options string)
cfg <- list()
addc <- function(name, data, fit, rargs, sopts) cfg[[length(cfg) + 1]] <<- list(name = name, data = data, fit = fit, rargs = rargs, sopts = sopts)
for (k in c("bartlett", "uniform")) for (d in c("haversine", "spherical", "chord")) {
  addc(sprintf("cs_%s_%s", k, d), "cs", "f_cs",
       list(kernel = k, dist_fn = d, dist_cutoff = 300),
       sprintf("absorb(region) lat(lat) lon(lon) cutoff(300) kernel(%s) dist(%s)", k, d))
}
addc("cs_bartlett_haversine_pixel25", "cs", "f_cs",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 300, pixel = 25),
     "absorb(region) lat(lat) lon(lon) cutoff(300) pixel(25)")
addc("cs_bartlett_haversine_cut3000_tile200", "cs", "f_cs",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 3000),
     "absorb(region) lat(lat) lon(lon) cutoff(3000) tile(200)")
addc("csw_bartlett_haversine_aw", "cs", "f_csw",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 300),
     "[aw=w] absorb(region) lat(lat) lon(lon) cutoff(300)")
addc("world_bartlett_haversine_500", "world", "f_world",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 500),
     "absorb(region) lat(lat) lon(lon) cutoff(500)")
addc("world_uniform_chord_1500", "world", "f_world",
     list(kernel = "uniform", dist_fn = "chord", dist_cutoff = 1500),
     "absorb(region) lat(lat) lon(lon) cutoff(1500) kernel(uniform) dist(chord)")
# raster engine configurations are meaningful for the plugin only; the Mata
# engine runs them through its pairwise path (same answer, slower).
addc("ras_grid_uniform_spherical", "ras", "f_ras",
     list(kernel = "uniform", dist_fn = "spherical", dist_cutoff = 250, method = "grid"),
     "absorb(region) lat(lat) lon(lon) cutoff(250) kernel(uniform) dist(spherical) method(auto)")
addc("ras_grid_bartlett_haversine", "ras", "f_ras",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 250, method = "grid"),
     "absorb(region) lat(lat) lon(lon) cutoff(250) method(auto)")
addc("wrap_grid_bartlett_spherical", "wrap", "f_wrap",
     list(kernel = "bartlett", dist_fn = "spherical", dist_cutoff = 1500, method = "grid"),
     "absorb(region) lat(lat) lon(lon) cutoff(1500) dist(spherical) method(auto)")
# IV: Stata orders coefficients (exog, endog); fixest puts fit_ first, so the
# reference is permuted to (x1, fit_x2) below (see `iv_order`).
addc("csiv_bartlett_haversine", "cs", "f_csiv",
     list(kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 300),
     "absorb(region) lat(lat) lon(lon) cutoff(300)")
addc("csiv_uniform_chord_aw", "cs", "f_csivw",
     list(kernel = "uniform", dist_fn = "chord", dist_cutoff = 300),
     "[aw=w] absorb(region) lat(lat) lon(lon) cutoff(300) kernel(uniform) dist(chord)")
addc("baliv_bartlett_spherical_lag2_balanced", "bal", "f_baliv",
     list(unit = "unit", time = "time", kernel = "bartlett", dist_fn = "spherical", dist_cutoff = 300,
          lag_cutoff = 2, balanced_pnl = TRUE),
     "absorb(unit time) lat(lat) lon(lon) cutoff(300) unit(unit) time(time) lag(2) balanced dist(spherical)")
addc("bal_bartlett_haversine_lag0", "bal", "f_bal",
     list(unit = "unit", time = "time", kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 300),
     "absorb(unit time) lat(lat) lon(lon) cutoff(300) unit(unit) time(time)")
addc("bal_bartlett_haversine_lag2_balanced", "bal", "f_bal",
     list(unit = "unit", time = "time", kernel = "bartlett", dist_fn = "haversine", dist_cutoff = 300,
          lag_cutoff = 2, balanced_pnl = TRUE),
     "absorb(unit time) lat(lat) lon(lon) cutoff(300) unit(unit) time(time) lag(2) balanced")
addc("bal_uniform_chord_lag1", "bal", "f_bal",
     list(unit = "unit", time = "time", kernel = "uniform", dist_fn = "chord", dist_cutoff = 300, lag_cutoff = 1),
     "absorb(unit time) lat(lat) lon(lon) cutoff(300) unit(unit) time(time) lag(1) kernel(uniform) dist(chord)")
addc("unbal_bartlett_spherical_lag2", "unbal", "f_unbal",
     list(unit = "unit", time = "time", kernel = "bartlett", dist_fn = "spherical", dist_cutoff = 300, lag_cutoff = 2),
     "absorb(unit time) lat(lat) lon(lon) cutoff(300) unit(unit) time(time) lag(2) dist(spherical)")

rows <- list()
for (cc in cfg) {
  fit <- get(cc$fit); d <- get(cc$data)
  V <- do.call(vcovSpHAC, c(list(fit, lat = "lat", lon = "lon", ssc = FALSE, psd_fix = FALSE,
                                 ncores = 2, data = d), cc$rargs))
  b <- coef(fit)
  if (any(grepl("^fit_", names(b)))) {   # IV: reorder to (exog, endog)
    iv_order <- c(which(!grepl("^fit_", names(b))), which(grepl("^fit_", names(b))))
    V <- V[iv_order, iv_order]; b <- b[iv_order]
  }
  V <- unname(V)
  write.table(format(V, digits = 17), file.path(out, paste0(cc$name, "_R.csv")),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(format(unname(b), digits = 17), file.path(out, paste0(cc$name, "_Rb.csv")),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  rows[[length(rows) + 1]] <- data.frame(name = cc$name, data = cc$data, sopts = cc$sopts)
}
write.csv(do.call(rbind, rows), file.path(out, "configs.csv"), row.names = FALSE)
cat(length(cfg), "configurations written to", out, "\n")
