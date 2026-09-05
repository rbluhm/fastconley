#!/usr/bin/env python3
"""Render stata/PERFORMANCE.md from the Stata benchmark results.

    python3 stata/bench/render_performance.py [results_dir] [out_md]

Reads <results_dir>/bench.csv and session.txt (written by run_bench.sh) and,
for the same-machine comparison, the R package's shipped benchmark files
inst/benchmarks/fastconley-benchmark-results.csv and
inst/benchmarks/fastconley-large-data-results.csv.
"""
import csv
import os
import sys
from collections import OrderedDict

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
res_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "stata", "bench", "results")
out_md = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "stata", "PERFORMANCE.md")


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt_s(x):
    v = num(x)
    if v is None:
        return "n/a"
    if v >= 100:
        return f"{v:,.0f}"
    if v >= 10:
        return f"{v:.1f}"
    return f"{v:.3f}" if v < 1 else f"{v:.2f}"


def fmt_e(x):
    v = num(x)
    return "n/a" if v is None else (f"{v:.1e}" if v != 0 else "0")


def fmt_x(a, b):
    a, b = num(a), num(b)
    return "n/a" if a is None or b is None or b == 0 else f"{a / b:.0f}x" if a / b >= 10 else f"{a / b:.1f}x"


rows = read_csv(os.path.join(res_dir, "bench.csv"))
r_rows = read_csv(os.path.join(root, "inst", "benchmarks", "fastconley-benchmark-results.csv"))
r_large = read_csv(os.path.join(root, "inst", "benchmarks", "fastconley-large-data-results.csv"))
session = open(os.path.join(res_dir, "session.txt")).read() if os.path.exists(os.path.join(res_dir, "session.txt")) else ""


def sess(key, default=""):
    for line in session.splitlines():
        if line.strip().startswith(key):
            rest = line.strip()[len(key):]
            return rest.lstrip(":= ").strip()
    return default


def first(rs, **kw):
    for r in rs:
        if all(str(r.get(k, "")) == str(v) for k, v in kw.items()):
            return r
    return None


def stata_time(r):
    """Covariance-step seconds for fastconley rows, whole-command seconds for acreg."""
    if r is None:
        return None
    return r["total_seconds"] if r["method"] == "acreg" else r["vce_seconds"]


def r_time(section, method, **kw):
    r = first([x for x in r_rows if x["section"] == section and x["method"] == method], **kw)
    return None if r is None else r["seconds"]


meta = rows[0] if rows else {}
threads_all = sorted({int(r["threads"]) for r in rows if r["engine"] == "plugin"})
acreg_version = next((r["acreg_version"] for r in rows if r["method"] == "acreg" and r["acreg_version"]), "")

L = []
L.append("# fastconley for Stata: performance\n")
L.append("This is the Stata counterpart of the R package's performance vignette. The same six benchmark "
         "sections are run through the `fastconley` command with the compiled plugin at several thread "
         "counts and with the pure-Mata fallback, and, where its O(n²) loop is feasible, through "
         "acreg (Colella, Lalive, Sakalli and Thoenig, 2019, \"Inference with Arbitrary Clustering\", IZA "
         "Discussion Paper 12584; on SSC), the established Stata implementation of Conley standard errors. The tables were produced by "
         "`stata/bench/run_bench.sh` on the machine below and rendered by `stata/bench/render_performance.py`; "
         "they are not regenerated automatically.\n")
L.append("## What is timed\n")
L.append("- **fastconley (plugin, Mata)**: `e(vce_seconds)`, the covariance step alone, after reghdfe has "
         "partialled out the fixed effects and solved the regression. This is the same quantity the R vignette "
         "reports for `vcovSpHAC()` (post-estimation only). The whole-command time including reghdfe is in "
         "the results file (`total_seconds`) and is typically 0.03 to 0.05 s longer on the small cases and a "
         "few seconds longer at one million rows.")
L.append("- **acreg**: the whole command, because acreg estimates and corrects in one pass and exposes no "
         "separate timer. On these configurations its regression is a negligible part of the total.")
L.append("- **R (same machine)**: the numbers shipped with the R package (`inst/benchmarks/`), run on the same "
         "CPU on 2026-06-10 with fastconley 0.8.0; the engine has since become faster (chord-form weights), so "
         "today's R package is at least as fast as the R column shows.")
L.append("- Every fastconley call uses `nossc nopsdfix` so that its covariance is comparable to acreg, which "
         "applies no small-sample correction. The relative-difference columns compare the slope block of `e(V)` "
         "with the plugin result on the same data (`mreldif`).\n")
L.append("## Machine\n")
L.append("| item | value |\n|---|---|")
L.append(f"| run date | {sess('run_date')} |")
L.append(f"| CPU | {sess('Model name')} |")
L.append(f"| logical CPUs | {sess('CPU(s)')} |")
L.append(f"| memory | {' '.join(session.split('Mem:')[1].split()[:1]) + ' GB' if 'Mem:' in session else ''} |")
L.append(f"| OS | {sess('PRETTY_NAME').strip(chr(34))} |")
L.append(f"| Stata | {meta.get('stata_version', '')} |")
L.append(f"| fastconley ado / engine build | {meta.get('ado_version', '')} / {meta.get('engine_build', '')} |")
L.append(f"| acreg | {acreg_version} |")
L.append(f"| plugin thread counts | {', '.join(str(t) for t in threads_all)} |")
L.append("")

# ---------------------------------------------------------------- dense
L.append("## Dense baseline\n")
L.append("Small cross-sections (5 regressors, 500 km, uniform kernel, spherical distance) where the R vignette "
         "checks the engine against an explicit dense weight matrix. acreg runs comfortably here.\n")
L.append("| observations | plugin (8 thr.) | Mata | acreg | R fastconley (8 thr.) | R dense | Mata vs plugin | acreg vs plugin |")
L.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
for n in (1000, 2000, 4000):
    p = first(rows, section="dense", engine="plugin", n_obs=str(n))
    m = first(rows, section="dense", engine="mata", n_obs=str(n))
    a = first(rows, section="dense", method="acreg", n_obs=str(n))
    L.append(f"| {n:,} | {fmt_s(stata_time(p))} | {fmt_s(stata_time(m))} | {fmt_s(stata_time(a))} | "
             f"{fmt_s(r_time('dense', 'fastconley', n_obs=str(n)))} | {fmt_s(r_time('dense', 'dense', n_obs=str(n)))} | "
             f"{fmt_e(m and m['rel_diff_vs_plugin'])} | {fmt_e(a and a['rel_diff_vs_plugin'])} |")
L.append("\nThe acreg column differs from fastconley at the 1e-5 level because acreg measures distance on an "
         "equirectangular plane (111 km per degree of latitude, scaled by the cosine of the latitude for "
         "longitude) rather than on the sphere; a few pairs near the cutoff boundary change status.\n")

# ---------------------------------------------------------------- cross-section
L.append("## Scattered cross-sections\n")
L.append("10 regressors, uniform kernel, spherical distance, points drawn uniformly over the contiguous United "
         "States. The R columns are the vignette's `vcovSpHAC()` and `fixest::vcov_conley()` times at the same "
         "thread count where the vignette ran one (1 and 8 threads).\n")
L.append("| observations | cutoff km | threads | plugin | Mata (1 thr.) | acreg | R fastconley | R fixest | plugin vs Mata |")
L.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for n, c in ((50000, 100), (50000, 500), (100000, 100), (100000, 500)):
    m = first(rows, section="cross_section", engine="mata", n_obs=str(n), cutoff_km=str(c))
    a = first(rows, section="cross_section", method="acreg", n_obs=str(n), cutoff_km=str(c))
    for t in threads_all:
        p = first(rows, section="cross_section", engine="plugin", n_obs=str(n), cutoff_km=str(c), threads=str(t))
        if p is None:
            continue
        rf = r_time("cross_section", "fastconley", n_obs=str(n), cutoff_km=str(c), threads=str(t))
        rx = r_time("cross_section", "fixest", n_obs=str(n), cutoff_km=str(c), threads=str(t))
        L.append(f"| {n:,} | {c} | {t} | {fmt_s(stata_time(p))} | {fmt_s(stata_time(m)) if t == threads_all[0] else ''} | "
                 f"{(fmt_s(stata_time(a)) if a else 'not run') if t == threads_all[0] else ''} | {fmt_s(rf)} | {fmt_s(rx)} | "
                 f"{fmt_x(stata_time(m), stata_time(p))} |")
acreg_notes = sorted({r["notes"] for r in rows if r["method"] == "acreg" and r["notes"].startswith("aborted")})
L.append("")
if acreg_notes:
    L.append("acreg configurations marked *not run* were not attempted or were stopped: " + "; ".join(acreg_notes) + ".")
L.append("acreg's Mata loop touches every pair of observations once per observation, so its time grows with "
         "n² regardless of the cutoff, while fastconley's cell grid only visits candidate pairs within the "
         "cutoff. The Mata fallback of fastconley uses the same cell grid, so it stays proportional to the "
         "number of pairs but runs single-threaded in interpreted Mata.\n")

# ---------------------------------------------------------------- panel
L.append("## Balanced panel with serial HAC\n")
L.append("10,000 units observed in 4 periods (40,000 rows), unit and time fixed effects absorbed, 5 regressors, "
         "500 km, uniform kernel, and a one-lag serial Bartlett term (`lag(1) balanced`). fastconley follows the "
         "Hsiang (2010) convention used by the R package: contemporaneous spatial correlation within each period "
         "plus own-unit serial correlation. acreg's `hac` option instead applies the product of the temporal "
         "Bartlett weight and the spatial indicator to every cross-unit pair, a different estimator, so the two "
         "are timed but not compared numerically.\n")
L.append("| method | threads | seconds | R (8 thr.) |")
L.append("|---|---:|---:|---:|")
for t in threads_all:
    p = first(rows, section="panel", engine="plugin", threads=str(t))
    if p:
        L.append(f"| plugin | {t} | {fmt_s(stata_time(p))} | {fmt_s(r_time('panel', 'fastconley')) if t == 8 else ''} |")
m = first(rows, section="panel", engine="mata")
if m:
    L.append(f"| Mata | 1 | {fmt_s(stata_time(m))} | |")
a = first(rows, section="panel", method="acreg")
L.append(f"| acreg (`hac lag(1) pfe1 pfe2`) | 1 | {fmt_s(stata_time(a)) if a else 'not run'} | |")
L.append(f"| R fixest composition (conley + NW - hetero) | 8 | | {fmt_s(r_time('panel', 'fixest composition'))} |")
L.append("")

# ---------------------------------------------------------------- raster
L.append("## Regular raster\n")
L.append("A 180 × 180 latitude/longitude lattice (32,400 cells at 0.05°), 3 regressors, 250 km. The plugin has "
         "a dedicated grid engine (prefix sums for the uniform kernel, per-ring FFT convolutions for Bartlett) "
         "that the Mata fallback does not have; both are exact, so the grid and pairwise results agree to "
         "roundoff.\n")
L.append("| kernel | plugin grid (8 thr.) | plugin pairwise (8 thr.) | Mata pairwise | acreg | R grid | R pairwise | grid vs pairwise |")
L.append("|---|---:|---:|---:|---:|---:|---:|---:|")
for kern in ("uniform", "bartlett"):
    g = first(rows, section="grid", method="grid", engine="plugin", kernel=kern)
    pw = first(rows, section="grid", method="pairwise", engine="plugin", kernel=kern)
    m = first(rows, section="grid", engine="mata", kernel=kern)
    a = first(rows, section="grid", method="acreg", kernel=kern)
    L.append(f"| {kern} | {fmt_s(stata_time(g))} | {fmt_s(stata_time(pw))} | {fmt_s(stata_time(m))} | "
             f"{fmt_s(stata_time(a)) if a else 'not run'} | {fmt_s(r_time('grid', 'grid', kernel=kern))} | "
             f"{fmt_s(r_time('grid', 'pairwise', kernel=kern))} | {fmt_e(pw and pw['rel_diff_vs_plugin'])} |")
L.append("")

# ---------------------------------------------------------------- pixel
L.append("## Repeated locations and pixel aggregation\n")
L.append("100,000 rows on 20,000 distinct locations (5 rows each), 5 regressors, 250 km, Bartlett kernel. "
         "`pixel(0)` merges rows with identical coordinates exactly; `pixel(10)` and `pixel(25)` snap "
         "coordinates to a 10 km or 25 km lattice first, which is approximate. acreg has no aggregation and "
         "is not attempted at this size.\n")
L.append("| pixel km | plugin (8 thr.) | Mata | R fastconley (8 thr.) | vs exact (pixel 0) |")
L.append("|---:|---:|---:|---:|---:|")
for px in ("0", "10", "25"):
    p = first(rows, section="pixel", method=f"pixel={px}", engine="plugin")
    m = first(rows, section="pixel", method=f"pixel={px}", engine="mata")
    L.append(f"| {px} | {fmt_s(stata_time(p))} | {fmt_s(stata_time(m))} | "
             f"{fmt_s(r_time('pixel', f'pixel={px}'))} | {fmt_e(p and p['rel_diff_vs_plugin'])} |")
L.append("")

# ---------------------------------------------------------------- large
L.append("## One million observations\n")
large_n = next((r["n_obs"] for r in rows if r["section"] == "large_xsection"), "1000000")
L.append(f"A global cross-section of {int(float(large_n)):,} points (latitude -55 to 70, all longitudes), 10 "
         "regressors, 100 km, uniform kernel. A dense weight matrix would be about 7.3 TiB.\n")
L.append("| threads | plugin | Mata | R fastconley |")
L.append("|---:|---:|---:|---:|")
m = first(rows, section="large_xsection", engine="mata")
for t in threads_all:
    p = first(rows, section="large_xsection", engine="plugin", threads=str(t))
    if p is None:
        continue
    rl = first([x for x in r_large if x["section"] == "large_xsection"], threads=str(t))
    L.append(f"| {t} | {fmt_s(stata_time(p))} | {fmt_s(stata_time(m)) if t == threads_all[0] else ''} | "
             f"{fmt_s(rl['seconds']) if rl else 'n/a'} |")
if m and num(m["vce_seconds"]) is None:
    L.append(f"\nThe Mata fallback did not finish this case: {m['notes']}. Its cell grid keeps the pair count "
             "proportional to the data, but interpreted Mata pays a fixed cost per cell pair, and at one million "
             "points over the whole globe there are hundreds of thousands of occupied cells.")
elif m:
    L.append(f"\nThe Mata fallback agrees with the plugin to {fmt_e(m['rel_diff_vs_plugin'])} relative on this case.")
L.append("acreg is not attempted here: its cost grows with n², and its whole-command time at 100,000 points "
         "is already about 700 s.\n")

# ---------------------------------------------------------------- overhead
L.append("## Fixed preparation cost\n")
L.append("`cutoff(-1)` keeps only the diagonal of the meat, so `e(vce_seconds)` then measures everything except "
         "the pair work: reading the sample into Mata, sorting by period, merging identical coordinates, "
         "marshalling scores into temporary variables for the plugin, and assembling the sandwich.\n")
L.append("| observations | plugin (8 thr.) | Mata |")
L.append("|---:|---:|---:|")
for n in sorted({int(float(r["n_obs"])) for r in rows if r["section"] == "overhead"}):
    p = first(rows, section="overhead", engine="plugin", n_obs=str(n))
    m2 = first(rows, section="overhead", engine="mata", n_obs=str(n))
    L.append(f"| {n:,} | {fmt_s(stata_time(p))} | {fmt_s(stata_time(m2))} |")
L.append("\nAt one million rows this preparation is most of the covariance time reported above, which is why "
         "the plugin's thread scaling looks flat there: the engine itself needs well under a second at 16 "
         "threads, as the R vignette's 1.1 s on the same machine shows. Trimming this Mata-side work (skipping "
         "the coordinate merge when locations are unique, streaming scores to the plugin without temporary "
         "variables) is the obvious next optimisation for very large samples.\n")

# ---------------------------------------------------------------- how to read / reproduce
L.append("## Reading the numbers\n")
L.append("- The plugin is the same C++ engine as the R package, so plugin and R times differ only by the "
         "front-end (Stata tempvars versus R memory aliasing) and by the engine improvements since the R "
         "numbers were recorded.")
L.append("- Thread scaling flattens beyond 8 threads on this 12-core, 16-thread laptop (hybrid P/E cores and "
         "memory bandwidth), and at one million rows the fixed preparation cost dominates (see the previous section). "
         "Stata's own licence (MP with 4 cores here) does not limit the plugin's `threads()`.")
L.append("- The Mata fallback is 5 to 15 times slower than the single-threaded plugin but has the same "
         "complexity, so it remains usable at a few hundred thousand observations. It is what `engine(auto)` "
         "uses when no plugin is available for the platform.")
L.append("- acreg and fastconley agree to about 1e-5 on the uniform kernel and 1e-6 on Bartlett at these "
         "cutoffs; the residual is acreg's planar distance approximation, not a difference in the estimator.\n")
L.append("## Reproducing\n")
L.append("```bash\n# from the repository root; needs reghdfe, ftools, require, and acreg on the adopath\n"
         "bash stata/bench/run_bench.sh stata/bench/results     # ~1-3 hours, mostly acreg\n"
         "python3 stata/bench/render_performance.py stata/bench/results stata/PERFORMANCE.md\n```\n")
L.append("`run_bench.sh` accepts `THREADS`, `ACREG_MAX_N`, `ACREG_TIMEOUT`, `LARGE_N`, and `REPS` in the "
         "environment; the driver `stata/bench/bench_vignette.do` can also be run directly for one section "
         "with the `BENCH_*` globals it documents. Results are one CSV row per timed call with the machine and "
         "software versions attached.\n")

with open(out_md, "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"wrote {out_md} ({len(L)} lines) from {len(rows)} result rows")
