#!/usr/bin/env bash
# run_bench.sh -- run the Stata counterpart of the R vignette benchmarks.
#
#   bash stata/bench/run_bench.sh [out_dir]        (from the repository root)
#
# Each section runs in its own Stata process so that a timeout on one acreg
# configuration cannot take the rest down. acreg's O(n^2) Mata loop is only
# attempted up to ACREG_MAX_N observations (default 50000) and, when the
# 50,000-point run finishes within ACREG_PROMOTE_SECS, once more at 100,000.
# Results: <out_dir>/bench.csv (one row per timed call), per-section logs,
# and session.txt (machine and software versions). Render with
#   python3 stata/bench/render_performance.py <out_dir>
set -u
STATA=${STATA:-stata-mp}
OUT=${1:-stata/bench/results}
THREADS=${THREADS:-"1 4 8 16"}
ACREG_MAX_N=${ACREG_MAX_N:-50000}
ACREG_TIMEOUT=${ACREG_TIMEOUT:-2700}
ACREG_PROMOTE_SECS=${ACREG_PROMOTE_SECS:-900}
LARGE_N=${LARGE_N:-1000000}
LARGE_TIMEOUT=${LARGE_TIMEOUT:-5400}
REPS=${REPS:-2}
mkdir -p "$OUT"
CSV="$OUT/bench.csv"
rm -f "$CSV"

{
  echo "run_date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "threads: $THREADS"; echo "reps: $REPS"; echo "acreg_max_n: $ACREG_MAX_N"; echo "large_n: $LARGE_N"
  echo; echo "# CPU"; lscpu 2>/dev/null | grep -E "Model name|^CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|L3"
  echo; echo "# Memory"; free -g 2>/dev/null | head -2
  echo; echo "# OS"; uname -srm; grep PRETTY_NAME /etc/os-release 2>/dev/null
  echo; echo "# Stata"
  printf 'di c(stata_version) " " c(flavor) " licensed=" c(processors) " machine=" c(processors_mach)\n' > "$OUT/_stata_version.do"
  ( cd "$OUT" && "$STATA" -b do _stata_version.do >/dev/null 2>&1; grep -E "^[0-9]+" _stata_version.log | tail -1; rm -f _stata_version.do _stata_version.log )
} > "$OUT/session.txt"

run_section () {  # name methods threads timeout [extra global assignments...]
  local name=$1 methods=$2 threads=$3 tmo=$4; shift 4
  local do="$OUT/section_${name}_$(echo "$methods" | tr ' ' '-').do"
  {
    echo "global BENCH_SECTION $name"
    echo "global BENCH_METHODS \"$methods\""
    echo "global BENCH_THREADS \"$threads\""
    echo "global BENCH_OUT \"$CSV\""
    echo "global BENCH_REPS $REPS"
    echo "global BENCH_LARGE_N $LARGE_N"
    for g in "$@"; do echo "global $g"; done
    echo "cd \"$PWD\""
    echo "do stata/bench/bench_vignette.do"
  } > "$do"
  local t0=$(date +%s)
  echo "[$(date +%H:%M:%S)] $name ($methods) start"
  ( cd "$OUT" && timeout "$tmo" "$STATA" -b do "$(basename "$do")" )
  local rc=$? t1=$(date +%s)
  echo "[$(date +%H:%M:%S)] $name ($methods) rc=$rc elapsed=$((t1-t0))s"
  if [ $rc -eq 124 ]; then
    # the last method of the list is the one that was running when the cap hit
    local last=${methods##* }
    echo "$name,$name,${last},${last},1,NA,NA,NA,NA,NA,NA,NA,0,0,NA,NA,NA,\"aborted: exceeded ${tmo}s timeout ($*)\",\"$(date -u +%FT%TZ)\",\"\",\"\",\"\",\"\"" >> "$CSV"
  fi
  LAST_RC=$rc; LAST_ELAPSED=$((t1-t0))
}

# fastconley (plugin at every thread count + Mata) on every section first.
run_section dense    "plugin mata" "8"       1800
run_section xsection "plugin mata" "$THREADS" 3600
run_section panel    "plugin mata" "$THREADS" 3600
run_section raster   "plugin mata" "8"       3600
run_section pixel    "plugin mata" "8"       3600
run_section large    "plugin mata" "$THREADS" "$LARGE_TIMEOUT"
run_section overhead "plugin mata" "8"       1800

# acreg, smallest to largest, each under its own timeout.
run_section dense  "acreg" "8" 1800
run_section raster "acreg" "8" "$ACREG_TIMEOUT"
run_section panel  "acreg" "8" "$ACREG_TIMEOUT"
promote=1
for case in "50000:100" "50000:500"; do
  n=${case%%:*}
  [ "$n" -le "$ACREG_MAX_N" ] || continue
  run_section xsection "acreg" "8" "$ACREG_TIMEOUT" "BENCH_XS_CASES \"$case\""
  [ $LAST_RC -eq 0 ] && [ $LAST_ELAPSED -le "$ACREG_PROMOTE_SECS" ] || promote=0
done
if [ $promote -eq 1 ]; then
  for case in "100000:100" "100000:500"; do
    run_section xsection "acreg" "8" "$ACREG_TIMEOUT" "BENCH_XS_CASES \"$case\""
  done
else
  echo "acreg at 100,000 observations not attempted (50,000-point runs exceeded ${ACREG_PROMOTE_SECS}s or failed)"
fi
echo "[$(date +%H:%M:%S)] all sections done -> $CSV"
