#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)
out_dir=${1:-"$repo_root/.scratch/stata-parity"}
stata_log="$repo_root/run_stata.log"

cleanup() {
    rm -f "$stata_log"
}
trap cleanup EXIT

mkdir -p "$out_dir"
cd "$repo_root"
Rscript stata/test/parity/gen_reference.R "$out_dir"

for engine in mata plugin; do
    rm -f "$stata_log"
    FASTCONLEY_PARITY_DIR="$out_dir" FASTCONLEY_PARITY_ENGINE="$engine" \
        stata-mp -b do stata/test/parity/run_stata.do
    grep -F "run_stata.do: 21 configurations written (engine $engine)" "$stata_log"
    Rscript stata/test/parity/compare.R "$out_dir"
done

echo "run_parity.sh: Mata and plugin parity checks passed"
