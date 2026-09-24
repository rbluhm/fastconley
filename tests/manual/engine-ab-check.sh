#!/usr/bin/env bash
# Bit-identity A/B check of the working-tree engine header against the header
# at a git revision (default HEAD). See tests/manual/engine_ab_check.cpp.
#
#   tests/manual/engine-ab-check.sh [REV] [ITERATIONS]
#
# Run from the package root. Needs git, a C++14 compiler (CXX, default g++),
# and the Armadillo headers (ARMA_INC, default RcppArmadillo's include dir).
set -euo pipefail
rev="${1:-HEAD}"
iterations="${2:-300}"
cxx="${CXX:-g++}"
arma_inc="${ARMA_INC:-$(Rscript -e 'cat(system.file("include", package = "RcppArmadillo"))')}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git show "$rev:src/conley_core.h" | sed \
  -e 's/FASTCONLEY_CONLEY_CORE_H/FASTCONLEY_CONLEY_CORE_OLD_H/g' \
  -e 's/^namespace conley {/namespace conley_old {/' \
  -e 's/^} \/\/ namespace conley$/} \/\/ namespace conley_old/' \
  -e 's/CONLEY_CORE_VERSION/CONLEY_CORE_VERSION_OLD/g' > "$work/conley_core_old.h"
grep -q '^namespace conley_old {' "$work/conley_core_old.h"

"$cxx" -std=c++14 -O2 -pthread -DNDEBUG -DARMA_64BIT_WORD -DARMA_DONT_USE_WRAPPER \
  -DARMA_DONT_USE_BLAS -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU \
  -Isrc -I"$work" -I"$arma_inc" tests/manual/engine_ab_check.cpp -o "$work/engine_ab_check"
echo "engine_ab_check: working tree vs $rev ($(git rev-parse --short "$rev")), $iterations iterations, $cxx"
"$work/engine_ab_check" "$iterations"
