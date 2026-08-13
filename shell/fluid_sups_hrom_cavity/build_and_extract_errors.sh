#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}/solvers/fluid_sups"

log="${TMPDIR:-/tmp}/fluid_st_dev_build.$$.log"

set +e
make -f Makefile_ST_cavity clean
make -f Makefile_ST_cavity dev_build 2>&1 | tee "${log}"
rc=${PIPESTATUS[0]}
set -e

echo
echo "============================================================"
echo "Compiler/linker failure summary"
echo "log: ${log}"
echo "make exit: ${rc}"
echo "============================================================"

grep -n -C 4 -E \
'(^|[[:space:]])(fatal error:|error:)|undefined symbols|Undefined symbols|ld:|collect2: error|clang: error|g\+\+: error|cc1plus:.*error' \
"${log}" || true

echo
echo "Warnings about old Accelerate CBLAS interface:"
grep -n -E 'deprecated: An updated CBLAS interface|deprecated-declarations' "${log}" | head -n 30 || true

echo
echo "Actual core_NR compile command:"
grep -n 'core_NR.c -o fluid_sups_core/core_NR.o' "${log}" | head -n 5 || true

exit "${rc}"
