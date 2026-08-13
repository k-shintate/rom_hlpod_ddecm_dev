#!/bin/bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
hdr="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_jold.h"

echo "Analytic J_old markers:"
grep -n 'BBFE_elemmat_fluid_sups_mat_NR_old' "${hdr}" "${src}"
grep -n 'FLUID_ST_JOLD_MODE' "${src}"
grep -n 'analytic J_old verification PASS' "${src}"
echo
echo "Old dedicated dG0 FD block should be absent:"
grep -n 'fluid_add_window_previous_jacobian_block_fd' "${src}" || true
