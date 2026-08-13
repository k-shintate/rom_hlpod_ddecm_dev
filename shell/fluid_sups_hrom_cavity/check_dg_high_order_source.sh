#!/bin/bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
main="${repo_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c"

echo "Higher-order dG implementation markers:"
grep -n 'fluid_sups_st_solve_window_all_at_once_dg(' "${src}"
grep -n 'fluid_time_basis_derivative' "${src}"
grep -n 'psi_beta J_current' "${src}"
grep -n 'ALL-AT-ONCE dG' "${src}"
grep -n 'FLUID_ST_TIME_DEGREE' "${main}"
echo
echo "J_inv_ip occurrences (must be none):"
grep -n 'J_inv_ip' "${src}" || true
