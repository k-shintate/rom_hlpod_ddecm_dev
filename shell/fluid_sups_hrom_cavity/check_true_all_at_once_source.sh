#!/bin/bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
main="${repo_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c"

echo "All-at-once implementation markers:"
grep -n 'fluid_sups_st_solve_window_all_at_once_dg0' "${src}"
grep -n 'fluid_add_window_previous_jacobian_block_fd' "${src}"
grep -n 'ONE MONOLIS solve' "${src}"
grep -n 'FLUID_ST_SOLVE_ALL_AT_ONCE' "${main}"
