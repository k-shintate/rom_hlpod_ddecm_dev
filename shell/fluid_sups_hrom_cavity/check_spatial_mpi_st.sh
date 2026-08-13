#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
model="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
header="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h"
common="${repo_root}/solvers/fluid_sups/fluid_sups_core/st_fom_common.c"

echo "[1] spatial graph commonization"
grep -n 'ST_fom_spatial_graph_build_from_connectivity' "${common}"

echo
echo "[2] window pattern is lifted from the spatial graph"
grep -n 'monolis_get_nonzero_pattern_by_nodal_graph_R' "${model}"
grep -n 'monolis_window.mat.N = owned' "${model}"
grep -n 'monolis_window.mat.NP = local_nodes' "${model}"

echo
echo "[3] window must reuse base spatial communicator"
if grep -R -n 'mono_com_window' "${model}" "${header}"; then
    echo "ERROR: stale mono_com_window reference remains." >&2
    exit 3
else
    echo "OK: no mono_com_window references."
fi

echo
echo "[4] solve/update communicator"
grep -n -B2 -A4 'BBFE_sys_monowrap_solve' "${model}" | \
    grep -E 'monolis_window|mono_com' | head -n 20 || true
grep -n -A3 'monolis_mpi_update_R(' "${model}" | \
    grep -E 'mono_com|window_dof' | head -n 30 || true

echo
echo "STATUS: source-level spatial-MPI wiring looks present."
