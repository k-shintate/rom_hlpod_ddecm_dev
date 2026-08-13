#!/bin/bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2' ERR

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 CASE_DIRECTORY [NP]" >&2
    exit 2
fi

case_dir="$(cd "$1" && pwd -P)"
np="${2:-1}"
slabs="${CONVDIFF_STFOM_SLABS_PER_WINDOW:-4}"
tol="${CONVDIFF_STFOM_COMPARE_TOL:-1.0e-10}"
write_vtk="${CONVDIFF_STFOM_WRITE_COMPARE_VTK:-1}"

if [ ! -x ./convdiff_stfom_dg0_validation ]; then
    echo "ERROR: ./convdiff_stfom_dg0_validation does not exist." >&2
    echo "Add the validation target to Makefile_ST before using this script." >&2
    exit 1
fi

env \
  CONVDIFF_STFOM_SLABS_PER_WINDOW="${slabs}" \
  CONVDIFF_STFOM_COMPARE_TOL="${tol}" \
  CONVDIFF_STFOM_WRITE_COMPARE_VTK="${write_vtk}" \
  mpirun -np "${np}" ./convdiff_stfom_dg0_validation "${case_dir}" \
  2>&1 | tee convdiff_stfom_dg0_validation.log
