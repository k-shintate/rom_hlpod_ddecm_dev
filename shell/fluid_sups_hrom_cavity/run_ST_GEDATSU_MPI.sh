#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
cd "$project_root"

e="${FLUID_ST_MESH_DIVISIONS:-20}"
ep="${FLUID_ST_DOMAIN_SIZE:-1}"
nm="${FLUID_ST_RESULT_MODES_TAG:-10}"
nd="${FLUID_ST_RESULT_DD_TAG:-1}"
np="${FLUID_ST_MPI_RANKS:-4}"
pa="${FLUID_ST_PA_TAG:-0}"
st="${FLUID_ST_STAGE_TAG:-3}"

exec bash shell/fluid_sups_hrom_cavity/execution_ST_GEDATSU_MPI.sh \
    "$e" "$ep" "$nm" "$nd" "$np" "$pa" "$st"
