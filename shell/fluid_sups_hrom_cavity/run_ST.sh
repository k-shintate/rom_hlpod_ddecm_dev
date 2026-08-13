#!/bin/bash
set -Eeuo pipefail

trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2; exit "$rc"' ERR

# mesh
e=20
ep=1

# Existing result-directory tags
num_modes=(10)
num_1stdd=(1)
num_parallel=(1)
pa=0
st=3

export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"
export FLUID_ST_SLABS_PER_WINDOW="${FLUID_ST_SLABS_PER_WINDOW:-4}"
export FLUID_ST_WRITE_WINDOW_BINARY="${FLUID_ST_WRITE_WINDOW_BINARY:-1}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
cd "${project_root}"

for nm in "${num_modes[@]}"; do
    for nd in "${num_1stdd[@]}"; do
        for np in "${num_parallel[@]}"; do
            echo
            echo "============================================================"
            echo "Fluid ST-FOM case"
            echo "  mesh divisions        : ${e}"
            echo "  domain size           : ${ep}"
            echo "  MPI ranks             : ${np}"
            echo "  slabs/window          : ${FLUID_ST_SLABS_PER_WINDOW}"
            echo "  initial condition     : ${FLUID_ST_INITIAL_CONDITION}"
            echo "  result tag nm,np,nd   : ${nm},${np},${nd}"
            echo "============================================================"

            # Execute child scripts. Do not source them.
            bash shell/fluid_sups_hrom_cavity/meshgen.sh \
                "${e}" "${ep}" "${nm}" "${nd}" "${np}" "${pa}"

            bash shell/fluid_sups_hrom_cavity/execution_ST.sh \
                "${e}" "${ep}" "${nm}" "${nd}" "${np}" "${pa}" "${st}"
        done
    done
done

echo
echo "Fluid ST-FOM driver completed."
