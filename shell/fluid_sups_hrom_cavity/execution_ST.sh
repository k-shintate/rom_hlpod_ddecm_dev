#!/bin/bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2; exit "$rc"' ERR

if [ "$#" -lt 7 ]; then
    echo "Usage: $0 e ep nm nd np pa st" >&2
    exit 2
fi

e="$1"; ep="$2"; nm="$3"; nd="$4"; np="$5"; pa="$6"; st="$7"
: "${e}" "${ep}" "${pa}" "${st}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
directory="${project_root}/result_fluid_sups_cavity/${nm}-${np}-${nd}"
solver_dir="${project_root}/solvers/fluid_sups"

makefile="${FLUID_ST_MAKEFILE:-Makefile_ST_cavity}"
target="${FLUID_ST_TARGET:-fluid_sups_st_fom}"
executable="${FLUID_ST_EXECUTABLE:-fluid_sups_st_fom}"
slabs_per_window="${FLUID_ST_SLABS_PER_WINDOW:-4}"
initial_condition="${FLUID_ST_INITIAL_CONDITION:-cavity}"
write_window_binary="${FLUID_ST_WRITE_WINDOW_BINARY:-1}"
openblas_threads="${OPENBLAS_NUM_THREADS:-4}"

source "${project_root}/shell/install.sh"
mkdir -p "${directory}"

cd "${solver_dir}"

detector="${project_root}/shell/fluid_sups_hrom_cavity/detect_fluid_nr_api_linkage.sh"
nr_api_linkage="$(bash "${detector}" --mode ../../../test_thermal)"

if [ "${nr_api_linkage}" = "missing" ] || [ "${nr_api_linkage}" = "mixed" ]; then
    # The current fluid ST-FOM links the verified NR solve-path kernels from
    # fluid_sups_core/core_FOM_st_provider.o, not from libBBFE_elemmat.
    if [[ "${FLUID_ST_REQUIRE_BBFE_NR_API:-0}" == "1" ]]; then
        bash "${detector}" ../../../test_thermal
        exit 1
    else
        echo "Skipping legacy libBBFE_elemmat NR API detector; using local core_FOM ST provider."
        nr_api_linkage="local_core_FOM_cpp"
    fi
fi

echo
echo "Building fluid ST-FOM"
echo "  makefile       : ${makefile}"
echo "  target         : ${target}"
echo "  NR API linkage : ${nr_api_linkage}"

make -f "${makefile}" clean
make -f "${makefile}" NR_API_LINKAGE="${nr_api_linkage}" "${target}"

cp -f "./${executable}" "${directory}/${executable}"
cd "${directory}"

mkdir -p pod_modes_vtk pod_modes_v pod_modes_p fem_solver_prm pod_solver_prm \
    hr_solver_prm calctime DDECM hr_prm hot_start
for ((i=0; i<nd; i++)); do
    mkdir -p "pod_modes_v/subdomain${i}" "pod_modes_p/subdomain${i}"
done

window_dir="${FLUID_ST_WINDOW_DIR:-${directory}/fluid_st_windows_s${slabs_per_window}}"
mkdir -p "${window_dir}"
logfile="fluid_st_fom_np${np}_s${slabs_per_window}.log"

env \
    OPENBLAS_NUM_THREADS="${openblas_threads}" \
    OMP_NUM_THREADS=1 \
    FLUID_ST_INITIAL_CONDITION="${initial_condition}" \
    FLUID_ST_SLABS_PER_WINDOW="${slabs_per_window}" \
    FLUID_ST_WRITE_WINDOW_BINARY="${write_window_binary}" \
    FLUID_ST_WINDOW_DIR="${window_dir}" \
    mpirun -np "${np}" "./${executable}" ./ \
    2>&1 | tee "${logfile}"

echo
echo "Fluid ST-FOM completed successfully."
echo "Log: ${directory}/${logfile}"
cd "${project_root}"
