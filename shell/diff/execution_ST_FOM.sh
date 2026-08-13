#!/bin/bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2' ERR

if [ "$#" -lt 7 ]; then
    echo "Usage: $0 e ep nm nd np pa st" >&2
    exit 2
fi

e="$1"
ep="$2"
nm="$3"
nd="$4"
np="$5"
pa="$6"
st="$7"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd "${script_dir}/../.." && pwd -P)"
directory="${project_root}/result_diff/${nm}-${np}-${nd}"

fom_mpi_np="${ST_FOM_MPI_NP:-1}"
write_output="${ST_DDROM_WRITE_OUTPUT:-${ST_ROM_WRITE_OUTPUT:-1}}"

# install.shの副作用はこの子bash内だけに閉じ込める。
source "${project_root}/shell/install.sh"

mkdir -p "${directory}"
mkdir -p \
    "${directory}/pod_modes_vtk" \
    "${directory}/pod_modes" \
    "${directory}/fem_solver_prm" \
    "${directory}/pod_solver_prm" \
    "${directory}/calctime" \
    "${directory}/DDECM"
for ((i=0; i<nd; i++)); do
    mkdir -p "${directory}/pod_modes/subdomain${i}"
done

cd "${project_root}/solvers/diff"
make -f Makefile_ST clean
make -f Makefile_ST

if [ ! -x ./diff_ST_MPI ]; then
    echo "ERROR: diff_ST_MPI was not created." >&2
    exit 1
fi

cp -f ./diff_ST_MPI "${directory}/"
cd "${directory}"

cat > stfom_run_info.txt <<INFO
mode=fom
mpi_np=${fom_mpi_np}
e=${e}
ep=${ep}
nm_label=${nm}
nd=${nd}
np_label=${np}
pa=${pa}
st=${st}
ST_DDROM_WRITE_OUTPUT=${write_output}
ST_ROM_WRITE_OUTPUT=${write_output}
INFO

echo
echo "Running FOM-only windowed ST calculation"
echo "  ST_DDROM_MODE    : fom"
echo "  ST_ROM_MODE      : fom (compatibility)"
echo "  MPI processes    : ${fom_mpi_np}"
echo "  output           : ${write_output}"
echo "  executable       : ${directory}/diff_ST_MPI"
echo

# IMPORTANT:
# ST_DDROM_MODE=fom / ST_ROM_MODE=fom は、どちらのdriverがMakefile_STで
# 選択されていてもsnapshot/POD/offline/onlineを実行しない純粋FOM経路。
# dG degree / slabs per window はdriverに応じて
# ST_DDROM_DEGREE / ST_DDROM_WINDOW_SLABS または
# ST_ROM_DEGREE / ST_ROM_WINDOW_SLABS のcompile-time macroで決まる。
env \
    ST_DDROM_MODE=fom \
    ST_DDROM_WRITE_OUTPUT="${write_output}" \
    ST_ROM_MODE=fom \
    ST_ROM_WRITE_OUTPUT="${write_output}" \
    mpirun -np "${fom_mpi_np}" ./diff_ST_MPI \
    2>&1 | tee stfom_only.log

echo
echo "ST-FOM only run completed successfully."
echo "Log: ${directory}/stfom_only.log"
