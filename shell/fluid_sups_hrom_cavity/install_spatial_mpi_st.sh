#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

install_one() {
    local rel="$1"
    local src="${package_root}/${rel}"
    local dst="${repo_root}/${rel}"

    if [ ! -f "${src}" ]; then
        echo "ERROR: package source not found: ${src}" >&2
        exit 2
    fi

    mkdir -p "$(dirname "${dst}")"

    if [ -f "${dst}" ] && \
       [ ! -f "${dst}.before_spatial_mpi_st" ]; then
        cp "${dst}" "${dst}.before_spatial_mpi_st"
    fi

    cp "${src}" "${dst}"
    echo "installed: ${dst}"
}

install_one solvers/fluid_sups/fluid_sups_core/st_fom_common.c
install_one solvers/fluid_sups/fluid_sups_core/st_fom_common.h
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h
install_one solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c
install_one shell/fluid_sups_hrom_cavity/run_ST_parallel.sh

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/st_fom_common.o"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o"
echo "  rm -f fluid_sups_main_cavity/main_ST_FOM.o"
echo "  rm -f fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
