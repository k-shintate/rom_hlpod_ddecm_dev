#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

install_one() {
    local rel="$1"
    local src="${package_root}/${rel}"
    local dst="${repo_root}/${rel}"

    [[ -f "$src" ]] || {
        echo "ERROR: package source missing: $src" >&2
        exit 2
    }

    mkdir -p "$(dirname "$dst")"

    if [[ -f "$dst" && ! -f "${dst}.before_gedatsu_st_partition" ]]; then
        cp "$dst" "${dst}.before_gedatsu_st_partition"
    fi

    cp "$src" "$dst"
    echo "installed: $dst"
}

install_one solvers/fluid_sups/fluid_sups_core/st_fom_common.c
install_one solvers/fluid_sups/fluid_sups_core/st_fom_common.h
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h
install_one solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c

install_one shell/fluid_sups_hrom_cavity/execution_ST_GEDATSU_MPI.sh
install_one shell/fluid_sups_hrom_cavity/run_ST_GEDATSU_MPI.sh
install_one shell/fluid_sups_hrom_cavity/run_ST_parallel.sh
install_one shell/fluid_sups_hrom_cavity/check_GEDATSU_partition.sh

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/st_fom_common.o"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o"
echo "  rm -f fluid_sups_main_cavity/main_ST_FOM.o"
echo "  rm -f fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
