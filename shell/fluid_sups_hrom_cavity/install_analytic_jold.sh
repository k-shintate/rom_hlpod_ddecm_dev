#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

install_one() {
    local rel="$1"
    local src="${package_root}/${rel}"
    local dst="${repo_root}/${rel}"

    if [ ! -f "${src}" ]; then
        echo "ERROR: package source missing: ${src}" >&2
        exit 2
    fi

    mkdir -p "$(dirname "${dst}")"

    if [ -f "${dst}" ] && \
       [ ! -f "${dst}.before_analytic_jold" ]; then
        cp "${dst}" "${dst}.before_analytic_jold"
    fi

    cp "${src}" "${dst}"
    echo "installed: ${dst}"
}

install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_jold.h
install_one solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
