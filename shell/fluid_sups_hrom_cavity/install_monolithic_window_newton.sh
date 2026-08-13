#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

copy_one() {
    local src="$1"
    local dst="$2"

    if [ ! -f "${src}" ]; then
        echo "ERROR: package source not found: ${src}" >&2
        exit 2
    fi

    mkdir -p "$(dirname "${dst}")"
    if [ -f "${dst}" ] && [ ! -f "${dst}.before_monolithic_window_newton" ]; then
        cp "${dst}" "${dst}.before_monolithic_window_newton"
    fi
    cp "${src}" "${dst}"
    echo "installed: ${dst}"
}

copy_one \
  "${package_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c" \
  "${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"

copy_one \
  "${package_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h" \
  "${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h"

copy_one \
  "${package_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c" \
  "${repo_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_FOM.c"

echo
echo "Rebuild from scratch:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  make -f Makefile_ST_cavity clean"
echo "  make -f Makefile_ST_cavity dev_build"
