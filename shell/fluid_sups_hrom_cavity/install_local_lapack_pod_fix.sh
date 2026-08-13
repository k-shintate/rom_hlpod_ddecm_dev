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
    [[ -e "$dst" ]] || mkdir -p "$(dirname "$dst")"

    if [[ -f "$dst" && ! -f "${dst}.before_local_lapack_pod" ]]; then
        cp "$dst" "${dst}.before_local_lapack_pod"
    fi
    cp "$src" "$dst"
    echo "installed: $dst"
}

install_one solvers/fluid_sups/fluid_sups_core/rom_std_stddrom.c
install_one solvers/fluid_sups/Makefile_ST_cavity
install_one shell/fluid_sups_hrom_cavity/check_local_lapack_pod.sh

echo
echo "Rebuild only the affected DDROM objects first:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/rom_std_stddrom.o fluid_sups_st_ddrom"
echo "  make -f Makefile_ST_cavity st_ddrom"
