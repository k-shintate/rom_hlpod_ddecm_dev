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

    if [[ -f "$dst" && ! -f "${dst}.before_fluid_st_ddrom" ]]; then
        cp "$dst" "${dst}.before_fluid_st_ddrom"
    fi

    cp "$src" "$dst"
    echo "installed: $dst"
}

install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c
install_one solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.h
install_one solvers/fluid_sups/fluid_sups_core/rom_std_stddrom_hlpod_bridge.c
install_one solvers/fluid_sups/fluid_sups_core/rom_std_stddrom_hlpod_bridge.h
install_one solvers/fluid_sups/fluid_sups_core/rom_std_stddrom_fluid.c
install_one solvers/fluid_sups/fluid_sups_core/rom_std_stddrom_fluid.h
install_one solvers/fluid_sups/fluid_sups_main_cavity/main_ST_DDROM.c
install_one solvers/fluid_sups/Makefile_ST_DDROM.fragment
install_one shell/fluid_sups_hrom_cavity/patch_makefile_ST_DDROM.py
install_one shell/fluid_sups_hrom_cavity/run_ST_DDROM.sh

python3 \
    "${repo_root}/shell/fluid_sups_hrom_cavity/patch_makefile_ST_DDROM.py" \
    "${repo_root}"

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o"
echo "  rm -f fluid_sups_core/rom_std_stddrom_hlpod_bridge.o"
echo "  rm -f fluid_sups_core/rom_std_stddrom_fluid.o"
echo "  rm -f fluid_sups_main_cavity/main_ST_DDROM.o"
echo "  rm -f fluid_sups_st_ddrom"
echo "  make -f Makefile_ST_cavity dev_build"
