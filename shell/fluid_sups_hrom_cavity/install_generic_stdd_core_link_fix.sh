#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

src_core="${package_root}/solvers/fluid_sups/fluid_sups_core/rom_std_stddrom.c"
src_make="${package_root}/solvers/fluid_sups/Makefile_ST_cavity"

dst_core="${repo_root}/solvers/fluid_sups/fluid_sups_core/rom_std_stddrom.c"
dst_make="${repo_root}/solvers/fluid_sups/Makefile_ST_cavity"

for f in "$src_core" "$src_make"; do
    [[ -f "$f" ]] || {
        echo "ERROR: package file missing: $f" >&2
        exit 2
    }
done

[[ -f "$dst_make" ]] || {
    echo "ERROR: repository Makefile not found: $dst_make" >&2
    exit 2
}

if [[ ! -f "${dst_make}.before_generic_stdd_core_link" ]]; then
    cp "$dst_make" "${dst_make}.before_generic_stdd_core_link"
fi

if [[ -f "$dst_core" && ! -f "${dst_core}.before_generic_stdd_core_link" ]]; then
    cp "$dst_core" "${dst_core}.before_generic_stdd_core_link"
fi

cp "$src_core" "$dst_core"
cp "$src_make" "$dst_make"

echo "installed generic ST-DDROM implementation:"
echo "  $dst_core"
echo "installed Makefile:"
echo "  $dst_make"
echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/rom_std_stddrom.o"
echo "  rm -f fluid_sups_core/rom_std_stddrom_hlpod_bridge.o"
echo "  rm -f fluid_sups_core/rom_std_stddrom_fluid.o"
echo "  rm -f fluid_sups_main_cavity/main_ST_DDROM.o"
echo "  rm -f fluid_sups_st_ddrom"
echo "  make -f Makefile_ST_cavity st_ddrom"
