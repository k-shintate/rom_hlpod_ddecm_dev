#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

src="${package_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_DDROM.c"
dst="${repo_root}/solvers/fluid_sups/fluid_sups_main_cavity/main_ST_DDROM.c"

[[ -f "$src" ]] || {
    echo "ERROR: package source missing: $src" >&2
    exit 2
}
[[ -f "$dst" ]] || {
    echo "ERROR: repository target missing: $dst" >&2
    exit 2
}

if [[ ! -f "${dst}.before_training_window_fix" ]]; then
    cp "$dst" "${dst}.before_training_window_fix"
fi

cp "$src" "$dst"

echo "installed: $dst"
echo
echo "Rebuild only the changed driver and ST-DDROM executable:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_main_cavity/main_ST_DDROM.o fluid_sups_st_ddrom"
echo "  make -f Makefile_ST_cavity st_ddrom"
