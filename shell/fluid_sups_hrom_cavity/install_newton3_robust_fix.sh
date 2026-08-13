#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"
src="${package_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
dst="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"

[[ -f "$src" ]] || { echo "ERROR: package source missing: $src" >&2; exit 2; }
[[ -f "$dst" ]] || { echo "ERROR: target missing: $dst" >&2; exit 2; }

if [[ ! -f "${dst}.before_newton3_robust_fix" ]]; then
    cp "$dst" "${dst}.before_newton3_robust_fix"
fi

cp "$src" "$dst"
echo "installed: $dst"

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
