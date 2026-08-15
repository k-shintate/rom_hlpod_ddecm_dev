#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

src="${package_root}/solvers/fluid_sups/Makefile_ST_cavity"
dst="${repo_root}/solvers/fluid_sups/Makefile_ST_cavity"

[[ -f "$src" ]] || { echo "ERROR: missing package Makefile: $src" >&2; exit 2; }
[[ -f "$dst" ]] || { echo "ERROR: missing repository Makefile: $dst" >&2; exit 2; }

if [[ ! -f "${dst}.before_linux_link_fix" ]]; then
    cp "$dst" "${dst}.before_linux_link_fix"
fi

cp "$src" "$dst"

echo "installed: $dst"
echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  make -f Makefile_ST_cavity clean_st"
echo "  make -f Makefile_ST_cavity verify_refactor"
echo "  make -f Makefile_ST_cavity dev_build"
