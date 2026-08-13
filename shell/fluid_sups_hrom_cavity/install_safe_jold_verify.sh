#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
pkg_root="$(cd "$(dirname "$0")/../.." && pwd)"
src="${pkg_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"
dst="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_st_model.c"

if [ ! -f "${src}" ] || [ ! -f "${dst}" ]; then
    echo "ERROR: source or target fluid_sups_st_model.c not found" >&2
    exit 2
fi

if [ ! -f "${dst}.before_safe_jold_verify" ]; then
    cp "${dst}" "${dst}.before_safe_jold_verify"
fi
cp "${src}" "${dst}"

echo "installed safe J_old verification model: ${dst}"
echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
