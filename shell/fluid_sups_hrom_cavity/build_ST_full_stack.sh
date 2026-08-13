#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"

cd "${repo_root}"

bash shell/fluid_sups_hrom_cavity/patch_core_HROM_compile.sh "${repo_root}"

cd solvers/fluid_sups

echo
echo "=== Compile full future fluid stack + ST-FOM ==="
make -f Makefile_ST_cavity clean
make -f Makefile_ST_cavity dev_build

echo
echo "Build succeeded."
echo "ST-FOM executable:"
ls -lh fluid_sups_st_fom
