#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}/solvers/fluid_sups"

src="fluid_sups_core/rom_std_stddrom.c"
obj="fluid_sups_core/rom_std_stddrom.o"
exe="fluid_sups_st_ddrom"

if grep -q 'monolis_scalapack_gesvd_R' "$src"; then
    echo "ERROR: stale ScaLAPACK POD call remains in $src" >&2
    exit 1
fi

echo "OK: generic ST-DDROM POD uses local LAPACK dgesvd."

echo "Checking unresolved ScaLAPACK/BLACS references in generic object..."
if [[ -f "$obj" ]]; then
    if nm -u "$obj" 2>/dev/null | grep -Eiq 'scalapack|blacs'; then
        nm -u "$obj" 2>/dev/null | grep -Ei 'scalapack|blacs' || true
        echo "ERROR: generic object still depends on ScaLAPACK/BLACS." >&2
        exit 1
    fi
    echo "OK: no ScaLAPACK/BLACS dependency in $obj."
else
    echo "NOTE: $obj has not been rebuilt yet."
fi

if [[ -x "$exe" ]]; then
    echo "ST-DDROM executable exists:"
    ls -lh "$exe"
fi
