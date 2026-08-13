#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}/solvers/fluid_sups"

obj="./fluid_sups_core/rom_std_stddrom.o"
exe="./fluid_sups_st_ddrom"

[[ -f "$obj" ]] || {
    echo "ERROR: missing object: $obj" >&2
    exit 1
}

symbols=(
    ROM_std_stdd_system_initialize
    ROM_std_stdd_system_finalize
    ROM_std_stdd_compute_shared_local_pod
    ROM_std_stdd_write_basis_rank
    ROM_std_stdd_read_basis_rank
    ROM_std_stdd_project_operators_rank_sweep
    ROM_std_stdd_project_local_rhs
    ROM_std_stdd_solve_window_marching
)

echo "Checking generic ST-DDROM object symbols..."
for s in "${symbols[@]}"; do
    if nm -g "$obj" 2>/dev/null | grep -q "$s"; then
        echo "  OK: $s"
    else
        echo "ERROR: missing symbol in $obj: $s" >&2
        exit 1
    fi
done

if [[ -x "$exe" ]]; then
    echo
    echo "ST-DDROM executable exists:"
    ls -lh "$exe"
else
    echo
    echo "WARNING: $exe not present yet; build st_ddrom first."
fi
