#!/usr/bin/env bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
pkg="$(cd "$(dirname "$0")" && pwd)"
stamp="$(date +%Y%m%d_%H%M%S)"
backup_dir="${repo_root}/refactor_backups/verify_v3_${stamp}"
mkdir -p "$backup_dir"

solver="${repo_root}/solvers/fluid_sups"
cp "$solver/Makefile_ST_cavity" "$backup_dir/Makefile_ST_cavity"
cp "$solver/st_refactor_tools/check_refactor.py" "$backup_dir/check_refactor.py"

cp "$pkg/solvers/fluid_sups/Makefile_ST_cavity" "$solver/Makefile_ST_cavity"
cp "$pkg/solvers/fluid_sups/st_refactor_tools/check_refactor.py" \
   "$solver/st_refactor_tools/check_refactor.py"
chmod +x "$solver/st_refactor_tools/check_refactor.py"

legacy="$solver/Makefile_ST_cavity.before_linux_link_fix"
if [[ -f "$legacy" ]]; then
    mv "$legacy" "$backup_dir/Makefile_ST_cavity.before_linux_link_fix"
fi

echo "installed v3 verification fix"
echo "backup: $backup_dir"
