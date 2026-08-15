#!/usr/bin/env bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
pkg="$(cd "$(dirname "$0")" && pwd)"
dst="${repo_root}/shell/fluid_sups_hrom_cavity"
stamp="$(date +%Y%m%d_%H%M%S)"
backup="${repo_root}/refactor_backups/case_workflow_v4_${stamp}"
mkdir -p "$backup"

for name in st_pipeline.sh meshgen.sh run_validated_case.sh; do
    [[ -f "$dst/$name" ]] && cp "$dst/$name" "$backup/$name"
done

for name in st_pipeline.sh meshgen.sh run_validated_case.sh run_validated_existing.sh check_ST_case.sh; do
    cp "$pkg/shell/fluid_sups_hrom_cavity/$name" "$dst/$name"
    chmod +x "$dst/$name"
done

echo "installed ST validated-case workflow v4"
echo "backup: $backup"
