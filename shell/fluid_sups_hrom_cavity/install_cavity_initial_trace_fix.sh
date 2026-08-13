#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
package_root="$(cd "$(dirname "$0")/../.." && pwd)"

src="${package_root}/shell/fluid_sups_hrom_cavity/run_ST_DDROM.sh"
dst="${repo_root}/shell/fluid_sups_hrom_cavity/run_ST_DDROM.sh"

[[ -f "$src" ]] || { echo "ERROR: missing package runner: $src" >&2; exit 2; }
[[ -f "$dst" ]] || { echo "ERROR: missing repository runner: $dst" >&2; exit 2; }

if [[ ! -f "${dst}.before_cavity_initial_trace_fix" ]]; then
    cp "$dst" "${dst}.before_cavity_initial_trace_fix"
fi

cp "$src" "$dst"
chmod +x "$dst"

echo "installed: $dst"
echo "default FLUID_ST_INITIAL_CONDITION is now cavity"
