#!/usr/bin/env bash
# Rank-controlled diffusion ST-DDROM sweep.
# Runs offline + online AOA + online window-marching for ranks 1,2,4,8 by default.

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
    echo "Usage: $0 <e> <ep> <nm> <nd> <np> <pa> <st>" >&2
    echo "Optional: DIFF_RANK_SWEEP=\"1 2 4 8\"" >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
compare="${ST_DDROM_COMPARE_LAUNCHER:-$script_dir/execution_ST_compare_AOA_vs_marching.sh}"
if [[ ! -x "$compare" ]]; then
    echo "ERROR: comparison launcher not found/executable: $compare" >&2
    exit 1
fi

find_project_root() {
    local candidate
    for candidate in \
        "${ST_PROJECT_ROOT:-}" \
        "$(pwd)" \
        "${script_dir}/../.." \
        "${script_dir}/../../.." \
        "${script_dir}/../../../.."; do
        [[ -n "$candidate" ]] || continue
        candidate="$(cd "$candidate" 2>/dev/null && pwd || true)"
        [[ -n "$candidate" ]] || continue
        if [[ -f "$candidate/shell/install.sh" && \
              -f "$candidate/solvers/diff/Makefile_ST" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

project_root="$(find_project_root || true)"
if [[ -z "$project_root" ]]; then
    echo "ERROR: project root not found; set ST_PROJECT_ROOT." >&2
    exit 1
fi

ranks="${DIFF_RANK_SWEEP:-1 2 4 8}"
base_nm="$3"
nd="$4"
np="$5"
user_pa="$6"
summary="${DIFF_RANK_SWEEP_SUMMARY:-$project_root/result_diff/stddrom_rank_sweep_${base_nm}-${np}-${nd}.csv}"
rm -f "$summary"
first=1
rebuild="${ST_FOM_REBUILD:-1}"

for rank in $ranks; do
    if ! [[ "$rank" =~ ^[0-9]+$ ]] || (( rank <= 0 )); then
        echo "ERROR: invalid rank in DIFF_RANK_SWEEP: $rank" >&2
        exit 2
    fi

    max_modes="$user_pa"
    if (( max_modes < rank )); then
        max_modes="$rank"
    fi

    result_dir="$project_root/result_diff/rank${rank}_${base_nm}-${np}-${nd}"

    echo
    echo "============================================================"
    echo "Rank-controlled benchmark: target rank = $rank"
    echo "Result directory: $result_dir"
    echo "============================================================"

    env \
        ST_FOM_REBUILD="$rebuild" \
        ST_RESULT_DIR="$result_dir" \
        DIFF_MANUFACTURED_RANK="$rank" \
        DIFF_ST_DDROM_MAX_MODES="$max_modes" \
        DIFF_ST_DDROM_WRITE_OUTPUT="${DIFF_ST_DDROM_WRITE_OUTPUT:-0}" \
        "$compare" "$@"

    rebuild=0

    timing="$result_dir/stddrom_timing.csv"
    if [[ ! -f "$timing" ]]; then
        echo "ERROR: missing timing file: $timing" >&2
        exit 1
    fi

    if (( first == 1 )); then
        head -n 1 "$timing" > "$summary"
        first=0
    fi
    awk -F, '$1 == "online" {print}' "$timing" >> "$summary"
done

echo
echo "Rank sweep completed."
echo "Summary: $summary"
if command -v column >/dev/null 2>&1; then
    column -s, -t "$summary" || cat "$summary"
else
    cat "$summary"
fi
