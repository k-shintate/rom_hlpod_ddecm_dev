#!/usr/bin/env bash

# Reference space-time spatial-block FOM launcher.
#
# The current Makefile_ST builds one executable, hlpod_diff_st_ddrom.  The executable
# selects the verified reference-FOM path with ST_DDROM_MODE=fom.  Therefore
# this wrapper reuses execution_ST_DDROM_MPI.sh instead of requiring a second
# executable or a second Makefile target.

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  execution_ST_spaceblock_MPI.sh <e> <ep> <nm> <nd> <np> <pa> <st>

This runs the verified spatial-block FOM path only.  It does not collect
snapshots or run the reduced model.

Environment:
  STSB_CHECK_EXTERNAL_RESIDUAL=0|1  default: 1
  STSB_WRITE_OUTPUT=0|1             default: 1
  STSB_REFERENCE_LOG=<path>         default:
                                     result_diff/<nm>-<np>-<nd>/reference_fom_np<np>.log
  ST_CLEAN_RESULT=0|1               default: 0
  ST_FOM_REBUILD=0|1                default: 1
  ST_PROJECT_ROOT=<project root>    optional
USAGE
}

if [[ "$#" -ne 7 ]]; then
    usage
    exit 2
fi

script_path="${BASH_SOURCE[0]:-$0}"
script_dir="$(cd "$(dirname "$script_path")" && pwd)"
ddrom_launcher="$script_dir/execution_ST_DDROM_MPI.sh"

if [[ ! -f "$ddrom_launcher" ]]; then
    echo "ERROR: required launcher was not found: $ddrom_launcher" >&2
    exit 1
fi

nm="$3"
nd="$4"
np="$5"

for pair in "nm:$nm" "nd:$nd" "np:$np"; do
    name="${pair%%:*}"
    value="${pair#*:}"
    if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: ${name} must be a positive integer: ${value}" >&2
        exit 2
    fi
done

find_project_root() {
    local candidate

    for candidate in \
        "${ST_PROJECT_ROOT:-}" \
        "$(pwd)" \
        "${script_dir}/../.." \
        "${script_dir}/../../.."; do
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
    echo "ERROR: project root was not found." >&2
    echo "Set ST_PROJECT_ROOT to rom_hlpod_ddecm_dev." >&2
    exit 1
fi

check_residual="${STSB_CHECK_EXTERNAL_RESIDUAL:-1}"
write_output="${STSB_WRITE_OUTPUT:-1}"
clean_result="${ST_CLEAN_RESULT:-0}"

for pair in \
    "STSB_CHECK_EXTERNAL_RESIDUAL:$check_residual" \
    "STSB_WRITE_OUTPUT:$write_output" \
    "ST_CLEAN_RESULT:$clean_result"; do
    name="${pair%%:*}"
    value="${pair#*:}"
    if [[ "$value" != 0 && "$value" != 1 ]]; then
        echo "ERROR: ${name} must be 0 or 1: ${value}" >&2
        exit 2
    fi
done

# Keep existing snapshots/bases by default.  Only the FOM branch is selected.
env \
    ST_PROJECT_ROOT="$project_root" \
    ST_DDROM_ACTION=fom \
    ST_DDROM_WRITE_OUTPUT="$write_output" \
    ST_DDROM_CHECK_FOM_RESIDUAL="$check_residual" \
    ST_DDROM_STRICT_FOM_RESIDUAL=0 \
    ST_CLEAN_RESULT="$clean_result" \
    bash "$ddrom_launcher" "$@"

result_dir="$project_root/result_diff/${nm}-${np}-${nd}"
source_log="$result_dir/stddrom_fom_np${np}.log"
target_log="${STSB_REFERENCE_LOG:-$result_dir/reference_fom_np${np}.log}"

if [[ ! -f "$source_log" ]]; then
    echo "ERROR: reference FOM log was not created: $source_log" >&2
    exit 1
fi

if [[ "$target_log" != /* ]]; then
    target_log="$(pwd)/$target_log"
fi

mkdir -p "$(dirname "$target_log")"
if [[ "$source_log" != "$target_log" ]]; then
    cp -f "$source_log" "$target_log"
fi

echo
echo "Reference spatial-block FOM completed."
echo "Reference log: $target_log"
