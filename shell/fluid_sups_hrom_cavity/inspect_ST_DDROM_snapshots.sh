#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "$#" -lt 1 || "$#" -gt 2 ]]; then
    echo "usage: $0 <snapshot-directory> [np]" >&2
    exit 2
fi

directory="$1"
np="${2:-${FLUID_ST_MPI_RANKS:-4}}"

[[ -d "$directory" ]] || {
    echo "ERROR: directory not found: $directory" >&2
    exit 1
}

echo "Snapshot inventory: $directory"
for ((rank=0; rank<np; rank++)); do
    count=0
    while printf -v file "%s/fluid_st_window_%04d_rank_%06d.bin" \
        "$directory" "$count" "$rank" && [[ -f "$file" ]]; do
        ((count+=1))
    done
    echo "  rank $rank : consecutive windows 0..$((count-1))  count=$count"
done
