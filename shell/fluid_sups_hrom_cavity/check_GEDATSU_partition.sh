#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "usage: $0 <result-dir> <np>" >&2
    exit 2
fi

result_dir="$1"
np="$2"

if ! [[ "$np" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: np must be positive: $np" >&2
    exit 2
fi

for global_file in node.dat elem.dat graph.dat D_bc_v.dat; do
    [[ -f "$result_dir/$global_file" ]] || {
        echo "ERROR: missing global file: $result_dir/$global_file" >&2
        exit 1
    }
done

for ((rank=0; rank<np; rank++)); do
    for local_file in \
        "node.dat.${rank}" \
        "elem.dat.${rank}" \
        "graph.dat.${rank}" \
        "D_bc_v.dat.${rank}"; do
        [[ -f "$result_dir/parted.0/$local_file" ]] || {
            echo "ERROR: missing rank file: $result_dir/parted.0/$local_file" >&2
            exit 1
        }
    done
done

echo "GEDATSU spatial partition OK"
echo "  result_dir : $result_dir"
echo "  np         : $np"
echo "  rank files : node / elem / graph / D_bc_v"
echo "  time mesh  : not partitioned"
