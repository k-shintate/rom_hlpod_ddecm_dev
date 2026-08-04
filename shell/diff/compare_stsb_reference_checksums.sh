#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 <reference_fom.log> <stddrom_collect.log>" >&2
    exit 2
fi

reference_log="$1"
collect_log="$2"

for file in "$reference_log" "$collect_log"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: log file not found: $file" >&2
        exit 1
    fi
done

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

grep 'STSB reference checksum window' "$reference_log" \
    > "$tmp_dir/reference.txt" || true
grep 'STSB reference checksum window' "$collect_log" \
    > "$tmp_dir/collect.txt" || true

if [[ ! -s "$tmp_dir/reference.txt" ]]; then
    echo "ERROR: no checksum records in reference log." >&2
    exit 1
fi
if [[ ! -s "$tmp_dir/collect.txt" ]]; then
    echo "ERROR: no checksum records in collect log." >&2
    exit 1
fi

if diff -u "$tmp_dir/reference.txt" "$tmp_dir/collect.txt"; then
    echo "PASS: reference FOM and ST-DDROM collect checksums are identical."
else
    echo "FAIL: checksum records differ." >&2
    exit 1
fi
