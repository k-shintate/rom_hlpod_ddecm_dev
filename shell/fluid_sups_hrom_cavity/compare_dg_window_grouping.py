#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit(
        "usage: compare_window_grouping.py <slabs1_step_checksums.csv> "
        "<slabs4_step_checksums.csv>"
    )

def read(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))

a = read(sys.argv[1])
b = read(sys.argv[2])

if len(a) != len(b):
    raise SystemExit(f"row count mismatch: {len(a)} vs {len(b)}")

fields = ["state_norm", "velocity_norm", "pressure_norm"]
max_rel = {name: 0.0 for name in fields}

for ra, rb in zip(a, b):
    if ra["step"] != rb["step"]:
        raise SystemExit(
            f"step mismatch: {ra['step']} vs {rb['step']}")
    for name in fields:
        xa = float(ra[name])
        xb = float(rb[name])
        rel = abs(xa - xb) / max(abs(xa), 1.0e-300)
        max_rel[name] = max(max_rel[name], rel)

print("dG window-grouping checksum comparison")
print("(diagnostic only; owned-vector comparison is the stronger test)")
for name, value in max_rel.items():
    print(f"  max relative {name:14s}: {value:.15e}")
