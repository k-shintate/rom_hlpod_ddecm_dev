#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit(
        "usage: compare_step_checksums.py "
        "<causal_fluid_st_step_checksums.csv> "
        "<monolithic_fluid_st_step_checksums.csv>"
    )

def read(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))

a = read(sys.argv[1])
b = read(sys.argv[2])

if len(a) != len(b):
    raise SystemExit(f"row-count mismatch: {len(a)} vs {len(b)}")

fields = ["state_norm", "velocity_norm", "pressure_norm"]
max_rel = {k: 0.0 for k in fields}

for ra, rb in zip(a, b):
    if ra["step"] != rb["step"]:
        raise SystemExit(f"step mismatch: {ra['step']} vs {rb['step']}")
    for key in fields:
        va = float(ra[key])
        vb = float(rb[key])
        rel = abs(vb - va) / max(abs(va), 1.0e-300)
        max_rel[key] = max(max_rel[key], rel)

print("checksum-norm comparison (diagnostic only; not a vector identity proof)")
for key in fields:
    print(f"  max relative {key:14s}: {max_rel[key]:.15e}")
