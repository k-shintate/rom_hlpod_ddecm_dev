#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit(
        "usage: compare_mpi_checksums.py "
        "<np1-result-dir> <npN-result-dir>"
    )

a_dir = Path(sys.argv[1])
b_dir = Path(sys.argv[2])

def load(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

def rel(a, b):
    return abs(b - a) / max(abs(a), 1.0e-300)

checks = [
    (
        "step",
        "fluid_st_step_checksums.csv",
        ["state_norm", "velocity_norm", "pressure_norm"],
    ),
    (
        "window",
        "fluid_st_window_checksums.csv",
        ["window_norm", "velocity_norm", "pressure_norm", "right_trace_norm"],
    ),
]

for label, name, fields in checks:
    a = load(a_dir / name)
    b = load(b_dir / name)

    if len(a) != len(b):
        raise SystemExit(
            f"{label}: row-count mismatch: {len(a)} vs {len(b)}")

    maxima = {field: 0.0 for field in fields}

    for ra, rb in zip(a, b):
        key = "step" if label == "step" else "window"
        if ra[key] != rb[key]:
            raise SystemExit(
                f"{label}: index mismatch: {ra[key]} vs {rb[key]}")

        for field in fields:
            maxima[field] = max(
                maxima[field],
                rel(float(ra[field]), float(rb[field])))

    print(f"{label} checksum comparison")
    for field in fields:
        print(f"  max relative {field:18s}: {maxima[field]:.15e}")
    print()

print(
    "NOTE: checksum/norm agreement is a first MPI diagnostic only; "
    "owned velocity/pressure vector comparison is the stronger validation."
)
