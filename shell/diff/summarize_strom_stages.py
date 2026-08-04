#!/usr/bin/env python3
from __future__ import annotations

import csv
import pathlib
import sys


def read_stage(path: pathlib.Path) -> dict[str, str]:
    rows: list[dict[str, str]] = []
    summary: dict[str, str] = {}

    with path.open(newline="") as fp:
        for line in fp:
            if line.startswith("# summary,"):
                for item in line.strip().split(",")[1:]:
                    if "=" in item:
                        key, value = item.split("=", 1)
                        summary[key] = value
                continue
            if line.startswith("#"):
                continue
            if not rows:
                header = next(csv.reader([line]))
                rows.append({"__header__": "\0".join(header)})
            else:
                header = rows[0]["__header__"].split("\0")
                values = next(csv.reader([line]))
                if len(values) == len(header):
                    rows.append(dict(zip(header, values)))

    data_rows = [row for row in rows[1:] if "window" in row]
    for key in (
        "g_relative",
        "a_relative",
        "q_relative",
        "q_right_relative",
    ):
        if data_rows:
            summary.setdefault(
                f"max_{key}",
                f"{max(float(row[key]) for row in data_rows):.15e}",
            )
    return summary


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_strom_stages.py <stage-root>", file=sys.stderr)
        return 2

    root = pathlib.Path(sys.argv[1])
    output = root / "stage_summary.csv"
    records: list[dict[str, str]] = []

    for stage_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        comparison = stage_dir / "stddrom_vs_sequential_strom.csv"
        if not comparison.exists():
            continue
        record = read_stage(comparison)
        record["stage"] = stage_dir.name
        records.append(record)

    fields = [
        "stage",
        "degree",
        "slabs_per_window",
        "global_modes",
        "D",
        "L",
        "max_effective_rhs",
        "max_g",
        "max_a",
        "max_q",
        "max_q_right",
        "tolerance",
        "status",
    ]

    with output.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
