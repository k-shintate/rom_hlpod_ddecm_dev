#!/usr/bin/env python3
"""Build D_bc_v.dat directly from Physical-group node lists.

BC convention preserved from the old shell:
  Inlet         : (u,v,w) = (1,0,0)
  Cylinder_wall : (0,0,0)
  Left/Right    : slip -> constrain y component only (dof=1)
  Top/Bottom    : slip -> constrain z component only (dof=2)
  Outlet        : no Dirichlet velocity BC

Later groups overwrite earlier groups at corner/shared nodes, matching merge_bc.py.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path


def read_ids(path: Path) -> list[int]:
    with path.open(encoding="utf-8") as f:
        lines = [x.strip() for x in f if x.strip()]
    if not lines:
        raise ValueError(f"empty node-id file: {path}")
    n = int(lines[0])
    ids = [int(x) for x in lines[1:]]
    if len(ids) != n:
        raise ValueError(f"{path}: header={n}, parsed={len(ids)}")
    return ids


def rows(ids: list[int], values: dict[int, float]):
    return {n: [(n, dof, val) for dof, val in sorted(values.items())] for n in ids}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh_dir", type=Path)
    ap.add_argument("--output", default="D_bc_v.dat")
    ap.add_argument("--inlet-u", type=float, default=1.0)
    args = ap.parse_args()
    d = args.mesh_dir

    # Same priority as the old merge call: later groups win for shared nodes.
    specs = [
        ("Right_wall", {1: 0.0}),
        ("Cylinder_wall", {0: 0.0, 1: 0.0, 2: 0.0}),
        ("Left_wall", {1: 0.0}),
        ("Bottom_wall", {2: 0.0}),
        ("Top_wall", {2: 0.0}),
        ("Inlet", {0: args.inlet_u, 1: 0.0, 2: 0.0}),
    ]

    merged: OrderedDict[int, list[tuple[int, int, float]]] = OrderedDict()
    for name, vals in specs:
        p = d / f"{name}_node_ids.dat"
        if not p.exists():
            raise SystemExit(f"required Physical group node list is missing: {p}")
        for n, block in rows(read_ids(p), vals).items():
            merged[n] = block

    flat = [r for block in merged.values() for r in block]
    out = d / args.output
    with out.open("w", encoding="utf-8") as f:
        f.write(f"{len(flat)} 3\n")
        for n, dof, val in flat:
            f.write(f"{n} {dof} {val:.15e}\n")
    print(f"[bc] constrained_rows={len(flat)} output={out}")


if __name__ == "__main__":
    main()
