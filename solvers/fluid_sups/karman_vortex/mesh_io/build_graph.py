#!/usr/bin/env python3
"""Build one nodal graph from tetra + prism connectivity.

Output format expected by STEP_set_nzpattern_from_graph():
  Nnode
  node_id Nadj adj0 adj1 ...
  ... one line for every node, including isolated nodes.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def read_node_count(path: Path) -> int:
    with path.open(encoding="utf-8") as f:
        line = f.readline().split()
    if not line:
        raise ValueError(f"empty node file: {path}")
    return int(line[0])


def read_elem(path: Path, expected_nen: int) -> list[list[int]]:
    with path.open(encoding="utf-8") as f:
        header = f.readline().split()
        if len(header) < 2:
            raise ValueError(f"invalid element header: {path}")
        nelem, nen = map(int, header[:2])
        if nen != expected_nen:
            raise ValueError(f"{path}: expected NEN={expected_nen}, got {nen}")
        elems = []
        for line_no, raw in enumerate(f, 2):
            if not raw.strip():
                continue
            e = list(map(int, raw.split()))
            if len(e) != nen:
                raise ValueError(f"{path}:{line_no}: expected {nen} node ids, got {len(e)}")
            elems.append(e)
    if len(elems) != nelem:
        raise ValueError(f"{path}: header nelem={nelem}, parsed={len(elems)}")
    return elems


def add_elements(adj: list[set[int]], elems: list[list[int]], label: str) -> None:
    nnode = len(adj)
    for ie, elem in enumerate(elems):
        for n in elem:
            if not 0 <= n < nnode:
                raise ValueError(f"{label} element {ie}: node id {n} outside [0,{nnode})")
        # A finite element creates a complete nodal clique in its element matrix.
        for i in elem:
            adj[i].update(j for j in elem if j != i)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--node", type=Path, default=Path("node.dat"))
    ap.add_argument("--tet", type=Path, default=Path("elem_tet.dat"))
    ap.add_argument("--prism", type=Path, default=Path("elem_prism.dat"))
    ap.add_argument("--output", type=Path, default=Path("graph.dat"))
    args = ap.parse_args()

    nnode = read_node_count(args.node)
    tet = read_elem(args.tet, 4)
    pri = read_elem(args.prism, 6)
    adj = [set() for _ in range(nnode)]
    add_elements(adj, tet, "tetra")
    add_elements(adj, pri, "prism")

    with args.output.open("w", encoding="utf-8") as f:
        f.write(f"{nnode}\n")
        for i, ngh in enumerate(adj):
            vals = sorted(ngh)
            f.write(f"{i} {len(vals)}")
            if vals:
                f.write(" " + " ".join(map(str, vals)))
            f.write("\n")

    nedge = sum(len(a) for a in adj) // 2
    print(f"[graph] nodes={nnode} undirected_edges={nedge} output={args.output}")


if __name__ == "__main__":
    main()
