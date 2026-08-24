#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_graph(path):
    with path.open() as f:
        n = int(f.readline().split()[0])
        adj = [set() for _ in range(n)]
        for row in range(n):
            vals = f.readline().split()
            if len(vals) < 2:
                raise RuntimeError(f"truncated graph at row {row}")
            node = int(vals[0])
            deg = int(vals[1])
            ns = [int(x) for x in vals[2:]]
            while len(ns) < deg:
                ns.extend(int(x) for x in f.readline().split())
            if node < 0 or node >= n:
                raise RuntimeError(f"bad graph node id {node}")
            adj[node].update(ns[:deg])
    return adj

def read_prisms(path):
    with path.open() as f:
        ne = int(f.readline().split()[0])
        elems = []
        for e in range(ne):
            row = f.readline().split()
            if not row:
                raise RuntimeError(f"truncated prism file at element {e}")
            if len(row) >= 8 and int(row[1]) == 6:
                c = [int(x) for x in row[2:8]]
            elif len(row) >= 6:
                c = [int(x) for x in row[:6]]
            else:
                raise RuntimeError(f"invalid prism row {e}: {row}")
            elems.append(c)
    return elems

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part_dir", type=Path)
    ap.add_argument("--rank", type=int, required=True)
    ap.add_argument("--graph-base", default="graph.dat")
    ap.add_argument("--prism-base", default="graph_elem_prism.dat")
    ap.add_argument("--show", type=int, default=30)
    args = ap.parse_args()

    gp = args.part_dir / f"{args.graph_base}.{args.rank}"
    ep = args.part_dir / f"{args.prism_base}.{args.rank}"

    adj = read_graph(gp)
    elems = read_prisms(ep)

    missing = []
    total_pairs = 0
    prism_nodes = set()

    for e, c in enumerate(elems):
        prism_nodes.update(c)
        for i in c:
            for j in c:
                if i == j:
                    continue
                total_pairs += 1
                if j not in adj[i]:
                    missing.append((e, i, j, c))

    zero_degree_prism = [i for i in sorted(prism_nodes) if len(adj[i]) == 0]

    print(f"rank={args.rank}")
    print(f"graph_nodes={len(adj)}")
    print(f"prisms={len(elems)} unique_prism_nodes={len(prism_nodes)}")
    print(f"required_directed_pairs={total_pairs}")
    print(f"missing_directed_pairs={len(missing)}")
    print(f"prism_nodes_with_zero_graph_degree={len(zero_degree_prism)}")

    for rec in missing[:args.show]:
        e, i, j, c = rec
        print(f"  missing e={e} {i}->{j} conn={c}")

    if zero_degree_prism:
        print("zero-degree prism nodes:", zero_degree_prism[:args.show])

    if missing:
        raise SystemExit(1)

    print("OK: every PRI6 node pair exists in graph sparsity.")

if __name__ == "__main__":
    main()

