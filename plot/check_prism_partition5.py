#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
from collections import Counter

def read_nodes(path):
    with path.open() as f:
        n = int(f.readline().split()[0])
        xyz = []
        for i in range(n):
            row = f.readline().split()
            if len(row) < 3:
                raise RuntimeError(f"truncated node file at {i}")
            xyz.append([float(row[0]), float(row[1]), float(row[2])])
    return np.asarray(xyz)

def read_elem(path, nen):
    with path.open() as f:
        ne = int(f.readline().split()[0])
        out = []
        for e in range(ne):
            row = f.readline().split()
            if not row:
                raise RuntimeError(f"truncated {path} e={e}")
            if len(row) >= 2 + nen and int(row[1]) == nen:
                c = [int(x) for x in row[2:2+nen]]
            elif len(row) >= nen:
                c = [int(x) for x in row[:nen]]
            else:
                raise RuntimeError(f"bad row in {path}: {row}")
            out.append(c)
    return out

def tet_faces(c):
    # all 4 TRI3 faces
    ids = ((0,1,2),(0,1,3),(0,2,3),(1,2,3))
    return [tuple(sorted((c[a],c[b],c[d]))) for a,b,d in ids]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part_dir", type=Path)
    ap.add_argument("--rank", type=int, required=True)
    ap.add_argument("--node-base", default="node.dat")
    ap.add_argument("--tet-base", default="graph_elem_tet.dat")
    ap.add_argument("--prism-base", default="graph_elem_prism.dat")
    ap.add_argument("--interface-radius", type=float, default=0.65)
    ap.add_argument("--tol", type=float, default=5e-4)
    ap.add_argument("--show", type=int, default=20)
    args = ap.parse_args()

    r = args.rank
    xyz = read_nodes(args.part_dir / f"{args.node_base}.{r}")
    tet = read_elem(args.part_dir / f"{args.tet_base}.{r}", 4)
    pri = read_elem(args.part_dir / f"{args.prism_base}.{r}", 6)

    tet_nodes = set(n for c in tet for n in c)
    pri_nodes = set(n for c in pri for n in c)
    shared_nodes = tet_nodes & pri_nodes

    print(f"rank={r}")
    print(f"tet={len(tet)} pri={len(pri)}")
    print(f"tet_unique_nodes={len(tet_nodes)} pri_unique_nodes={len(pri_nodes)}")
    print(f"shared_tet_pri_nodes={len(shared_nodes)}")

    # Boundary faces of tet mesh: face appearing once.
    tf = Counter()
    for c in tet:
        tf.update(tet_faces(c))
    tet_boundary_faces = {f for f,n in tf.items() if n == 1}

    # For each PRI6, decide which triangular face is the outer one by radius.
    pri_outer_faces = []
    pri_inner_faces = []
    for e,c in enumerate(pri):
        f0 = c[:3]
        f1 = c[3:]
        r0 = np.mean(np.linalg.norm(xyz[f0,:2], axis=1))
        r1 = np.mean(np.linalg.norm(xyz[f1,:2], axis=1))
        if r0 >= r1:
            outer, inner = f0, f1
        else:
            outer, inner = f1, f0
        pri_outer_faces.append(tuple(sorted(outer)))
        pri_inner_faces.append(tuple(sorted(inner)))

    matched = [f for f in pri_outer_faces if f in tet_boundary_faces]
    missing = [f for f in pri_outer_faces if f not in tet_boundary_faces]

    # Unique outer faces (layered prism file has only outermost-layer faces near r=R)
    outer_near_R = []
    for f in pri_outer_faces:
        rr = np.linalg.norm(xyz[list(f),:2], axis=1)
        if np.all(np.abs(rr - args.interface_radius) <= args.tol):
            outer_near_R.append(f)

    outer_near_R_unique = set(outer_near_R)
    matched_R = outer_near_R_unique & tet_boundary_faces
    missing_R = outer_near_R_unique - tet_boundary_faces

    interface_nodes_pri = set(n for f in outer_near_R_unique for n in f)
    interface_nodes_tet = set(n for f in matched_R for n in f)

    print(f"tet_boundary_TRI3={len(tet_boundary_faces)}")
    print(f"all_prism_outer_face_records={len(pri_outer_faces)}")
    print(f"all_prism_outer_face_records_matching_tet_boundary={len(matched)}")
    print(f"interface_PRI_faces_unique_near_r={len(outer_near_R_unique)}")
    print(f"interface_faces_matched_exactly={len(matched_R)}")
    print(f"interface_faces_missing={len(missing_R)}")
    print(f"interface_pri_nodes={len(interface_nodes_pri)} matched_interface_nodes={len(interface_nodes_tet)}")

    if missing_R:
        print("Missing interface face examples:")
        for f in list(sorted(missing_R))[:args.show]:
            rr = np.linalg.norm(xyz[list(f),:2], axis=1)
            print(f"  face={f} radii={rr.tolist()}")
        raise SystemExit(1)

    if not outer_near_R_unique:
        print("ERROR: no PRI interface faces detected near requested radius.")
        raise SystemExit(2)

    print("OK: PRI6 outer interface TRI3 faces are shared exactly with TET4 boundary faces.")

if __name__ == "__main__":
    main()

