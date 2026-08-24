#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np


def read_nodes(path: Path):
    with path.open() as f:
        first = f.readline().split()
        if not first:
            raise RuntimeError(f"empty node file: {path}")
        n = int(first[0])
        xyz = []
        for i in range(n):
            row = f.readline().split()
            if len(row) < 3:
                raise RuntimeError(f"truncated node file {path} at node {i}")
            xyz.append([float(row[0]), float(row[1]), float(row[2])])
    return np.asarray(xyz, dtype=float)


def read_prisms(path: Path):
    with path.open() as f:
        first = f.readline().split()
        if not first:
            raise RuntimeError(f"empty element file: {path}")
        ne = int(first[0])
        conn = []
        for e in range(ne):
            row = f.readline().split()
            if len(row) < 2:
                raise RuntimeError(f"truncated element file {path} at element {e}")

            # Supported forms:
            #   elem_id 6 n0 n1 n2 n3 n4 n5
            #   n0 n1 n2 n3 n4 n5   (when header is "Nelem 6")
            if len(row) >= 8 and int(row[1]) == 6:
                nodes = list(map(int, row[2:8]))
            elif len(row) >= 6:
                nodes = list(map(int, row[:6]))
            else:
                raise RuntimeError(f"invalid PRI6 row in {path}: {' '.join(row)}")
            conn.append(nodes)
    return np.asarray(conn, dtype=int)


def dshape(r, s, z):
    c = 0.5
    dr = np.array([
        -c*(1-z), c*(1-z), 0.0,
        -c*(1+z), c*(1+z), 0.0
    ])
    ds = np.array([
        -c*(1-z), 0.0, c*(1-z),
        -c*(1+z), 0.0, c*(1+z)
    ])
    dz = np.array([
        -c*(1-r-s), -c*r, -c*s,
         c*(1-r-s),  c*r,  c*s
    ])
    return np.column_stack((dr, ds, dz))


def prism_gauss_points():
    a = 1.0 / math.sqrt(3.0)
    out = []
    # 2x2x2 Gauss on cube, mapped with Duffy transform.
    for u in (-a, a):
        for v in (-a, a):
            for w in (-a, a):
                r = 0.5*(1.0 + u)
                s = 0.25*(1.0-u)*(1.0+v)
                z = w
                out.append((r, s, z))
    # Also test the reference center.
    out.append((1.0/3.0, 1.0/3.0, 0.0))
    return out


def check_rank(part_dir: Path, rank: int, node_base: str, prism_base: str):
    node_path = part_dir / f"{node_base}.{rank}"
    pri_path  = part_dir / f"{prism_base}.{rank}"

    if not node_path.exists():
        return {"rank": rank, "error": f"missing {node_path}"}
    if not pri_path.exists():
        return {"rank": rank, "error": f"missing {pri_path}"}

    xyz = read_nodes(node_path)
    conn = read_prisms(pri_path)

    result = {
        "rank": rank,
        "nodes": len(xyz),
        "prisms": len(conn),
    }

    if len(conn) == 0:
        result.update({
            "status": "EMPTY",
            "max_node_id": -1,
            "bad_index": 0,
            "neg_det": 0,
            "zero_det": 0,
        })
        return result

    max_node = int(conn.max())
    min_node = int(conn.min())
    result["max_node_id"] = max_node
    result["min_node_id"] = min_node

    bad_index = int(np.count_nonzero((conn < 0) | (conn >= len(xyz))))
    result["bad_index"] = bad_index
    if bad_index:
        result["status"] = "BAD_NODE_NUMBERING"
        return result

    det_min = math.inf
    det_max = -math.inf
    neg_det = 0
    zero_det = 0
    worst_e = -1
    worst_p = None

    radial_edge_min = math.inf
    radial_edge_max = 0.0
    radial_cos_min = 1.0
    inner_r_min = math.inf
    inner_r_max = 0.0
    outer_r_min = math.inf
    outer_r_max = 0.0

    qps = prism_gauss_points()

    for e, c in enumerate(conn):
        X = xyz[c]

        # Determine inner/outer triangular faces from mean xy radius.
        r0 = np.linalg.norm(X[:3, :2], axis=1)
        r1 = np.linalg.norm(X[3:, :2], axis=1)
        if float(r0.mean()) <= float(r1.mean()):
            rin, rout = r0, r1
        else:
            rin, rout = r1, r0

        inner_r_min = min(inner_r_min, float(rin.min()))
        inner_r_max = max(inner_r_max, float(rin.max()))
        outer_r_min = min(outer_r_min, float(rout.min()))
        outer_r_max = max(outer_r_max, float(rout.max()))

        # Canonical through-thickness PRI6 edges: 0-3, 1-4, 2-5.
        for a, b in ((0, 3), (1, 4), (2, 5)):
            d = X[b] - X[a]
            L = float(np.linalg.norm(d))
            radial_edge_min = min(radial_edge_min, L)
            radial_edge_max = max(radial_edge_max, L)

            mid = 0.5*(X[a] + X[b])
            rr = np.array([mid[0], mid[1], 0.0])
            nr = float(np.linalg.norm(rr))
            if L > 0.0 and nr > 0.0:
                cosv = abs(float(np.dot(d, rr) / (L * nr)))
                radial_cos_min = min(radial_cos_min, cosv)

        for qp in qps:
            D = dshape(*qp)              # 6 x 3
            J = X.T @ D                  # 3 x 3
            det = float(np.linalg.det(J))

            if det < det_min:
                det_min = det
                worst_e = e
                worst_p = qp
            det_max = max(det_max, det)

            if det < -1.0e-12:
                neg_det += 1
            elif abs(det) <= 1.0e-12:
                zero_det += 1

    result.update({
        "status": "OK" if neg_det == 0 and zero_det == 0 else "BAD_JACOBIAN",
        "det_min": det_min,
        "det_max": det_max,
        "neg_det": neg_det,
        "zero_det": zero_det,
        "worst_elem": worst_e,
        "worst_qp": worst_p,
        "radial_edge_min": radial_edge_min,
        "radial_edge_max": radial_edge_max,
        "radial_cos_min": radial_cos_min,
        "inner_r_min": inner_r_min,
        "inner_r_max": inner_r_max,
        "outer_r_min": outer_r_min,
        "outer_r_max": outer_r_max,
    })
    return result


def main():
    ap = argparse.ArgumentParser(
        description="Check partitioned PRI6 connectivity and signed Jacobians."
    )
    ap.add_argument("part_dir", type=Path, help="directory such as mesh_karman_vortex/parted.0")
    ap.add_argument("--nranks", type=int, required=True)
    ap.add_argument("--node-base", default="node.dat")
    ap.add_argument("--prism-base", default="graph_elem_prism.dat")
    args = ap.parse_args()

    global_bad = False

    for rank in range(args.nranks):
        r = check_rank(args.part_dir, rank, args.node_base, args.prism_base)

        if "error" in r:
            global_bad = True
            print(f"[rank {rank}] ERROR {r['error']}")
            continue

        if r["prisms"] == 0:
            print(f"[rank {rank}] EMPTY nodes={r['nodes']} prisms=0")
            continue

        print(
            f"[rank {rank}] {r['status']} "
            f"nodes={r['nodes']} prisms={r['prisms']} "
            f"node_id=[{r['min_node_id']},{r['max_node_id']}] bad_index={r['bad_index']}"
        )

        if r["bad_index"]:
            global_bad = True
            continue

        print(
            f"           detJ=[{r['det_min']:.6e},{r['det_max']:.6e}] "
            f"neg={r['neg_det']} zero={r['zero_det']} "
            f"worst_elem={r['worst_elem']} qp={r['worst_qp']}"
        )
        print(
            f"           through-edge=[{r['radial_edge_min']:.6e},{r['radial_edge_max']:.6e}] "
            f"min|cos(radial)|={r['radial_cos_min']:.6f}"
        )
        print(
            f"           face radius inner=[{r['inner_r_min']:.6e},{r['inner_r_max']:.6e}] "
            f"outer=[{r['outer_r_min']:.6e},{r['outer_r_max']:.6e}]"
        )

        if r["status"] != "OK":
            global_bad = True

    raise SystemExit(1 if global_bad else 0)


if __name__ == "__main__":
    main()

