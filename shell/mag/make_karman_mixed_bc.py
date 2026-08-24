#!/usr/bin/env python3
"""
Generate D_bc_v.dat for the mixed TET4+PRI6 Karman mesh.

Strong velocity BCs:
  inlet (x=xmin)        : ux=Uin, uy=0, uz=0
  y-min / y-max walls   : uy=0 only       (slip)
  z-min / z-max walls   : uz=0 only       (slip)
  outlet (x=xmax)       : no strong BC
  cylinder wall         : no strong BC    (Nitsche)

Input node.dat format:
  N
  x y z
  ...

Output D_bc_v.dat format:
  Nbc 3
  node_id dof value
  ...
"""
import argparse
from pathlib import Path

def read_nodes(path):
    with path.open() as f:
        n = int(f.readline().split()[0])
        xyz = []
        for i in range(n):
            row = f.readline().split()
            if len(row) < 3:
                raise RuntimeError(f"node.dat truncated at node {i}")
            xyz.append(tuple(float(v) for v in row[:3]))
    return xyz

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("node_dat", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--uin", type=float, default=1.0)
    ap.add_argument("--tol", type=float, default=None,
                    help="absolute geometry tolerance; default = 1e-9 * domain length")
    args = ap.parse_args()

    xyz = read_nodes(args.node_dat)

    xs = [p[0] for p in xyz]
    ys = [p[1] for p in xyz]
    zs = [p[2] for p in xyz]

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)

    L = max(xmax-xmin, ymax-ymin, zmax-zmin, 1.0)
    tol = args.tol if args.tol is not None else 1.0e-9 * L

    # (node,dof) -> value
    bc = {}

    def set_bc(i, d, value):
        key = (i, d)
        if key in bc and abs(bc[key] - value) > 1.0e-12:
            raise RuntimeError(
                f"conflicting BC: node={i} dof={d}: {bc[key]} vs {value}"
            )
        bc[key] = value

    n_inlet = n_ywall = n_zwall = 0

    for i, (x,y,z) in enumerate(xyz):
        on_inlet = abs(x-xmin) <= tol
        on_ywall = abs(y-ymin) <= tol or abs(y-ymax) <= tol
        on_zwall = abs(z-zmin) <= tol or abs(z-zmax) <= tol

        # Slip walls first.
        if on_ywall:
            set_bc(i, 1, 0.0)  # uy only
            n_ywall += 1

        if on_zwall:
            set_bc(i, 2, 0.0)  # uz only
            n_zwall += 1

        # Inlet supersedes/intersects walls consistently:
        # uy=uz=0 are the same imposed values.
        if on_inlet:
            set_bc(i, 0, args.uin)
            set_bc(i, 1, 0.0)
            set_bc(i, 2, 0.0)
            n_inlet += 1

    rows = sorted((i,d,v) for (i,d),v in bc.items())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        f.write(f"{len(rows)} 3\n")
        for i,d,v in rows:
            f.write(f"{i} {d} {v:.16e}\n")

    counts = [0,0,0]
    for _,d,_ in rows:
        if 0 <= d < 3:
            counts[d] += 1

    print(f"nodes={len(xyz)}")
    print(f"bounds x=[{xmin:.16e},{xmax:.16e}]")
    print(f"       y=[{ymin:.16e},{ymax:.16e}]")
    print(f"       z=[{zmin:.16e},{zmax:.16e}]")
    print(f"tol={tol:.3e}")
    print(f"inlet_nodes={n_inlet} ywall_nodes={n_ywall} zwall_nodes={n_zwall}")
    print(f"BC rows={len(rows)}")
    print(f"ux rows={counts[0]} uy rows={counts[1]} uz rows={counts[2]}")
    print(f"wrote {args.output}")

if __name__ == "__main__":
    main()
