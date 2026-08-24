#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_elements(path, nen_expected=None):
    with path.open() as f:
        ne = int(f.readline().split()[0])
        elems = []
        for e in range(ne):
            row = f.readline().split()
            if not row:
                raise RuntimeError(f"truncated {path} at e={e}")
            if len(row) >= 2 and row[1].lstrip("+-").isdigit():
                nen = int(row[1])
                if nen_expected is not None and nen != nen_expected:
                    raise RuntimeError(f"{path}: e={e} NEN={nen}, expected {nen_expected}")
                if len(row) < 2 + nen:
                    raise RuntimeError(f"{path}: short row e={e}")
                elems.append([int(x) for x in row[2:2+nen]])
            else:
                nen = nen_expected
                if nen is None:
                    raise RuntimeError("nen_expected required for headerless rows")
                elems.append([int(x) for x in row[:nen]])
        return elems

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part_dir", type=Path)
    ap.add_argument("--nranks", type=int, required=True)
    ap.add_argument("--prism-base", default="graph_elem_prism.dat")
    ap.add_argument("--surf-base", default="surf_graph.dat")
    args = ap.parse_args()

    gsurf = gmap = gunmap = 0

    for r in range(args.nranks):
        pp = args.part_dir / f"{args.prism_base}.{r}"
        sp = args.part_dir / f"{args.surf_base}.{r}"

        if not pp.exists() or not sp.exists():
            print(f"[rank {r}] missing prism/surface file")
            continue

        pri = read_elements(pp, 6)
        surf = read_elements(sp, 3)

        faces = {}
        for e,c in enumerate(pri):
            for lf,loc in enumerate(((0,1,2),(3,4,5))):
                key = tuple(sorted(c[i] for i in loc))
                faces.setdefault(key, []).append((e,lf))

        mapped = 0
        unmapped = []
        for es,c in enumerate(surf):
            key = tuple(sorted(c))
            if key in faces:
                mapped += 1
            else:
                unmapped.append((es,c))

        print(
            f"[rank {r}] prisms={len(pri)} surf={len(surf)} "
            f"mapped={mapped} unmapped={len(unmapped)}"
        )
        for es,c in unmapped[:5]:
            print(f"  unmapped surf e={es} conn={c}")

        gsurf += len(surf)
        gmap += mapped
        gunmap += len(unmapped)

    print(f"LOCAL-COPY TOTALS surf={gsurf} mapped={gmap} unmapped={gunmap}")
    print("Note: overlap copies can make these totals larger than the physical wall-face count.")

if __name__ == "__main__":
    main()

