#!/usr/bin/env python3
"""
Detailed ROM-vs-HROM reduced-matrix diagnostic.

Reads from the current case directory:
  DDECM/rom_block_entries.*.csv
  DDECM/hrom_block_entries.*.csv

For every physical block it computes:
  ROM norm      = ||R||_F
  HROM norm     = ||H||_F
  exact error   = ||H-R||_F / ||R||_F
  gain alpha    = <H,R> / <R,R>
  cosine        = <H,R> / (||H|| ||R||)
  shape error   = ||H-alpha R||_F / ||R||_F

Interpretation:
  cosine ~ 1 and shape_error small, alpha far from 1
      -> mainly a scaling / duplicate-weight issue
  cosine substantially below 1
      -> block shape / modal composition is wrong

Usage:
  python3 compare_hrom_blocks_scaling.py CASE_DIR
"""

import csv
import glob
import math
import os
import sys
from collections import defaultdict


def load(pattern):
    files = sorted(glob.glob(pattern))
    rows = []
    for fn in files:
        with open(fn, newline="") as f:
            rows.extend(csv.DictReader(f))
    return files, rows


def to_i(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def to_f(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def entry_key(r, namespace):
    if namespace == "sid":
        row = to_i(r.get("row_subdomain_id"))
        col = to_i(r.get("col_subdomain_id"))
    else:
        row = to_i(r.get("row_global_id"))
        col = to_i(r.get("col_global_id"))

    return (
        r.get("type", ""),
        row,
        col,
        to_i(r.get("ir")),
        to_i(r.get("ic")),
    )


def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 2

    case = sys.argv[1]
    d = os.path.join(case, "DDECM")

    rom_files, rom_rows = load(os.path.join(d, "rom_block_entries.*.csv"))
    hrom_files, hrom_rows = load(os.path.join(d, "hrom_block_entries.*.csv"))

    print("=" * 96)
    print(case)
    print("=" * 96)
    print("ROM entry files :", len(rom_files))
    print("HROM entry files:", len(hrom_files))

    if not rom_rows or not hrom_rows:
        print("ERROR: scalar entry CSVs are missing.")
        return 1

    # sid/gid are already known to coincide in the current runs.
    namespace = "gid"

    rom = {entry_key(r, namespace): to_f(r.get("value")) for r in rom_rows}
    hrom = {entry_key(r, namespace): to_f(r.get("value")) for r in hrom_rows}

    rkeys = set(rom)
    hkeys = set(hrom)
    common = sorted(rkeys & hkeys)

    print("common scalar entries:", len(common))
    print("ROM-only entries     :", len(rkeys - hkeys))
    print("HROM-only entries    :", len(hkeys - rkeys))

    if not common:
        return 1

    # Per-block scalar vectors
    block_rom = defaultdict(dict)
    block_hrom = defaultdict(dict)

    for typ, row, col, ir, ic in common:
        bk = (typ, row, col)
        block_rom[bk][(ir, ic)] = rom[(typ, row, col, ir, ic)]
        block_hrom[bk][(ir, ic)] = hrom[(typ, row, col, ir, ic)]

    records = []

    total_r2 = 0.0
    total_d2 = 0.0

    type_r2 = defaultdict(float)
    type_d2 = defaultdict(float)

    mode_row_r2 = defaultdict(float)
    mode_row_d2 = defaultdict(float)
    mode_col_r2 = defaultdict(float)
    mode_col_d2 = defaultdict(float)

    for bk in sorted(block_rom):
        rr = block_rom[bk]
        hh = block_hrom[bk]

        coords = sorted(set(rr) & set(hh))

        r2 = 0.0
        h2 = 0.0
        d2 = 0.0
        dot = 0.0

        for c in coords:
            rv = rr[c]
            hv = hh[c]
            dv = hv - rv
            ir, ic = c

            r2 += rv * rv
            h2 += hv * hv
            d2 += dv * dv
            dot += rv * hv

            mode_row_r2[ir] += rv * rv
            mode_row_d2[ir] += dv * dv
            mode_col_r2[ic] += rv * rv
            mode_col_d2[ic] += dv * dv

        rn = math.sqrt(r2)
        hn = math.sqrt(h2)
        dn = math.sqrt(d2)

        rel = dn / max(rn, 1.0e-30)
        ratio = hn / max(rn, 1.0e-30)
        alpha = dot / max(r2, 1.0e-300)

        if rn > 0.0 and hn > 0.0:
            cosine = dot / (rn * hn)
        else:
            cosine = 0.0

        shape2 = 0.0
        for c in coords:
            rv = rr[c]
            hv = hh[c]
            ev = hv - alpha * rv
            shape2 += ev * ev

        shape_rel = math.sqrt(shape2) / max(rn, 1.0e-30)

        typ, row, col = bk

        records.append({
            "type": typ,
            "row": row,
            "col": col,
            "rom": rn,
            "hrom": hn,
            "ratio": ratio,
            "rel": rel,
            "alpha": alpha,
            "cosine": cosine,
            "shape_rel": shape_rel,
        })

        total_r2 += r2
        total_d2 += d2
        type_r2[typ] += r2
        type_d2[typ] += d2

    print()
    print("[global exact matrix errors]")
    print("all blocks rel :", math.sqrt(total_d2) / max(math.sqrt(total_r2), 1e-30))
    for typ in sorted(type_r2):
        print(
            f"{typ:7s} rel : "
            f"{math.sqrt(type_d2[typ]) / max(math.sqrt(type_r2[typ]),1e-30):.8e}"
        )

    print()
    print("[mode-wise exact relative error]")
    for m in sorted(mode_row_r2):
        er = math.sqrt(mode_row_d2[m]) / max(math.sqrt(mode_row_r2[m]), 1e-30)
        ec = math.sqrt(mode_col_d2[m]) / max(math.sqrt(mode_col_r2[m]), 1e-30)
        print(f"mode {m}: row={er:.8e} col={ec:.8e}")

    print()
    print("[largest exact block errors]")
    for r in sorted(records, key=lambda x: x["rel"], reverse=True)[:40]:
        print(
            f"{r['type']:7s} ({r['row']:3d},{r['col']:3d}) "
            f"R={r['rom']:.6e} H={r['hrom']:.6e} "
            f"H/R={r['ratio']:.4f} rel={r['rel']:.4f} "
            f"alpha={r['alpha']:.4f} cos={r['cosine']:.4f} "
            f"shape={r['shape_rel']:.4f}"
        )

    print()
    print("[largest gain alpha]")
    for r in sorted(records, key=lambda x: x["alpha"], reverse=True)[:30]:
        print(
            f"{r['type']:7s} ({r['row']:3d},{r['col']:3d}) "
            f"alpha={r['alpha']:.4f} cos={r['cosine']:.4f} "
            f"H/R={r['ratio']:.4f} rel={r['rel']:.4f} shape={r['shape_rel']:.4f}"
        )

    print()
    print("[lowest cosine blocks: modal-shape distortion]")
    nonzero = [r for r in records if r["rom"] > 1e-14 and r["hrom"] > 1e-14]
    for r in sorted(nonzero, key=lambda x: x["cosine"])[:30]:
        print(
            f"{r['type']:7s} ({r['row']:3d},{r['col']:3d}) "
            f"cos={r['cosine']:.4f} alpha={r['alpha']:.4f} "
            f"H/R={r['ratio']:.4f} rel={r['rel']:.4f} shape={r['shape_rel']:.4f}"
        )

    print()
    print("[ROM nonzero -> HROM zero]")
    lost = [r for r in records if r["rom"] > 1e-14 and r["hrom"] <= 1e-14]
    print("lost blocks:", len(lost))
    for r in sorted(lost, key=lambda x: x["rom"], reverse=True)[:40]:
        print(
            f"{r['type']:7s} ({r['row']:3d},{r['col']:3d}) "
            f"ROM={r['rom']:.6e}"
        )

    # Symmetry check at block level
    recmap = {(r["type"], r["row"], r["col"]): r for r in records}
    pair_diff = []
    for r in records:
        if r["type"] != "offdiag":
            continue
        rk = ("offdiag", r["col"], r["row"])
        if rk not in recmap:
            continue
        s = recmap[rk]
        # Compare aggregate diagnostics; exact transpose comparison would need entry transpose.
        pair_diff.append((
            abs(r["rel"] - s["rel"]),
            r["row"], r["col"], r["rel"], s["rel"]
        ))

    print()
    print("[largest directional asymmetry in block relative error]")
    for drel, row, col, a, b in sorted(pair_diff, reverse=True)[:20]:
        print(
            f"({row:3d},{col:3d}) vs ({col:3d},{row:3d}) "
            f"rel={a:.6e}/{b:.6e} delta={drel:.3e}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
