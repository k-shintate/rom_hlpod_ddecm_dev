#!/usr/bin/env python3
"""
Compare ROM reference blocks and HROM blocks.

This script expects:
  DDECM/rom_block_reference.<rank>.csv
  DDECM/rom_block_entries.<rank>.csv

For exact scalar comparison, HROM must also export:
  DDECM/hrom_block_entries.<rank>.csv

The HROM entry CSV should contain at least:
  rank,type,row_subdomain_id,col_subdomain_id,
  row_global_id,col_global_id,ir,ic,value

If hrom_block_entries are not present, the script still summarizes the
ROM structural mapping and reports that scalar comparison is unavailable.

Usage:
  python3 compare_rom_hrom_blocks.py CASE_DIR
"""

import csv
import glob
import math
import os
import sys
from collections import defaultdict


def load(pattern):
    rows = []
    files = sorted(glob.glob(pattern))
    for fn in files:
        with open(fn, newline="") as f:
            rows.extend(csv.DictReader(f))
    return files, rows


def f(x):
    try:
        return float(x)
    except Exception:
        return 0.0


def i(x):
    try:
        return int(x)
    except Exception:
        return 0


def entry_key(r, namespace):
    if namespace == "sid":
        a = i(r["row_subdomain_id"])
        b = i(r["col_subdomain_id"])
    else:
        a = i(r["row_global_id"])
        b = i(r["col_global_id"])

    return (
        r["type"],
        a,
        b,
        i(r["ir"]),
        i(r["ic"]),
    )


def compare(case):
    d = os.path.join(case, "DDECM")

    rom_block_files, rom_blocks = load(
        os.path.join(d, "rom_block_reference.*.csv")
    )
    rom_entry_files, rom_entries = load(
        os.path.join(d, "rom_block_entries.*.csv")
    )
    hrom_entry_files, hrom_entries = load(
        os.path.join(d, "hrom_block_entries.*.csv")
    )

    print("=" * 78)
    print(case)
    print("=" * 78)
    print("ROM block files       :", len(rom_block_files))
    print("ROM entry files       :", len(rom_entry_files))
    print("HROM entry files      :", len(hrom_entry_files))

    missing = sum(i(r.get("match_count", "0")) == 0
                  for r in rom_blocks if r.get("type") == "offdiag")
    multi = sum(i(r.get("match_count", "0")) > 1
                for r in rom_blocks if r.get("type") == "offdiag")

    print("ROM missing mappings  :", missing)
    print("ROM multiple mappings :", multi)

    if not hrom_entries:
        print()
        print("HROM scalar entry CSV is not present.")
        print("Add scalar export at the HROM reduced_mat assembly and rerun.")
        return

    for namespace in ("sid", "gid"):
        rom = {}
        hrom = {}

        for r in rom_entries:
            rom[entry_key(r, namespace)] = f(r["value"])

        for r in hrom_entries:
            hrom[entry_key(r, namespace)] = f(r["value"])

        common = sorted(set(rom) & set(hrom))
        only_rom = set(rom) - set(hrom)
        only_hrom = set(hrom) - set(rom)

        diff_sq = defaultdict(float)
        ref_sq = defaultdict(float)
        hrom_sq = defaultdict(float)

        for key in common:
            typ, row, col, ir, ic = key
            bkey = (typ, row, col)

            vr = rom[key]
            vh = hrom[key]
            dval = vh - vr

            diff_sq[bkey] += dval * dval
            ref_sq[bkey] += vr * vr
            hrom_sq[bkey] += vh * vh

        bad = []
        for key in ref_sq:
            rn = math.sqrt(ref_sq[key])
            hn = math.sqrt(hrom_sq[key])
            dn = math.sqrt(diff_sq[key])
            rel = dn / max(rn, 1e-30)
            bad.append((rel, key, rn, hn, dn))

        bad.sort(reverse=True)

        print()
        print(f"[namespace={namespace}]")
        print("common scalar entries :", len(common))
        print("ROM-only entries      :", len(only_rom))
        print("HROM-only entries     :", len(only_hrom))

        print("worst block relative errors:")
        for rel, key, rn, hn, dn in bad[:30]:
            typ, row, col = key
            print(
                f"  {typ:7s} ({row},{col}) "
                f"ROM={rn:.6e} HROM={hn:.6e} "
                f"diff={dn:.6e} rel={rel:.6e}"
            )


def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 2

    compare(sys.argv[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
