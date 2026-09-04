#!/usr/bin/env python3
import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt

EXPECTED_DTS = (0.1, 0.05, 0.01)
EXPECTED_DEGREES = (0, 1, 2)


def read_rows(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            r = dict(row)
            r["degree"] = int(r["degree"])
            r["dt"] = float(r["dt"])
            r["final_step"] = int(r["final_step"])
            r["final_time"] = float(r["final_time"])
            r["relative_l2"] = float(r["relative_l2"])
            r["absolute_l2"] = float(r["absolute_l2"])
            r["linf"] = float(r["linf"])
            rows.append(r)
    return rows


def plot_by_degree(rows, metric, ylabel, filename):
    degrees = sorted({r["degree"] for r in rows})
    dts = sorted({r["dt"] for r in rows}, reverse=True)
    by = {(r["degree"], round(r["dt"], 12)): r for r in rows}

    plt.figure(figsize=(7.4, 5.4))
    for dt in dts:
        y = [by[(d, round(dt, 12))][metric] for d in degrees]
        plt.plot(degrees, y, marker="o", label=f"dt={dt:g}")

    plt.xticks(degrees, [f"dG({d})" for d in degrees])
    plt.xlabel("Time discretization")
    plt.ylabel(ylabel)
    plt.yscale("log")
    plt.title("ST-FOM final-time accuracy at fixed spatial resolution")
    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=220)
    plt.close()


def plot_vs_dt(rows, metric, ylabel, filename):
    by_degree = defaultdict(list)
    for r in rows:
        by_degree[r["degree"]].append(r)

    plt.figure(figsize=(7.4, 5.4))
    for degree in sorted(by_degree):
        pts = sorted(by_degree[degree], key=lambda r: r["dt"])
        x = [r["dt"] for r in pts]
        y = [r[metric] for r in pts]
        plt.loglog(x, y, marker="o", label=f"dG({degree})")

    plt.xlabel("Time step dt")
    plt.ylabel(ylabel)
    plt.title("ST-FOM final-time temporal accuracy")
    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=220)
    plt.close()


def pairwise_orders(rows):
    by_degree = defaultdict(dict)
    for r in rows:
        by_degree[r["degree"]][round(r["dt"], 12)] = r

    pairs = [(0.1, 0.05), (0.05, 0.01)]
    out = []
    for degree in sorted(by_degree):
        for coarse, fine in pairs:
            rc = by_degree[degree][round(coarse, 12)]
            rf = by_degree[degree][round(fine, 12)]
            ratio = coarse / fine
            row = {
                "degree": degree,
                "dt_coarse": coarse,
                "dt_fine": fine,
            }
            for metric in ("relative_l2", "absolute_l2", "linf"):
                ec = rc[metric]
                ef = rf[metric]
                row[metric] = (
                    math.log(ec / ef) / math.log(ratio)
                    if ec > 0.0 and ef > 0.0
                    else float("nan")
                )
            out.append(row)
    return out


def write_orders(orders, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "degree", "dt_coarse", "dt_fine",
            "relative_l2_order", "absolute_l2_order", "linf_order"
        ])
        for o in orders:
            w.writerow([
                o["degree"], f'{o["dt_coarse"]:.16g}', f'{o["dt_fine"]:.16g}',
                f'{o["relative_l2"]:.16g}',
                f'{o["absolute_l2"]:.16g}',
                f'{o["linf"]:.16g}',
            ])


def write_table(rows, orders, path):
    with open(path, "w") as f:
        f.write("# ST-FOM final-time accuracy comparison\n\n")
        f.write("Fixed spatial mesh and MPI count; only dG degree and dt are varied.\n\n")
        f.write("| degree | dt | final time | relative L2 | absolute L2 | Linf |\n")
        f.write("|---:|---:|---:|---:|---:|---:|\n")
        for r in sorted(rows, key=lambda r: (r["degree"], -r["dt"])):
            f.write(
                f'| dG({r["degree"]}) | {r["dt"]:g} | {r["final_time"]:g} | '
                f'{r["relative_l2"]:.10e} | {r["absolute_l2"]:.10e} | {r["linf"]:.10e} |\n'
            )

        f.write("\n## Pairwise observed temporal order\n\n")
        f.write("| degree | dt coarse | dt fine | p(relative L2) | p(absolute L2) | p(Linf) |\n")
        f.write("|---:|---:|---:|---:|---:|---:|\n")
        for o in orders:
            f.write(
                f'| dG({o["degree"]}) | {o["dt_coarse"]:g} | {o["dt_fine"]:g} | '
                f'{o["relative_l2"]:.6f} | {o["absolute_l2"]:.6f} | {o["linf"]:.6f} |\n'
            )


def print_summary(rows, orders):
    print("\nFinal-time accuracy:")
    print("degree   dt       relative_L2        absolute_L2        Linf")
    for r in sorted(rows, key=lambda r: (r["degree"], -r["dt"])):
        print(
            f'dG({r["degree"]})   {r["dt"]:<7g} '
            f'{r["relative_l2"]: .8e}  {r["absolute_l2"]: .8e}  {r["linf"]: .8e}'
        )

    print("\nPairwise observed temporal order:")
    for o in orders:
        print(
            f'dG({o["degree"]}) {o["dt_coarse"]:g}->{o["dt_fine"]:g}: '
            f'p_relL2={o["relative_l2"]:.6f}, '
            f'p_absL2={o["absolute_l2"]:.6f}, '
            f'p_Linf={o["linf"]:.6f}'
        )


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: plot_st_accuracy_compare_v10_dt3.py <summary_csv> <output_dir>",
            file=sys.stderr,
        )
        return 2

    summary_csv, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    rows = read_rows(summary_csv)

    expected = {
        (d, round(dt, 12))
        for d in EXPECTED_DEGREES
        for dt in EXPECTED_DTS
    }
    actual = {(r["degree"], round(r["dt"], 12)) for r in rows}
    missing = expected - actual
    if missing:
        raise SystemExit(f"missing comparison cases: {sorted(missing)}")

    plot_by_degree(
        rows, "relative_l2", "Final relative L2 error",
        os.path.join(outdir, "final_relative_l2_by_degree.png"))
    plot_by_degree(
        rows, "absolute_l2", "Final absolute L2 error",
        os.path.join(outdir, "final_absolute_l2_by_degree.png"))
    plot_by_degree(
        rows, "linf", "Final Linf error",
        os.path.join(outdir, "final_linf_by_degree.png"))

    plot_vs_dt(
        rows, "relative_l2", "Final relative L2 error",
        os.path.join(outdir, "final_relative_l2_vs_dt.png"))
    plot_vs_dt(
        rows, "absolute_l2", "Final absolute L2 error",
        os.path.join(outdir, "final_absolute_l2_vs_dt.png"))
    plot_vs_dt(
        rows, "linf", "Final Linf error",
        os.path.join(outdir, "final_linf_vs_dt.png"))

    orders = pairwise_orders(rows)
    write_orders(
        orders, os.path.join(outdir, "observed_time_order_pairwise.csv"))
    write_table(
        rows, orders, os.path.join(outdir, "final_accuracy_table.md"))
    print_summary(rows, orders)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
