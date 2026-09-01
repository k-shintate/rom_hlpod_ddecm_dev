#!/usr/bin/env python3
import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt


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


def value(rows_by_key, degree, dt, metric):
    return rows_by_key[(degree, dt)][metric]


def plot_by_degree(rows, metric, ylabel, filename, logy=False):
    degrees = sorted({r["degree"] for r in rows})
    dts = sorted({r["dt"] for r in rows})
    by = {(r["degree"], r["dt"]): r for r in rows}

    plt.figure(figsize=(7.2, 5.2))
    for dt in dts:
        y = [value(by, d, dt, metric) for d in degrees]
        if logy:
            plt.semilogy(degrees, y, marker="o", label=f"dt={dt:g}")
        else:
            plt.plot(degrees, y, marker="o", label=f"dt={dt:g}")

    plt.xticks(degrees, [f"dG({d})" for d in degrees])
    plt.xlabel("Time discretization")
    plt.ylabel(ylabel)
    plt.title("ST-FOM final-time accuracy at fixed spatial resolution")
    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=220)
    plt.close()


def plot_vs_dt(rows, metric, ylabel, filename):
    degrees = sorted({r["degree"] for r in rows})
    by_degree = defaultdict(list)
    for r in rows:
        by_degree[r["degree"]].append(r)

    plt.figure(figsize=(7.2, 5.2))
    for degree in degrees:
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


def observed_orders(rows):
    by_degree = defaultdict(dict)
    for r in rows:
        by_degree[r["degree"]][r["dt"]] = r

    results = []
    for degree in sorted(by_degree):
        dts = sorted(by_degree[degree])
        if len(dts) < 2:
            continue
        dt_fine = dts[0]
        dt_coarse = dts[-1]
        ratio = dt_coarse / dt_fine
        if ratio <= 1.0:
            continue
        out = {"degree": degree, "dt_coarse": dt_coarse, "dt_fine": dt_fine}
        for metric in ("relative_l2", "absolute_l2", "linf"):
            e_coarse = by_degree[degree][dt_coarse][metric]
            e_fine = by_degree[degree][dt_fine][metric]
            if e_coarse > 0.0 and e_fine > 0.0:
                out[metric] = math.log(e_coarse / e_fine) / math.log(ratio)
            else:
                out[metric] = float("nan")
        results.append(out)
    return results


def write_order_csv(orders, path):
    fields = [
        "degree", "dt_coarse", "dt_fine",
        "relative_l2_order", "absolute_l2_order", "linf_order"
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for o in orders:
            writer.writerow({
                "degree": o["degree"],
                "dt_coarse": f'{o["dt_coarse"]:.16g}',
                "dt_fine": f'{o["dt_fine"]:.16g}',
                "relative_l2_order": f'{o["relative_l2"]:.16g}',
                "absolute_l2_order": f'{o["absolute_l2"]:.16g}',
                "linf_order": f'{o["linf"]:.16g}',
            })


def write_table(rows, orders, path):
    order_map = {o["degree"]: o for o in orders}
    rows = sorted(rows, key=lambda r: (r["degree"], r["dt"]))
    with open(path, "w") as f:
        f.write("# ST-FOM final-time accuracy comparison\n\n")
        f.write("Fixed spatial resolution. Only dG degree and time step are varied.\n\n")
        f.write("| degree | dt | final time | relative L2 | absolute L2 | Linf |\n")
        f.write("|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(
                f'| dG({r["degree"]}) | {r["dt"]:g} | {r["final_time"]:g} | '
                f'{r["relative_l2"]:.10e} | {r["absolute_l2"]:.10e} | {r["linf"]:.10e} |\n'
            )
        f.write("\n## Observed time-convergence order (dt=0.05 -> 0.01)\n\n")
        f.write("| degree | relative L2 order | absolute L2 order | Linf order |\n")
        f.write("|---:|---:|---:|---:|\n")
        for degree in sorted(order_map):
            o = order_map[degree]
            f.write(
                f'| dG({degree}) | {o["relative_l2"]:.6f} | '
                f'{o["absolute_l2"]:.6f} | {o["linf"]:.6f} |\n'
            )


def print_summary(rows, orders):
    print("\nFinal-time accuracy:")
    print("degree   dt       relative_L2        absolute_L2        Linf")
    for r in sorted(rows, key=lambda r: (r["degree"], r["dt"])):
        print(
            f'dG({r["degree"]})   {r["dt"]:<7g} '
            f'{r["relative_l2"]: .8e}  {r["absolute_l2"]: .8e}  {r["linf"]: .8e}'
        )
    print("\nObserved order from dt=0.05 to 0.01:")
    for o in orders:
        print(
            f'dG({o["degree"]}): p_relL2={o["relative_l2"]:.6f}, '
            f'p_absL2={o["absolute_l2"]:.6f}, p_Linf={o["linf"]:.6f}'
        )


def main():
    if len(sys.argv) != 3:
        print("Usage: plot_st_accuracy_compare_v8.py <summary_csv> <output_dir>", file=sys.stderr)
        return 2

    summary_csv = sys.argv[1]
    outdir = sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    rows = read_rows(summary_csv)
    expected = {(d, dt) for d in (0, 1, 2) for dt in (0.01, 0.05)}
    actual = {(r["degree"], round(r["dt"], 12)) for r in rows}
    expected_rounded = {(d, round(dt, 12)) for d, dt in expected}
    missing = expected_rounded - actual
    if missing:
        raise SystemExit(f"missing comparison cases: {sorted(missing)}")

    plot_by_degree(
        rows, "relative_l2", "Final relative L2 error",
        os.path.join(outdir, "final_relative_l2_by_degree.png"), False)
    plot_by_degree(
        rows, "relative_l2", "Final relative L2 error",
        os.path.join(outdir, "final_relative_l2_by_degree_log.png"), True)
    plot_by_degree(
        rows, "absolute_l2", "Final absolute L2 error",
        os.path.join(outdir, "final_absolute_l2_by_degree.png"), True)
    plot_by_degree(
        rows, "linf", "Final Linf error",
        os.path.join(outdir, "final_linf_by_degree.png"), True)
    plot_vs_dt(
        rows, "relative_l2", "Final relative L2 error",
        os.path.join(outdir, "final_relative_l2_vs_dt.png"))

    orders = observed_orders(rows)
    write_order_csv(orders, os.path.join(outdir, "observed_time_order.csv"))
    write_table(rows, orders, os.path.join(outdir, "final_accuracy_table.md"))
    print_summary(rows, orders)

    print("\nGenerated plots:")
    for name in (
        "final_relative_l2_by_degree.png",
        "final_relative_l2_by_degree_log.png",
        "final_absolute_l2_by_degree.png",
        "final_linf_by_degree.png",
        "final_relative_l2_vs_dt.png",
    ):
        print(os.path.join(outdir, name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

