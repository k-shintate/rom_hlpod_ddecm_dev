#!/usr/bin/env python3
"""
plot_hrom_prm1_sweep.py

Collect and compare HROM Stage-3 runs stored as

    hrom_stage3_results/
      prm1_1.0e-7_prm2_0/
        DDECM/stage_reduction_summary.txt
        l2_error_fem_hrom.txt
        l2_error_rom.txt              # optional
      prm1_1.0e-7_prm2_4/
        ...

Typical use
-----------
python plot_hrom_prm1_sweep.py hrom_stage3_results \
    --modes 0 4 \
    --expected-final-time 0.2 \
    --out paper_prm1_sweep

Outputs
-------
summary_all.csv
paired_mode0_mode4.csv
paired_mode0_mode4.tex
prm1_vs_stage3_selected.png
prm1_vs_cumulative_retained.png
prm1_vs_normalized_residual.png
prm1_vs_l2_max.png
prm1_vs_l2_rms.png
nsel_vs_l2_max.png
nsel_vs_offline_residual.png
l2_error_history_all.png

Notes
-----
* prm1 is read from the directory name.
* prm2 is read from the directory name, but the actual coverage_mode saved in
  stage_reduction_summary.txt is treated as authoritative.
* Duplicate (time,error) rows in l2_error_*.txt are removed automatically.
* A run is marked "reached_final_time" only if the error history reaches
  --expected-final-time.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


DEFAULT_MODE_LABELS = {
    0: "Stage-3 NNLS",
    1: "Expanded candidates",
    2: "Physical interface coverage",
    3: "Interface operator coverage",
    4: "Interface + internal operator coverage",
    5: "Adaptive operator enrichment",
}


def parse_run_dir_name(name: str) -> Optional[Tuple[float, int]]:
    m = re.fullmatch(rf"prm1_({FLOAT_RE})_prm2_(-?\d+)", name)
    if not m:
        return None
    return float(m.group(1)), int(m.group(2))


def read_stage_summary(path: Path) -> Tuple[Optional[int], str, Dict[str, float]]:
    text = path.read_text(errors="replace")

    mode: Optional[int] = None
    desc = ""

    m = re.search(r"(?m)^coverage_mode\s+(\d+)\s*$", text)
    if m:
        mode = int(m.group(1))

    m = re.search(r"(?m)^coverage_mode_description\s+(.+?)\s*$", text)
    if m:
        desc = m.group(1).strip()

    metrics: Dict[str, float] = {}
    m = re.search(r"\[Global reduction\]\s*\n(.*)", text, flags=re.S)
    if m:
        for line in m.group(1).splitlines():
            parts = line.strip().split(maxsplit=1)
            if len(parts) != 2:
                continue
            key, value = parts
            try:
                metrics[key] = float(value)
            except ValueError:
                # Stop only if a new bracketed section is encountered.
                if line.lstrip().startswith("["):
                    break

    return mode, desc, metrics


def read_error_history(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["time", "error"])

    rows: List[Tuple[float, float]] = []
    pat = re.compile(rf"^\s*({FLOAT_RE})\s+({FLOAT_RE})\s*$")

    for line in path.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if not m:
            continue
        t = float(m.group(1))
        e = float(m.group(2))
        if math.isfinite(t) and math.isfinite(e):
            rows.append((t, e))

    if not rows:
        return pd.DataFrame(columns=["time", "error"])

    df = pd.DataFrame(rows, columns=["time", "error"])
    # The current output can contain exact duplicate rows.
    df = df.drop_duplicates(subset=["time", "error"])
    # If the same time exists with different values, keep the last written one.
    df = df.sort_values(["time"]).drop_duplicates(subset=["time"], keep="last")
    return df.reset_index(drop=True)


def error_stats(df: pd.DataFrame, prefix: str) -> Dict[str, float]:
    if df.empty:
        return {
            f"{prefix}_points": 0,
            f"{prefix}_t_final": np.nan,
            f"{prefix}_max": np.nan,
            f"{prefix}_rms": np.nan,
            f"{prefix}_final": np.nan,
        }

    e = df["error"].to_numpy(dtype=float)
    t = df["time"].to_numpy(dtype=float)
    return {
        f"{prefix}_points": int(len(df)),
        f"{prefix}_t_final": float(np.max(t)),
        f"{prefix}_max": float(np.max(e)),
        f"{prefix}_rms": float(np.sqrt(np.mean(e * e))),
        f"{prefix}_final": float(df.sort_values("time").iloc[-1]["error"]),
    }


def discover_cases(root: Path, modes: Optional[Iterable[int]]) -> Tuple[pd.DataFrame, Dict[Tuple[float, int], pd.DataFrame]]:
    requested = None if modes is None else set(int(x) for x in modes)

    rows: List[Dict[str, object]] = []
    histories: Dict[Tuple[float, int], pd.DataFrame] = {}

    for d in sorted(root.rglob("prm1_*_prm2_*")):
        if not d.is_dir():
            continue

        parsed = parse_run_dir_name(d.name)
        if parsed is None:
            continue
        prm1, prm2 = parsed

        summary_file = d / "DDECM" / "stage_reduction_summary.txt"
        if not summary_file.exists():
            continue

        coverage_mode, coverage_description, metrics = read_stage_summary(summary_file)
        mode = prm2 if coverage_mode is None else coverage_mode

        if requested is not None and mode not in requested:
            continue

        if coverage_mode is not None and coverage_mode != prm2:
            print(
                f"WARNING: directory {d.name} says prm2={prm2}, "
                f"but summary says coverage_mode={coverage_mode}. "
                "Using coverage_mode from the summary.",
                file=sys.stderr,
            )

        hrom = read_error_history(d / "l2_error_fem_hrom.txt")
        rom = read_error_history(d / "l2_error_rom.txt")
        histories[(prm1, mode)] = hrom

        row: Dict[str, object] = {
            "case_dir": str(d),
            "prm1": prm1,
            "prm2": prm2,
            "coverage_mode": mode,
            "coverage_description": coverage_description,
            "mode_label": DEFAULT_MODE_LABELS.get(mode, f"Mode {mode}"),
        }
        row.update(metrics)
        row.update(error_stats(hrom, "hrom_l2"))
        row.update(error_stats(rom, "rom_l2"))
        rows.append(row)

    if not rows:
        return pd.DataFrame(), histories

    df = pd.DataFrame(rows).sort_values(["coverage_mode", "prm1"]).reset_index(drop=True)
    return df, histories


def add_online_status(df: pd.DataFrame, expected_final_time: Optional[float]) -> pd.DataFrame:
    df = df.copy()
    if expected_final_time is None:
        df["reached_final_time"] = "unknown"
        return df

    tol = max(1.0e-12, abs(expected_final_time) * 1.0e-8)
    t = pd.to_numeric(df["hrom_l2_t_final"], errors="coerce")
    df["reached_final_time"] = t.notna() & (t >= expected_final_time - tol)
    return df


def mode_style(mode: int) -> Tuple[str, str]:
    # Keep color assignment automatic.  Distinguish methods by line/marker only.
    if mode == 0:
        return "--", "o"
    if mode == 4:
        return "-", "s"
    marker_cycle = ["^", "D", "v", "P", "X", "*"]
    return "-", marker_cycle[mode % len(marker_cycle)]


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=250, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_prm1_metric(
    df: pd.DataFrame,
    out: Path,
    ycol: str,
    ylabel: str,
    filename: str,
    logy: bool = False,
    prm1_label: str = "prm1",
) -> None:
    if ycol not in df.columns:
        return

    tmp = df.copy()
    tmp["__x"] = pd.to_numeric(tmp["prm1"], errors="coerce")
    tmp["__y"] = pd.to_numeric(tmp[ycol], errors="coerce")
    tmp = tmp.dropna(subset=["__x", "__y"])
    if tmp.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    for mode, g in tmp.groupby("coverage_mode"):
        g = g.sort_values("__x")
        linestyle, marker = mode_style(int(mode))
        label = str(g.iloc[0]["mode_label"])
        ax.plot(
            g["__x"],
            g["__y"],
            linestyle=linestyle,
            marker=marker,
            linewidth=1.8,
            markersize=6,
            label=label,
        )

    if np.all(tmp["__x"] > 0):
        ax.set_xscale("log")
    if logy and np.all(tmp["__y"] > 0):
        ax.set_yscale("log")

    ax.set_xlabel(prm1_label)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save_figure(fig, out / filename)


def plot_tradeoff(
    df: pd.DataFrame,
    out: Path,
    xcol: str,
    ycol: str,
    xlabel: str,
    ylabel: str,
    filename: str,
    logy: bool = True,
) -> None:
    if xcol not in df.columns or ycol not in df.columns:
        return

    tmp = df.copy()
    tmp["__x"] = pd.to_numeric(tmp[xcol], errors="coerce")
    tmp["__y"] = pd.to_numeric(tmp[ycol], errors="coerce")
    tmp = tmp.dropna(subset=["__x", "__y"])
    if tmp.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    for mode, g in tmp.groupby("coverage_mode"):
        g = g.sort_values("prm1")
        linestyle, marker = mode_style(int(mode))
        label = str(g.iloc[0]["mode_label"])
        ax.plot(
            g["__x"],
            g["__y"],
            linestyle=linestyle,
            marker=marker,
            linewidth=1.5,
            markersize=6,
            label=label,
        )

        for _, r in g.iterrows():
            ax.annotate(
                f"{float(r['prm1']):.0e}",
                (r["__x"], r["__y"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
            )

    if logy and np.all(tmp["__y"] > 0):
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save_figure(fig, out / filename)


def plot_all_error_histories(
    df: pd.DataFrame,
    histories: Dict[Tuple[float, int], pd.DataFrame],
    out: Path,
) -> None:
    if not histories:
        return

    prm_values = sorted(df["prm1"].dropna().unique())
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    prm_marker = {p: markers[i % len(markers)] for i, p in enumerate(prm_values)}

    fig, ax = plt.subplots(figsize=(8.4, 5.2))

    for _, r in df.sort_values(["coverage_mode", "prm1"]).iterrows():
        prm1 = float(r["prm1"])
        mode = int(r["coverage_mode"])
        hist = histories.get((prm1, mode))
        if hist is None or hist.empty:
            continue

        linestyle, _ = mode_style(mode)
        label = f"{r['mode_label']}, prm1={prm1:.0e}"

        ax.plot(
            hist["time"],
            hist["error"],
            linestyle=linestyle,
            marker=prm_marker[prm1],
            markevery=max(1, len(hist) // 8),
            markersize=4,
            linewidth=1.3,
            label=label,
        )

    ax.set_xlabel("Time")
    ax.set_ylabel("FEM-HROM relative L2 error")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    save_figure(fig, out / "l2_error_history_all.png")


def make_paired_table(df: pd.DataFrame, mode_a: int, mode_b: int) -> pd.DataFrame:
    keep = [
        "prm1",
        "coverage_mode",
        "stage3_selected",
        "cumulative_retained_percent",
        "normalized_residual",
        "hrom_l2_max",
        "hrom_l2_rms",
        "hrom_l2_final",
        "hrom_l2_t_final",
        "reached_final_time",
    ]
    work = df[[c for c in keep if c in df.columns]].copy()
    work = work[work["coverage_mode"].isin([mode_a, mode_b])]

    if work.empty:
        return pd.DataFrame()

    a = work[work["coverage_mode"] == mode_a].drop(columns=["coverage_mode"]).set_index("prm1")
    b = work[work["coverage_mode"] == mode_b].drop(columns=["coverage_mode"]).set_index("prm1")

    paired = a.add_prefix(f"m{mode_a}_").join(
        b.add_prefix(f"m{mode_b}_"),
        how="outer",
    ).sort_index()

    n0 = pd.to_numeric(paired.get(f"m{mode_a}_stage3_selected"), errors="coerce")
    n1 = pd.to_numeric(paired.get(f"m{mode_b}_stage3_selected"), errors="coerce")
    e0 = pd.to_numeric(paired.get(f"m{mode_a}_hrom_l2_max"), errors="coerce")
    e1 = pd.to_numeric(paired.get(f"m{mode_b}_hrom_l2_max"), errors="coerce")
    r0 = pd.to_numeric(paired.get(f"m{mode_a}_normalized_residual"), errors="coerce")
    r1 = pd.to_numeric(paired.get(f"m{mode_b}_normalized_residual"), errors="coerce")

    paired["delta_Nsel"] = n1 - n0
    paired["delta_Nsel_percent_vs_mode0"] = 100.0 * (n1 - n0) / n0
    paired["L2max_ratio_mode4_over_mode0"] = e1 / e0
    paired["offline_residual_ratio_mode4_over_mode0"] = r1 / r0

    return paired.reset_index()


def write_paired_latex(df: pd.DataFrame, path: Path, mode_a: int, mode_b: int) -> None:
    if df.empty:
        return

    def fmt(x, kind="sci") -> str:
        if x is None or pd.isna(x):
            return "--"
        if kind == "int":
            return str(int(round(float(x))))
        if kind == "pct":
            return f"{float(x):.2f}"
        return f"{float(x):.3e}"

    lines = [
        r"\begin{table}[tb]",
        r"\centering",
        rf"\caption{{Stage-3 parameter sweep: {DEFAULT_MODE_LABELS.get(mode_a, f'Mode {mode_a}')} versus {DEFAULT_MODE_LABELS.get(mode_b, f'Mode {mode_b}')}.}}",
        r"\label{tab:stage3_prm1_sweep}",
        r"\begin{tabular}{rrrrrrrr}",
        r"\toprule",
        rf"$\mathrm{{prm1}}$ & $N_{{\rm sel}}^{{({mode_a})}}$ & "
        rf"$N_{{\rm sel}}^{{({mode_b})}}$ & "
        rf"$E_{{\rm off}}^{{({mode_a})}}$ & $E_{{\rm off}}^{{({mode_b})}}$ & "
        rf"$E_{{L^2,\max}}^{{({mode_a})}}$ & $E_{{L^2,\max}}^{{({mode_b})}}$ & "
        rf"$\Delta N_{{\rm sel}}$ \\",
        r"\midrule",
    ]

    for _, r in df.iterrows():
        lines.append(
            f"{fmt(r.get('prm1'))} & "
            f"{fmt(r.get(f'm{mode_a}_stage3_selected'), 'int')} & "
            f"{fmt(r.get(f'm{mode_b}_stage3_selected'), 'int')} & "
            f"{fmt(r.get(f'm{mode_a}_normalized_residual'))} & "
            f"{fmt(r.get(f'm{mode_b}_normalized_residual'))} & "
            f"{fmt(r.get(f'm{mode_a}_hrom_l2_max'))} & "
            f"{fmt(r.get(f'm{mode_b}_hrom_l2_max'))} & "
            f"{fmt(r.get('delta_Nsel'), 'int')} \\\\"
        )

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compare multiple prm1 values and coverage modes in HROM Stage-3 results."
    )
    ap.add_argument("root", type=Path, help="root containing prm1_*_prm2_* result directories")
    ap.add_argument("--out", type=Path, default=Path("paper_prm1_sweep"))
    ap.add_argument("--modes", nargs="+", type=int, default=[0, 4])
    ap.add_argument("--expected-final-time", type=float, default=None)
    ap.add_argument(
        "--prm1-label",
        default="Stage-3 NNLS tolerance (prm1)",
        help="x-axis label used in plots",
    )
    ap.add_argument("--pair", nargs=2, type=int, default=[0, 4], metavar=("MODE_A", "MODE_B"))
    args = ap.parse_args()

    if not args.root.exists():
        ap.error(f"root does not exist: {args.root}")

    df, histories = discover_cases(args.root, args.modes)
    if df.empty:
        print(
            "No cases found. Expected directories like "
            "prm1_1.0e-7_prm2_0/DDECM/stage_reduction_summary.txt",
            file=sys.stderr,
        )
        return 2

    df = add_online_status(df, args.expected_final_time)

    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "summary_all.csv", index=False)

    plot_prm1_metric(
        df, args.out,
        "stage3_selected",
        "Selected Stage-3 elements",
        "prm1_vs_stage3_selected.png",
        logy=False,
        prm1_label=args.prm1_label,
    )
    plot_prm1_metric(
        df, args.out,
        "cumulative_retained_percent",
        "Retained physical elements [%]",
        "prm1_vs_cumulative_retained.png",
        logy=False,
        prm1_label=args.prm1_label,
    )
    plot_prm1_metric(
        df, args.out,
        "normalized_residual",
        "Normalized Stage-3 residual",
        "prm1_vs_normalized_residual.png",
        logy=True,
        prm1_label=args.prm1_label,
    )
    plot_prm1_metric(
        df, args.out,
        "hrom_l2_max",
        "Maximum FEM-HROM relative L2 error",
        "prm1_vs_l2_max.png",
        logy=True,
        prm1_label=args.prm1_label,
    )
    plot_prm1_metric(
        df, args.out,
        "hrom_l2_rms",
        "RMS FEM-HROM relative L2 error",
        "prm1_vs_l2_rms.png",
        logy=True,
        prm1_label=args.prm1_label,
    )

    plot_tradeoff(
        df, args.out,
        "stage3_selected",
        "hrom_l2_max",
        "Selected Stage-3 elements",
        "Maximum FEM-HROM relative L2 error",
        "nsel_vs_l2_max.png",
        logy=True,
    )
    plot_tradeoff(
        df, args.out,
        "stage3_selected",
        "normalized_residual",
        "Selected Stage-3 elements",
        "Normalized Stage-3 residual",
        "nsel_vs_offline_residual.png",
        logy=True,
    )
    plot_all_error_histories(df, histories, args.out)

    mode_a, mode_b = args.pair
    paired = make_paired_table(df, mode_a, mode_b)
    if not paired.empty:
        paired.to_csv(args.out / f"paired_mode{mode_a}_mode{mode_b}.csv", index=False)
        write_paired_latex(
            paired,
            args.out / f"paired_mode{mode_a}_mode{mode_b}.tex",
            mode_a,
            mode_b,
        )

    print(f"Found {len(df)} cases.")
    print(df[
        [
            "prm1",
            "coverage_mode",
            "stage3_selected",
            "cumulative_retained_percent",
            "normalized_residual",
            "hrom_l2_max",
            "hrom_l2_rms",
            "hrom_l2_t_final",
            "reached_final_time",
        ]
    ].to_string(index=False))
    print(f"\nOutput directory: {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
