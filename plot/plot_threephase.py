#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make TRANSMAG comparison table from:
  - terminal_voltage_log.csv
  - transmag_cycle_summary.csv
  - transmag_peak_summary.txt

Usage:
    python make_transmag_comparison_table.py \
        --log terminal_voltage_log.csv \
        --cycle transmag_cycle_summary.csv \
        --peak transmag_peak_summary.txt \
        --out transmag_comparison_table.csv
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Any

import pandas as pd


def read_peak_summary(path: Path) -> Dict[str, float]:
    """
    Read key-value style summary file:
        Run_Max_B = ...
        Run_Max_B_Time = ...
        ...
    """
    data: Dict[str, float] = {}

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            try:
                if "." in value or "e" in value.lower():
                    data[key] = float(value)
                else:
                    data[key] = int(value)
            except ValueError:
                # Keep raw string if not numeric
                data[key] = value

    return data


def safe_get(series: pd.Series, key: str, default: Any = float("nan")) -> Any:
    return series[key] if key in series.index else default


def make_comparison_table(
    log_path: Path,
    cycle_path: Path,
    peak_path: Path,
) -> pd.DataFrame:
    """
    Build a comparison table focused on TRANSMAG-style quantities:
      - final-cycle RMS currents
      - final-cycle average power/loss
      - running max B and quantities at that instant
      - final-step values as reference
    """
    df_log = pd.read_csv(log_path)
    df_cycle = pd.read_csv(cycle_path)
    peak = read_peak_summary(peak_path)

    if df_log.empty:
        raise ValueError(f"{log_path} is empty.")
    if df_cycle.empty:
        raise ValueError(f"{cycle_path} is empty.")

    last_log = df_log.iloc[-1]
    last_cycle = df_cycle.iloc[-1]

    comparison_rows = [
        {
            "Category": "Final cycle",
            "Quantity": "Irms_u",
            "Value": safe_get(last_cycle, "Irms_u"),
            "Unit": "A",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の U 相 RMS 電流",
        },
        {
            "Category": "Final cycle",
            "Quantity": "Irms_v",
            "Value": safe_get(last_cycle, "Irms_v"),
            "Unit": "A",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の V 相 RMS 電流",
        },
        {
            "Category": "Final cycle",
            "Quantity": "Irms_w",
            "Value": safe_get(last_cycle, "Irms_w"),
            "Unit": "A",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の W 相 RMS 電流",
        },
        {
            "Category": "Final cycle",
            "Quantity": "Pin_avg",
            "Value": safe_get(last_cycle, "Pin_avg"),
            "Unit": "W",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の平均入力電力",
        },
        {
            "Category": "Final cycle",
            "Quantity": "LossCore_avg",
            "Value": safe_get(last_cycle, "LossCore_avg"),
            "Unit": "W",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の平均コア損",
        },
        {
            "Category": "Final cycle",
            "Quantity": "LossTotal_avg",
            "Value": safe_get(last_cycle, "LossTotal_avg"),
            "Unit": "W",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の平均総損失",
        },
        {
            "Category": "Final cycle",
            "Quantity": "CoreSatRatio_avg",
            "Value": safe_get(last_cycle, "CoreSatRatio_avg"),
            "Unit": "-",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期の平均飽和体積率",
        },
        {
            "Category": "Final cycle",
            "Quantity": "CycleMax_B",
            "Value": safe_get(last_cycle, "CycleMax_B"),
            "Unit": "T",
            "Source": "transmag_cycle_summary.csv",
            "Comment": "最後の1周期内の最大 B",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_Max_B",
            "Value": peak.get("Run_Max_B", float("nan")),
            "Unit": "T",
            "Source": "transmag_peak_summary.txt",
            "Comment": "全計算中の最大 B",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_Max_B_Time",
            "Value": peak.get("Run_Max_B_Time", float("nan")),
            "Unit": "s",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B が出た時刻",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_Max_B_Step",
            "Value": peak.get("Run_Max_B_Step", float("nan")),
            "Unit": "-",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B が出た step",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_Iu",
            "Value": peak.get("Run_MaxB_Iu", float("nan")),
            "Unit": "A",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の U 相電流",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_Iv",
            "Value": peak.get("Run_MaxB_Iv", float("nan")),
            "Unit": "A",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の V 相電流",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_Iw",
            "Value": peak.get("Run_MaxB_Iw", float("nan")),
            "Unit": "A",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の W 相電流",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_P1_B",
            "Value": peak.get("Run_MaxB_P1_B", float("nan")),
            "Unit": "T",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の Probe 1 B",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_P2_B",
            "Value": peak.get("Run_MaxB_P2_B", float("nan")),
            "Unit": "T",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の Probe 2 B",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_P3_B",
            "Value": peak.get("Run_MaxB_P3_B", float("nan")),
            "Unit": "T",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の Probe 3 B",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_CoreSatRatio",
            "Value": peak.get("Run_MaxB_CoreSatRatio", float("nan")),
            "Unit": "-",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の飽和体積率",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_LossCore",
            "Value": peak.get("Run_MaxB_LossCore", float("nan")),
            "Unit": "W",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻のコア損",
        },
        {
            "Category": "Peak saturation instant",
            "Quantity": "Run_MaxB_LossTotal",
            "Value": peak.get("Run_MaxB_LossTotal", float("nan")),
            "Unit": "W",
            "Source": "transmag_peak_summary.txt",
            "Comment": "最大 B 時刻の総損失",
        },
        {
            "Category": "Final step reference",
            "Quantity": "Final_Time",
            "Value": safe_get(last_log, "Time"),
            "Unit": "s",
            "Source": "terminal_voltage_log.csv",
            "Comment": "最終 step の時刻",
        },
        {
            "Category": "Final step reference",
            "Quantity": "Final_Max_B",
            "Value": safe_get(last_log, "Max_B"),
            "Unit": "T",
            "Source": "terminal_voltage_log.csv",
            "Comment": "最終 step の最大 B",
        },
        {
            "Category": "Final step reference",
            "Quantity": "Final_Loss_Total",
            "Value": safe_get(last_log, "Loss_Total"),
            "Unit": "W",
            "Source": "terminal_voltage_log.csv",
            "Comment": "最終 step の総損失",
        },
        {
            "Category": "Final step reference",
            "Quantity": "Final_Core_SatRatio",
            "Value": safe_get(last_log, "Core_SatRatio"),
            "Unit": "-",
            "Source": "terminal_voltage_log.csv",
            "Comment": "最終 step の飽和体積率",
        },
    ]

    out_df = pd.DataFrame(comparison_rows)
    return out_df


def print_pretty_table(df: pd.DataFrame) -> None:
    """
    Print simple console table without external deps.
    """
    if df.empty:
        print("No data.")
        return

    display_df = df.copy()
    display_df["Value"] = display_df["Value"].map(
        lambda x: f"{x:.6e}" if isinstance(x, (int, float)) and not pd.isna(x) else str(x)
    )

    widths = {}
    for col in display_df.columns:
        widths[col] = max(len(col), display_df[col].astype(str).map(len).max())

    header = " | ".join(col.ljust(widths[col]) for col in display_df.columns)
    sep = "-+-".join("-" * widths[col] for col in display_df.columns)

    print(header)
    print(sep)
    for _, row in display_df.iterrows():
        print(" | ".join(str(row[col]).ljust(widths[col]) for col in display_df.columns))


def main() -> None:
    parser = argparse.ArgumentParser(description="Create TRANSMAG comparison table.")
    parser.add_argument("--log", type=Path, default=Path("terminal_voltage_log.csv"))
    parser.add_argument("--cycle", type=Path, default=Path("transmag_cycle_summary.csv"))
    parser.add_argument("--peak", type=Path, default=Path("transmag_peak_summary.txt"))
    parser.add_argument("--out", type=Path, default=Path("transmag_comparison_table.csv"))
    args = parser.parse_args()

    for p in [args.log, args.cycle, args.peak]:
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {p}")

    df = make_comparison_table(args.log, args.cycle, args.peak)
    df.to_csv(args.out, index=False, encoding="utf-8-sig")

    print(f"[OK] Wrote comparison table: {args.out}")
    print()
    print_pretty_table(df)


if __name__ == "__main__":
    main()