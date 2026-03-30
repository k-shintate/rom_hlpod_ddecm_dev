#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from typing import List


def parse_ints(line: str) -> List[int]:
    parts = line.strip().split()
    return [int(p) for p in parts]


def convert(in_path: str, out_path: str) -> None:
    with open(in_path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        raise ValueError("入力ファイルが空です。")

    # ---- header ----
    header = parse_ints(lines[0])
    if len(header) == 1:
        header_id = header[0]
        # group_id は最初のデータ行の2列目から取る（例: 8）
        if len(lines) < 2:
            raise ValueError("ヘッダ行の次にデータ行がありません。")
        first_data = parse_ints(lines[1])
        if len(first_data) < 2:
            raise ValueError("最初のデータ行に group_id（例: 8）がありません。")
        group_id = first_data[1]
    else:
        # 例: "1047 8" のように2つ以上あれば先頭2つを採用
        header_id, group_id = header[0], header[1]

    # ---- write ----
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(f"{header_id} {group_id}\n")

        # 各データ行: 先頭ID と group_id(8) を落として残りだけ出力
        for idx, line in enumerate(lines[1:], start=2):
            vals = parse_ints(line)

            # 想定: "1290 8 0 1 3 2 ..."
            if len(vals) < 3:
                raise ValueError(f"{idx}行目のデータが短すぎます: {line!r}")

            # 念のため2列目が group_id であることを軽くチェック（違っても処理は可能だが警告）
            if vals[1] != group_id:
                print(
                    f"Warning: {idx}行目の2列目({vals[1]})がヘッダのgroup_id({group_id})と異なります。",
                    file=sys.stderr,
                )

            converted = vals[2:]  # drop [ID, group_id]
            out.write(" ".join(map(str, converted)) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="変換: ヘッダを '<id> <group_id>' にし、各データ行の先頭2列(ID, group_id)を削除して出力"
    )
    parser.add_argument("-i", "--input", required=True, help="入力ファイルパス")
    parser.add_argument("-o", "--output", required=True, help="出力ファイルパス")
    args = parser.parse_args()

    try:
        convert(args.input, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
