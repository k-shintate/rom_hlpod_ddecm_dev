#!/usr/bin/env python3
# convert_scalarize.py
import sys

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} input.txt", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]

    with open(path, "r", encoding="utf-8") as f:
        # 空行は無視
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        print("Empty file", file=sys.stderr)
        sys.exit(1)

    # 1行目: ID (例: 24684)
    first_tokens = lines[0].split()
    if len(first_tokens) != 1:
        print("First line must contain a single integer (e.g., 24684).", file=sys.stderr)
        sys.exit(1)
    header_id = first_tokens[0]  # 文字列のまま保持

    # 2行目以降: N行 × dim列
    data_lines = lines[1:]
    if not data_lines:
        print("No data lines found.", file=sys.stderr)
        sys.exit(1)

    dim = len(data_lines[0].split())
    if dim <= 0:
        print("Invalid data line.", file=sys.stderr)
        sys.exit(1)

    # ヘッダ出力: "24684 3"
    print(f"{header_id} {dim}")

    # 各行 i の各成分 j を 1行ずつ出力: "i j value"
    for i, ln in enumerate(data_lines):
        cols = ln.split()
        if len(cols) != dim:
            print(f"Column count mismatch at line {i+2}: expected {dim}, got {len(cols)}",
                  file=sys.stderr)
            sys.exit(1)

        for j, v in enumerate(cols):
            print(f"{i} {j} {v}")

if __name__ == "__main__":
    main()
