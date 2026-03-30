#!/usr/bin/env python3
# reverse_convert_scalarize.py
import sys

def main():
    if len(sys.argv) not in (2, 3):
        print(f"Usage: {sys.argv[0]} input_scalar.txt [output.txt]", file=sys.stderr)
        print(f"Example: {sys.argv[0]} output.txt > restored.txt", file=sys.stderr)
        sys.exit(1)

    in_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) == 3 else None

    with open(in_path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        print("Empty file", file=sys.stderr)
        sys.exit(1)

    # 1行目: "<header_id> <dim>"
    header = lines[0].split()
    if len(header) != 2:
        print("First line must be: '<header_id> <dim>'", file=sys.stderr)
        sys.exit(1)

    header_id = header[0]
    try:
        dim = int(header[1])
    except ValueError:
        print("dim must be an integer", file=sys.stderr)
        sys.exit(1)

    if dim <= 0:
        print("dim must be > 0", file=sys.stderr)
        sys.exit(1)

    # 2行目以降: "i j value"
    data = {}  # data[i][j] = value(str)
    for line_no, ln in enumerate(lines[1:], start=2):
        # value に余計な空白があっても壊れないよう maxsplit=2
        parts = ln.split(maxsplit=2)
        if len(parts) != 3:
            print(f"Invalid line {line_no}: '{ln}' (need: i j value)", file=sys.stderr)
            sys.exit(1)

        try:
            i = int(parts[0])
            j = int(parts[1])
        except ValueError:
            print(f"Invalid i/j at line {line_no}: '{ln}'", file=sys.stderr)
            sys.exit(1)

        if j < 0 or j >= dim:
            print(f"Invalid component j={j} at line {line_no} (dim={dim})", file=sys.stderr)
            sys.exit(1)

        value = parts[2]  # 文字列のまま保持（情報落ち防止）

        if i not in data:
            data[i] = {}
        if j in data[i]:
            print(f"Duplicate entry for i={i}, j={j} at line {line_no}", file=sys.stderr)
            sys.exit(1)

        data[i][j] = value

    if not data:
        print("No data lines found.", file=sys.stderr)
        sys.exit(1)

    # i を昇順に出力（convert側は 0..N-1 の想定）
    indices = sorted(data.keys())

    # 各 i について j=0..dim-1 が揃っているかチェック
    for i in indices:
        missing = [j for j in range(dim) if j not in data[i]]
        if missing:
            print(f"Missing components for i={i}: missing j={missing}", file=sys.stderr)
            sys.exit(1)

    out_lines = []
    out_lines.append(f"{header_id}")
    for i in indices:
        row = [data[i][j] for j in range(dim)]
        out_lines.append(" ".join(row))

    text = "\n".join(out_lines) + "\n"

    if out_path is None:
        sys.stdout.write(text)
    else:
        with open(out_path, "w", encoding="utf-8") as f_out:
            f_out.write(text)

if __name__ == "__main__":
    main()
