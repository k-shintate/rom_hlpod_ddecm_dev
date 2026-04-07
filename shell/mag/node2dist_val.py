#!/usr/bin/env python3
# fill_dist_val_argN.py
import sys

def is_int_single_token(line: str) -> bool:
    parts = line.split()
    if len(parts) != 1:
        return False
    try:
        int(parts[0])
        return True
    except ValueError:
        return False

def read_n_from_file(path: str) -> int:
    try:
        with open(path, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()
    except OSError as e:
        print(f"Error: failed to read N file '{path}': {e}", file=sys.stderr)
        sys.exit(1)

    if not first_line:
        print(f"Error: N file '{path}' is empty", file=sys.stderr)
        sys.exit(1)

    parts = first_line.split()
    if len(parts) != 1:
        print(f"Error: first line of N file must contain exactly one integer: '{first_line}'",
              file=sys.stderr)
        sys.exit(1)

    try:
        n = int(parts[0])
    except ValueError:
        print(f"Error: first line of N file is not an integer: '{first_line}'",
              file=sys.stderr)
        sys.exit(1)

    if n < 0:
        print("Error: N must be non-negative", file=sys.stderr)
        sys.exit(1)

    return n

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.txt output.txt n_file.txt", file=sys.stderr)
        print(f"Example: {sys.argv[0]} data.txt dist_val.txt n.txt", file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    n_file_path = sys.argv[3]

    N = read_n_from_file(n_file_path)

    with open(input_path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        print("Error: empty input file", file=sys.stderr)
        sys.exit(1)

    # 旧形式互換: 1行目が整数のみ、かつ2行目が3列なら1行目をスキップ
    data_lines = lines
    if len(lines) >= 2 and is_int_single_token(lines[0]):
        if len(lines[1].split()) == 3:
            data_lines = lines[1:]

    # 3列チェックしつつ保持（値は文字列のまま：情報落ちなし）
    rows = []
    for line_no, ln in enumerate(data_lines, start=1):
        cols = ln.split()
        if len(cols) != 3:
            print(f"Error: expected 3 columns at data line {line_no}, got {len(cols)}: '{ln}'",
                  file=sys.stderr)
            sys.exit(1)
        rows.append(cols)

    M = len(rows)
    if M > N:
        print(f"Error: input has {M} data rows but N={N} (too many rows)", file=sys.stderr)
        sys.exit(1)

    with open(output_path, "w", encoding="utf-8") as out:
        out.write("# dist_val\n")
        out.write(f"{N} 3\n")

        for cols in rows:
            out.write(" ".join(cols) + "\n")

        zeros = "0.0 0.0 0.0"
        for _ in range(N - M):
            out.write(zeros + "\n")

if __name__ == "__main__":
    main()