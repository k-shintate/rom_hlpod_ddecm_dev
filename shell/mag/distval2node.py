#!/usr/bin/env python3
# remove_zero_block_anywhere.py
import sys

def to_float(tok: str) -> float:
    # Fortran の D 指数にも対応
    t = tok.strip().replace("D", "E").replace("d", "e")
    return float(t)

def is_all_zero(tokens) -> bool:
    try:
        return all(to_float(x) == -10000.0 for x in tokens)
    except ValueError:
        return False

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.txt output.txt", file=sys.stderr)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, "r", encoding="utf-8") as f:
        raw_lines = [ln.rstrip("\n") for ln in f]

    # ヘッダ(N dim)を探す：最初の「空行でも#でもない行」
    pre_lines = []
    header_line = None
    header_idx = None

    for idx, ln in enumerate(raw_lines):
        s = ln.strip()
        if not s:
            pre_lines.append(ln)
            continue
        if s.startswith("#"):
            pre_lines.append(ln)
            continue
        header_line = ln
        header_idx = idx
        break

    if header_line is None:
        print("Error: header line 'N dim' not found.", file=sys.stderr)
        sys.exit(1)

    h = header_line.split()
    if len(h) != 2:
        print("Error: header must be: 'N dim' (e.g., '23597 3')", file=sys.stderr)
        sys.exit(1)

    try:
        _N_in = int(h[0])
        dim = int(h[1])
    except ValueError:
        print("Error: header values must be integers: 'N dim'", file=sys.stderr)
        sys.exit(1)

    if dim <= 0:
        print("Error: dim must be > 0", file=sys.stderr)
        sys.exit(1)

    kept_data_lines = []
    post_header_nondata = []  # ヘッダ以降のコメント/空行（必要なら残す）

    # ヘッダ以降を処理：データ行は dim 列、ゼロ行は削除
    for ln in raw_lines[header_idx + 1:]:
        s = ln.strip()
        if not s:
            post_header_nondata.append(ln)
            continue
        if s.startswith("#"):
            post_header_nondata.append(ln)
            continue

        toks = s.split()
        if len(toks) != dim:
            print(f"Error: data line has {len(toks)} cols, expected {dim}: '{ln}'",
                  file=sys.stderr)
            sys.exit(1)

        if is_all_zero(toks):
            # ここで削除
            continue

        # 残す行は“元の行”をそのまま保持
        kept_data_lines.append(ln)

    N_out = len(kept_data_lines)

    with open(out_path, "w", encoding="utf-8") as out:
        # ヘッダは更新して出力
        out.write(f"{N_out}\n")

        # ヘッダ以降のコメント/空行を残したいなら出す（不要ならこのブロックを消してOK）
        for ln in post_header_nondata:
            out.write(ln + "\n")

        # データ出力
        for ln in kept_data_lines:
            out.write(ln + "\n")

if __name__ == "__main__":
    main()
