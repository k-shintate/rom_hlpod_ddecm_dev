#!/usr/bin/env python3
import sys


def rewrite_and_append(path: str) -> None:
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if not lines:
        raise ValueError(f"{path}: ファイルが空です")

    first_line = lines[0].strip()
    if not first_line:
        raise ValueError(f"{path}: 1行目が空です")

    parts = first_line.split()
    if not parts:
        raise ValueError(f"{path}: 1行目に数値がありません")

    num1 = int(parts[0])
    new_num1 = num1 + 1

    # 1行目を書き換え
    lines[0] = str(new_num1)

    # 末尾に追加する行を作成
    new_line = [str(new_num1), str(num1)] + [str(i) for i in range(num1 + 1)]
    lines.append(" ".join(new_line))

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main(argv):
    if len(argv) != 1:
        print(f"使い方: python {sys.argv[0]} <対象ファイル>", file=sys.stderr)
        sys.exit(1)

    rewrite_and_append(argv[0])


if __name__ == "__main__":
    main(sys.argv[1:])