#!/usr/bin/env python3
import sys


def read_num1(path: str) -> int:
    with open(path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()

    if not first_line:
        raise ValueError(f"{path}: 1行目が空です")

    return int(first_line.split()[0]) - 1


def rewrite_mesh_file(input_file: str, output_file: str, num1_file: str) -> None:
    num1 = read_num1(num1_file)

    with open(input_file, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if not lines:
        raise ValueError(f"{input_file}: ファイルが空です")

    new_lines = [lines[0]]

    for lineno, line in enumerate(lines[1:], start=2):
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()
        if len(parts) < 2:
            raise ValueError(f"{input_file}:{lineno}: 列数が不足しています: {line}")

        try:
            parts[1] = str(int(parts[1]) + 1)
        except ValueError as e:
            raise ValueError(
                f"{input_file}:{lineno}: 2列目が整数ではありません: {parts[1]}"
            ) from e

        parts.append(str(num1))
        new_lines.append(" ".join(parts))

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(new_lines) + "\n")


def main(argv):
    if len(argv) != 3:
        print(
            f"使い方: python {sys.argv[0]} <入力ファイル> <出力ファイル> <num1を読むファイル>",
            file=sys.stderr,
        )
        sys.exit(1)

    rewrite_mesh_file(argv[0], argv[1], argv[2])


if __name__ == "__main__":
    main(sys.argv[1:])