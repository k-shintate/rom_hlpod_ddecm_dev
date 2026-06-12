# extract_local_boundary.py

import argparse


def read_local_to_global(mapping_file):
    """
    local_to_global.txt から global_id -> local_id を作る。

    入力例:
        # local_id global_id
        0 12964
        1 12988
        2 10678
    """

    global_to_local = {}

    with open(mapping_file, "r") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            values = line.split()

            if len(values) < 2:
                raise ValueError(
                    f"{mapping_file}:{line_no}: local_id global_id の2列が必要です"
                )

            local_id = int(values[0])
            global_id = int(values[1])

            global_to_local[global_id] = local_id

    return global_to_local


def extract_local_boundary(
    local_to_global_file,
    global_boundary_file,
    local_boundary_file,
):
    global_to_local = read_local_to_global(local_to_global_file)

    with open(global_boundary_file, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("境界条件ファイルが空です")

    header = lines[0].split()

    if len(header) < 2:
        raise ValueError("境界条件ファイルの1行目は '境界条件数 成分数' を想定しています")

    num_components = header[1]

    local_boundary_lines = []

    for line_no, line in enumerate(lines[1:], start=2):
        line = line.strip()

        if not line:
            continue

        values = line.split()

        if len(values) < 3:
            raise ValueError(
                f"{global_boundary_file}:{line_no}: "
                "global_id dof value の形式を想定しています"
            )

        global_id = int(values[0])

        # この領域に含まれない節点IDは出力しない
        if global_id not in global_to_local:
            continue

        local_id = global_to_local[global_id]

        # 1列目だけ local_id に置換し、2列目以降はそのまま
        rest = values[1:]

        local_boundary_lines.append([str(local_id)] + rest)

    with open(local_boundary_file, "w") as f:
        f.write(f"{len(local_boundary_lines)} {num_components}\n")

        for values in local_boundary_lines:
            f.write(" ".join(values) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="グローバル境界条件ファイルを、領域内ローカルIDの境界条件ファイルへ変換する"
    )

    parser.add_argument(
        "local_to_global_file",
        help="local_id global_id の対応表",
    )

    parser.add_argument(
        "global_boundary_file",
        help="未領域分割のグローバル境界条件ファイル",
    )

    parser.add_argument(
        "local_boundary_file",
        help="出力するローカル境界条件ファイル",
    )

    args = parser.parse_args()

    extract_local_boundary(
        args.local_to_global_file,
        args.global_boundary_file,
        args.local_boundary_file,
    )


if __name__ == "__main__":
    main()