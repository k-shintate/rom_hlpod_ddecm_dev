# make_local_region.py

import argparse


def read_connectivity(connectivity_file):
    with open(connectivity_file, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("コネクティビティファイルが空です")

    header = lines[0].strip()
    elements = []

    for line_no, line in enumerate(lines[1:], start=2):
        line = line.strip()

        if not line:
            continue

        values = list(map(int, line.split()))
        num_nodes = values[0]
        global_ids = values[1:]

        if len(global_ids) != num_nodes:
            raise ValueError(
                f"{connectivity_file}:{line_no}: "
                f"節点数が一致しません "
                f"(先頭の値={num_nodes}, ID数={len(global_ids)})"
            )

        elements.append(global_ids)

    return header, elements


def read_coordinates(coord_file):
    with open(coord_file, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("座標ファイルが空です")

    num_coords = int(lines[0].strip())
    coords = []

    for line_no, line in enumerate(lines[1:], start=2):
        line = line.strip()

        if not line:
            continue

        xyz = list(map(float, line.split()))

        if len(xyz) != 3:
            raise ValueError(
                f"{coord_file}:{line_no}: 座標は x y z の3成分である必要があります"
            )

        coords.append(xyz)

    if len(coords) != num_coords:
        raise ValueError(
            f"座標数が一致しません: "
            f"ヘッダ={num_coords}, 実際の座標行数={len(coords)}"
        )

    return coords


def build_local_mapping(elements):
    global_to_local = {}
    local_to_global = {}
    next_local_id = 0

    for element in elements:
        for gid in element:
            if gid not in global_to_local:
                global_to_local[gid] = next_local_id
                local_to_global[next_local_id] = gid
                next_local_id += 1

    return global_to_local, local_to_global


def write_local_connectivity(output_file, header, elements, global_to_local):
    with open(output_file, "w") as f:
        f.write(header + "\n")

        for element in elements:
            local_ids = [global_to_local[gid] for gid in element]
            f.write(str(len(local_ids)) + " ")
            f.write(" ".join(map(str, local_ids)))
            f.write("\n")


def write_local_coordinates(output_file, coords, local_to_global, global_index_base):
    num_local_nodes = len(local_to_global)

    with open(output_file, "w") as f:
        f.write(str(num_local_nodes) + "\n")

        for local_id in range(num_local_nodes):
            gid = local_to_global[local_id]
            coord_index = gid - global_index_base

            if coord_index < 0 or coord_index >= len(coords):
                raise IndexError(
                    f"global_id={gid} に対応する座標がありません。"
                    f"coord_index={coord_index}, 座標数={len(coords)}"
                )

            x, y, z = coords[coord_index]
            f.write(f"{x} {y} {z}\n")


def write_local_to_global(output_file, local_to_global):
    with open(output_file, "w") as f:
        f.write("# local_id global_id\n")

        for local_id in sorted(local_to_global):
            f.write(f"{local_id} {local_to_global[local_id]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="FEMのグローバルIDコネクティビティと座標を領域内ローカルIDへ変換する"
    )

    parser.add_argument(
        "connectivity_file",
        help="グローバルIDで書かれたコネクティビティファイル",
    )

    parser.add_argument(
        "coord_file",
        help="グローバル座標ファイル",
    )

    parser.add_argument(
        "local_connectivity_file",
        help="ローカルID化したコネクティビティの出力ファイル",
    )

    parser.add_argument(
        "local_coord_file",
        help="領域分割後のローカル座標データの出力ファイル",
    )

    parser.add_argument(
        "local_to_global_file",
        help="local_id -> global_id 対応表の出力ファイル",
    )

    parser.add_argument(
        "--global-index-base",
        type=int,
        choices=[0, 1],
        default=0,
        help=(
            "座標ファイルにおけるグローバルIDの開始番号。"
            "0ならglobal_id=0が座標ファイル2行目、"
            "1ならglobal_id=1が座標ファイル2行目。"
            "デフォルトは0。"
        ),
    )

    args = parser.parse_args()

    header, elements = read_connectivity(args.connectivity_file)
    coords = read_coordinates(args.coord_file)

    global_to_local, local_to_global = build_local_mapping(elements)

    write_local_connectivity(
        args.local_connectivity_file,
        header,
        elements,
        global_to_local,
    )

    write_local_coordinates(
        args.local_coord_file,
        coords,
        local_to_global,
        args.global_index_base,
    )

    write_local_to_global(
        args.local_to_global_file,
        local_to_global,
    )


if __name__ == "__main__":
    main()