import argparse
import numpy as np


def read_non_comment_lines(path):
    lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue
            if line.startswith("#"):
                continue

            lines.append(line)

    return lines


def read_elem_file(path):
    """
    elem ファイルを読む。

    例:
        #elem_id
        78499 10
        778 4972 5049 ...
        ...
    """
    lines = read_non_comment_lines(path)

    n_elem, n_node_per_elem = map(int, lines[0].split())

    elem = np.array(
        [[int(v) for v in line.split()] for line in lines[1:]],
        dtype=int
    )

    if elem.shape != (n_elem, n_node_per_elem):
        raise ValueError(
            f"elem shape error: expected {(n_elem, n_node_per_elem)}, "
            f"got {elem.shape}"
        )

    return elem


def read_element_coord_file(path, n_node_per_elem):
    """
    要素ごとの node_coordinate ファイルを読む。

    例:
        # node_coordinate
        78499 30
        x1 y1 z1 x2 y2 z2 ... x10 y10 z10
        ...
    """
    lines = read_non_comment_lines(path)

    n_elem, n_coord = map(int, lines[0].split())
    expected_coord = n_node_per_elem * 3

    if n_coord != expected_coord:
        raise ValueError(
            f"coordinate column error: expected {expected_coord}, got {n_coord}"
        )

    coord_flat = np.array(
        [[float(v) for v in line.split()] for line in lines[1:]],
        dtype=float
    )

    if coord_flat.shape != (n_elem, n_coord):
        raise ValueError(
            f"coordinate shape error: expected {(n_elem, n_coord)}, "
            f"got {coord_flat.shape}"
        )

    coord = coord_flat.reshape(n_elem, n_node_per_elem, 3)

    return coord


def read_local_to_global_file(path):
    """
    local_id -> global_id の対応表を読む。

    形式1:
        # local_to_global
        14148
        778
        4972
        5049
        ...

    この場合:
        local_id 1 -> global_id 778
        local_id 2 -> global_id 4972

    形式2:
        # local_to_global
        14148
        1 778
        2 4972
        3 5049
        ...

    この場合:
        明示的に local_id global_id として読む。
    """
    lines = read_non_comment_lines(path)

    n_local = int(lines[0].split()[0])
    data_lines = lines[1:]

    if len(data_lines) != n_local:
        raise ValueError(
            f"local_to_global shape error: expected {n_local} lines, "
            f"got {len(data_lines)}"
        )

    local_to_global = {}

    first_values = data_lines[0].split()

    if len(first_values) == 1:
        for i, line in enumerate(data_lines):
            local_id = i + 1
            global_id = int(line.split()[0])
            local_to_global[local_id] = global_id

    elif len(first_values) == 2:
        for line in data_lines:
            local_id, global_id = map(int, line.split())
            local_to_global[local_id] = global_id

    else:
        raise ValueError(
            "local_to_global format error: each line must have 1 or 2 integers"
        )

    return local_to_global


def get_global_id(node_id_in_elem, local_to_global):
    """
    elem に入っている番号を global_id に変換する。

    local_to_global が None の場合:
        elem の番号をそのまま global_id とみなす。

    local_to_global がある場合:
        elem の番号を local_id とみなし、global_id に変換する。
    """
    if local_to_global is None:
        return node_id_in_elem

    if node_id_in_elem not in local_to_global:
        return None

    return local_to_global[node_id_in_elem]


def build_node_coordinate_map(elem, coord, local_to_global=None, tolerance=1.0e-10):
    """
    global_id -> xyz の辞書を作る。

    elem[e, j] に対応する座標は coord[e, j, :]。

    同じ global_id が複数回出てきた場合、
    座標が一致しているか確認する。

    0 0 0 は通常の座標として扱う。
    """
    node_coord = {}
    original_ids_for_global = {}

    conflict_records = []
    missing_mapping_records = []

    n_elem, n_node_per_elem = elem.shape

    for e in range(n_elem):
        for j in range(n_node_per_elem):
            node_id_in_elem = int(elem[e, j])
            xyz = coord[e, j, :]

            global_id = get_global_id(node_id_in_elem, local_to_global)

            if global_id is None:
                missing_mapping_records.append(
                    {
                        "element_index": e + 1,
                        "local_index": j + 1,
                        "node_id_in_elem": node_id_in_elem
                    }
                )
                continue

            if global_id in node_coord:
                old_xyz = node_coord[global_id]

                if not np.allclose(old_xyz, xyz, atol=tolerance, rtol=0.0):
                    conflict_records.append(
                        {
                            "global_id": global_id,
                            "node_id_in_elem": node_id_in_elem,
                            "old_coord": old_xyz.copy(),
                            "new_coord": xyz.copy(),
                            "element_index": e + 1,
                            "local_index": j + 1
                        }
                    )
            else:
                node_coord[global_id] = xyz.copy()

            if global_id not in original_ids_for_global:
                original_ids_for_global[global_id] = set()

            original_ids_for_global[global_id].add(node_id_in_elem)

    return (
        node_coord,
        original_ids_for_global,
        conflict_records,
        missing_mapping_records
    )


def build_global_to_new_local_id(node_coord):
    """
    global_id 昇順に並べたときの新しい local_id を作る。

    sorted_node_coordinate の 1行目 -> local_id 1
    sorted_node_coordinate の 2行目 -> local_id 2
    ...
    """
    sorted_global_ids = sorted(node_coord.keys())

    global_to_new_local = {
        global_id: i
        for i, global_id in enumerate(sorted_global_ids)
    }

    return global_to_new_local


def convert_elem_to_new_local_elem(elem, local_to_global, global_to_new_local):
    """
    elem を、ソート後座標ファイルに対応する new local_id の elem に変換する。

    local_to_global が None の場合:
        elem に入っている番号を global_id とみなす。

    local_to_global がある場合:
        elem に入っている番号を local_id とみなし、
        local_id -> global_id -> new_local_id に変換する。
    """
    n_elem, n_node_per_elem = elem.shape

    local_elem = np.empty_like(elem)
    missing_records = []

    for e in range(n_elem):
        for j in range(n_node_per_elem):
            node_id_in_elem = int(elem[e, j])

            global_id = get_global_id(node_id_in_elem, local_to_global)

            if global_id is None:
                local_elem[e, j] = -1
                missing_records.append(
                    {
                        "element_index": e + 1,
                        "local_index": j + 1,
                        "node_id_in_elem": node_id_in_elem,
                        "global_id": None
                    }
                )
                continue

            if global_id not in global_to_new_local:
                local_elem[e, j] = -1
                missing_records.append(
                    {
                        "element_index": e + 1,
                        "local_index": j + 1,
                        "node_id_in_elem": node_id_in_elem,
                        "global_id": global_id
                    }
                )
                continue

            local_elem[e, j] = global_to_new_local[global_id]

    return local_elem, missing_records


def write_sorted_coordinates(path, node_coord, with_global_id=False):
    """
    global_id 昇順に node coordinate を出力する。

    with_global_id=False:
        14148
        x y z
        x y z
        ...

    with_global_id=True:
        14148
        global_id x y z
        global_id x y z
        ...
    """
    sorted_global_ids = sorted(node_coord.keys())

    with open(path, "w") as f:
        f.write(f"{len(sorted_global_ids)}\n")

        for global_id in sorted_global_ids:
            x, y, z = node_coord[global_id]

            if with_global_id:
                f.write(f"{global_id:d} {x:.15g} {y:.15g} {z:.15g}\n")
            else:
                f.write(f"{x:.15g} {y:.15g} {z:.15g}\n")


def write_local_elem_file(path, local_elem):
    """
    new local_id に変換した elem を出力する。

    ここで出力される local_id は、
    sorted_node_coordinate の行番号に対応する。
    """
    n_elem, n_node_per_elem = local_elem.shape

    with open(path, "w") as f:
        f.write("#elem_id\n")
        f.write(f"{n_elem} {n_node_per_elem}\n")

        for e in range(n_elem):
            f.write(" ".join(str(v) for v in local_elem[e]) + "\n")


def write_global_local_map(path, node_coord, original_ids_for_global, global_to_new_local):
    """
    global_id, new_local_id, 元の elem に出てきた番号, 座標を出力する。

    確認用。
    """
    sorted_global_ids = sorted(node_coord.keys())

    with open(path, "w") as f:
        f.write("# new_local_id global_id original_id_list x y z\n")

        for global_id in sorted_global_ids:
            new_local_id = global_to_new_local[global_id]
            original_ids = sorted(original_ids_for_global.get(global_id, []))
            original_ids_text = ",".join(str(v) for v in original_ids)

            x, y, z = node_coord[global_id]

            f.write(
                f"{new_local_id:d} "
                f"{global_id:d} "
                f"{original_ids_text} "
                f"{x:.15g} {y:.15g} {z:.15g}\n"
            )


def write_check_report(
    path,
    conflict_records,
    missing_mapping_records,
    missing_local_elem_records
):
    with open(path, "w") as f:
        f.write("# check report\n\n")

        f.write(f"coordinate conflict count: {len(conflict_records)}\n")
        f.write(f"missing local-to-global mapping count: {len(missing_mapping_records)}\n")
        f.write(f"missing local elem conversion count: {len(missing_local_elem_records)}\n")
        f.write("\n")

        if missing_mapping_records:
            f.write("# missing local-to-global mappings\n")
            f.write("element_index local_index node_id_in_elem\n")

            for r in missing_mapping_records:
                f.write(
                    f'{r["element_index"]} '
                    f'{r["local_index"]} '
                    f'{r["node_id_in_elem"]}\n'
                )

            f.write("\n")

        if missing_local_elem_records:
            f.write("# missing local elem conversions\n")
            f.write("element_index local_index node_id_in_elem global_id\n")

            for r in missing_local_elem_records:
                f.write(
                    f'{r["element_index"]} '
                    f'{r["local_index"]} '
                    f'{r["node_id_in_elem"]} '
                    f'{r["global_id"]}\n'
                )

            f.write("\n")

        if conflict_records:
            f.write("# coordinate conflicts\n")
            f.write(
                "global_id node_id_in_elem element_index local_index "
                "old_x old_y old_z new_x new_y new_z\n"
            )

            for r in conflict_records:
                ox, oy, oz = r["old_coord"]
                nx, ny, nz = r["new_coord"]

                f.write(
                    f'{r["global_id"]} '
                    f'{r["node_id_in_elem"]} '
                    f'{r["element_index"]} '
                    f'{r["local_index"]} '
                    f'{ox:.15g} {oy:.15g} {oz:.15g} '
                    f'{nx:.15g} {ny:.15g} {nz:.15g}\n'
                )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create sorted node coordinate table and local-id element connectivity "
            "from domain-decomposed FEM connectivity."
        )
    )

    parser.add_argument(
        "-e", "--elem",
        required=True,
        help="input element connectivity file"
    )

    parser.add_argument(
        "-c", "--coord",
        required=True,
        help="input element-wise node_coordinate file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="output sorted node coordinate file"
    )

    parser.add_argument(
        "--output-local-elem",
        required=True,
        help="output elem file using new local IDs corresponding to sorted coordinates"
    )

    parser.add_argument(
        "-m", "--local-to-global",
        default=None,
        help="optional local_id to global_id mapping file"
    )

    parser.add_argument(
        "--map-output",
        default="global_local_map.txt",
        help="output new_local_id-global_id correspondence file"
    )

    parser.add_argument(
        "--report",
        default="check_report.txt",
        help="output check report file"
    )

    parser.add_argument(
        "--with-global-id",
        action="store_true",
        help="output global_id x y z instead of x y z only"
    )

    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-10,
        help="tolerance for checking duplicated node coordinates"
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="stop if conflicts or missing mappings are found"
    )

    args = parser.parse_args()

    elem = read_elem_file(args.elem)

    coord = read_element_coord_file(
        args.coord,
        n_node_per_elem=elem.shape[1]
    )

    if elem.shape[0] != coord.shape[0]:
        raise ValueError(
            f"number of elements mismatch: "
            f"elem={elem.shape[0]}, coord={coord.shape[0]}"
        )

    if args.local_to_global is not None:
        local_to_global = read_local_to_global_file(args.local_to_global)
    else:
        local_to_global = None

    (
        node_coord,
        original_ids_for_global,
        conflict_records,
        missing_mapping_records
    ) = build_node_coordinate_map(
        elem,
        coord,
        local_to_global=local_to_global,
        tolerance=args.tolerance
    )

    global_to_new_local = build_global_to_new_local_id(node_coord)

    local_elem, missing_local_elem_records = convert_elem_to_new_local_elem(
        elem,
        local_to_global,
        global_to_new_local
    )

    write_check_report(
        args.report,
        conflict_records,
        missing_mapping_records,
        missing_local_elem_records
    )

    if args.strict and (
        conflict_records or
        missing_mapping_records or
        missing_local_elem_records
    ):
        raise ValueError(
            f"check failed: "
            f"conflicts={len(conflict_records)}, "
            f"missing_mappings={len(missing_mapping_records)}, "
            f"missing_local_elem={len(missing_local_elem_records)}. "
            f"See {args.report}"
        )

    write_sorted_coordinates(
        args.output,
        node_coord,
        with_global_id=args.with_global_id
    )

    write_local_elem_file(
        args.output_local_elem,
        local_elem
    )

    write_global_local_map(
        args.map_output,
        node_coord,
        original_ids_for_global,
        global_to_new_local
    )


if __name__ == "__main__":
    main()