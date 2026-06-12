#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
領域分割後の要素IDリストを使って、
elem.dat 形式の要素コネクティビティを抽出する。

elem.dat は以下の形式を想定する。

    要素数 1要素あたりの節点数
    node1 node2 node3 node4
    node1 node2 node3 node4
    ...

注意:
    elem.dat には要素ID列がないため、
    要素IDは行番号から決める。

    --id-base 0 の場合:
        elem.dat の1行目の要素ID = 0

    --id-base 1 の場合:
        elem.dat の1行目の要素ID = 1

使い方:
    python3 nedelecelem2nodalelem.py id_file elem.dat output_elem.dat

例:
    python3 shell/mag/nedelecelem2nodalelem.py \
      result_mag/2-8-8/parted.0/graph_nedelec_conn_all.dat.id.0 \
      result_mag/2-8-8/elem.dat \
      result_mag/2-8-8/parted.0/elem.dat.0

1始まりの要素IDとして扱う場合:
    python3 shell/mag/nedelecelem2nodalelem.py \
      result_mag/2-8-8/parted.0/graph_nedelec_conn_all.dat.id.0 \
      result_mag/2-8-8/elem.dat \
      result_mag/2-8-8/parted.0/elem.dat.0 \
      --id-base 1
"""

import argparse
import sys


def read_partition_ids(filename):
    """
    領域分割後の要素IDを読み込む。

    # で始まる行は無視する。
    各行の先頭の整数だけを要素IDとして読む。

    例:
        #id
        10861 1
        10454
        10473
    """

    partition_ids = set()

    with open(filename, "r") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            values = line.split()

            try:
                elem_id = int(values[0])
            except ValueError:
                raise ValueError(
                    f"{filename}:{line_no}: "
                    f"要素IDとして読めない行があります: {line}"
                )

            partition_ids.add(elem_id)

    if not partition_ids:
        raise ValueError(
            f"{filename}: 有効な要素IDが1つも読み込めませんでした。"
        )

    return partition_ids


def read_elem_dat(filename):
    """
    elem.dat を読み込む。

    形式:
        num_elem nodes_per_elem
        node1 node2 node3 node4
        node1 node2 node3 node4
        ...

    返り値:
        num_elem_header
        nodes_per_elem
        elements
    """

    elements = []

    with open(filename, "r") as f:
        header_line = f.readline().strip()

        if not header_line:
            raise ValueError(f"{filename}: ファイルが空です。")

        header = header_line.split()

        if len(header) < 2:
            raise ValueError(
                f"{filename}: ヘッダが不正です: {header_line}"
            )

        try:
            num_elem_header = int(header[0])
            nodes_per_elem = int(header[1])
        except ValueError:
            raise ValueError(
                f"{filename}: ヘッダを整数として読めません: {header_line}"
            )

        for line_no, line in enumerate(f, start=2):
            line = line.strip()

            if not line:
                continue

            try:
                values = list(map(int, line.split()))
            except ValueError:
                raise ValueError(
                    f"{filename}:{line_no}: "
                    f"整数として読めない値があります: {line}"
                )

            if len(values) != nodes_per_elem:
                raise ValueError(
                    f"{filename}:{line_no}: "
                    f"節点数が不正です。"
                    f"期待値={nodes_per_elem}, 実際={len(values)}, 行={line}"
                )

            elements.append(values)

    if num_elem_header != len(elements):
        raise ValueError(
            f"{filename}: ヘッダの要素数と実際の行数が一致しません。"
            f"ヘッダ={num_elem_header}, 実際={len(elements)}"
        )

    return num_elem_header, nodes_per_elem, elements


def extract_partitioned_elements(elements, partition_ids, id_base):
    """
    要素IDに対応する行だけを抽出する。

    見つからない ID がある場合はエラーにする。
    抽出結果が 0 件の場合もエラーにする。
    """

    extracted = []
    found_ids = set()

    num_elements = len(elements)

    for elem_id in sorted(partition_ids):
        index = elem_id - id_base

        if index < 0 or index >= num_elements:
            continue

        extracted.append(elements[index])
        found_ids.add(elem_id)

    missing_ids = partition_ids - found_ids

    if missing_ids:
        preview = sorted(missing_ids)[:20]
        raise ValueError(
            "partition_id_file に elem.dat 内で見つからない要素IDがあります。\n"
            f"見つからないID数: {len(missing_ids)}\n"
            f"例: {preview}\n"
            f"elem.dat の有効ID範囲: "
            f"{id_base} 〜 {id_base + num_elements - 1}\n"
            f"id_base: {id_base}"
        )

    if not extracted:
        raise ValueError(
            "抽出結果が 0 件です。partition_id_file と elem.dat の対応、"
            "または --id-base を確認してください。"
        )

    return extracted


def write_elem_dat(filename, elements, nodes_per_elem):
    """
    分割後の elem.dat を書き出す。
    """

    with open(filename, "w") as f:
        f.write(f"{len(elements)}\n")

        for elem_nodes in elements:
            f.write(f"{nodes_per_elem} ")
            f.write(" ".join(map(str, elem_nodes)) + "\n")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="領域分割後の要素IDを使って elem.dat を抽出します。"
    )

    parser.add_argument(
        "partition_id_file",
        help="領域分割後の要素IDファイル"
    )

    parser.add_argument(
        "input_elem_file",
        help="入力 elem.dat"
    )

    parser.add_argument(
        "output_elem_file",
        help="出力 elem.dat"
    )

    parser.add_argument(
        "--id-base",
        type=int,
        choices=[0, 1],
        default=0,
        help=(
            "elem.dat の最初の要素行を ID 0 とみなすか ID 1 とみなすか。"
            "default: 0"
        )
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    try:
        partition_ids = read_partition_ids(args.partition_id_file)

        num_elem_header, nodes_per_elem, elements = read_elem_dat(
            args.input_elem_file
        )

        extracted_elements = extract_partitioned_elements(
            elements,
            partition_ids,
            args.id_base
        )

        write_elem_dat(
            args.output_elem_file,
            extracted_elements,
            nodes_per_elem
        )

    except Exception as e:
        print("ERROR:", e, file=sys.stderr)
        sys.exit(1)

    print("変換完了")
    print(f"入力 elem.dat ヘッダ要素数: {num_elem_header}")
    print(f"入力 elem.dat 実要素数: {len(elements)}")
    print(f"分割後 ID 数: {len(partition_ids)}")
    print(f"抽出された要素数: {len(extracted_elements)}")
    print(f"id_base: {args.id_base}")
    print(f"出力ファイル: {args.output_elem_file}")


if __name__ == "__main__":
    main()