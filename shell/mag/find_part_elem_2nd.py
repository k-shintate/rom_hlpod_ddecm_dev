import argparse


def make_key(nodes, ignore_order=True):
    if ignore_order:
        return tuple(sorted(nodes))
    return tuple(nodes)


def read_graph2(filename, ignore_order=True):
    """
    グラフ2を読む。
    4節点要素を key にして、元の節点並びを保存する。
    """
    graph2_map = {}

    with open(filename, "r") as f:
        n_elem = int(f.readline().strip())

        for line in f:
            cols = list(map(int, line.split()))
            if not cols:
                continue

            elem_id = cols[0]
            elem_type = cols[1]

            if elem_type == 4:
                nodes = cols[2:6]
                key = make_key(nodes, ignore_order)

                graph2_map[key] = nodes

    return graph2_map


def convert_graph1(graph1_file, graph2_map, ignore_order=True):
    """
    グラフ1を読み、指定形式の出力行を作る。

    出力形式:
      {graph1の要素数}
      {graph1_id} 4 n1 n2 n3 n4

    マッチしない場合:
      {graph1_id} 4 0 0 0 0
    """
    output_lines = []

    with open(graph1_file, "r") as f:
        n_elem = int(f.readline().strip())
        output_lines.append(str(n_elem))

        for line in f:
            cols = list(map(int, line.split()))
            if not cols:
                continue

            graph1_id = cols[0]
            elem_type = cols[1]

            matched_nodes = None

            if elem_type == 10:
                # 10節点要素の先頭4節点を使って照合
                corner_nodes = cols[2:6]
                key = make_key(corner_nodes, ignore_order)

                if key in graph2_map:
                    matched_nodes = graph2_map[key]

            if matched_nodes is None:
                output_lines.append(f"{graph1_id} 4 -1 -1 -1 -1")
            else:
                n1, n2, n3, n4 = matched_nodes
                output_lines.append(f"{graph1_id} 4 {n1} {n2} {n3} {n4}")

    return output_lines


def main():
    parser = argparse.ArgumentParser(
        description="グラフ1のIDを使って、グラフ2と一致する4節点要素を出力する"
    )

    parser.add_argument(
        "graph1",
        help="要素コネクティビティグラフ1のファイル"
    )

    parser.add_argument(
        "graph2",
        help="要素コネクティビティグラフ2のファイル"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="出力ファイル"
    )

    parser.add_argument(
        "--keep-order",
        action="store_true",
        help="節点順も含めて完全一致で比較する"
    )

    args = parser.parse_args()

    ignore_order = not args.keep_order

    graph2_map = read_graph2(args.graph2, ignore_order)
    output_lines = convert_graph1(args.graph1, graph2_map, ignore_order)

    with open(args.output, "w") as f:
        f.write("\n".join(output_lines))
        f.write("\n")


if __name__ == "__main__":
    main()