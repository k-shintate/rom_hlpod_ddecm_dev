import sys

def elem_to_bc(connectivity_file, output_file, block_length, values):
    """
    要素コネクティビティファイルから、全要素を構成するユニークな節点を抽出し、
    それらの節点に境界条件を付与したファイルを出力する。
    """

    try:
        with open(connectivity_file, "r") as f:
            lines = f.readlines()

        if not lines:
            raise ValueError("Input file is empty.")

        # ヘッダー読み込み
        header = lines[0].strip().split()
        total_elements, nodes_per_element = map(int, header)

        nodes = set()

        # 各要素の節点番号を収集
        for i, line in enumerate(lines[1:], start=1):
            parts = line.strip().split()

            if not parts:
                continue

            node_ids = list(map(int, parts))

            if len(node_ids) != nodes_per_element:
                raise ValueError(
                    f"Line {i + 1}: expected {nodes_per_element} nodes, "
                    f"but got {len(node_ids)}."
                )

            nodes.update(node_ids)

        unique_nodes = sorted(nodes)

        # 各節点につき block_length 個のBCを書くので、
        # ヘッダーにはBCデータの総行数を書く
        boundary_condition_dofs = len(unique_nodes) * block_length

        with open(output_file, "w") as f_out:
            f_out.write(f"{boundary_condition_dofs} {block_length}\n")

            for node_id in unique_nodes:
                for block_index in range(block_length):
                    value = values[block_index]
                    f_out.write(
                        f"{node_id} {block_index} {value:.6f}\n"
                    )

        print(f"Boundary condition file saved to: {output_file}")
        print(f"Total unique nodes: {len(unique_nodes)}")
        print(f"Block length: {block_length}")
        print(f"Total BCs: {boundary_condition_dofs}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage:")
        print(
            "  python elem_to_bc.py "
            "<connectivity_file> <output_file> "
            "<block_length> <values>"
        )
        print()
        print("Example:")
        print(
            "  python elem_to_bc.py "
            "elem.dat bc.dat 3 0.0 1.0 1.0"
        )
        sys.exit(1)

    connectivity_file = sys.argv[1]
    output_file = sys.argv[2]
    block_length = int(sys.argv[3])

    try:
        values = list(map(float, sys.argv[4:]))

        if len(values) != block_length:
            raise ValueError(
                f"Number of values ({len(values)}) does not match "
                f"block length ({block_length})."
            )

    except ValueError as e:
        print(f"Value Error: {e}")
        sys.exit(1)

    elem_to_bc(
        connectivity_file,
        output_file,
        block_length,
        values,
    )
