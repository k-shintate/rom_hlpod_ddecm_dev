from collections import defaultdict
import sys

def convert_element_to_node_connectivity(input_file, output_file, n_nodes=None, unique=True):
    node_conn = defaultdict(list)
    max_node_id = -1

    with open(input_file, "r") as f:
        n_elem = int(f.readline().strip())

        for line_no in range(1, n_elem + 1):
            line = f.readline()
            if not line:
                raise ValueError(f"{line_no} 行目でファイルが終了しました")

            parts = list(map(int, line.split()))

            if len(parts) < 3:
                continue

            elem_id = parts[0]
            n = parts[1]
            nodes = parts[2:]

            if len(nodes) != n:
                raise ValueError(
                    f"要素 {elem_id}: 節点数が一致しません "
                    f"expected {n}, got {len(nodes)}"
                )

            max_node_id = max(max_node_id, max(nodes))

            # 同じ要素内の節点どうしをすべて連結
            for i, ni in enumerate(nodes):
                for j, nj in enumerate(nodes):
                    if i != j:
                        node_conn[ni].append(nj)

    if n_nodes is None:
        n_nodes = max_node_id + 1

    with open(output_file, "w") as f:
        f.write(f"{n_nodes}\n")

        for node_id in range(n_nodes):
            neighbors = node_conn.get(node_id, [])

            if unique:
                neighbors = sorted(set(neighbors))

            f.write(
                f"{node_id} {len(neighbors)}"
                + ("" if not neighbors else " " + " ".join(map(str, neighbors)))
                + "\n"
            )

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    convert_element_to_node_connectivity(
        input_file,
        output_file,
        n_nodes=None,
        unique=True
    )