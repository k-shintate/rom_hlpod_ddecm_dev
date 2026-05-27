from collections import defaultdict
import sys

def convert_connectivity(input_file, output_file, n_nodes=None, unique=True):
    conn = defaultdict(list)

    with open(input_file, "r") as f:
        n_elem = int(f.readline().strip())

        for _ in range(n_elem):
            parts = list(map(int, f.readline().split()))
            if not parts:
                continue

            n = parts[0]
            nodes = parts[1:]

            if len(nodes) != n:
                raise ValueError(f"節点数が一致しません: expected {n}, got {len(nodes)}")

            for i, ni in enumerate(nodes):
                for j, nj in enumerate(nodes):
                    if i != j:
                        conn[ni].append(nj)

    if n_nodes is None:
        max_node = max(conn.keys())
        n_nodes = max_node + 1

    with open(output_file, "w") as f:
        f.write(f"{n_nodes}\n")

        for node_id in range(n_nodes):
            neighbors = conn.get(node_id, [])

            if unique:
                neighbors = sorted(set(neighbors))
            else:
                neighbors = neighbors

            f.write(
                f"{node_id} {len(neighbors)} "
                + " ".join(map(str, neighbors))
                + "\n"
            )

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    convert_connectivity(
        input_file,
        output_file,
        n_nodes=None,
        unique=True
    )
