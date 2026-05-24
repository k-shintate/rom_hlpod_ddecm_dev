import numpy as np
import argparse


def read_connectivity(path):
    """
    Connectivity file format:

        num_elements nodes_per_element
        n1 n2 n3 n4
        n1 n2 n3 n4
        ...

    Example:

        49987 4
        15568 51649 15552 15566
        51651 15565 15563 15535
        ...
    """
    with open(path, "r") as f:
        nelem, nnode_per_elem = map(int, f.readline().split())
        conn = np.loadtxt(f, dtype=np.int64)

    conn = conn.reshape(nelem, nnode_per_elem)

    return conn


def write_connectivity(path, conn):
    """
    Write local connectivity file.
    """
    nelem, nnode_per_elem = conn.shape

    with open(path, "w") as f:
        f.write(f"{nelem} {nnode_per_elem}\n")

        for elem in conn:
            f.write(" ".join(map(str, elem)) + "\n")


def write_mapping(path, global_nodes, local_start):
    """
    Write node mapping file.

    Output format:

        number_of_local_nodes
        local_node_id global_node_id
        local_node_id global_node_id
        ...

    Example:

        12345
        1 1505
        2 1511
        3 15523
        ...
    """
    num_local_nodes = len(global_nodes)

    with open(path, "w") as f:
        # 1行目: 変換後のローカル節点数
        f.write(f"{num_local_nodes}\n")

        # 2行目以降: local_node_id global_node_id
        for i, gid in enumerate(global_nodes):
            local_id = i + local_start
            f.write(f"{local_id} {gid}\n")


def renumber_connectivity(conn_global, local_start):
    """
    Convert global node IDs in connectivity to local node IDs.
    """

    # サブメッシュに出現するグローバル節点番号を取得
    # np.unique は昇順に並べる
    global_nodes = np.unique(conn_global)

    # global_node_id -> local_node_id の辞書を作る
    global_to_local = {
        int(gid): i + local_start
        for i, gid in enumerate(global_nodes)
    }

    # コネクティビティをローカル節点番号へ変換
    conn_local = np.vectorize(global_to_local.get)(conn_global)

    return conn_local.astype(np.int64), global_nodes


def main():
    parser = argparse.ArgumentParser(
        description="Renumber submesh connectivity from global node IDs to local node IDs."
    )

    parser.add_argument(
        "input_conn",
        help="input connectivity file with global node IDs",
    )

    parser.add_argument(
        "output_conn",
        help="output connectivity file with local node IDs",
    )

    parser.add_argument(
        "--mapping",
        default="node_mapping.txt",
        help="output mapping file: first line is number of local nodes",
    )

    parser.add_argument(
        "--local-start",
        type=int,
        default=1,
        choices=[0, 1],
        help="local node index start: 0 or 1",
    )

    args = parser.parse_args()

    # 入力コネクティビティを読む
    conn_global = read_connectivity(args.input_conn)

    # グローバル節点番号からローカル節点番号へ変換
    conn_local, global_nodes = renumber_connectivity(
        conn_global,
        args.local_start,
    )

    # ローカル番号化したコネクティビティを出力
    write_connectivity(args.output_conn, conn_local)

    # 対応表を出力
    write_mapping(args.mapping, global_nodes, args.local_start)

    print("Done.")
    print(f"number of elements    = {conn_local.shape[0]}")
    print(f"nodes per element     = {conn_local.shape[1]}")
    print(f"number of local nodes = {len(global_nodes)}")
    print(f"local connectivity    = {args.output_conn}")
    print(f"node mapping          = {args.mapping}")


if __name__ == "__main__":
    main()