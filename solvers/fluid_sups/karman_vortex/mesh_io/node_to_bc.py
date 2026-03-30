import sys

def process_nodes(input_file, output_file, block_length, values):
    """
    境界条件が付与された節点自由度数、ブロック長さ、値を指定してフォーマットされたファイルを出力

    Parameters:
        input_file (str): 読み込むファイルのパス
        output_file (str): 出力ファイルのパス
        block_length (int): ブロック長さ (b)
        values (list of float): ブロック長さに対応する値のリスト
    """
    try:
        # ファイルを読み込む
        with open(input_file, "r") as f:
            lines = f.readlines()

        # 最初の行はノード数
        node_count = int(lines[0].strip())
        # 残りの行は節点ID
        node_ids = [int(line.strip()) for line in lines[1:]]

        # 境界条件ファイルに記述した行数（節点数）
        boundary_condition_dofs = len(node_ids)

        # 出力ファイルにフォーマットを保存
        with open(output_file, "w") as f_out:
            # ファイルヘッダーを記述
            f_out.write(f"{boundary_condition_dofs * block_length} {block_length}\n")

            # 各節点IDに対してデータを記述
            for node_id in node_ids:
                for block_index in range(block_length):
                    # 値を浮動小数点数としてフォーマット
                    value = values[block_index]
                    f_out.write(f"{node_id} {block_index} {value:.6f}\n")

        print(f"Processed file saved to: {output_file}")
        print(f"Total nodes processed: {boundary_condition_dofs}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    # コマンドライン引数からファイル名を取得
    if len(sys.argv) < 5:
        print("Usage: python script.py <input_file> <output_file> <block_length> <values>")
        print("Example: python script.py input.dat output.dat 3 0.0 1.0 1.0")
    else:
        input_file = sys.argv[1]      # 入力ファイル
        output_file = sys.argv[2]     # 出力ファイル
        block_length = int(sys.argv[3])  # ブロック長さ

        # 値のリストを浮動小数点数として取得
        try:
            values = list(map(float, sys.argv[4:]))
            if len(values) != block_length:
                raise ValueError(f"Number of values ({len(values)}) does not match block length ({block_length}).")
        except ValueError as ve:
            print(f"Value Error: {ve}")
            sys.exit(1)

        # 関数を実行
        process_nodes(input_file, output_file, block_length, values)
