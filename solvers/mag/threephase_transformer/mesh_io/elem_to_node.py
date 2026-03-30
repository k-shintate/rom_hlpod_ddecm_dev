import sys

def extract_unique_nodes(connectivity_file, output_file):
    """
    要素コネクティビティフォーマットを読み込み、全要素を構成する節点番号を出力する。

    Parameters:
        connectivity_file (str): 入力コネクティビティファイルのパス。
        output_file (str): 節点番号を保存するファイルのパス。
    """
    try:
        # ファイルを読み込む
        with open(connectivity_file, "r") as f:
            lines = f.readlines()

        # 最初の行から総要素数と要素内節点数を取得
        header = lines[0].strip().split()
        total_elements, nodes_per_element = map(int, header)

        # 残りの行を読み込み、節点番号を取得
        nodes = set()
        for line in lines[1:]:
            nodes.update(map(int, line.strip().split()))

        # ユニークな節点番号をリスト化してソート
        unique_nodes = sorted(nodes)

        # 結果を保存
        with open(output_file, "w") as f_out:
            # 総ノード数を最初に書き込む
            f_out.write(f"{len(unique_nodes)}\n")
            
            # 各ノードを改行して書き込む
            for node in unique_nodes:
                f_out.write(f"{node}\n")

        print(f"Unique nodes saved to: {output_file}")
        print(f"Total unique nodes: {len(unique_nodes)}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    # コマンドライン引数からファイル名を取得
    if len(sys.argv) != 3:
        print("Usage: python script.py <connectivity_file> <output_file>")
    else:
        connectivity_file = sys.argv[1]  # 入力ファイル
        output_file = sys.argv[2]       # 出力ファイル

        # 関数を実行
        extract_unique_nodes(connectivity_file, output_file)
