import sys

def filter_rows_by_column(input_file, output_file, column_index, exclude_values, header_divisor):
    """
    指定された列の値に基づいて行をフィルタリングし、ヘッダーを処理する

    Parameters:
    - input_file (str): 入力ファイル名
    - output_file (str): 出力ファイル名
    - column_index (int): フィルタリングに使用する列のインデックス（0から始まる）
    - exclude_values (list of int): 除外する値のリスト
    - header_divisor (int): ヘッダー行の最初の値をこの値で割る
    """
    with open(input_file, "r") as infile:
        lines = infile.readlines()

    # ヘッダー行を処理
    header = lines[0].strip()
    try:
        header_parts = header.split()
        original_value = int(header_parts[0])
        processed_value = original_value / header_divisor
        header_parts[0] = str(int(processed_value))
        new_header = " ".join(header_parts) + "\n"
    except (IndexError, ValueError):
        print(f"Error processing header: {header}")
        sys.exit(1)

    # データ行をフィルタリング
    filtered_lines = [new_header]
    for line in lines[1:]:
        # 指定された列の値を取得
        try:
            value = int(line.split()[column_index])
            if value not in exclude_values:
                filtered_lines.append(line)
        except (IndexError, ValueError):
            print(f"Skipping line due to parsing error: {line.strip()}")

    # フィルタリング後のデータを保存
    with open(output_file, "w") as outfile:
        outfile.writelines(filtered_lines)

    print(f"Filtered data saved to {output_file}")

# コマンドライン引数から入力ファイル、出力ファイル、列インデックス、除外する値を取得
if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python script.py <input_file> <output_file> <column_index> <header_divisor> <exclude_values>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    column_index = int(sys.argv[3])  # フィルタリングに使用する列のインデックス（0始まり）
    header_divisor = int(sys.argv[4])  # ヘッダー行を処理するための値
    exclude_values = list(map(int, sys.argv[5:]))  # 5番目以降の引数をリストに変換

    filter_rows_by_column(input_file, output_file, column_index, exclude_values, header_divisor)
