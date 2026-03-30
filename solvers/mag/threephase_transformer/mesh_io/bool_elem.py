#!/usr/bin/env python3
import sys
from pathlib import Path

def read_cells(path):
    """ファイルを読み込み、ヘッダーと要素(タプルのリスト)を返す"""
    try:
        content = Path(path).read_text(encoding="utf-8")
    except FileNotFoundError:
        sys.stderr.write(f"Error: File not found '{path}'\n")
        sys.exit(1)

    lines = [ln.strip() for ln in content.splitlines() if ln.strip()]
    if not lines:
        return [], []
    
    # 1行目はヘッダーと仮定
    header = [int(x) for x in lines[0].split()]
    # 2行目以降が要素定義
    cells = [tuple(int(x) for x in ln.split()) for ln in lines[1:]]
    return header, cells

def main(args):
    # 最低でも [入力1, 統合, 出力] の3つの引数が必要
    if len(args) < 3:
        print("Usage: python3 bool_elem.py <file1> [<file2> ...] <merged_file> <out_file>")
        sys.exit(1)

    # 引数のパース
    # 最後尾が出力ファイル
    out_file = args[-1]
    # 後ろから2番目が統合ファイル(チェック対象)
    merged_file = args[-2]
    # それ以外(先頭～後ろから3番目まで)が入力ファイル群
    input_files = args[:-2]

    # 入力ファイルをすべて読み込み、setに変換してリストに格納
    # input_sets[0] -> file1の要素セット(Label 1用)
    # input_sets[1] -> file2の要素セット(Label 2用) ...
    input_sets = []
    for f in input_files:
        _, cells = read_cells(f)
        input_sets.append(set(cells))

    # 統合ファイルを読み込み
    h_merged, c_merged = read_cells(merged_file)

    labels = []
    unknown = []

    # 統合ファイルの各要素に対してラベル判定
    for i, cell in enumerate(c_merged):
        found_label = 0
        # 入力ファイル順にチェック (file1優先)
        for label_idx, cell_set in enumerate(input_sets):
            if cell in cell_set:
                found_label = label_idx + 1 # 1始まりのラベルにする
                break
        
        labels.append(found_label)
        if found_label == 0:
            unknown.append(i)

    # 出力処理
    # 1行目: #bool_elem
    # 2行目: "<全要素数> 1"
    # 以降: 各ラベルを1行ずつ
    out_lines = ["#bool_elem", f"{len(labels)} 1"] + [str(x) for x in labels]
    Path(out_file).write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    # ログ出力
    sys.stderr.write(f"Merged header: {h_merged}\n")
    sys.stderr.write(f"Total elements: {len(labels)}\n")
    sys.stderr.write(f"Input files processed: {len(input_files)}\n")
    
    for idx, f in enumerate(input_files):
        sys.stderr.write(f"  Label {idx+1}: {f} ({len(input_sets[idx])} elements)\n")

    if unknown:
        sys.stderr.write(f"WARNING: {len(unknown)} cells in merged not found in any input file. First few indices: {unknown[:10]}\n")

if __name__ == "__main__":
    # sys.argv[0] はスクリプト名なので除外して渡す
    main(sys.argv[1:])