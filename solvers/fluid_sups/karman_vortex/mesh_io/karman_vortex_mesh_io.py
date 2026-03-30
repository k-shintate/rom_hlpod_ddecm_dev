import os
import sys
import numpy as np
import meshio

# コマンドライン引数で入力ファイルを取得
if len(sys.argv) < 2:
    print("Usage: python script_name.py <input_mesh_file>")
    sys.exit(1)

input_file = sys.argv[1]  # コマンドライン引数から入力ファイルを取得

# 出力ディレクトリを指定
output_dir = "./mesh_tmp/"  # 出力先ディレクトリ
combined_output_file = os.path.join(output_dir, "combined_mesh_connectivity.dat")  # 全体メッシュデータの出力ファイル

# 出力ディレクトリが存在しない場合は作成
os.makedirs(output_dir, exist_ok=True)

# メッシュを読み込む
try:
    mesh = meshio.read(input_file)
except Exception as e:
    print(f"Error reading mesh file '{input_file}': {e}")
    sys.exit(1)

# 2Dデータに3Dの0成分を追加（必要な場合）
if mesh.points.shape[1] == 2:
    mesh.points = np.hstack([mesh.points, np.zeros((mesh.points.shape[0], 1))])

# ラベル番号と名前のマッピング（`field_data`から取得）
if mesh.field_data:
    label_names = {v[0]: k for k, v in mesh.field_data.items()}
else:
    print("No field data found. Ensure Physical Names are properly defined.")
    label_names = {}

# 全体メッシュデータを統合して保存する準備
combined_data = []

# メッシュデータをラベルごとに分離して保存
if "gmsh:physical" in mesh.cell_data:
    for cell_block, labels in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
        cell_type = cell_block.type  # セルタイプ（例: triangle, tetra）

        # ラベルごとにコネクティビティデータを出力
        for label in np.unique(labels):
            if label in label_names:
                # 該当ラベルの要素をフィルタ
                element_indices = np.where(labels == label)[0]
                elements = cell_block.data[element_indices]

                # 出力ファイル名を設定
                connectivity_filename = f"{output_dir}{label_names[label]}_{cell_type}_connectivity.dat"

                # ラベルごとのコネクティビティデータを書き込む
                with open(connectivity_filename, "w") as f:
                    # 総要素数と要素内節点数
                    total_elements = len(elements)
                    nodes_per_element = elements.shape[1]
                    f.write(f"{total_elements} {nodes_per_element}\n")

                    # 要素コネクティビティ情報
                    for element in elements:
                        f.write(" ".join(map(str, element)) + "\n")

                print(f"Saved connectivity file: {connectivity_filename}")

                # Fluidの場合は節点座標も出力
                if label_names[label] == "Fluid":
                    nodes_in_elements = np.unique(elements)
                    node_coords = mesh.points[nodes_in_elements]

                    # 節点座標ファイル
                    node_coords_filename = f"{output_dir}{label_names[label]}_node_coordinates.dat"
                    with open(node_coords_filename, "w") as f:
                        f.write(f"{len(nodes_in_elements)}\n")  # 節点数と次元数
                        # 節点座標情報（IDなしで座標のみ）
                        for coord in node_coords:
                            f.write(" ".join(map(str, coord)) + "\n")

                    print(f"Saved node coordinates file: {node_coords_filename}")

                # 全体データに追加
                combined_data.append((label_names[label], cell_type, elements))
else:
    print("No 'gmsh:physical' data found in cell_data. Verify the mesh file.")

# 全体メッシュデータを1つのファイルに保存
with open(combined_output_file, "w") as f:
    for label, cell_type, elements in combined_data:
        f.write(f"# Label: {label}, Type: {cell_type}\n")
        for element in elements:
            f.write(" ".join(map(str, element)) + "\n")
print(f"Saved combined mesh data to {combined_output_file}")
