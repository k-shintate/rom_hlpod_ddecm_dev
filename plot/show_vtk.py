import vtk
import pyvista as pv
import glob
import os
import re
from collections import defaultdict
import numpy as np
import imageio
import argparse

#コマンドの例
#ffmpeg -framerate 30 -i ./../../data_karmanvortex/data12/visualize_data.%04d.png -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p ./../../data_karmanvortex/data12/output.mp4
#vlc ./../../data_karmanvortex/data10/output.mp4 

class VTKVisualizer:
    def __init__(self, prefix="result", mode="timeseries"):
        """
        Args:
            prefix (str): ファイル名のプレフィックス
            mode (str): 'timeseries' または 'single'
        """
        self.time_steps = []
        self.prefix = prefix
        self.mode = mode

    def load_files(self, folder_path, pattern=None):
        """
        VTKファイルを読み込む（タイムステップ有無に関わらず）

        Args:
            folder_path (str): VTKファイルが存在するディレクトリのパス
            pattern (str, optional): グロブパターン。指定がない場合はprefixを使用 (例: "result_*.vtk*")

        Returns:
            list: 読み込んだファイルのリスト
        """
        if pattern is None:
            pattern = f"{self.prefix}*.vtk*"

        matched_files = glob.glob(os.path.join(folder_path, pattern))
        print(f"Glob pattern '{pattern}' matched {len(matched_files)} files.")
        print("Matched files:")
        for f in matched_files:
            print(f" - {f}")

        if self.mode == "timeseries":
            # タイムステップモード: ファイルをタイムステップごとにグループ化
            regex_pattern = rf"^{re.escape(self.prefix)}[_.](\d+)\.vtk(?:\.\d+)?$"
            file_groups = defaultdict(list)

            for filepath in matched_files:
                base_name = os.path.basename(filepath)
                print(f"Checking file: {base_name}")
                match = re.match(regex_pattern, base_name)
                if match:
                    time_step = int(match.group(1))
                    file_groups[time_step].append(filepath)
                    print(f"Grouped file: {filepath} (Time Step: {time_step})")
                else:
                    print(f"Ignored non-matching or incomplete file: {filepath}")

            self.time_steps = sorted(file_groups.items())
            print(f"Total valid time steps loaded: {len(self.time_steps)}")
            return len(self.time_steps)

        elif self.mode == "single":
            # シングルモード: 全ファイルを一つのタイムステップとして扱う
            all_files = sorted(matched_files)
            self.time_steps = [(0, all_files)]  # タイムステップ0にすべてのファイルをまとめる
            print(f"Total files loaded and grouped into one time step: {len(all_files)}")
            return len(all_files)
        else:
            print(f"Unsupported mode: {self.mode}")
            return 0

    def read_and_merge_files_vtk(self, file_paths, time_step):
        """
        複数のメタファイルを読み込み、1つのメッシュにマージする

        Args:
            file_paths (list of str): 読み込むメタファイルのパス一覧
            time_step (int): タイムステップ番号（マージ後のファイル名に使用）

        Returns:
            pyvista.UnstructuredGrid または None: マージされたVTKデータ
        """
        if not file_paths:
            print("No files provided for merging.")
            return None

        if len(file_paths) == 1:
            # ファイルが1つの場合、そのまま読み込む
            return self.read_vtk_file(file_paths[0])

        # vtkAppendFilterを使用して複数のメッシュをマージ
        append_filter = vtk.vtkAppendFilter()
        for file_path in file_paths:
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(file_path)
            reader.Update()
            mesh = reader.GetOutput()
            if mesh.GetNumberOfPoints() == 0:
                print(f"Warning: Mesh from {file_path} has no points.")
                continue
            append_filter.AddInputData(mesh)
            print(f"Added mesh from {file_path} to append filter.")

        if append_filter.GetNumberOfInputConnections(0) == 0:
            print("Error: No valid meshes to append.")
            return None

        append_filter.Update()
        combined_mesh = pv.wrap(append_filter.GetOutput())
        print(f"Merged {append_filter.GetNumberOfInputConnections(0)} meshes into a single mesh.")

        # マージ後のメッシュを一時ファイルとして保存
        temp_merged_file = os.path.join(os.path.dirname(file_paths[0]), f"{self.prefix}_merged_{time_step}.vtk")
        combined_mesh.save(temp_merged_file)
        print(f"Merged mesh saved to {temp_merged_file}")

        return combined_mesh

    def read_vtk_file(self, file_path):
        """
        単一のVTKファイルを読み込む

        Args:
            file_path (str): 読み込むVTKファイルのパス

        Returns:
            pyvista.UnstructuredGrid または None: 読み込んだVTKデータ
        """
        try:
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(file_path)
            reader.Update()
            mesh = reader.GetOutput()
            if mesh.GetNumberOfPoints() == 0:
                print(f"Warning: Mesh from {file_path} has no points.")
                return None
            pyvista_mesh = pv.wrap(mesh)
            print(f"Successfully read file: {file_path}")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return None

        # デバッグ情報の表示
        print(f"Mesh type: {type(pyvista_mesh)}")
        if pyvista_mesh.point_arrays:
            print(f"Point Arrays: {list(pyvista_mesh.point_arrays.keys())}")
        else:
            print("No point arrays found.")
        if pyvista_mesh.cell_arrays:
            print(f"Cell Arrays: {list(pyvista_mesh.cell_arrays.keys())}")
        else:
            print("No cell arrays found.")

        return pyvista_mesh

    def create_animation(self, output_path, camera_position=None):
        """
        アニメーションを作成する

        Args:
            output_path (str): 出力動画ファイルのパス
            camera_position (tuple, optional): カメラの位置を指定 (デフォルト: 自動)
        """
        # 各タイムステップのメッシュを読み込み
        meshes = []
        for time_step, files in self.time_steps:
            print(f"Processing time step {time_step} with files: {files}")
            # メタファイルを読み込み、必要に応じてマージ
            mesh = self.read_and_merge_files_vtk(files, time_step)
            if mesh is None:
                print(f"Skipping time step {time_step} due to read or merge errors.")
                continue
            meshes.append((time_step, mesh))
            print(f"Time step {time_step} loaded and merged.")

        if not meshes:
            print("No valid meshes to visualize.")
            return

        # アニメーション用のフレームを保存するリスト
        frames = []

        # カラーマップの一貫性を保つためにスカラー値の範囲を決定
        scalar_values = []
        scalar_name = None  # スカラー配列の名前を保持
        for _, mesh in meshes:
            if mesh.point_arrays:
                current_scalar = next(iter(mesh.point_arrays))
                scalar_values.append(mesh.point_arrays[current_scalar])
                scalar_name = current_scalar
            elif mesh.cell_arrays:
                current_scalar = next(iter(mesh.cell_arrays))
                scalar_values.append(mesh.cell_arrays[current_scalar])
                scalar_name = current_scalar

        if scalar_values and scalar_name:
            all_scalars = np.concatenate([arr for arr in scalar_values if arr is not None])
            vmin = np.min(all_scalars)
            vmax = np.max(all_scalars)
            print(f"Scalar range across all meshes: {vmin} to {vmax}")
        else:
            vmin = None
            vmax = None
            print("No scalar data found across all meshes.")

        # Plotterの初期化
        plotter = pv.Plotter(off_screen=True)
        plotter.add_axes()
        plotter.set_background("white")

        # 各タイムステップごとにフレームを作成
        for idx, (time_step, mesh) in enumerate(meshes):
            plotter.clear()
            plotter.add_axes()

            # メッシュの属性をデバッグ出力
            print(f"Time Step {time_step}:")
            print(f" - Mesh type: {type(mesh)}")
            print(f" - Available point arrays: {list(mesh.point_arrays.keys()) if mesh.point_arrays else 'None'}")
            print(f" - Available cell arrays: {list(mesh.cell_arrays.keys()) if mesh.cell_arrays else 'None'}")

            # スカラー配列が存在するか確認
            if mesh.point_arrays and scalar_name in mesh.point_arrays:
                mesh_scalar = mesh.point_arrays.get(scalar_name)
                plotter.add_mesh(mesh, scalars=mesh_scalar, cmap="viridis", clim=[vmin, vmax])
                print(f" - Using point scalar '{scalar_name}'.")
            elif mesh.cell_arrays and scalar_name in mesh.cell_arrays:
                mesh_scalar = mesh.cell_arrays.get(scalar_name)
                plotter.add_mesh(mesh, scalars=mesh_scalar, cmap="viridis", preference='cells', clim=[vmin, vmax])
                print(f" - Using cell scalar '{scalar_name}'.")
            else:
                # スカラー配列が存在しない場合は単色で表示
                plotter.add_mesh(mesh, color="lightblue")
                print(" - No scalar data. Displaying in light blue.")

            plotter.add_text(f"Step: {time_step}", position="upper_left", font_size=12)

            if camera_position:
                plotter.camera_position = camera_position

            # 画像を取得してフレームに追加
            try:
                img = plotter.screenshot(None, return_img=True)
                frames.append(img)
                print(f" - Frame for time step {time_step} captured.")
            except Exception as e:
                print(f" - Error capturing frame for time step {time_step}: {e}")

        plotter.close()

        # アニメーションの保存
        if output_path.lower().endswith('.gif'):
            try:
                imageio.mimsave(output_path, frames, fps=10)
                print(f"アニメーションをGIF形式で {output_path} に保存しました。")
            except Exception as e:
                print(f"Error saving GIF animation: {e}")
        elif output_path.lower().endswith(('.mp4', '.avi')):
            try:
                writer = imageio.get_writer(output_path, fps=10)
                for frame in frames:
                    writer.append_data(frame)
                writer.close()
                print(f"アニメーションを動画形式で {output_path} に保存しました。")
            except Exception as e:
                print(f"Error saving video animation: {e}")
        else:
            print("エラー: サポートされていないファイル形式です。拡張子を .gif, .mp4, または .avi にしてください。")

def main():
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description='VTKファイルの可視化とアニメーションを作成')
    parser.add_argument('input_dir', help='VTKファイルが存在するディレクトリのパス')
    parser.add_argument('--output', '-o', default='animation.mp4',
                        help='出力動画ファイルのパス (デフォルト: animation.mp4)')
    parser.add_argument('--pattern', type=str, default=None,
                        help='VTKファイルのグロブパターン。指定しない場合はprefixを使用 (例: "result_*.vtk*")')
    parser.add_argument('--prefix', type=str, default='result',
                        help='VTKファイルのプレフィックス (デフォルト: "result")')
    parser.add_argument('--camera', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                        help='カメラの位置を指定 (例: --camera 1 1 1)')
    parser.add_argument('--mode', type=str, choices=['timeseries', 'single'], default='timeseries',
                        help='データのモードを指定: "timeseries"（デフォルト）または "single"')

    args = parser.parse_args()

    # 入力ディレクトリの存在確認
    if not os.path.exists(args.input_dir):
        print(f"エラー: 指定されたディレクトリ '{args.input_dir}' が見つかりません。")
        return

    # 出力ディレクトリの作成（必要な場合）
    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"出力ディレクトリ '{output_dir}' を作成しました。")
        except Exception as e:
            print(f"出力ディレクトリの作成中にエラーが発生しました: {e}")
            return

    # 可視化の実行
    visualizer = VTKVisualizer(prefix=args.prefix, mode=args.mode)
    num_steps = visualizer.load_files(args.input_dir, args.pattern)

    if num_steps > 0:
        print(f"処理するエントリー数: {num_steps}")
        print(f"入力ディレクトリ: {args.input_dir}")
        print(f"出力ファイル: {args.output}")

        camera_position = tuple(args.camera) if args.camera else None

        visualizer.create_animation(
            output_path=args.output,
            camera_position=camera_position
        )
    else:
        print("エラー: 指定されたディレクトリに処理対象のVTKファイルが見つかりません。")

if __name__ == "__main__":
    main()
