import argparse
import glob
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pyvista as pv


def natural_key(path):
    nums = re.findall(r"\d+", Path(path).name)
    return [int(x) for x in nums] if nums else [0]


def read_any_vtk(pattern):
    files = sorted(glob.glob(pattern), key=natural_key)

    if not files:
        raise FileNotFoundError(f"No files matched: {pattern}")

    # result_000001.vtk.0 のような拡張子だと読めないことがあるので一時的に .vtk にする
    tmp = tempfile.TemporaryDirectory()
    tmpdir = Path(tmp.name)

    blocks = pv.MultiBlock()

    for i, f in enumerate(files):
        src = Path(f)
        if src.suffix.lower() == ".vtk":
            read_path = src
        else:
            read_path = tmpdir / f"piece_{i:04d}.vtk"
            shutil.copy(src, read_path)

        blocks.append(pv.read(read_path))

    if len(blocks) == 1:
        mesh = blocks[0]
    else:
        mesh = blocks.combine()

    mesh._tmpdir_ref = tmp  # 一時ファイルが消えないよう保持
    return mesh


def print_arrays(mesh):
    print("\nPoint data:")
    for name in mesh.point_data:
        arr = mesh.point_data[name]
        print(f"  {name}: shape={arr.shape}")

    print("\nCell data:")
    for name in mesh.cell_data:
        arr = mesh.cell_data[name]
        print(f"  {name}: shape={arr.shape}")


def find_vector_array(mesh):
    preferred = ["Velocity", "velocity", "Vel", "vel", "U", "u"]

    for name in preferred:
        if name in mesh.point_data and mesh.point_data[name].ndim == 2 and mesh.point_data[name].shape[1] == 3:
            return name, "point"
        if name in mesh.cell_data and mesh.cell_data[name].ndim == 2 and mesh.cell_data[name].shape[1] == 3:
            return name, "cell"

    for name in mesh.point_data:
        arr = mesh.point_data[name]
        if arr.ndim == 2 and arr.shape[1] == 3:
            return name, "point"

    for name in mesh.cell_data:
        arr = mesh.cell_data[name]
        if arr.ndim == 2 and arr.shape[1] == 3:
            return name, "cell"

    raise ValueError("3成分の速度ベクトルらしい配列が見つかりません。--velocity で指定してください。")


def ensure_point_vector(mesh, velocity_name, association):
    if association == "point":
        return mesh, velocity_name

    # Cell data の速度を point data に変換
    mesh2 = mesh.cell_data_to_point_data()
    if velocity_name not in mesh2.point_data:
        raise ValueError(f"Cell data から point data への変換後に {velocity_name} が見つかりません。")

    return mesh2, velocity_name


def add_vorticity(mesh, velocity_name):
    # PyVistaでは vector でも scalars= に配列名を渡す
    out = mesh.compute_derivative(
        scalars=velocity_name,
        gradient=False,
        divergence=False,
        vorticity=True,
        qcriterion=False,
        preference="point",
    )

    # 環境によって名前が Vorticity / vorticity になる可能性を吸収
    vort_name = None
    for candidate in ["vorticity", "Vorticity"]:
        if candidate in out.point_data:
            vort_name = candidate
            break

    if vort_name is None:
        raise RuntimeError("渦度配列が生成されませんでした。速度場が3D volume上にあるか確認してください。")

    vort = out.point_data[vort_name]
    out.point_data["vorticity_mag"] = np.linalg.norm(vort, axis=1)

    return out


def choose_iso_value(mesh, percentile=None, value=None):
    w = np.asarray(mesh.point_data["vorticity_mag"])
    w = w[np.isfinite(w)]
    w = w[w > 0]

    if len(w) == 0:
        raise ValueError("vorticity_mag が全て0または不正値です。")

    if value is not None:
        return value

    if percentile is None:
        percentile = 97.0

    return float(np.percentile(w, percentile))


def find_color_array(mesh, requested):
    if requested:
        if requested in mesh.point_data:
            return requested
        if requested in mesh.cell_data:
            return requested
        print(f"Warning: color array '{requested}' が見つかりません。vorticity_mag で着色します。")
        return "vorticity_mag"

    preferred = ["Mach", "mach", "MachNumber", "mach_number", "Pressure", "pressure", "p"]
    for name in preferred:
        if name in mesh.point_data or name in mesh.cell_data:
            return name

    return "vorticity_mag"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help='例: "result_000001.vtk.*" または "flow.vtk"')
    parser.add_argument("--velocity", default=None, help="速度ベクトル配列名")
    parser.add_argument("--color", default=None, help="着色に使う配列名。例: Mach")
    parser.add_argument("--iso", type=float, default=None, help="渦度等値面の値を直接指定")
    parser.add_argument("--percentile", type=float, default=97.0, help="渦度の上位何%を等値面にするか")
    parser.add_argument("--out", default="vorticity_iso.png", help="出力画像")
    parser.add_argument("--save-iso", default=None, help="等値面を .vtp などで保存")
    parser.add_argument("--list-arrays", action="store_true", help="配列名だけ表示して終了")
    args = parser.parse_args()

    mesh = read_any_vtk(args.input)

    print(mesh)
    print_arrays(mesh)

    if args.list_arrays:
        return

    if args.velocity:
        if args.velocity in mesh.point_data:
            velocity_name, association = args.velocity, "point"
        elif args.velocity in mesh.cell_data:
            velocity_name, association = args.velocity, "cell"
        else:
            raise ValueError(f"--velocity {args.velocity} が見つかりません。")
    else:
        velocity_name, association = find_vector_array(mesh)

    print(f"\nUsing velocity array: {velocity_name} ({association} data)")

    mesh, velocity_name = ensure_point_vector(mesh, velocity_name, association)
    mesh = add_vorticity(mesh, velocity_name)

    iso_value = choose_iso_value(mesh, percentile=args.percentile, value=args.iso)
    print(f"Using iso value: {iso_value}")

    iso = mesh.contour(isosurfaces=[iso_value], scalars="vorticity_mag")

    if iso.n_points == 0:
        raise RuntimeError("等値面が空です。--iso を下げるか --percentile を 90 などにしてください。")

    color_name = find_color_array(iso, args.color)
    print(f"Coloring by: {color_name}")

    if args.save_iso:
        iso.save(args.save_iso)
        print(f"Saved iso surface: {args.save_iso}")

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(
        iso,
        scalars=color_name,
        smooth_shading=True,
        show_scalar_bar=True,
    )
    plotter.add_axes()
    plotter.show_bounds()
    plotter.camera_position = "iso"
    plotter.screenshot(args.out, window_size=(1600, 1000))

    print(f"Saved screenshot: {args.out}")


if __name__ == "__main__":
    main()