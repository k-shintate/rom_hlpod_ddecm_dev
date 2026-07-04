import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


def get_array_names(point_data):
    return [
        point_data.GetArrayName(i)
        for i in range(point_data.GetNumberOfArrays())
    ]


def read_velocity_x_sorted_by_y(
    vtp_path: str,
    array_name: str = "Velocity",
):
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()

    polydata = reader.GetOutput()
    if polydata is None:
        raise RuntimeError(f"Failed to read VTP file: {vtp_path}")

    points_vtk = polydata.GetPoints()
    if points_vtk is None:
        raise RuntimeError(f"No points found in VTP file: {vtp_path}")

    coords = vtk_to_numpy(points_vtk.GetData())

    point_data = polydata.GetPointData()
    if point_data is None:
        raise RuntimeError(f"No point data found in VTP file: {vtp_path}")

    velocity_vtk = point_data.GetArray(array_name)
    if velocity_vtk is None:
        available = get_array_names(point_data)
        raise RuntimeError(
            f'Array "{array_name}" not found in {vtp_path}\n'
            f"Available PointData arrays: {available}"
        )

    velocity = vtk_to_numpy(velocity_vtk)

    if velocity.ndim != 2 or velocity.shape[1] < 3:
        raise RuntimeError(
            f'"{array_name}" must be vector data with 3 components. '
            f"Shape = {velocity.shape}"
        )

    y = coords[:, 1]
    velocity_x = velocity[:, 0]

    sort_idx = np.argsort(y)

    y_sorted = y[sort_idx]
    velocity_x_sorted = velocity_x[sort_idx]

    return y_sorted, velocity_x_sorted


def plot_velocity_x_vs_y(vtp_paths, array_name: str = "Velocity"):
    plt.figure(figsize=(7, 9))

    for i, vtp_path in enumerate(vtp_paths, start=1):
        y_sorted, vx_sorted = read_velocity_x_sorted_by_y(
            vtp_path,
            array_name=array_name,
        )

        plt.plot(
            vx_sorted,
            y_sorted,
            linestyle="-",
            marker="o",
            markersize=3,
            label=f"File {i}: {Path(vtp_path).name}",
        )

    plt.xlabel("x-direction velocity [m/s]")
    plt.ylabel("y [m]")
    plt.title("y-coordinate vs x-direction velocity")

    plt.ylim(0.0, 3.0)

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage:")
        print("  python plot_velocity_x_vs_y.py file1.vtp file2.vtp file3.vtp file4.vtp")
        sys.exit(1)

    vtp_files = sys.argv[1:]

    plot_velocity_x_vs_y(
        vtp_files,
        array_name="Velocity",
    )