import sys
import numpy as np
import matplotlib.pyplot as plt

from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


# ------------------------------------------------------------
# Problem 7 Table 2(a)
# Bz along line A1-B1
# y = 72 mm, z = 34 mm
# unit: Bz = 1e-3 T
# ------------------------------------------------------------
X_REF_MM = np.array([
    0, 18, 36, 54, 72, 90, 108, 126, 144,
    162, 180, 198, 216, 234, 252, 270, 288
], dtype=float)

# f = 50 Hz, wt = 0 deg
BZ_A1B1_50_WT0 = np.array([
    -4.9, -17.8, -22.1, -20.1, -15.67, 0.36,
    43.64, 78.11, 71.5, 60.44, 53.82, 58.81,
    56.91, 59.24, 52.78, 27.61, 2.36
], dtype=float) * 1.0e-3


def get_reference_a1b1():
    x_ref = X_REF_MM * 1.0e-3  # mm -> m
    return x_ref, BZ_A1B1_50_WT0


def read_bnode_z_sorted_by_x(vtp_path: str, array_name: str = "B_node"):
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

    x = coords[:, 0]

    point_data = polydata.GetPointData()
    if point_data is None:
        raise RuntimeError(f"No point data found in VTP file: {vtp_path}")

    bnode_vtk = point_data.GetArray(array_name)
    if bnode_vtk is None:
        available = [
            point_data.GetArrayName(i)
            for i in range(point_data.GetNumberOfArrays())
        ]
        raise RuntimeError(
            f'Array "{array_name}" not found.\n'
            f"Available arrays: {available}"
        )

    bnode = vtk_to_numpy(bnode_vtk)

    if bnode.ndim != 2 or bnode.shape[1] < 3:
        raise RuntimeError(
            f'"{array_name}" must be vector data with 3 components. '
            f"Shape = {bnode.shape}"
        )

    # Bz component
    bnode_z = bnode[:, 2]

    sort_idx = np.argsort(x)

    x_sorted = x[sort_idx]
    bnode_z_sorted = bnode_z[sort_idx]

    return x_sorted, bnode_z_sorted


def plot_bz_a1b1_comparison(vtp_paths, array_name: str = "B_node"):
    plt.figure(figsize=(9, 6))

    # VTP simulation result
    for i, vtp_path in enumerate(vtp_paths, start=1):
        x_sorted, bz_sorted = read_bnode_z_sorted_by_x(
            vtp_path,
            array_name=array_name
        )

        plt.plot(
            x_sorted,
            bz_sorted,
            linestyle="-",
            marker="o",
            markersize=3,
            label=f"Simulation {i}: {vtp_path}"
        )

    # Reference data
    x_ref, bz_ref = get_reference_a1b1()

    plt.plot(x_ref,bz_ref,linestyle="None",marker="o",markersize=6,label="Reference A1-B1, f=50 Hz, wt=0°")

    plt.xlim(0.0, 0.288)
    plt.xlabel("x [m]")
    plt.ylabel("Bz [T]")
    plt.title("Comparison of Bz along A1-B1: y=72 mm, z=34 mm")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python plot_bz_a1b1.py input1.vtp")
        print("  python plot_bz_a1b1.py input1.vtp input2.vtp")
        sys.exit(1)

    if len(sys.argv) > 3:
        print("Error: please specify at most two VTP files.")
        sys.exit(1)

    vtp_files = sys.argv[1:]
    plot_bz_a1b1_comparison(vtp_files)