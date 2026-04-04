import sys
import numpy as np
import matplotlib.pyplot as plt

# VTK imports
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


def plot_bnode_x_sorted_by_z(vtp_path: str, array_name: str = "B_node") -> None:
    # --- read .vtp ---
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()

    polydata = reader.GetOutput()
    if polydata is None:
        raise RuntimeError("Failed to read the VTP file.")

    # --- coordinates ---
    points_vtk = polydata.GetPoints()
    if points_vtk is None:
        raise RuntimeError("No points found in the VTP file.")

    coords = vtk_to_numpy(points_vtk.GetData())  # shape: (N, 3)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise RuntimeError(f"Unexpected point coordinate shape: {coords.shape}")

    z = coords[:, 2]

    # --- point data array ---
    point_data = polydata.GetPointData()
    if point_data is None:
        raise RuntimeError("No point data found in the VTP file.")

    bnode_vtk = point_data.GetArray(array_name)
    if bnode_vtk is None:
        available = [
            point_data.GetArrayName(i)
            for i in range(point_data.GetNumberOfArrays())
        ]
        raise RuntimeError(
            f'Array "{array_name}" not found.\n'
            f"Available point-data arrays: {available}"
        )

    bnode = vtk_to_numpy(bnode_vtk)

    # Expect vector data like (N, 3)
    if bnode.ndim != 2:
        raise RuntimeError(
            f'"{array_name}" is not a vector array. Shape = {bnode.shape}'
        )
    if bnode.shape[1] < 1:
        raise RuntimeError(
            f'"{array_name}" has no x component. Shape = {bnode.shape}'
        )
    if bnode.shape[0] != z.shape[0]:
        raise RuntimeError(
            f"Point count mismatch: coords={z.shape[0]}, {array_name}={bnode.shape[0]}"
        )

    bnode_x = bnode[:, 0]

    # --- sort by z ascending ---
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    bnode_x_sorted = bnode_x[sort_idx]

    # --- plot ---
    plt.figure(figsize=(8, 5))
    plt.plot(z_sorted, bnode_x_sorted, marker="o", linestyle="-", markersize=3)
    plt.xlim(0.0, 0.4)
    plt.xlabel("z")
    plt.ylabel(f"{array_name}_x")
    plt.title(f"{array_name} x-component sorted by z")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_bnode_x.py input.vtp")
        sys.exit(1)

    vtp_file = sys.argv[1]
    plot_bnode_x_sorted_by_z(vtp_file)