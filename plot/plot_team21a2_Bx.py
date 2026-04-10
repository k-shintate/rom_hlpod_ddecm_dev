import sys
import numpy as np
import matplotlib.pyplot as plt

# VTK imports
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy


# ---------------------------------------------------------------------
# Table A3-2: RESULTS OF BX FOR PROBLEM 21a-2
# (Presented at TEAM-Yichang, China, 1993)
#
# Columns:
# Z(mm),
# X=+5.76mm Measured (×10^-4 T),
# X=+5.76mm Calculated (×10^-4 T),
# X=-5.76mm Measured (×10^-4 T),
# X=-5.76mm Calculated (×10^-4 T)
# ---------------------------------------------------------------------
TABLE_A3_2 = np.array([
    [  6.0,   81.25,   81.70,   64.60,   65.60],
    [ 30.0,   63.80,   63.10,   53.40,   54.28],
    [ 66.0,   34.70,   37.50,   31.60,   34.63],
    [102.0,   16.20,   20.20,   14.90,   17.79],
    [139.0,    1.00,    1.96,    1.30,    0.43],
    [175.0,  -15.50,  -15.04,  -12.50,  -13.02],
    [211.0,  -35.10,  -30.73,  -27.50,  -25.25],
    [230.0,  -42.00,  -41.72,  -32.20,  -31.66],
    [246.0,  -37.20,  -36.87,  -29.90,  -29.53],
    [280.0,  -23.70,  -24.26,  -20.70,  -21.90],
    [313.0,  -15.30,  -15.10,  -13.60,  -14.36],
    [344.0,  -11.20,  -10.18,   -9.40,   -9.80],
    [372.0,   -8.70,   -7.17,   -6.80,   -7.07],
    [398.0,   -7.60,   -5.27,   -5.20,   -5.25],
], dtype=float)


def get_table_a3_2_in_si():
    """
    Convert:
      Z: mm -> m
      Bx: (×10^-4 T) -> T
    """
    z_m = TABLE_A3_2[:, 0] * 1.0e-3
    xp_meas_T = TABLE_A3_2[:, 1] * 1.0e-4
    xp_calc_T = TABLE_A3_2[:, 2] * 1.0e-4
    xm_meas_T = TABLE_A3_2[:, 3] * 1.0e-4
    xm_calc_T = TABLE_A3_2[:, 4] * 1.0e-4
    return z_m, xp_meas_T, xp_calc_T, xm_meas_T, xm_calc_T


def read_bnode_x_sorted_by_z(vtp_path: str, array_name: str = "B_node"):
    # --- read .vtp ---
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()

    polydata = reader.GetOutput()
    if polydata is None:
        raise RuntimeError(f"Failed to read the VTP file: {vtp_path}")

    # --- coordinates ---
    points_vtk = polydata.GetPoints()
    if points_vtk is None:
        raise RuntimeError(f"No points found in the VTP file: {vtp_path}")

    coords = vtk_to_numpy(points_vtk.GetData())  # shape: (N, 3)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise RuntimeError(f"Unexpected point coordinate shape in {vtp_path}: {coords.shape}")

    z = coords[:, 2]

    # --- point data array ---
    point_data = polydata.GetPointData()
    if point_data is None:
        raise RuntimeError(f"No point data found in the VTP file: {vtp_path}")

    bnode_vtk = point_data.GetArray(array_name)
    if bnode_vtk is None:
        available = [
            point_data.GetArrayName(i)
            for i in range(point_data.GetNumberOfArrays())
        ]
        raise RuntimeError(
            f'Array "{array_name}" not found in {vtp_path}.\n'
            f"Available point-data arrays: {available}"
        )

    bnode = vtk_to_numpy(bnode_vtk)

    # Expect vector data like (N, 3)
    if bnode.ndim != 2:
        raise RuntimeError(
            f'"{array_name}" is not a vector array in {vtp_path}. Shape = {bnode.shape}'
        )
    if bnode.shape[1] < 1:
        raise RuntimeError(
            f'"{array_name}" has no x component in {vtp_path}. Shape = {bnode.shape}'
        )
    if bnode.shape[0] != z.shape[0]:
        raise RuntimeError(
            f"Point count mismatch in {vtp_path}: coords={z.shape[0]}, {array_name}={bnode.shape[0]}"
        )

    bnode_x = bnode[:, 0]

    # --- sort by z ascending ---
    sort_idx = np.argsort(z)
    z_sorted = z[sort_idx]
    bnode_x_sorted = bnode_x[sort_idx]

    return z_sorted, bnode_x_sorted


def plot_bnode_x_sorted_by_z(vtp_paths, array_name: str = "B_node") -> None:
    plt.figure(figsize=(9, 6))

    # --- plot VTP results ---
    for i, vtp_path in enumerate(vtp_paths, start=1):
        z_sorted, bnode_x_sorted = read_bnode_x_sorted_by_z(vtp_path, array_name=array_name)
        plt.plot(
            z_sorted,
            bnode_x_sorted,
            marker="o",
            linestyle="-",
            markersize=3,
            label=f"{array_name}_x from VTP #{i}: {vtp_path}"
        )

    # --- reference data from Table A3-2 ---
    z_ref, xp_meas, xp_calc, xm_meas, xm_calc = get_table_a3_2_in_si()



    plt.xlim(0.0, 0.4)
    plt.xlabel("z [m]")
    plt.ylabel("Bx [T]")
    plt.title(f"{array_name} x-component sorted by z + Table A3-2 reference")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python plot_bnode_x.py input1.vtp")
        print("  python plot_bnode_x.py input1.vtp input2.vtp")
        sys.exit(1)

    if len(sys.argv) > 3:
        print("Error: please specify at most two VTP files.")
        print("Usage:")
        print("  python plot_bnode_x.py input1.vtp")
        print("  python plot_bnode_x.py input1.vtp input2.vtp")
        sys.exit(1)

    vtp_files = sys.argv[1:]
    plot_bnode_x_sorted_by_z(vtp_files)