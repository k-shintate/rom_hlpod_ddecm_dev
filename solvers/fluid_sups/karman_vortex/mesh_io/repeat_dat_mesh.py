#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
repeat_dat_mesh.py

1組の node.dat / elem.dat を「1コンポーネント」とみなし、
同じコンポーネントを一定間隔で複製して、
1つの node.dat / elem.dat にまとめる。

さらに、各コピーを1つの subdomain として、
gedatsu_partitioner_simple_mesh の

    --given_subdomain

に渡すための given_subdomain.dat も生成する。

接触・node weld は行わない。
各コンポーネントは独立した disconnected mesh のまま、
同じメッシュファイル内に格納される。

Input
-----
node.dat
    Nnode
    x y z
    ...

elem.dat
    Nelem NEN
    n1 n2 ... nNEN
    ...

Output
------
<out-dir>/node.dat
<out-dir>/elem.dat
<out-dir>/given_subdomain.dat
<out-dir>/repeat_report.txt

given_subdomain.dat
-------------------
GEDATSU の graph format:

    {num_points}
    {point id} {num of contains} {id 1} {id 2} ... {id m}
    ...

ここでは

    1 copy = 1 point = 1 subdomain

として出力する。

GEDATSU 内部ではファイルの

    1行目 -> domain 0
    2行目 -> domain 1
    3行目 -> domain 2
    ...

として扱われる。

Example
-------
python repeat_dat_mesh.py \
    node.dat elem.dat \
    --copies 10 \
    --axis y \
    --gap 0.4 \
    --out-dir merged_mesh

この場合、各コンポーネントの y 方向 bounding-box 長さを Ly として

    pitch = Ly + 0.4

で等間隔配置する。

固定 pitch を直接指定する場合:

python repeat_dat_mesh.py \
    node.dat elem.dat \
    --copies 10 \
    --axis x \
    --pitch 3.0 \
    --out-dir merged_mesh

その後 GEDATSU を例えば次のように実行する:

gedatsu_partitioner_simple_mesh \
    -n 10 \
    -in merged_mesh/node.dat \
    -ie merged_mesh/elem.dat \
    --given_subdomain merged_mesh/given_subdomain.dat \
    -d merged_mesh/parted.0
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def read_nodes(path: Path) -> np.ndarray:
    """
    Read node.dat.

    Format:
        Nnode
        x y z
        ...
    """

    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split()

    if not header:
        raise ValueError(f"{path}: empty file")

    nnode = int(header[0])

    nodes = np.loadtxt(path, skiprows=1, dtype=float)

    if nodes.ndim == 1:
        nodes = nodes.reshape(1, -1)

    if nodes.shape != (nnode, 3):
        raise ValueError(
            f"{path}: header={nnode} nodes, actual shape={nodes.shape}"
        )

    return nodes


def read_elements(path: Path):
    """
    Read elem.dat.

    Format:
        Nelem NEN
        n1 n2 ... nNEN
        ...
    """

    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split()

    if len(header) < 2:
        raise ValueError(
            f"{path}: first line must be 'Nelem NEN'"
        )

    nelem = int(header[0])
    nen = int(header[1])

    elems = np.loadtxt(
        path,
        skiprows=1,
        dtype=np.int64,
    )

    if elems.ndim == 1:
        elems = elems.reshape(1, -1)

    if elems.shape != (nelem, nen):
        raise ValueError(
            f"{path}: header=({nelem}, {nen}), "
            f"actual shape={elems.shape}"
        )

    return elems, nen


def detect_base(
    elems: np.ndarray,
    nnode: int,
) -> int:
    """
    Detect whether connectivity is 0-based or 1-based.

    Returns
    -------
    0
        0-based connectivity

    1
        1-based connectivity
    """

    mn = int(elems.min())
    mx = int(elems.max())

    # 0-based
    if mn == 0:
        if mx >= nnode:
            raise ValueError(
                "Invalid 0-based connectivity: "
                f"max={mx}, nnode={nnode}"
            )

        return 0

    # 1-based
    if mx == nnode:
        if mn < 1:
            raise ValueError(
                "Invalid 1-based connectivity"
            )

        return 1

    raise ValueError(
        "Index base could not be determined automatically. "
        f"Connectivity range=[{mn}, {mx}], nnode={nnode}. "
        "Use a dataset containing node 0 (0-based) "
        "or node N (1-based), "
        "or modify INPUT_BASE in the script."
    )


def write_nodes(
    path: Path,
    nodes: np.ndarray,
):
    """
    Write node.dat.
    """

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(nodes)}\n")

        for x, y, z in nodes:
            f.write(
                f"{x:.16g} "
                f"{y:.16g} "
                f"{z:.16g}\n"
            )


def write_elements(
    path: Path,
    elems: np.ndarray,
    nen: int,
):
    """
    Write elem.dat.
    """

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"{len(elems)} {nen}\n"
        )

        for e in elems:
            f.write(
                " ".join(
                    str(int(v))
                    for v in e
                )
                + "\n"
            )


def write_given_subdomain(
    path: Path,
    nnode_per_component: int,
    copies: int,
    base: int,
):
    """
    Write GEDATSU --given_subdomain graph file.

    Each repeated mesh component is assigned to
    exactly one subdomain.

    Graph format
    ------------
        num_points
        point_id num_contains node_id_1 node_id_2 ...
        ...

    Example: 4 nodes/component, 3 copies, 1-based

        3
        1 4 1 2 3 4
        2 4 5 6 7 8
        3 4 9 10 11 12

    This corresponds to

        domain 0 -> nodes 1,2,3,4
        domain 1 -> nodes 5,6,7,8
        domain 2 -> nodes 9,10,11,12

    GEDATSU determines the domain ID by line order,
    not by point_id:

        line 1 -> domain 0
        line 2 -> domain 1
        ...

    Node IDs use the same index base as elem.dat.
    """

    with path.open("w", encoding="utf-8") as f:

        # Number of graph points = number of domains
        f.write(f"{copies}\n")

        for k in range(copies):

            # Graph point ID.
            #
            # 1-based mesh:
            #   1, 2, 3, ...
            #
            # 0-based mesh:
            #   0, 1, 2, ...
            point_id = k + base

            # Global node ID range belonging to copy k.
            #
            # Example:
            #
            # nnode_per_component = 4
            #
            # 1-based:
            #   copy 0 -> 1 2 3 4
            #   copy 1 -> 5 6 7 8
            #
            # 0-based:
            #   copy 0 -> 0 1 2 3
            #   copy 1 -> 4 5 6 7
            start = (
                k * nnode_per_component
                + base
            )

            end = (
                start
                + nnode_per_component
            )

            node_ids = range(
                start,
                end,
            )

            f.write(
                f"{point_id} "
                f"{nnode_per_component} "
            )

            f.write(
                " ".join(
                    str(node_id)
                    for node_id in node_ids
                )
            )

            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Repeat one node.dat/elem.dat component "
            "at equal spacing without contact or "
            "node welding, and generate a GEDATSU "
            "--given_subdomain file."
        )
    )

    parser.add_argument(
        "node_dat",
        type=Path,
    )

    parser.add_argument(
        "elem_dat",
        type=Path,
    )

    parser.add_argument(
        "--copies",
        type=int,
        required=True,
        help=(
            "Number of component copies. "
            "This is also the number of GEDATSU "
            "subdomains."
        ),
    )

    parser.add_argument(
        "--axis",
        choices=(
            "x",
            "y",
            "z",
        ),
        default="x",
        help=(
            "Direction in which components "
            "are arranged"
        ),
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "--gap",
        type=float,
        help=(
            "Empty distance between neighboring "
            "component bounding boxes. "
            "pitch = component length + gap"
        ),
    )

    group.add_argument(
        "--pitch",
        type=float,
        help=(
            "Constant origin-to-origin translation "
            "between copies"
        ),
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("merged_mesh"),
    )

    parser.add_argument(
        "--start-offset",
        type=float,
        default=0.0,
        help=(
            "Translation of the first copy "
            "along the selected axis"
        ),
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Validate arguments
    # ------------------------------------------------------------

    if args.copies < 1:
        raise ValueError(
            "--copies must be >= 1"
        )

    # ------------------------------------------------------------
    # Read input mesh
    # ------------------------------------------------------------

    nodes = read_nodes(
        args.node_dat
    )

    elems, nen = read_elements(
        args.elem_dat
    )

    nnode = len(nodes)
    nelem = len(elems)

    # ------------------------------------------------------------
    # Detect connectivity index base
    # ------------------------------------------------------------

    base = detect_base(
        elems,
        nnode,
    )

    print(
        f"Detected connectivity base: "
        f"{base}-based"
    )

    # ------------------------------------------------------------
    # Work internally as 0-based
    # ------------------------------------------------------------

    elems0 = elems - base

    # ------------------------------------------------------------
    # Arrangement direction
    # ------------------------------------------------------------

    axis_id = {
        "x": 0,
        "y": 1,
        "z": 2,
    }[args.axis]

    # ------------------------------------------------------------
    # Original component bounding box
    # ------------------------------------------------------------

    xyz_min = nodes.min(axis=0)
    xyz_max = nodes.max(axis=0)

    length = float(
        xyz_max[axis_id]
        - xyz_min[axis_id]
    )

    # ------------------------------------------------------------
    # Determine pitch / gap
    # ------------------------------------------------------------

    if args.pitch is not None:

        pitch = float(
            args.pitch
        )

        gap = (
            pitch
            - length
        )

    else:

        gap = float(
            args.gap
        )

        pitch = (
            length
            + gap
        )

    if pitch <= 0:
        raise ValueError(
            "pitch must be positive. "
            f"component length={length}, "
            f"pitch={pitch}"
        )

    if gap < 0:
        print(
            "WARNING: "
            f"gap={gap} < 0, "
            "so neighboring bounding boxes overlap."
        )

    # ------------------------------------------------------------
    # Repeat mesh
    # ------------------------------------------------------------

    all_nodes = []
    all_elems = []

    for k in range(args.copies):

        # Translation vector
        shift = np.zeros(
            3,
            dtype=float,
        )

        shift[axis_id] = (
            args.start_offset
            + k * pitch
        )

        # --------------------------------------------------------
        # Coordinates
        # --------------------------------------------------------

        nodes_k = (
            nodes
            + shift
        )

        # --------------------------------------------------------
        # Connectivity
        #
        # Each copy gets its own global node IDs.
        # No node welding is performed.
        # --------------------------------------------------------

        node_offset = (
            k * nnode
        )

        elems_k0 = (
            elems0
            + node_offset
        )

        all_nodes.append(
            nodes_k
        )

        all_elems.append(
            elems_k0
        )

    # ------------------------------------------------------------
    # Merge copies
    # ------------------------------------------------------------

    merged_nodes = np.vstack(
        all_nodes
    )

    merged_elems0 = np.vstack(
        all_elems
    )

    # Restore the same index base as the input mesh
    merged_elems = (
        merged_elems0
        + base
    )

    # ------------------------------------------------------------
    # Output directory
    # ------------------------------------------------------------

    args.out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    node_out = (
        args.out_dir
        / "node.dat"
    )

    elem_out = (
        args.out_dir
        / "elem.dat"
    )

    subdomain_out = (
        args.out_dir
        / "given_subdomain.dat"
    )

    report_out = (
        args.out_dir
        / "repeat_report.txt"
    )

    # ------------------------------------------------------------
    # Write mesh
    # ------------------------------------------------------------

    write_nodes(
        node_out,
        merged_nodes,
    )

    write_elements(
        elem_out,
        merged_elems,
        nen,
    )

    # ------------------------------------------------------------
    # Write GEDATSU given-subdomain definition
    #
    # copy 0 -> domain 0
    # copy 1 -> domain 1
    # ...
    # ------------------------------------------------------------

    write_given_subdomain(
        subdomain_out,
        nnode_per_component=nnode,
        copies=args.copies,
        base=base,
    )

    # ------------------------------------------------------------
    # Write report
    # ------------------------------------------------------------

    with report_out.open(
        "w",
        encoding="utf-8",
    ) as f:

        f.write(
            "Repeated component mesh report\n"
        )

        f.write(
            "==============================\n"
        )

        f.write(
            f"input node file     : "
            f"{args.node_dat}\n"
        )

        f.write(
            f"input elem file     : "
            f"{args.elem_dat}\n"
        )

        f.write(
            f"input nodes         : "
            f"{nnode}\n"
        )

        f.write(
            f"input elements      : "
            f"{nelem}\n"
        )

        f.write(
            f"nodes/element       : "
            f"{nen}\n"
        )

        f.write(
            f"index base          : "
            f"{base}\n"
        )

        f.write(
            f"copies              : "
            f"{args.copies}\n"
        )

        f.write(
            f"GEDATSU domains     : "
            f"{args.copies}\n"
        )

        f.write(
            f"axis                : "
            f"{args.axis}\n"
        )

        f.write(
            f"component bbox min  : "
            f"{xyz_min.tolist()}\n"
        )

        f.write(
            f"component bbox max  : "
            f"{xyz_max.tolist()}\n"
        )

        f.write(
            f"component length    : "
            f"{length:.16g}\n"
        )

        f.write(
            f"gap                 : "
            f"{gap:.16g}\n"
        )

        f.write(
            f"pitch               : "
            f"{pitch:.16g}\n"
        )

        f.write(
            f"start offset        : "
            f"{args.start_offset:.16g}\n"
        )

        f.write(
            f"output nodes        : "
            f"{len(merged_nodes)}\n"
        )

        f.write(
            f"output elements     : "
            f"{len(merged_elems)}\n"
        )

        f.write(
            "node welding        : "
            "disabled\n"
        )

        f.write(
            "subdomain strategy  : "
            "one repeated component "
            "per subdomain\n"
        )

        f.write(
            f"given subdomain file: "
            f"{subdomain_out}\n"
        )

        f.write("\n")

        f.write(
            "GEDATSU example command\n"
        )

        f.write(
            "-----------------------\n"
        )

        f.write(
            "gedatsu_partitioner_simple_mesh "
            f"-n {args.copies} "
            f"-in {node_out} "
            f"-ie {elem_out} "
            f"--given_subdomain {subdomain_out} "
            f"-d {args.out_dir / 'parted.0'}\n"
        )

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print()
    print("Done.")
    print()

    print(
        f"component nodes     : "
        f"{nnode}"
    )

    print(
        f"component elements  : "
        f"{nelem}"
    )

    print(
        f"copies              : "
        f"{args.copies}"
    )

    print(
        f"GEDATSU domains     : "
        f"{args.copies}"
    )

    print(
        f"axis                : "
        f"{args.axis}"
    )

    print(
        f"component length    : "
        f"{length}"
    )

    print(
        f"gap                 : "
        f"{gap}"
    )

    print(
        f"pitch               : "
        f"{pitch}"
    )

    print(
        f"output nodes        : "
        f"{len(merged_nodes)}"
    )

    print(
        f"output elements     : "
        f"{len(merged_elems)}"
    )

    print()

    print(
        f"node.dat            : "
        f"{node_out}"
    )

    print(
        f"elem.dat            : "
        f"{elem_out}"
    )

    print(
        f"given_subdomain.dat : "
        f"{subdomain_out}"
    )

    print(
        f"report              : "
        f"{report_out}"
    )

    print()

    print("GEDATSU command:")
    print()

    print(
        "gedatsu_partitioner_simple_mesh "
        f"-n {args.copies} "
        f"-in {node_out} "
        f"-ie {elem_out} "
        f"--given_subdomain {subdomain_out} "
        f"-d {args.out_dir / 'parted.0'}"
    )


if __name__ == "__main__":
    main()