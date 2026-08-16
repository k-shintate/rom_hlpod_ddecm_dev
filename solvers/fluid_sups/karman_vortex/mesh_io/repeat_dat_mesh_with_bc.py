#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
repeat_dat_mesh_with_bc.py

1組の node.dat / elem.dat を1コンポーネントとして等間隔に複製し、
以下を生成する。

- node.dat
- elem.dat
- given_subdomain.dat
- repeat_report.txt
- D_bc.dat              (--bc 指定時)

通常は各コピーを disconnected mesh として複製する。
--continuous を指定した場合は、bounding-box 面ではなく、ユーザーが指定した
「左共有面 connectivity」と「右共有面 connectivity」から共有節点を取得して
隣接コピー間を node welding する。

連続モードの処理
----------------
1. --interface-left / --interface-right の surface connectivity を読む。
2. 各 surface connectivity に現れる node ID の unique 集合を取得する。
3. copy k の right surface と copy k+1 の left surface が一致するよう、
   right 側座標と「pitch 分だけ平行移動した left 側座標」を比較する。
4. 座標が weld_tol 内で1対1対応することを確認する。
5. copy k+1 の left interface node は、copy k の right interface node の
   global node ID を再利用する。
6. copy k+1 の全 volume element connectivity を local->global map で張り替える。

したがって、共有面を含む両側の体積要素は同一の interface node ID を参照し、
新しい bridging element を追加せずに conformal な連続メッシュになる。

注意
----
- --continuous では --interface-left と --interface-right が必須。
- 左右 surface は同じ節点離散化を持つ必要がある。
- surface connectivity の node ID base は node.dat/elem.dat と同じものを使う。
- --gap から計算した pitch で左右 surface が一致しない場合は --pitch を直接指定する。
- given_subdomain.dat では shared node を二重所属させず、前の copy/domain を owner とする。
- --continuous で bounding-box BC を使う場合、repeat-axis の内部 xmin/xmax 等は BC から除外する。

Surface connectivity format
---------------------------
volume elem.dat と同様に、例えば quad 面なら:

    Nface 4
    n1 n2 n3 n4
    ...

Example
-------
python repeat_dat_mesh_with_bc.py \
    node.dat elem.dat \
    --copies 4 \
    --axis x \
    --gap 0.0 \
    --continuous \
    --interface-left left_quad_connectivity.dat \
    --interface-right right_quad_connectivity.dat \
    --out-dir merged_mesh

BC も生成する場合:

python repeat_dat_mesh_with_bc.py \
    node.dat elem.dat \
    --copies 4 \
    --axis x \
    --gap 0.0 \
    --continuous \
    --interface-left left_quad_connectivity.dat \
    --interface-right right_quad_connectivity.dat \
    --block-length 3 \
    --bc xmin:1,2,3 \
    --bc xmax:1,2,3 \
    --bc ymin:1,2,3 \
    --bc ymax:1,2,3 \
    --bc zmin:1,2,3 \
    --bc zmax:1,2,3 \
    --bc-value 0.0 \
    --out-dir merged_mesh
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


# -----------------------------------------------------------------------------
# I/O
# -----------------------------------------------------------------------------


def read_nodes(path: Path) -> np.ndarray:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split()

    if not header:
        raise ValueError(f"{path}: empty file")

    nnode = int(header[0])
    nodes = np.loadtxt(path, skiprows=1, dtype=float)

    if nodes.size == 0:
        nodes = np.empty((0, 3), dtype=float)
    elif nodes.ndim == 1:
        nodes = nodes.reshape(1, -1)

    if nodes.shape != (nnode, 3):
        raise ValueError(
            f"{path}: header={nnode} nodes, actual shape={nodes.shape}"
        )

    return nodes


def read_elements(path: Path) -> Tuple[np.ndarray, int]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split()

    if len(header) < 2:
        raise ValueError(f"{path}: first line must be 'Nelem NEN'")

    nelem = int(header[0])
    nen = int(header[1])
    elems = np.loadtxt(path, skiprows=1, dtype=np.int64)

    if elems.size == 0:
        elems = np.empty((0, nen), dtype=np.int64)
    elif elems.ndim == 1:
        elems = elems.reshape(1, -1)

    if elems.shape != (nelem, nen):
        raise ValueError(
            f"{path}: header=({nelem}, {nen}), actual shape={elems.shape}"
        )

    return elems, nen


def detect_base(elems: np.ndarray, nnode: int, requested: str) -> int:
    if requested in ("0", "1"):
        base = int(requested)

        if elems.size:
            mn = int(elems.min())
            mx = int(elems.max())

            if base == 0 and (mn < 0 or mx >= nnode):
                raise ValueError(
                    f"Connectivity is incompatible with --input-base 0: "
                    f"range=[{mn}, {mx}], nnode={nnode}"
                )

            if base == 1 and (mn < 1 or mx > nnode):
                raise ValueError(
                    f"Connectivity is incompatible with --input-base 1: "
                    f"range=[{mn}, {mx}], nnode={nnode}"
                )

        return base

    if elems.size == 0:
        return 0

    mn = int(elems.min())
    mx = int(elems.max())

    if mn == 0:
        if mx >= nnode:
            raise ValueError(
                f"Invalid 0-based connectivity: max={mx}, nnode={nnode}"
            )
        return 0

    if mx == nnode:
        if mn < 1:
            raise ValueError("Invalid 1-based connectivity")
        return 1

    if 0 < mn and mx < nnode:
        raise ValueError(
            "Index base is ambiguous. Specify --input-base 0 or --input-base 1. "
            f"Connectivity range=[{mn}, {mx}], nnode={nnode}"
        )

    raise ValueError(
        f"Invalid connectivity range [{mn}, {mx}] for nnode={nnode}"
    )


def write_nodes(path: Path, nodes: np.ndarray) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(nodes)}\n")
        for x, y, z in nodes:
            f.write(f"{x:.16g} {y:.16g} {z:.16g}\n")


def write_elements(
    path: Path,
    elems_zero: np.ndarray,
    nen: int,
    base: int,
) -> None:
    elems_out = elems_zero + base

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(elems_out)} {nen}\n")
        for elem in elems_out:
            f.write(" ".join(str(int(v)) for v in elem) + "\n")


def write_given_subdomain(
    path: Path,
    domain_nodes_zero: Sequence[Sequence[int]],
    base: int,
) -> None:
    """
    Write GEDATSU --given_subdomain graph file.

    domain_nodes_zero[k] contains the zero-based global node IDs OWNED by domain k.
    Shared interface nodes must occur in exactly one domain.
    """

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(domain_nodes_zero)}\n")

        for k, ids_zero in enumerate(domain_nodes_zero):
            ids_zero = sorted(int(v) for v in ids_zero)

            if not ids_zero:
                raise ValueError(
                    f"Subdomain {k} owns zero nodes. "
                    "GEDATSU requires every domain to contain at least one node."
                )

            point_id = k + base
            ids_out = [v + base for v in ids_zero]

            f.write(f"{point_id} {len(ids_out)}")
            if ids_out:
                f.write(" " + " ".join(str(v) for v in ids_out))
            f.write("\n")


def write_bc(
    path: Path,
    bc_rows: Sequence[Tuple[int, int]],
    block_length: int,
    value: float = 0.0,
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(bc_rows)} {block_length}\n")
        for node_id, block_id in bc_rows:
            f.write(f"{node_id} {block_id} {value:.16g}\n")


# -----------------------------------------------------------------------------
# BC helpers
# -----------------------------------------------------------------------------


def parse_bc(spec: str, block_length: int) -> Tuple[str, List[int]]:
    if ":" not in spec:
        raise ValueError(
            f"Invalid --bc '{spec}'. Example: --bc ymin:1,2,3"
        )

    face, dof_text = spec.split(":", 1)
    face = face.strip().lower()

    allowed = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}
    if face not in allowed:
        raise ValueError(
            f"Unknown face '{face}'. Choose from {sorted(allowed)}"
        )

    dofs = [int(s.strip()) for s in dof_text.split(",") if s.strip()]
    if not dofs:
        raise ValueError(f"No DOFs specified in '{spec}'")

    for dof in dofs:
        if dof < 1 or dof > block_length:
            raise ValueError(
                f"DOF {dof} is outside 1..{block_length} in '{spec}'"
            )

    return face, sorted(set(dofs))


def local_face_nodes(nodes: np.ndarray, face: str, tol: float) -> np.ndarray:
    if len(nodes) == 0:
        return np.empty(0, dtype=np.int64)

    axis = {"x": 0, "y": 1, "z": 2}[face[0]]
    use_min = face.endswith("min")
    target = nodes[:, axis].min() if use_min else nodes[:, axis].max()

    return np.where(np.abs(nodes[:, axis] - target) <= tol)[0]


# -----------------------------------------------------------------------------
# Interface welding helpers
# -----------------------------------------------------------------------------


def extract_surface_node_ids(
    surface_elems: np.ndarray,
    base: int,
    nnode: int,
    surface_name: str,
) -> np.ndarray:
    """
    Extract unique zero-based local node IDs from a surface connectivity file.
    """

    if surface_elems.size == 0:
        raise ValueError(f"{surface_name}: surface connectivity is empty")

    surface0 = surface_elems - base

    mn = int(surface0.min())
    mx = int(surface0.max())

    if mn < 0 or mx >= nnode:
        raise ValueError(
            f"{surface_name}: surface node IDs are incompatible with the volume mesh. "
            f"zero-based range=[{mn}, {mx}], volume nnode={nnode}, base={base}"
        )

    return np.unique(surface0.reshape(-1)).astype(np.int64)


def _lexicographic_order(coords: np.ndarray) -> np.ndarray:
    """Return deterministic x,y,z lexicographic order for N x 3 coordinates."""

    return np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))


def build_interface_pairing_from_surfaces(
    nodes: np.ndarray,
    left_surface_elems: np.ndarray,
    right_surface_elems: np.ndarray,
    base: int,
    axis_id: int,
    pitch: float,
    tol: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[int, int], float]:
    """
    Build a conformal node pairing using explicitly supplied surface
    connectivity files.

    For identical repeated components:

        copy k right surface
            ==
        copy k+1 left surface

    Since copy k+1 is translated by +pitch along the repeat axis, reference
    component coordinates must satisfy:

        xyz(right_local) ~= xyz(left_local) + pitch * e_axis

    Returns
    -------
    left_ids
        Unique zero-based local node IDs on the left interface surface.
    right_ids
        Unique zero-based local node IDs on the right interface surface.
    right_to_left
        Mapping right local node ID -> corresponding left local node ID.
    max_pair_error
        Maximum absolute XYZ mismatch after shifting the left surface by pitch.
    """

    nnode = len(nodes)

    left_ids = extract_surface_node_ids(
        left_surface_elems,
        base,
        nnode,
        "left interface",
    )
    right_ids = extract_surface_node_ids(
        right_surface_elems,
        base,
        nnode,
        "right interface",
    )

    if len(left_ids) != len(right_ids):
        raise ValueError(
            "Cannot create a conformal continuous interface: "
            f"left-interface unique nodes={len(left_ids)}, "
            f"right-interface unique nodes={len(right_ids)}. "
            "The two supplied interface surfaces must have the same nodal discretization."
        )

    shift = np.zeros(3, dtype=float)
    shift[axis_id] = pitch

    # Compare copy-0 right coordinates with copy-1 left coordinates.
    left_shifted_xyz = nodes[left_ids] + shift
    right_xyz = nodes[right_ids]

    left_order = _lexicographic_order(left_shifted_xyz)
    right_order = _lexicographic_order(right_xyz)

    left_sorted_ids = left_ids[left_order]
    right_sorted_ids = right_ids[right_order]
    left_sorted_xyz = left_shifted_xyz[left_order]
    right_sorted_xyz = right_xyz[right_order]

    errors_xyz = np.abs(left_sorted_xyz - right_sorted_xyz)
    row_errors = np.max(errors_xyz, axis=1)
    max_error = float(row_errors.max()) if len(row_errors) else 0.0

    if np.any(row_errors > tol):
        bad = int(np.argmax(row_errors))
        li = int(left_sorted_ids[bad])
        ri = int(right_sorted_ids[bad])
        raise ValueError(
            "Cannot weld supplied interface surfaces: coordinate pairing failed. "
            f"Maximum XYZ mismatch={max_error:.16g} > weld_tol={tol:.16g}. "
            f"Example: right local node={ri}, xyz={nodes[ri].tolist()}, "
            f"left local node={li}, shifted xyz={left_sorted_xyz[bad].tolist()}. "
            "Check --interface-left/--interface-right and the repetition --pitch."
        )

    # Guard against duplicate coordinates on either interface, because a
    # coordinate-only pairing would then be ambiguous.
    if len(left_sorted_xyz) > 1:
        left_neighbor_error = np.max(
            np.abs(left_sorted_xyz[1:] - left_sorted_xyz[:-1]), axis=1
        )
        right_neighbor_error = np.max(
            np.abs(right_sorted_xyz[1:] - right_sorted_xyz[:-1]), axis=1
        )

        if np.any(left_neighbor_error <= tol):
            raise ValueError(
                "Left interface contains duplicate/indistinguishable node coordinates "
                f"within weld_tol={tol:.16g}; coordinate pairing is ambiguous."
            )

        if np.any(right_neighbor_error <= tol):
            raise ValueError(
                "Right interface contains duplicate/indistinguishable node coordinates "
                f"within weld_tol={tol:.16g}; coordinate pairing is ambiguous."
            )

    right_to_left = {
        int(right_id): int(left_id)
        for right_id, left_id in zip(right_sorted_ids, left_sorted_ids)
    }

    return left_ids, right_ids, right_to_left, max_error


def repeat_disconnected(
    nodes: np.ndarray,
    elems0: np.ndarray,
    copies: int,
    axis_id: int,
    pitch: float,
    start_offset: float,
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray], List[List[int]]]:
    """
    Original disconnected behavior.

    Returns:
      merged_nodes
      merged_elems0
      copy_local_to_global[k][local_node_id] -> global zero-based node ID
      domain_nodes_zero (owner nodes for given_subdomain)
    """

    nnode = len(nodes)
    all_nodes: List[np.ndarray] = []
    all_elems0: List[np.ndarray] = []
    copy_maps: List[np.ndarray] = []
    domain_nodes_zero: List[List[int]] = []

    for k in range(copies):
        shift = np.zeros(3, dtype=float)
        shift[axis_id] = start_offset + k * pitch

        nodes_k = nodes + shift
        node_offset = k * nnode
        local_to_global = np.arange(nnode, dtype=np.int64) + node_offset
        elems_k0 = local_to_global[elems0]

        all_nodes.append(nodes_k)
        all_elems0.append(elems_k0)
        copy_maps.append(local_to_global)
        domain_nodes_zero.append(local_to_global.tolist())

    return (
        np.vstack(all_nodes),
        np.vstack(all_elems0),
        copy_maps,
        domain_nodes_zero,
    )


def repeat_continuous(
    nodes: np.ndarray,
    elems0: np.ndarray,
    copies: int,
    axis_id: int,
    pitch: float,
    start_offset: float,
    left_face_ids: np.ndarray,
    right_face_ids: np.ndarray,
    right_to_left: Dict[int, int],
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray], List[List[int]]]:
    """
    Repeat the component and weld each copy's explicit left interface to
    the previous copy's explicit right interface.

    No bridging elements are created. Instead, the element connectivity of
    the next copy is rewritten so its left-interface nodes reuse the previous
    copy's right-interface global node IDs.

    Shared interface node ownership for given_subdomain:
      - shared nodes remain owned by the previous/lower-index domain;
      - only newly created nodes are owned by the current domain.
    """

    nnode = len(nodes)

    merged_nodes_list: List[np.ndarray] = []
    all_elems0: List[np.ndarray] = []
    copy_maps: List[np.ndarray] = []
    domain_nodes_zero: List[List[int]] = [[] for _ in range(copies)]

    left_face_set = set(int(v) for v in left_face_ids)

    previous_map: np.ndarray | None = None

    for k in range(copies):
        shift = np.zeros(3, dtype=float)
        shift[axis_id] = start_offset + k * pitch
        nodes_k = nodes + shift

        local_to_global = np.empty(nnode, dtype=np.int64)

        if k == 0:
            for local_id in range(nnode):
                gid = len(merged_nodes_list)
                merged_nodes_list.append(nodes_k[local_id].copy())
                local_to_global[local_id] = gid
                domain_nodes_zero[k].append(gid)
        else:
            assert previous_map is not None

            # Reuse previous copy's right-interface global IDs for this
            # copy's left-interface nodes.
            for prev_right_local, this_left_local in right_to_left.items():
                local_to_global[this_left_local] = previous_map[prev_right_local]

            # All remaining local nodes are genuinely new nodes and are owned
            # by the current domain.
            for local_id in range(nnode):
                if local_id in left_face_set:
                    continue

                gid = len(merged_nodes_list)
                merged_nodes_list.append(nodes_k[local_id].copy())
                local_to_global[local_id] = gid
                domain_nodes_zero[k].append(gid)

        elems_k0 = local_to_global[elems0]

        all_elems0.append(elems_k0)
        copy_maps.append(local_to_global)
        previous_map = local_to_global

    merged_nodes = np.asarray(merged_nodes_list, dtype=float)
    merged_elems0 = np.vstack(all_elems0)

    return merged_nodes, merged_elems0, copy_maps, domain_nodes_zero


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Repeat one node.dat/elem.dat component. Optionally weld adjacent "
            "copy interfaces to create one continuous conformal mesh."
        )
    )

    parser.add_argument("node_dat", type=Path)
    parser.add_argument("elem_dat", type=Path)

    parser.add_argument(
        "--copies",
        type=int,
        required=True,
        help="Number of component copies / GEDATSU subdomains",
    )

    parser.add_argument(
        "--axis",
        choices=("x", "y", "z"),
        default="x",
        help="Direction in which components are arranged",
    )

    spacing = parser.add_mutually_exclusive_group(required=True)
    spacing.add_argument(
        "--gap",
        type=float,
        help=(
            "Empty distance between neighboring component bounding boxes. "
            "pitch = component length + gap"
        ),
    )
    spacing.add_argument(
        "--pitch",
        type=float,
        help="Constant origin-to-origin translation between copies",
    )

    parser.add_argument(
        "--start-offset",
        type=float,
        default=0.0,
        help="Translation of the first copy along the selected axis",
    )

    parser.add_argument(
        "--input-base",
        choices=("auto", "0", "1"),
        default="auto",
        help="Connectivity/node ID base. Default: auto",
    )

    parser.add_argument(
        "--continuous",
        "--weld-interfaces",
        dest="continuous",
        action="store_true",
        help=(
            "Weld each copy's explicit left interface surface to the previous "
            "copy's explicit right interface surface. Requires "
            "--interface-left and --interface-right."
        ),
    )

    parser.add_argument(
        "--weld-tol",
        type=float,
        default=1.0e-10,
        help=(
            "Tolerance used to detect and validate conformal interface nodes. "
            "Default: 1e-10"
        ),
    )

    parser.add_argument(
        "--interface-left",
        "--interface-min",
        dest="interface_left",
        type=Path,
        default=None,
        help=(
            "Surface connectivity file for the LEFT interface of the original "
            "component. Required with --continuous."
        ),
    )

    parser.add_argument(
        "--interface-right",
        "--interface-max",
        dest="interface_right",
        type=Path,
        default=None,
        help=(
            "Surface connectivity file for the RIGHT interface of the original "
            "component. Required with --continuous."
        ),
    )

    # BC
    parser.add_argument(
        "--block-length",
        type=int,
        default=3,
        help="Block length b in D_bc.dat. Default: 3",
    )

    parser.add_argument(
        "--bc",
        action="append",
        default=[],
        metavar="FACE:DOFS",
        help=(
            "Boundary condition. Example: --bc ymin:1,2,3. "
            "Can be repeated. In continuous mode, internal repeat-axis "
            "interfaces are automatically excluded."
        ),
    )

    parser.add_argument(
        "--bc-tol",
        type=float,
        default=1.0e-10,
        help="Tolerance for detecting local bounding-box face nodes",
    )

    parser.add_argument(
        "--bc-value",
        type=float,
        default=0.0,
        help="Value written to all D_bc.dat rows. Default: 0.0",
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("merged_mesh"),
    )

    args = parser.parse_args()

    if args.copies < 1:
        raise ValueError("--copies must be >= 1")
    if args.block_length < 1:
        raise ValueError("--block-length must be >= 1")
    if args.bc_tol < 0:
        raise ValueError("--bc-tol must be >= 0")
    if args.weld_tol <= 0:
        raise ValueError("--weld-tol must be > 0")

    if args.continuous:
        if args.interface_left is None or args.interface_right is None:
            raise ValueError(
                "--continuous requires both --interface-left and --interface-right."
            )

    nodes = read_nodes(args.node_dat)
    elems, nen = read_elements(args.elem_dat)

    nnode = len(nodes)
    nelem = len(elems)

    if nnode == 0:
        raise ValueError("node.dat contains no nodes")
    if nelem == 0:
        raise ValueError("elem.dat contains no elements")

    base = detect_base(elems, nnode, args.input_base)
    elems0 = elems - base

    print(f"Detected connectivity base: {base}-based")

    axis_id = {"x": 0, "y": 1, "z": 2}[args.axis]
    axis_name = args.axis
    min_axis_face = f"{axis_name}min"
    max_axis_face = f"{axis_name}max"

    xyz_min = nodes.min(axis=0)
    xyz_max = nodes.max(axis=0)
    component_length = float(xyz_max[axis_id] - xyz_min[axis_id])

    if component_length <= 0:
        raise ValueError(
            f"Component has zero/negative length along axis {args.axis}: "
            f"{component_length}"
        )

    if args.pitch is not None:
        pitch = float(args.pitch)
        gap = pitch - component_length
    else:
        gap = float(args.gap)
        pitch = component_length + gap

    if pitch <= 0:
        raise ValueError(
            f"pitch must be positive: length={component_length}, pitch={pitch}"
        )

    if gap < 0 and not args.continuous:
        print(
            f"WARNING: gap={gap} < 0, so neighboring bounding boxes overlap."
        )

    # BC specs and original local face node lists.
    bc_specs = [parse_bc(spec, args.block_length) for spec in args.bc]
    face_cache: Dict[str, np.ndarray] = {}
    for face, _ in bc_specs:
        if face not in face_cache:
            face_cache[face] = local_face_nodes(nodes, face, args.bc_tol)

    interface_node_count = 0
    max_pair_error = 0.0

    if args.continuous:
        left_surface_elems, _ = read_elements(args.interface_left)
        right_surface_elems, _ = read_elements(args.interface_right)

        (
            left_face_ids,
            right_face_ids,
            right_to_left,
            max_pair_error,
        ) = build_interface_pairing_from_surfaces(
            nodes,
            left_surface_elems,
            right_surface_elems,
            base,
            axis_id,
            pitch,
            args.weld_tol,
        )
        interface_node_count = len(left_face_ids)

        (
            merged_nodes,
            merged_elems0,
            copy_maps,
            domain_nodes_zero,
        ) = repeat_continuous(
            nodes,
            elems0,
            args.copies,
            axis_id,
            pitch,
            args.start_offset,
            left_face_ids,
            right_face_ids,
            right_to_left,
        )
    else:
        (
            merged_nodes,
            merged_elems0,
            copy_maps,
            domain_nodes_zero,
        ) = repeat_disconnected(
            nodes,
            elems0,
            args.copies,
            axis_id,
            pitch,
            args.start_offset,
        )

    # Validate connectivity after welding/repetition.
    if merged_elems0.min() < 0 or merged_elems0.max() >= len(merged_nodes):
        raise RuntimeError(
            "Generated element connectivity contains an invalid node ID: "
            f"range=[{int(merged_elems0.min())}, {int(merged_elems0.max())}], "
            f"nnode={len(merged_nodes)}"
        )

    # BC construction using local->global maps.
    bc_set = set()

    for face, dofs in bc_specs:
        local_ids = face_cache[face]
        face_axis = {"x": 0, "y": 1, "z": 2}[face[0]]

        for k in range(args.copies):
            if args.continuous and face_axis == axis_id:
                # Internal welded interfaces are not external BC surfaces.
                if face == min_axis_face and k != 0:
                    continue
                if face == max_axis_face and k != args.copies - 1:
                    continue

            local_to_global = copy_maps[k]

            for local_id in local_ids:
                global_zero_id = int(local_to_global[int(local_id)])
                output_node_id = global_zero_id + base

                for dof in dofs:
                    bc_set.add((output_node_id, dof))

    bc_rows = sorted(bc_set, key=lambda x: (x[0], x[1]))

    # Validate node ownership in given_subdomain.
    ownership = np.full(len(merged_nodes), -1, dtype=np.int64)
    for domain_id, ids in enumerate(domain_nodes_zero):
        for gid in ids:
            if ownership[gid] != -1:
                raise RuntimeError(
                    f"Node {gid} is owned by more than one subdomain: "
                    f"{int(ownership[gid])} and {domain_id}"
                )
            ownership[gid] = domain_id

    missing = np.where(ownership < 0)[0]
    if len(missing):
        raise RuntimeError(
            f"{len(missing)} generated nodes have no given_subdomain owner. "
            f"First missing node={int(missing[0])}"
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)

    node_out = args.out_dir / "node.dat"
    elem_out = args.out_dir / "elem.dat"
    subdomain_out = args.out_dir / "given_subdomain.dat"
    bc_out = args.out_dir / "D_bc.dat"
    report_out = args.out_dir / "repeat_report.txt"

    write_nodes(node_out, merged_nodes)
    write_elements(elem_out, merged_elems0, nen, base)
    write_given_subdomain(subdomain_out, domain_nodes_zero, base)

    if bc_specs:
        write_bc(
            bc_out,
            bc_rows,
            args.block_length,
            value=args.bc_value,
        )

    expected_disconnected_nodes = nnode * args.copies
    welded_duplicates = expected_disconnected_nodes - len(merged_nodes)

    with report_out.open("w", encoding="utf-8") as f:
        f.write("Repeated component mesh report\n")
        f.write("==============================\n")
        f.write(f"input node file       : {args.node_dat}\n")
        f.write(f"input elem file       : {args.elem_dat}\n")
        f.write(f"input nodes           : {nnode}\n")
        f.write(f"input elements        : {nelem}\n")
        f.write(f"nodes/element         : {nen}\n")
        f.write(f"index base            : {base}\n")
        f.write(f"copies                : {args.copies}\n")
        f.write(f"GEDATSU domains       : {args.copies}\n")
        f.write(f"axis                  : {args.axis}\n")
        f.write(f"component bbox min    : {xyz_min.tolist()}\n")
        f.write(f"component bbox max    : {xyz_max.tolist()}\n")
        f.write(f"component length      : {component_length:.16g}\n")
        f.write(f"gap                   : {gap:.16g}\n")
        f.write(f"pitch                 : {pitch:.16g}\n")
        f.write(f"start offset          : {args.start_offset:.16g}\n")
        f.write(f"continuous            : {args.continuous}\n")
        if args.continuous:
            f.write(f"interface left file   : {args.interface_left}\n")
            f.write(f"interface right file  : {args.interface_right}\n")
        f.write(f"node welding          : {'enabled' if args.continuous else 'disabled'}\n")
        f.write(f"interface unique nodes : {interface_node_count}\n")
        f.write(f"interface pair error   : {max_pair_error:.16g}\n")
        f.write(f"welded duplicate nodes: {welded_duplicates}\n")
        f.write(f"output nodes          : {len(merged_nodes)}\n")
        f.write(f"output elements       : {len(merged_elems0)}\n")
        f.write("subdomain strategy    : one copy per domain, shared nodes owned by previous domain\n")
        f.write(f"given subdomain file  : {subdomain_out}\n")

        f.write("\nowned nodes / domain\n")
        for k, ids in enumerate(domain_nodes_zero):
            f.write(f"  domain {k}: {len(ids)}\n")

        if bc_specs:
            f.write(f"\nBC block length       : {args.block_length}\n")
            f.write(f"BC tolerance          : {args.bc_tol:.16g}\n")
            f.write(f"BC value              : {args.bc_value:.16g}\n")
            f.write(f"BC specs              : {args.bc}\n")
            f.write(f"D_bc rows             : {len(bc_rows)}\n")
            f.write(f"D_bc file             : {bc_out}\n")
            if args.continuous:
                f.write("internal interface BC : excluded on repeat-axis internal faces\n")
        else:
            f.write("\nD_bc                  : not generated (--bc not specified)\n")

        f.write("\nGEDATSU example command\n")
        f.write("-----------------------\n")
        f.write(f"cd {args.out_dir}\n")
        f.write("mkdir -p parted.0\n")
        f.write(
            "gedatsu_simple_mesh_partitioner "
            f"-n {args.copies} "
            "-in node.dat "
            "-ie elem.dat "
            "--given_subdomain given_subdomain.dat "
            "-d parted.0\n"
        )

    print()
    print("Done.")
    print()
    print(f"component nodes      : {nnode}")
    print(f"component elements   : {nelem}")
    print(f"copies               : {args.copies}")
    print(f"GEDATSU domains      : {args.copies}")
    print(f"axis                 : {args.axis}")
    print(f"component length     : {component_length}")
    print(f"gap                  : {gap}")
    print(f"pitch                : {pitch}")
    print(f"continuous           : {args.continuous}")
    print(f"interface unique nodes: {interface_node_count}")
    print(f"interface pair error : {max_pair_error}")
    print(f"welded duplicates    : {welded_duplicates}")
    print(f"output nodes         : {len(merged_nodes)}")
    print(f"output elements      : {len(merged_elems0)}")
    print()
    print(f"node.dat             : {node_out}")
    print(f"elem.dat             : {elem_out}")
    print(f"given_subdomain.dat  : {subdomain_out}")

    if bc_specs:
        print(f"D_bc.dat             : {bc_out}")
        print(f"D_bc rows            : {len(bc_rows)}")
        if args.continuous:
            print("internal interface BC: excluded")
    else:
        print("D_bc.dat             : not generated")

    print(f"report               : {report_out}")


if __name__ == "__main__":
    main()