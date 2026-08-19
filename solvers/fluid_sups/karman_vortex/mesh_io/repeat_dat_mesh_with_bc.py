#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
repeat_dat_mesh_with_bc_symmetric.py

Gmsh の .geo から生成した 1 個の半コンポーネントを基準として、

    x 方向:
        C, Mirror(C), C, Mirror(C), ...

    z 方向:
        同じ x 列を奥側（負の z 方向）へ連続的に複製

した conformal continuous mesh を生成する。

重要
====
1 GEO component = 1 GEDATSU domain

したがって、

    --copies 8
    --depth-copies 2

なら

    x方向       : 8 components
    z方向       : 2 layers
    GEDATSU領域 : 8 * 2 = 16 domains

となる。

Domain numbering
----------------

layer-major:

    domain_id = depth_layer * x_components + x_component

例: --copies 8 --depth-copies 2

    z layer 0:
        domain 0  1  2  3  4  5  6  7

    z layer 1:
        domain 8  9 10 11 12 13 14 15


x-direction structure
---------------------

--copies は x 方向の GEO component 数。

左右端を fin 側に残すため、symmetric mode では偶数を要求する。

    C == M == C == M == C == M ...

C = original GEO component
M = x-direction mirror of the same GEO component

x interface:
    C -> M:
        Gamma03 == Gamma03

    M -> next C:
        (Gamma11 + Gamma12) == (Gamma11 + Gamma12)


depth/z-direction structure
---------------------------

thermal_fin_fig3.geo は

    0 >= z >= -depth

へ extrusion しているため、追加 layer は負の z 方向へ配置する。

手前側 z=max:
    Gamma01 + Gamma07

奥側 z=min:
    Gamma05 + Gamma09

depth layer 0:
    original z range

depth layer 1:
    original range - depth

depth layer 2:
    original range - 2*depth

...

layer 間では

    previous back
        ==
    current front

となるよう node weld する。


Default GEO physical groups
---------------------------

Volume:
    Fluid

x center:
    Gamma11
    Gamma12

x gap:
    Gamma03

z front:
    Gamma01
    Gamma07

z back:
    Gamma05
    Gamma09

BC:
    All_external_boundaries


BC
--

--bc-surface を使用する場合、surface element 単位で内部面を除外する。

x 内部面:
    Gamma03
    Gamma11 + Gamma12（最外端を除く）

z 内部面:
    Gamma01 + Gamma07（最前面を除く）
    Gamma05 + Gamma09（最奥面を除く）

そのため edge/corner node が別の真の外部面に属する場合は、
その node は D_bc.dat に正しく残る。

convdiff のような scalar problem では

    --block-length 1
    --bc-surface-dofs 1

を使用する。


Recommended GEO usage
---------------------

python3 repeat_dat_mesh_with_bc_symmetric.py \
    --geo gmsh/thermal_fin_fig3.geo \
    --physical-group-script mesh_io/save_physical_groups.py \
    --copies 8 \
    --depth-copies 2 \
    --axis x \
    --gap 0.0 \
    --symmetric \
    --block-length 1 \
    --bc-surface-dofs 1 \
    --bc-value 0.0 \
    --out-dir merged_mesh

これで

    x components = 8
    z layers     = 2
    total domains= 16

となる。

GEDATSU:

    cd merged_mesh
    mkdir -p parted.0

    gedatsu_simple_mesh_partitioner \
        -n 16 \
        -in node.dat \
        -ie elem.dat \
        --given_subdomain given_subdomain.dat \
        -d parted.0

    gedatsu_bc_partitioner_R \
        -n 16 \
        -i D_bc.dat \
        -ig node.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
import shutil
import subprocess
import sys

import numpy as np


# =============================================================================
# GEO generation
# =============================================================================


def require_file(path: Path, label: str) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"{label} was not found: {path}"
        )

    return path


def require_executable(executable: str) -> str:
    p = Path(executable)

    if p.is_absolute() or p.parent != Path("."):
        p = p.resolve()

        if not p.is_file():
            raise FileNotFoundError(
                f"Executable was not found: {p}"
            )

        return str(p)

    found = shutil.which(executable)

    if found is None:
        raise FileNotFoundError(
            f"Executable '{executable}' was not found in PATH."
        )

    return found


def physical_surface_file(
    mesh_tmp_dir: Path,
    group_name: str,
) -> Path:
    return (
        mesh_tmp_dir
        / f"{group_name}_quad_connectivity.dat"
    )


def generate_component_from_geo(
    geo_path: Path,
    gmsh_executable: str,
    physical_group_script: Path,
    component_work_dir: Path,
    volume_group: str,
    center_groups: Sequence[str],
    gap_groups: Sequence[str],
    depth_front_groups: Sequence[str],
    depth_back_groups: Sequence[str],
    bc_group: str | None,
    clean_work: bool,
) -> Dict[str, object]:

    geo_path = require_file(
        geo_path,
        "Gmsh .geo file",
    )

    physical_group_script = require_file(
        physical_group_script,
        "physical-group extraction script",
    )

    gmsh_executable = require_executable(
        gmsh_executable
    )

    component_work_dir = (
        component_work_dir.resolve()
    )

    component_work_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    mesh_tmp_dir = (
        component_work_dir
        / "mesh_tmp"
    )

    if clean_work and mesh_tmp_dir.exists():
        shutil.rmtree(
            mesh_tmp_dir
        )

    msh_out = (
        component_work_dir
        / "component.msh"
    )

    if clean_work and msh_out.exists():
        msh_out.unlink()

    print()
    print("Generating component from GEO...")
    print(f"geo                 : {geo_path}")
    print(f"component work dir  : {component_work_dir}")
    print(f"gmsh                : {gmsh_executable}")
    print()

    subprocess.run(
        [
            gmsh_executable,
            "-3",
            "-format",
            "msh2",
            "-0",
            str(geo_path),
            "-o",
            str(msh_out),
        ],
        cwd=component_work_dir,
        check=True,
    )

    if not msh_out.is_file():
        raise RuntimeError(
            f"Gmsh did not create {msh_out}"
        )

    subprocess.run(
        [
            sys.executable,
            str(physical_group_script),
            str(msh_out),
        ],
        cwd=component_work_dir,
        check=True,
    )

    node_dat = require_file(
        mesh_tmp_dir
        / f"{volume_group}_node_coordinates.dat",
        f"{volume_group} node coordinate file",
    )

    elem_dat = require_file(
        mesh_tmp_dir
        / f"{volume_group}_hexahedron_connectivity.dat",
        f"{volume_group} hexahedron connectivity file",
    )

    def group_files(
        names: Sequence[str],
        label: str,
    ) -> List[Path]:
        return [
            require_file(
                physical_surface_file(
                    mesh_tmp_dir,
                    name,
                ),
                f"{label} physical surface '{name}'",
            )
            for name in names
        ]

    center_files = group_files(
        center_groups,
        "x-center",
    )

    gap_files = group_files(
        gap_groups,
        "x-gap",
    )

    depth_front_files = group_files(
        depth_front_groups,
        "depth-front",
    )

    depth_back_files = group_files(
        depth_back_groups,
        "depth-back",
    )

    if bc_group is None:
        bc_files: List[Path] = []
    else:
        bc_files = [
            require_file(
                physical_surface_file(
                    mesh_tmp_dir,
                    bc_group,
                ),
                f"BC physical surface '{bc_group}'",
            )
        ]

    print()
    print("GEO component generation completed.")
    print(f"node input          : {node_dat}")
    print(f"elem input          : {elem_dat}")
    print(f"x center surfaces   : {[str(p) for p in center_files]}")
    print(f"x gap surfaces      : {[str(p) for p in gap_files]}")
    print(f"z front surfaces    : {[str(p) for p in depth_front_files]}")
    print(f"z back surfaces     : {[str(p) for p in depth_back_files]}")

    if bc_files:
        print(f"BC surface          : {bc_files[0]}")

    return {
        "geo": geo_path,
        "msh": msh_out,
        "mesh_tmp": mesh_tmp_dir,
        "node_dat": node_dat,
        "elem_dat": elem_dat,
        "center_interface": center_files,
        "gap_interface": gap_files,
        "depth_front_interface": depth_front_files,
        "depth_back_interface": depth_back_files,
        "bc_surface": bc_files,
        "component_work_dir": component_work_dir,
    }


# =============================================================================
# Mesh I/O
# =============================================================================


def read_nodes(path: Path) -> np.ndarray:
    with path.open(
        "r",
        encoding="utf-8",
    ) as f:
        header = (
            f.readline()
            .strip()
            .split()
        )

    if not header:
        raise ValueError(
            f"{path}: empty file"
        )

    nnode = int(
        header[0]
    )

    nodes = np.loadtxt(
        path,
        skiprows=1,
        dtype=float,
    )

    if nodes.size == 0:
        nodes = np.empty(
            (0, 3),
            dtype=float,
        )
    elif nodes.ndim == 1:
        nodes = nodes.reshape(
            1,
            -1,
        )

    if nodes.shape != (
        nnode,
        3,
    ):
        raise ValueError(
            f"{path}: header={nnode}, "
            f"actual shape={nodes.shape}"
        )

    return nodes


def read_elements(
    path: Path,
) -> Tuple[np.ndarray, int]:
    with path.open(
        "r",
        encoding="utf-8",
    ) as f:
        header = (
            f.readline()
            .strip()
            .split()
        )

    if len(header) < 2:
        raise ValueError(
            f"{path}: first line must be 'Nelem NEN'"
        )

    nelem = int(
        header[0]
    )

    nen = int(
        header[1]
    )

    elems = np.loadtxt(
        path,
        skiprows=1,
        dtype=np.int64,
    )

    if elems.size == 0:
        elems = np.empty(
            (0, nen),
            dtype=np.int64,
        )
    elif elems.ndim == 1:
        elems = elems.reshape(
            1,
            -1,
        )

    if elems.shape != (
        nelem,
        nen,
    ):
        raise ValueError(
            f"{path}: header=({nelem},{nen}), "
            f"actual={elems.shape}"
        )

    return (
        elems,
        nen,
    )


def detect_base(
    elems: np.ndarray,
    nnode: int,
    requested: str,
) -> int:

    if requested in (
        "0",
        "1",
    ):
        base = int(
            requested
        )

        if elems.size:
            mn = int(
                elems.min()
            )

            mx = int(
                elems.max()
            )

            if (
                base == 0
                and (
                    mn < 0
                    or mx >= nnode
                )
            ):
                raise ValueError(
                    f"Invalid --input-base 0: "
                    f"range=[{mn},{mx}], nnode={nnode}"
                )

            if (
                base == 1
                and (
                    mn < 1
                    or mx > nnode
                )
            ):
                raise ValueError(
                    f"Invalid --input-base 1: "
                    f"range=[{mn},{mx}], nnode={nnode}"
                )

        return base

    if elems.size == 0:
        return 0

    mn = int(
        elems.min()
    )

    mx = int(
        elems.max()
    )

    if mn == 0:
        if mx >= nnode:
            raise ValueError(
                "Invalid 0-based connectivity"
            )

        return 0

    if mx == nnode:
        if mn < 1:
            raise ValueError(
                "Invalid 1-based connectivity"
            )

        return 1

    if (
        0 < mn
        and mx < nnode
    ):
        raise ValueError(
            "Index base is ambiguous. "
            "Specify --input-base 0 or 1."
        )

    raise ValueError(
        f"Invalid connectivity range "
        f"[{mn},{mx}] for nnode={nnode}"
    )


def write_nodes(
    path: Path,
    nodes: np.ndarray,
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            f"{len(nodes)}\n"
        )

        for x, y, z in nodes:
            f.write(
                f"{x:.16g} "
                f"{y:.16g} "
                f"{z:.16g}\n"
            )


def write_elements(
    path: Path,
    elems0: np.ndarray,
    nen: int,
    base: int,
) -> None:

    elems = (
        elems0
        + base
    )

    with path.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            f"{len(elems)} {nen}\n"
        )

        for row in elems:
            f.write(
                " ".join(
                    str(int(v))
                    for v in row
                )
                + "\n"
            )


def write_given_subdomain(
    path: Path,
    domain_nodes_zero: Sequence[
        Sequence[int]
    ],
    base: int,
) -> None:

    with path.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            f"{len(domain_nodes_zero)}\n"
        )

        for domain_id, ids in enumerate(
            domain_nodes_zero
        ):
            ids0 = sorted(
                set(
                    int(v)
                    for v in ids
                )
            )

            if not ids0:
                raise ValueError(
                    f"Domain {domain_id} "
                    "owns zero nodes."
                )

            point_id = (
                domain_id
                + base
            )

            ids_out = [
                gid
                + base
                for gid in ids0
            ]

            f.write(
                f"{point_id} "
                f"{len(ids_out)}"
            )

            if ids_out:
                f.write(
                    " "
                    + " ".join(
                        str(v)
                        for v in ids_out
                    )
                )

            f.write(
                "\n"
            )


def write_bc(
    path: Path,
    rows: Sequence[
        Tuple[int, int]
    ],
    block_length: int,
    value: float,
) -> None:

    with path.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            f"{len(rows)} "
            f"{block_length}\n"
        )

        for node_id, dof in rows:
            f.write(
                f"{node_id} "
                f"{dof} "
                f"{value:.16g}\n"
            )


# =============================================================================
# Surface helpers
# =============================================================================


def read_surface_faces_union(
    paths: Sequence[Path],
    base: int,
    nnode: int,
    label: str,
) -> np.ndarray:

    if not paths:
        return np.empty(
            (0, 0),
            dtype=np.int64,
        )

    unique_rows: Dict[
        Tuple[int, ...],
        np.ndarray,
    ] = {}

    nen_ref = None

    for path in paths:
        surf, nen = read_elements(
            path
        )

        if nen_ref is None:
            nen_ref = nen
        elif nen != nen_ref:
            raise ValueError(
                f"{label}: mixed surface NEN "
                f"{nen_ref} and {nen}"
            )

        surf0 = (
            surf
            - base
        )

        if surf0.size:
            mn = int(
                surf0.min()
            )

            mx = int(
                surf0.max()
            )

            if (
                mn < 0
                or mx >= nnode
            ):
                raise ValueError(
                    f"{label}: {path}: "
                    f"surface node range "
                    f"[{mn},{mx}] is invalid "
                    f"for nnode={nnode}"
                )

        for row in surf0:
            key = tuple(
                sorted(
                    int(v)
                    for v in row
                )
            )

            if key not in unique_rows:
                unique_rows[
                    key
                ] = (
                    np.asarray(
                        row,
                        dtype=np.int64,
                    )
                    .copy()
                )

    if not unique_rows:
        return np.empty(
            (
                0,
                nen_ref
                if nen_ref is not None
                else 0,
            ),
            dtype=np.int64,
        )

    return np.vstack(
        list(
            unique_rows.values()
        )
    )


def surface_node_ids(
    faces0: np.ndarray,
) -> np.ndarray:

    if faces0.size == 0:
        return np.empty(
            0,
            dtype=np.int64,
        )

    return np.unique(
        faces0.reshape(
            -1
        )
    ).astype(
        np.int64
    )


def face_key_set(
    faces0: np.ndarray,
) -> set:
    return {
        tuple(
            sorted(
                int(v)
                for v in row
            )
        )
        for row in faces0
    }


def bbox_face_nodes(
    nodes: np.ndarray,
    axis_id: int,
    use_min: bool,
    tol: float,
) -> np.ndarray:

    target = (
        nodes[:, axis_id].min()
        if use_min
        else nodes[:, axis_id].max()
    )

    return np.where(
        np.abs(
            nodes[:, axis_id]
            - target
        )
        <= tol
    )[0].astype(
        np.int64
    )


def compare_node_sets(
    actual: np.ndarray,
    expected: np.ndarray,
    actual_name: str,
    expected_name: str,
) -> None:

    a = set(
        int(v)
        for v in actual
    )

    e = set(
        int(v)
        for v in expected
    )

    if a == e:
        return

    missing = sorted(
        e - a
    )

    extra = sorted(
        a - e
    )

    raise ValueError(
        f"{actual_name} does not match "
        f"{expected_name}. "
        f"actual={len(a)}, "
        f"expected={len(e)}, "
        f"missing(first10)={missing[:10]}, "
        f"extra(first10)={extra[:10]}"
    )


def lexicographic_order(
    coords: np.ndarray,
) -> np.ndarray:
    return np.lexsort(
        (
            coords[:, 2],
            coords[:, 1],
            coords[:, 0],
        )
    )


# =============================================================================
# Hex8 mirror
# =============================================================================


HEX8_REF_COORDS = np.array(
    [
        [-1.0, -1.0, -1.0],
        [ 1.0, -1.0, -1.0],
        [ 1.0,  1.0, -1.0],
        [-1.0,  1.0, -1.0],
        [-1.0, -1.0,  1.0],
        [ 1.0, -1.0,  1.0],
        [ 1.0,  1.0,  1.0],
        [-1.0,  1.0,  1.0],
    ],
    dtype=float,
)

HEX8_DNDXI_CENTER = (
    HEX8_REF_COORDS
    / 8.0
)

# x-like reflection in the local Hex8 ordering.
HEX8_MIRROR_PERM = np.array(
    [
        1, 0, 3, 2,
        5, 4, 7, 6,
    ],
    dtype=np.int64,
)


def hex8_center_jacobian_det(
    nodes: np.ndarray,
    elems0: np.ndarray,
) -> np.ndarray:

    coords = nodes[
        elems0
    ]

    jac = np.einsum(
        "eia,ib->eab",
        coords,
        HEX8_DNDXI_CENTER,
    )

    return np.linalg.det(
        jac
    )


def mirrored_elements(
    nodes: np.ndarray,
    elems0: np.ndarray,
    axis_id: int,
) -> Tuple[np.ndarray, float]:

    if elems0.shape[1] != 8:
        raise ValueError(
            "--symmetric currently supports "
            "8-node hexahedra only."
        )

    axis_min = float(
        nodes[
            :,
            axis_id,
        ].min()
    )

    mirrored_nodes = (
        nodes.copy()
    )

    mirrored_nodes[
        :,
        axis_id,
    ] = (
        2.0
        * axis_min
        - nodes[
            :,
            axis_id,
        ]
    )

    mirrored_elems0 = (
        elems0[
            :,
            HEX8_MIRROR_PERM,
        ]
    )

    det0 = hex8_center_jacobian_det(
        nodes,
        elems0,
    )

    detm = hex8_center_jacobian_det(
        mirrored_nodes,
        mirrored_elems0,
    )

    tiny = 1.0e-14

    if np.any(
        np.abs(
            det0
        )
        <= tiny
    ):
        idx = int(
            np.argmin(
                np.abs(
                    det0
                )
            )
        )

        raise ValueError(
            "Original Hex8 contains "
            "near-zero center Jacobian. "
            f"element={idx}, "
            f"det={det0[idx]}"
        )

    if np.any(
        np.abs(
            detm
        )
        <= tiny
    ):
        idx = int(
            np.argmin(
                np.abs(
                    detm
                )
            )
        )

        raise ValueError(
            "Mirrored Hex8 contains "
            "near-zero center Jacobian. "
            f"element={idx}, "
            f"det={detm[idx]}"
        )

    if np.any(
        np.sign(
            det0
        )
        != np.sign(
            detm
        )
    ):
        idx = int(
            np.where(
                np.sign(det0)
                != np.sign(detm)
            )[0][0]
        )

        raise ValueError(
            "Mirrored Hex8 orientation "
            "does not match original. "
            f"element={idx}, "
            f"original={det0[idx]}, "
            f"mirror={detm[idx]}"
        )

    orientation_error = float(
        np.max(
            np.abs(
                detm
                / det0
                - 1.0
            )
        )
    )

    return (
        mirrored_elems0,
        orientation_error,
    )


# =============================================================================
# Geometry transforms
# =============================================================================


def transform_original_component(
    nodes: np.ndarray,
    axis_id: int,
    axis_min: float,
    left_center_position: float,
    depth_axis_id: int,
    depth_shift: float,
) -> np.ndarray:

    out = nodes.copy()

    u = (
        nodes[:, axis_id]
        - axis_min
    )

    out[:, axis_id] = (
        left_center_position
        + u
    )

    out[:, depth_axis_id] = (
        nodes[:, depth_axis_id]
        + depth_shift
    )

    return out


def transform_mirror_component(
    nodes: np.ndarray,
    axis_id: int,
    axis_min: float,
    right_center_position: float,
    depth_axis_id: int,
    depth_shift: float,
) -> np.ndarray:

    out = nodes.copy()

    u = (
        nodes[:, axis_id]
        - axis_min
    )

    out[:, axis_id] = (
        right_center_position
        - u
    )

    out[:, depth_axis_id] = (
        nodes[:, depth_axis_id]
        + depth_shift
    )

    return out


# =============================================================================
# Interface validation / pairing
# =============================================================================


def validate_x_interfaces(
    nodes: np.ndarray,
    axis_id: int,
    center_ids: np.ndarray,
    gap_ids: np.ndarray,
    tol: float,
) -> None:

    compare_node_sets(
        center_ids,
        bbox_face_nodes(
            nodes,
            axis_id,
            True,
            tol,
        ),
        "x-center surface union",
        "repeat-axis minimum face",
    )

    compare_node_sets(
        gap_ids,
        bbox_face_nodes(
            nodes,
            axis_id,
            False,
            tol,
        ),
        "x-gap surface union",
        "repeat-axis maximum face",
    )


def build_depth_pairing(
    nodes: np.ndarray,
    depth_axis_id: int,
    front_ids: np.ndarray,
    back_ids: np.ndarray,
    tol: float,
) -> Tuple[
    Dict[int, int],
    float,
    float,
]:

    compare_node_sets(
        front_ids,
        bbox_face_nodes(
            nodes,
            depth_axis_id,
            False,
            tol,
        ),
        "depth-front surface union",
        "z/depth maximum face",
    )

    compare_node_sets(
        back_ids,
        bbox_face_nodes(
            nodes,
            depth_axis_id,
            True,
            tol,
        ),
        "depth-back surface union",
        "z/depth minimum face",
    )

    if len(
        front_ids
    ) != len(
        back_ids
    ):
        raise ValueError(
            "Depth interface is not conformal: "
            f"front nodes={len(front_ids)}, "
            f"back nodes={len(back_ids)}"
        )

    zmin = float(
        nodes[
            :,
            depth_axis_id,
        ].min()
    )

    zmax = float(
        nodes[
            :,
            depth_axis_id,
        ].max()
    )

    depth_length = (
        zmax
        - zmin
    )

    if depth_length <= 0:
        raise ValueError(
            f"Invalid depth length: {depth_length}"
        )

    shift = np.zeros(
        3,
        dtype=float,
    )

    # A new layer is placed deeper in the negative z direction.
    # Its front face is therefore shifted by -depth_length
    # and must coincide with the previous layer's back face.
    shift[
        depth_axis_id
    ] = (
        -depth_length
    )

    shifted_front = (
        nodes[
            front_ids
        ]
        + shift
    )

    back_xyz = (
        nodes[
            back_ids
        ]
    )

    fo = lexicographic_order(
        shifted_front
    )

    bo = lexicographic_order(
        back_xyz
    )

    error = np.max(
        np.abs(
            shifted_front[
                fo
            ]
            - back_xyz[
                bo
            ]
        ),
        axis=1,
    )

    max_error = (
        float(
            error.max()
        )
        if len(error)
        else 0.0
    )

    if np.any(
        error > tol
    ):
        idx = int(
            np.where(
                error > tol
            )[0][0]
        )

        raise ValueError(
            "Depth front/back coordinate "
            "pairing failed. "
            f"max error={max_error:.16g}, "
            f"tol={tol:.16g}, "
            f"front local id="
            f"{int(front_ids[fo[idx]])}, "
            f"back local id="
            f"{int(back_ids[bo[idx]])}"
        )

    # current layer front local ID
    #       ->
    # previous layer back local ID
    front_to_back = {
        int(
            front_ids[
                fi
            ]
        ):
        int(
            back_ids[
                bi
            ]
        )
        for fi, bi in zip(
            fo,
            bo,
        )
    }

    return (
        front_to_back,
        max_error,
        depth_length,
    )


# =============================================================================
# 2D component-grid construction: x sequence x z layers
# =============================================================================


def assign_shared_gid(
    local_map: np.ndarray,
    lid: int,
    gid: int,
    label: str,
) -> None:

    current = int(
        local_map[
            lid
        ]
    )

    if current == -1:
        local_map[
            lid
        ] = gid
        return

    if current != gid:
        raise RuntimeError(
            "Inconsistent shared-node mapping "
            f"at {label}: "
            f"local node={lid}, "
            f"existing global={current}, "
            f"new global={gid}"
        )


def build_component_grid(
    nodes: np.ndarray,
    elems0: np.ndarray,
    x_components: int,
    depth_copies: int,
    axis_id: int,
    depth_axis_id: int,
    start_offset: float,
    center_ids: np.ndarray,
    gap_ids: np.ndarray,
    depth_front_to_back: Dict[int, int],
    depth_length: float,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    List[List[np.ndarray]],
    List[List[int]],
    float,
]:

    if x_components < 2:
        raise ValueError(
            "--copies must be >= 2 "
            "in symmetric mode"
        )

    if x_components % 2 != 0:
        raise ValueError(
            "--copies must be even in "
            "symmetric mode so that both "
            "global x ends are fin/center faces."
        )

    if depth_copies < 1:
        raise ValueError(
            "--depth-copies must be >= 1"
        )

    if (
        depth_copies > 1
        and axis_id == depth_axis_id
    ):
        raise ValueError(
            "The symmetric repeat axis and "
            "depth axis cannot both be z."
        )

    nnode = len(
        nodes
    )

    axis_min = float(
        nodes[
            :,
            axis_id,
        ].min()
    )

    axis_max = float(
        nodes[
            :,
            axis_id,
        ].max()
    )

    half_length = (
        axis_max
        - axis_min
    )

    if half_length <= 0:
        raise ValueError(
            "Invalid repeat-axis component length"
        )

    (
        mirror_elems0,
        orientation_error,
    ) = mirrored_elements(
        nodes,
        elems0,
        axis_id,
    )

    center_set = set(
        int(v)
        for v in center_ids
    )

    gap_set = set(
        int(v)
        for v in gap_ids
    )

    merged_nodes: List[
        np.ndarray
    ] = []

    all_elems: List[
        np.ndarray
    ] = []

    # [depth_layer][x_component]
    component_maps: List[
        List[np.ndarray]
    ] = []

    total_domains = (
        x_components
        * depth_copies
    )

    domain_nodes: List[
        List[int]
    ] = [
        []
        for _ in range(
            total_domains
        )
    ]

    for iz in range(
        depth_copies
    ):
        layer_maps: List[
            np.ndarray
        ] = []

        depth_shift = (
            -iz
            * depth_length
        )

        for ix in range(
            x_components
        ):
            domain_id = (
                iz
                * x_components
                + ix
            )

            is_original = (
                ix % 2 == 0
            )

            if is_original:
                left_center = (
                    start_offset
                    + ix
                    * half_length
                )

                xyz = (
                    transform_original_component(
                        nodes,
                        axis_id,
                        axis_min,
                        left_center,
                        depth_axis_id,
                        depth_shift,
                    )
                )

                elem_template = (
                    elems0
                )

            else:
                right_center = (
                    start_offset
                    + (
                        ix + 1
                    )
                    * half_length
                )

                xyz = (
                    transform_mirror_component(
                        nodes,
                        axis_id,
                        axis_min,
                        right_center,
                        depth_axis_id,
                        depth_shift,
                    )
                )

                elem_template = (
                    mirror_elems0
                )

            local_map = np.full(
                nnode,
                -1,
                dtype=np.int64,
            )

            # ---------------------------------------------------------
            # z/depth interface with previous layer
            # ---------------------------------------------------------
            if iz > 0:
                previous_depth_map = (
                    component_maps[
                        iz - 1
                    ][
                        ix
                    ]
                )

                for (
                    current_front_lid,
                    previous_back_lid,
                ) in (
                    depth_front_to_back.items()
                ):
                    assign_shared_gid(
                        local_map,
                        current_front_lid,
                        int(
                            previous_depth_map[
                                previous_back_lid
                            ]
                        ),
                        (
                            f"depth interface "
                            f"(z layer {iz}, "
                            f"x component {ix})"
                        ),
                    )

            # ---------------------------------------------------------
            # x interface with previous component in same layer
            # ---------------------------------------------------------
            if ix > 0:
                previous_x_map = (
                    layer_maps[
                        ix - 1
                    ]
                )

                if is_original:
                    # Previous = Mirror(C)
                    # center(min) == center(min)
                    shared_ids = (
                        center_ids
                    )

                    label = (
                        "x center interface"
                    )

                else:
                    # Previous = C
                    # gap(max) == gap(max)
                    shared_ids = (
                        gap_ids
                    )

                    label = (
                        "x gap interface"
                    )

                for lid in shared_ids:
                    lid_i = int(
                        lid
                    )

                    assign_shared_gid(
                        local_map,
                        lid_i,
                        int(
                            previous_x_map[
                                lid_i
                            ]
                        ),
                        (
                            f"{label} "
                            f"(z layer {iz}, "
                            f"x component {ix})"
                        ),
                    )

            # ---------------------------------------------------------
            # Create all remaining nodes.
            #
            # Newly created node is owned by this component/domain.
            # Shared nodes remain owned by an earlier domain.
            # ---------------------------------------------------------
            for lid in range(
                nnode
            ):
                if (
                    local_map[
                        lid
                    ]
                    != -1
                ):
                    continue

                gid = len(
                    merged_nodes
                )

                merged_nodes.append(
                    xyz[
                        lid
                    ].copy()
                )

                local_map[
                    lid
                ] = gid

                domain_nodes[
                    domain_id
                ].append(
                    gid
                )

            all_elems.append(
                local_map[
                    elem_template
                ]
            )

            layer_maps.append(
                local_map
            )

        component_maps.append(
            layer_maps
        )

    return (
        np.asarray(
            merged_nodes,
            dtype=float,
        ),
        np.vstack(
            all_elems
        ),
        component_maps,
        domain_nodes,
        orientation_error,
    )


# =============================================================================
# BC generation
# =============================================================================


def parse_dof_list(
    text: str,
    block_length: int,
    label: str,
) -> List[int]:

    dofs = sorted(
        set(
            int(s.strip())
            for s in text.split(
                ","
            )
            if s.strip()
        )
    )

    if not dofs:
        raise ValueError(
            f"{label}: no DOFs specified"
        )

    for dof in dofs:
        if (
            dof < 1
            or dof > block_length
        ):
            raise ValueError(
                f"{label}: DOF {dof} "
                f"is outside "
                f"1..{block_length}"
            )

    return dofs


def parse_bbox_bc(
    spec: str,
    block_length: int,
) -> Tuple[
    str,
    List[int],
]:

    if ":" not in spec:
        raise ValueError(
            f"Invalid --bc '{spec}'. "
            "Example: --bc ymin:1"
        )

    face, dofs = spec.split(
        ":",
        1,
    )

    face = (
        face.strip()
        .lower()
    )

    allowed = {
        "xmin",
        "xmax",
        "ymin",
        "ymax",
        "zmin",
        "zmax",
    }

    if face not in allowed:
        raise ValueError(
            f"Unknown BC face '{face}'"
        )

    return (
        face,
        parse_dof_list(
            dofs,
            block_length,
            f"--bc {face}",
        ),
    )


def build_global_bbox_bc(
    merged_nodes: np.ndarray,
    specs: Sequence[
        Tuple[
            str,
            Sequence[int],
        ]
    ],
    base: int,
    tol: float,
) -> List[
    Tuple[int, int]
]:

    result = set()

    for face, dofs in specs:
        axis_id = {
            "x": 0,
            "y": 1,
            "z": 2,
        }[
            face[0]
        ]

        use_min = (
            face.endswith(
                "min"
            )
        )

        ids = bbox_face_nodes(
            merged_nodes,
            axis_id,
            use_min,
            tol,
        )

        for gid in ids:
            for dof in dofs:
                # GEDATSU D_bc.dat node IDs are written in the global
                # node.dat ordering (0-based), independently of the
                # connectivity base used by elem.dat.  Using ``+ base``
                # here shifts every BC by one when the source connectivity
                # is 1-based.
                result.add(
                    (
                        int(gid),
                        int(dof),
                    )
                )

    return sorted(
        result,
        key=lambda x: (
            x[0],
            x[1],
        ),
    )


def build_surface_bc_grid(
    component_maps: Sequence[
        Sequence[np.ndarray]
    ],
    boundary_faces0: np.ndarray,
    x_center_faces0: np.ndarray,
    x_gap_faces0: np.ndarray,
    depth_front_faces0: np.ndarray,
    depth_back_faces0: np.ndarray,
    base: int,
    dofs: Sequence[int],
) -> List[
    Tuple[int, int]
]:

    center_keys = face_key_set(
        x_center_faces0
    )

    gap_keys = face_key_set(
        x_gap_faces0
    )

    front_keys = face_key_set(
        depth_front_faces0
    )

    back_keys = face_key_set(
        depth_back_faces0
    )

    nz = len(
        component_maps
    )

    nx = len(
        component_maps[0]
    )

    bc_nodes = set()

    for iz in range(
        nz
    ):
        for ix in range(
            nx
        ):
            local_map = (
                component_maps[
                    iz
                ][
                    ix
                ]
            )

            is_original = (
                ix % 2 == 0
            )

            for face in boundary_faces0:
                key = tuple(
                    sorted(
                        int(v)
                        for v in face
                    )
                )

                # -------------------------------------------------
                # x direction
                # -------------------------------------------------

                # Gamma03 is always an internal C<->M gap interface.
                if key in gap_keys:
                    continue

                if key in center_keys:
                    # Original C center/min:
                    # external only at global x-min.
                    if is_original:
                        if ix != 0:
                            continue

                    # Mirror(C) center/min:
                    # external only at global x-max.
                    else:
                        if ix != nx - 1:
                            continue

                # -------------------------------------------------
                # z/depth direction
                # -------------------------------------------------

                # Front z=max is external only on first depth layer.
                if (
                    key in front_keys
                    and iz != 0
                ):
                    continue

                # Back z=min is external only on last depth layer.
                if (
                    key in back_keys
                    and iz != nz - 1
                ):
                    continue

                # -------------------------------------------------
                # This face remains external.
                # Add its nodes.
                # -------------------------------------------------
                for lid in face:
                    bc_nodes.add(
                        int(
                            local_map[
                                int(lid)
                            ]
                        )
                    )

    rows: List[
        Tuple[int, int]
    ] = []

    for gid in sorted(
        bc_nodes
    ):
        for dof in dofs:
            # Same convention as build_global_bbox_bc(): D_bc.dat uses
            # the 0-based global node.dat ordering.
            rows.append(
                (
                    gid,
                    int(dof),
                )
            )

    return rows


# =============================================================================
# Validation
# =============================================================================


def validate_connectivity(
    elems0: np.ndarray,
    nnode: int,
) -> None:

    if elems0.size == 0:
        raise RuntimeError(
            "Generated mesh has zero elements"
        )

    mn = int(
        elems0.min()
    )

    mx = int(
        elems0.max()
    )

    if (
        mn < 0
        or mx >= nnode
    ):
        raise RuntimeError(
            "Generated element connectivity "
            "contains invalid node IDs: "
            f"range=[{mn},{mx}], "
            f"nnode={nnode}"
        )


def validate_domain_ownership(
    domain_nodes: Sequence[
        Sequence[int]
    ],
    nnode: int,
) -> None:

    owner = np.full(
        nnode,
        -1,
        dtype=np.int64,
    )

    for domain_id, ids in enumerate(
        domain_nodes
    ):
        if not ids:
            raise RuntimeError(
                f"Domain {domain_id} "
                "owns zero nodes."
            )

        for gid in ids:
            gid_i = int(
                gid
            )

            if owner[
                gid_i
            ] != -1:
                raise RuntimeError(
                    f"Node {gid_i} has "
                    "multiple owners: "
                    f"{int(owner[gid_i])} "
                    f"and {domain_id}"
                )

            owner[
                gid_i
            ] = domain_id

    missing = np.where(
        owner < 0
    )[0]

    if len(
        missing
    ):
        raise RuntimeError(
            f"{len(missing)} nodes have "
            "no domain owner. "
            f"first={int(missing[0])}"
        )


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Build a continuous symmetric thermal-fin component array "
            "in x and extend it continuously in the negative-z depth direction."
        )
    )

    # -----------------------------------------------------------------
    # Component input
    # -----------------------------------------------------------------

    p.add_argument(
        "node_dat",
        nargs="?",
        type=Path,
    )

    p.add_argument(
        "elem_dat",
        nargs="?",
        type=Path,
    )

    p.add_argument(
        "--geo",
        type=Path,
        default=None,
    )

    p.add_argument(
        "--gmsh",
        default="gmsh",
    )

    p.add_argument(
        "--physical-group-script",
        type=Path,
        default=Path(
            "mesh_io/save_physical_groups.py"
        ),
    )

    p.add_argument(
        "--component-work-dir",
        type=Path,
        default=Path(
            "generated_component"
        ),
    )

    p.add_argument(
        "--volume-group",
        default="Fluid",
    )

    # GEO groups
    p.add_argument(
        "--center-group",
        action="append",
        default=[],
    )

    p.add_argument(
        "--gap-group",
        action="append",
        default=[],
    )

    p.add_argument(
        "--depth-front-group",
        action="append",
        default=[],
    )

    p.add_argument(
        "--depth-back-group",
        action="append",
        default=[],
    )

    p.add_argument(
        "--bc-group",
        default="All_external_boundaries",
    )

    p.add_argument(
        "--no-auto-bc",
        action="store_true",
    )

    p.add_argument(
        "--no-clean-component-work",
        action="store_true",
    )

    # Direct surface files
    p.add_argument(
        "--center-interface",
        action="append",
        type=Path,
        default=[],
    )

    p.add_argument(
        "--gap-interface",
        action="append",
        type=Path,
        default=[],
    )

    p.add_argument(
        "--depth-front-interface",
        action="append",
        type=Path,
        default=[],
    )

    p.add_argument(
        "--depth-back-interface",
        action="append",
        type=Path,
        default=[],
    )

    # -----------------------------------------------------------------
    # Repetition
    # -----------------------------------------------------------------

    p.add_argument(
        "--copies",
        type=int,
        required=True,
        help=(
            "Number of GEO components in the symmetric x direction. "
            "Must be even."
        ),
    )

    p.add_argument(
        "--depth-copies",
        type=int,
        default=1,
        help=(
            "Number of continuous component layers in the negative-z "
            "depth direction. Total GEDATSU domains = "
            "copies * depth-copies. Default: 1"
        ),
    )

    p.add_argument(
        "--axis",
        choices=(
            "x",
            "y",
            "z",
        ),
        default="x",
    )

    p.add_argument(
        "--gap",
        type=float,
        default=0.0,
        help=(
            "Must be 0.0 for conformal continuous welding."
        ),
    )

    p.add_argument(
        "--symmetric",
        action="store_true",
        help=(
            "Use C,Mirror(C),C,Mirror(C),... in the repeat direction."
        ),
    )

    p.add_argument(
        "--start-offset",
        type=float,
        default=0.0,
    )

    p.add_argument(
        "--input-base",
        choices=(
            "auto",
            "0",
            "1",
        ),
        default="auto",
    )

    p.add_argument(
        "--weld-tol",
        type=float,
        default=1.0e-10,
    )

    # -----------------------------------------------------------------
    # BC
    # -----------------------------------------------------------------

    p.add_argument(
        "--bc-surface",
        action="append",
        type=Path,
        default=[],
    )

    p.add_argument(
        "--bc-surface-dofs",
        default="1",
    )

    p.add_argument(
        "--bc",
        action="append",
        default=[],
        metavar="FACE:DOFS",
        help=(
            "Optional global bounding-box BC, e.g. --bc ymin:1"
        ),
    )

    p.add_argument(
        "--block-length",
        type=int,
        default=1,
    )

    p.add_argument(
        "--bc-value",
        type=float,
        default=0.0,
    )

    p.add_argument(
        "--bc-tol",
        type=float,
        default=1.0e-10,
    )

    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path(
            "merged_mesh"
        ),
    )

    args = p.parse_args()

    # Remember whether the user explicitly requested a surface-based BC.
    # For this generator, bbox BCs (--bc ...) and surface BCs are mutually
    # exclusive unless the code is intentionally extended later.  In
    # particular, edge-heating with --bc xmin:1 must NEVER inherit the
    # default All_external_boundaries surface group.
    explicit_bc_surface = bool(args.bc_surface)

    if args.bc and explicit_bc_surface:
        raise ValueError(
            "Do not combine --bc with --bc-surface. "
            "For edge heating use bbox BC only, e.g. --bc xmin:1."
        )

    # -----------------------------------------------------------------
    # Basic validation
    # -----------------------------------------------------------------

    if not args.symmetric:
        raise ValueError(
            "This version is intended for the symmetric component array. "
            "Specify --symmetric."
        )

    if args.copies < 2:
        raise ValueError(
            "--copies must be >= 2"
        )

    if (
        args.copies
        % 2
        != 0
    ):
        raise ValueError(
            "--copies must be even."
        )

    if args.depth_copies < 1:
        raise ValueError(
            "--depth-copies must be >= 1"
        )

    if abs(
        args.gap
    ) > args.weld_tol:
        raise ValueError(
            "Continuous mode requires --gap 0.0"
        )

    if args.block_length < 1:
        raise ValueError(
            "--block-length must be >= 1"
        )

    if args.weld_tol <= 0:
        raise ValueError(
            "--weld-tol must be > 0"
        )

    axis_id = {
        "x": 0,
        "y": 1,
        "z": 2,
    }[
        args.axis
    ]

    depth_axis_id = 2

    if (
        args.depth_copies > 1
        and axis_id == depth_axis_id
    ):
        raise ValueError(
            "--axis z cannot be used together "
            "with --depth-copies > 1. "
            "The depth extension itself uses z."
        )

    # -----------------------------------------------------------------
    # Resolve component source
    # -----------------------------------------------------------------

    geo_info = None

    if args.geo is not None:
        if (
            args.node_dat is not None
            or args.elem_dat is not None
        ):
            raise ValueError(
                "Use either --geo OR positional "
                "node_dat/elem_dat."
            )

        center_groups = (
            args.center_group
            if args.center_group
            else [
                "Gamma11",
                "Gamma12",
            ]
        )

        gap_groups = (
            args.gap_group
            if args.gap_group
            else [
                "Gamma03",
            ]
        )

        depth_front_groups = (
            args.depth_front_group
            if args.depth_front_group
            else [
                "Gamma01",
                "Gamma07",
            ]
        )

        depth_back_groups = (
            args.depth_back_group
            if args.depth_back_group
            else [
                "Gamma05",
                "Gamma09",
            ]
        )

        # STRICT BC precedence:
        #   1) Any explicit bbox BC (--bc ...) disables the automatic GEO
        #      surface BC unconditionally.
        #   2) --bc-surface is not allowed together with --bc (checked above).
        # This guarantees that --bc xmin:1 produces ONLY the global xmin face.
        if args.bc:
            bc_group = None
            args.bc_surface = []
            print(
                "[BC-GENERATION] explicit bbox BC supplied; "
                "automatic GEO surface BC is disabled and ignored."
            )
        elif args.no_auto_bc:
            bc_group = None
        else:
            bc_group = args.bc_group

        geo_info = generate_component_from_geo(
            geo_path=args.geo,
            gmsh_executable=args.gmsh,
            physical_group_script=(
                args.physical_group_script
            ),
            component_work_dir=(
                args.component_work_dir
            ),
            volume_group=(
                args.volume_group
            ),
            center_groups=(
                center_groups
            ),
            gap_groups=(
                gap_groups
            ),
            depth_front_groups=(
                depth_front_groups
            ),
            depth_back_groups=(
                depth_back_groups
            ),
            bc_group=(
                bc_group
            ),
            clean_work=(
                not args.no_clean_component_work
            ),
        )

        args.node_dat = (
            geo_info[
                "node_dat"
            ]
        )

        args.elem_dat = (
            geo_info[
                "elem_dat"
            ]
        )

        if not args.center_interface:
            args.center_interface = list(
                geo_info[
                    "center_interface"
                ]
            )

        if not args.gap_interface:
            args.gap_interface = list(
                geo_info[
                    "gap_interface"
                ]
            )

        if not args.depth_front_interface:
            args.depth_front_interface = list(
                geo_info[
                    "depth_front_interface"
                ]
            )

        if not args.depth_back_interface:
            args.depth_back_interface = list(
                geo_info[
                    "depth_back_interface"
                ]
            )

        if (
            not args.bc
            and not args.bc_surface
            and geo_info[
                "bc_surface"
            ]
        ):
            args.bc_surface = list(
                geo_info[
                    "bc_surface"
                ]
            )

    else:
        if (
            args.node_dat is None
            or args.elem_dat is None
        ):
            raise ValueError(
                "Specify --geo or positional "
                "node.dat elem.dat."
            )

        args.node_dat = require_file(
            args.node_dat,
            "node.dat",
        )

        args.elem_dat = require_file(
            args.elem_dat,
            "elem.dat",
        )

    if not args.center_interface:
        raise ValueError(
            "x center interface is required"
        )

    if not args.gap_interface:
        raise ValueError(
            "x gap interface is required"
        )

    if args.depth_copies > 1:
        if not args.depth_front_interface:
            raise ValueError(
                "depth front interface is required "
                "when --depth-copies > 1"
            )

        if not args.depth_back_interface:
            raise ValueError(
                "depth back interface is required "
                "when --depth-copies > 1"
            )

    # -----------------------------------------------------------------
    # Read component mesh
    # -----------------------------------------------------------------

    nodes = read_nodes(
        args.node_dat
    )

    elems, nen = read_elements(
        args.elem_dat
    )

    nnode = len(
        nodes
    )

    nelem = len(
        elems
    )

    base = detect_base(
        elems,
        nnode,
        args.input_base,
    )

    elems0 = (
        elems
        - base
    )

    print(
        f"Detected connectivity base: "
        f"{base}-based"
    )

    # -----------------------------------------------------------------
    # Read x interfaces
    # -----------------------------------------------------------------

    center_faces0 = read_surface_faces_union(
        args.center_interface,
        base,
        nnode,
        "x center interface",
    )

    gap_faces0 = read_surface_faces_union(
        args.gap_interface,
        base,
        nnode,
        "x gap interface",
    )

    center_ids = surface_node_ids(
        center_faces0
    )

    gap_ids = surface_node_ids(
        gap_faces0
    )

    validate_x_interfaces(
        nodes,
        axis_id,
        center_ids,
        gap_ids,
        args.weld_tol,
    )

    # -----------------------------------------------------------------
    # Read z/depth interfaces
    # -----------------------------------------------------------------

    if args.depth_copies > 1:
        depth_front_faces0 = (
            read_surface_faces_union(
                args.depth_front_interface,
                base,
                nnode,
                "z front interface",
            )
        )

        depth_back_faces0 = (
            read_surface_faces_union(
                args.depth_back_interface,
                base,
                nnode,
                "z back interface",
            )
        )

        depth_front_ids = surface_node_ids(
            depth_front_faces0
        )

        depth_back_ids = surface_node_ids(
            depth_back_faces0
        )

        (
            depth_front_to_back,
            depth_pair_error,
            depth_length,
        ) = build_depth_pairing(
            nodes,
            depth_axis_id,
            depth_front_ids,
            depth_back_ids,
            args.weld_tol,
        )

    else:
        # Still read these surfaces if available, because BC face filtering
        # uses them. If absent, no depth face is removed for a single layer.
        if args.depth_front_interface:
            depth_front_faces0 = (
                read_surface_faces_union(
                    args.depth_front_interface,
                    base,
                    nnode,
                    "z front interface",
                )
            )
        else:
            depth_front_faces0 = np.empty(
                (0, 0),
                dtype=np.int64,
            )

        if args.depth_back_interface:
            depth_back_faces0 = (
                read_surface_faces_union(
                    args.depth_back_interface,
                    base,
                    nnode,
                    "z back interface",
                )
            )
        else:
            depth_back_faces0 = np.empty(
                (0, 0),
                dtype=np.int64,
            )

        depth_front_to_back = {}

        zmin = float(
            nodes[:, depth_axis_id].min()
        )

        zmax = float(
            nodes[:, depth_axis_id].max()
        )

        depth_length = (
            zmax - zmin
        )

        depth_pair_error = 0.0

    # -----------------------------------------------------------------
    # Build x-z component grid
    # -----------------------------------------------------------------

    (
        merged_nodes,
        merged_elems0,
        component_maps,
        domain_nodes,
        orientation_error,
    ) = build_component_grid(
        nodes=nodes,
        elems0=elems0,
        x_components=args.copies,
        depth_copies=args.depth_copies,
        axis_id=axis_id,
        depth_axis_id=depth_axis_id,
        start_offset=args.start_offset,
        center_ids=center_ids,
        gap_ids=gap_ids,
        depth_front_to_back=(
            depth_front_to_back
        ),
        depth_length=(
            depth_length
        ),
    )

    total_domains = (
        args.copies
        * args.depth_copies
    )

    validate_connectivity(
        merged_elems0,
        len(
            merged_nodes
        ),
    )

    validate_domain_ownership(
        domain_nodes,
        len(
            merged_nodes
        ),
    )

    expected_elements = (
        total_domains
        * nelem
    )

    if (
        len(
            merged_elems0
        )
        != expected_elements
    ):
        raise RuntimeError(
            "Element-count check failed: "
            f"actual={len(merged_elems0)}, "
            f"expected={expected_elements}"
        )

    # -----------------------------------------------------------------
    # BC
    # -----------------------------------------------------------------

    bc_rows_set = set()

    bbox_specs = [
        parse_bbox_bc(
            spec,
            args.block_length,
        )
        for spec in args.bc
    ]

    if bbox_specs and args.bc_surface:
        raise RuntimeError(
            "Internal BC configuration error: bbox BC and surface BC "
            "are both active. This generator requires bbox-only BC when "
            "--bc is supplied."
        )

    bbox_rows = build_global_bbox_bc(
        merged_nodes,
        bbox_specs,
        base,
        args.bc_tol,
    )

    for row in bbox_rows:
        bc_rows_set.add(row)

    surface_dofs = parse_dof_list(
        args.bc_surface_dofs,
        args.block_length,
        "--bc-surface-dofs",
    )

    if args.bc_surface:
        boundary_faces0 = (
            read_surface_faces_union(
                args.bc_surface,
                base,
                nnode,
                "BC surface",
            )
        )

        for row in build_surface_bc_grid(
            component_maps=component_maps,
            boundary_faces0=boundary_faces0,
            x_center_faces0=center_faces0,
            x_gap_faces0=gap_faces0,
            depth_front_faces0=(
                depth_front_faces0
            ),
            depth_back_faces0=(
                depth_back_faces0
            ),
            base=base,
            dofs=surface_dofs,
        ):
            bc_rows_set.add(
                row
            )

    # Hard invariant for bbox-only BC generation.  This catches accidental
    # re-introduction of All_external_boundaries or any other surface group.
    if bbox_specs:
        expected_bbox_rows = set(bbox_rows)
        if bc_rows_set != expected_bbox_rows:
            extra = sorted(bc_rows_set - expected_bbox_rows)
            missing = sorted(expected_bbox_rows - bc_rows_set)
            raise RuntimeError(
                "bbox-only BC invariant failed: "
                f"extra={len(extra)} missing={len(missing)} "
                f"extra(first10)={extra[:10]} missing(first10)={missing[:10]}"
            )
        print(
            f"[BC-GENERATION] mode=bbox-only specs={args.bc} "
            f"rows={len(expected_bbox_rows)}"
        )

    bc_rows = sorted(
        bc_rows_set,
        key=lambda x: (
            x[0],
            x[1],
        ),
    )

    generate_bc = bool(
        args.bc_surface
        or bbox_specs
    )

    # -----------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------

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

    bc_out = (
        args.out_dir
        / "D_bc.dat"
    )

    report_out = (
        args.out_dir
        / "repeat_report.txt"
    )

    write_nodes(
        node_out,
        merged_nodes,
    )

    write_elements(
        elem_out,
        merged_elems0,
        nen,
        base,
    )

    write_given_subdomain(
        subdomain_out,
        domain_nodes,
        base,
    )

    if generate_bc:
        write_bc(
            bc_out,
            bc_rows,
            args.block_length,
            args.bc_value,
        )
    elif bc_out.exists():
        bc_out.unlink()

    # -----------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------

    disconnected_nodes = (
        total_domains
        * nnode
    )

    welded_duplicates = (
        disconnected_nodes
        - len(
            merged_nodes
        )
    )

    xyz_min_out = (
        merged_nodes.min(
            axis=0
        )
    )

    xyz_max_out = (
        merged_nodes.max(
            axis=0
        )
    )

    with report_out.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "Symmetric x-z component array report\n"
        )
        f.write(
            "====================================\n"
        )

        if geo_info is not None:
            f.write(
                f"component source       : GEO\n"
            )
            f.write(
                f"component geo          : "
                f"{geo_info['geo']}\n"
            )
        else:
            f.write(
                "component source       : "
                "node.dat / elem.dat\n"
            )

        f.write(
            f"component nodes        : {nnode}\n"
        )
        f.write(
            f"component elements     : {nelem}\n"
        )
        f.write(
            f"nodes/element          : {nen}\n"
        )
        f.write(
            f"index base             : {base}\n"
        )
        f.write(
            "D_bc node id base       : 0 (GEDATSU global node ordering)\n"
        )
        f.write(
            f"x components           : {args.copies}\n"
        )
        f.write(
            f"depth copies           : {args.depth_copies}\n"
        )
        f.write(
            f"total GEDATSU domains  : {total_domains}\n"
        )
        f.write(
            f"domain numbering       : "
            "id = iz*x_components + ix\n"
        )
        f.write(
            f"x center nodes         : {len(center_ids)}\n"
        )
        f.write(
            f"x gap nodes            : {len(gap_ids)}\n"
        )

        if args.depth_copies > 1:
            f.write(
                f"depth front nodes      : "
                f"{len(surface_node_ids(depth_front_faces0))}\n"
            )
            f.write(
                f"depth back nodes       : "
                f"{len(surface_node_ids(depth_back_faces0))}\n"
            )

        f.write(
            f"depth length/component : {depth_length:.16g}\n"
        )
        f.write(
            f"depth pair error       : {depth_pair_error:.16g}\n"
        )
        f.write(
            f"mirror orient. error   : {orientation_error:.16g}\n"
        )
        f.write(
            f"output bbox min        : {xyz_min_out.tolist()}\n"
        )
        f.write(
            f"output bbox max        : {xyz_max_out.tolist()}\n"
        )
        f.write(
            f"welded duplicates      : {welded_duplicates}\n"
        )
        f.write(
            f"output nodes           : {len(merged_nodes)}\n"
        )
        f.write(
            f"output elements        : {len(merged_elems0)}\n"
        )

        f.write(
            "\nOwned nodes/domain\n"
        )

        for domain_id, ids in enumerate(
            domain_nodes
        ):
            iz = (
                domain_id
                // args.copies
            )

            ix = (
                domain_id
                % args.copies
            )

            f.write(
                f"  domain {domain_id}: "
                f"iz={iz}, ix={ix}, "
                f"owned={len(ids)}\n"
            )

        if generate_bc:
            f.write(
                "\nBC\n"
            )
            f.write(
                f"block length           : {args.block_length}\n"
            )
            f.write(
                f"surface files          : "
                f"{[str(p) for p in args.bc_surface]}\n"
            )
            f.write(
                f"surface dofs           : {surface_dofs}\n"
            )
            f.write(
                f"bbox specs             : {args.bc}\n"
            )
            f.write(
                f"BC value               : {args.bc_value:.16g}\n"
            )
            f.write(
                f"D_bc rows               : {len(bc_rows)}\n"
            )

        f.write(
            "\nGEDATSU\n"
        )
        f.write(
            f"cd {args.out_dir}\n"
        )
        f.write(
            "mkdir -p parted.0\n"
        )
        f.write(
            "gedatsu_simple_mesh_partitioner "
            f"-n {total_domains} "
            "-in node.dat "
            "-ie elem.dat "
            "--given_subdomain given_subdomain.dat "
            "-d parted.0\n"
        )
        if generate_bc:
            f.write(
                "gedatsu_bc_partitioner_R "
                f"-n {total_domains} "
                "-i D_bc.dat "
                "-ig node.dat\n"
            )

    # -----------------------------------------------------------------
    # Console
    # -----------------------------------------------------------------

    print()
    print("Done.")
    print()
    print(f"x components        : {args.copies}")
    print(f"depth copies        : {args.depth_copies}")
    print(f"GEDATSU domains     : {total_domains}")
    print(f"component nodes     : {nnode}")
    print(f"component elements  : {nelem}")
    print(f"depth length        : {depth_length}")
    print(f"depth pair error    : {depth_pair_error}")
    print(f"welded duplicates  : {welded_duplicates}")
    print(f"output nodes        : {len(merged_nodes)}")
    print(f"output elements     : {len(merged_elems0)}")
    print(f"output bbox min     : {xyz_min_out.tolist()}")
    print(f"output bbox max     : {xyz_max_out.tolist()}")
    print()
    print(f"node.dat            : {node_out}")
    print(f"elem.dat            : {elem_out}")
    print(f"given_subdomain.dat : {subdomain_out}")

    if generate_bc:
        print(f"D_bc.dat            : {bc_out}")
        print(f"D_bc rows           : {len(bc_rows)}")
        print(f"BC block length     : {args.block_length}")
    else:
        print("D_bc.dat            : not generated")

    print(f"report              : {report_out}")
    print()
    print("GEDATSU command:")
    print()
    print(f"cd {args.out_dir}")
    print("mkdir -p parted.0")
    print(
        "gedatsu_simple_mesh_partitioner "
        f"-n {total_domains} "
        "-in node.dat "
        "-ie elem.dat "
        "--given_subdomain given_subdomain.dat "
        "-d parted.0"
    )

    if generate_bc:
        print(
            "gedatsu_bc_partitioner_R "
            f"-n {total_domains} "
            "-i D_bc.dat "
            "-ig node.dat"
        )


if __name__ == "__main__":
    main()
