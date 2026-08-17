#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
repeat_dat_mesh_with_bc.py

node.dat / elem.dat で与えられた 1 つの「半コンポーネント」から、
通常の複製、または左右対称セルの連続メッシュを生成する。

主な出力
--------
- node.dat
- elem.dat
- given_subdomain.dat
- repeat_report.txt
- D_bc.dat              (--bc 指定時)


============================================================
GEO component mode (recommended for thermal_fin_fig3.geo)
============================================================

The half component can be generated directly from the .geo file.
The script runs:

    gmsh -3 -format msh2 -0
        ->
    component.msh
        ->
    mesh_io/save_physical_groups.py
        ->
    generated_component/mesh_tmp/

Default Physical Groups for the supplied thermal-fin GEO:

    volume:
        Fluid

    fin-center side:
        Gamma11
        Gamma12

    gap side:
        Gamma03

    BC:
        All_external_boundaries

Example:

python3 repeat_dat_mesh_with_bc_symmetric.py \
    --geo gmsh/thermal_fin_fig3.geo \
    --physical-group-script mesh_io/save_physical_groups.py \
    --copies 8 \
    --axis x \
    --gap 0.0 \
    --symmetric \
    --block-length 1 \
    --bc-surface-dofs 1 \
    --bc-value 0.0 \
    --out-dir merged_mesh

In this mode it is NOT necessary to manually provide:
    node.dat
    elem.dat
    Gamma11_quad_connectivity.dat
    Gamma12_quad_connectivity.dat
    Gamma03_quad_connectivity.dat
    All_external_boundaries_quad_connectivity.dat

They are taken from the GEO-generated Physical Groups automatically.

============================================================
1. 通常モード
============================================================

従来どおり 1 コンポーネントを平行移動して複製する。

disconnected:
    C | C | C | C

continuous:
    C == C == C == C

continuous の場合は

    --interface-left
    --interface-right

で明示した surface connectivity を使って node weld する。

============================================================
2. 左右対称セルモード: --symmetric
============================================================

今回の thermal-fin 半コンポーネントを想定したモード。

元コンポーネント C:

    repeat-axis min 側:
        Gamma11 + Gamma12
        （フィン側 + spreader 側）

    repeat-axis max 側:
        Gamma03
        （spreader 側）

1 個の左右対称ピッチを

    C == Mirror(C)

として作る。

この順序にすることで、全体の左端と右端に fin 側
(center interface 側) が残る。

1 ピッチ内の接続:
    C の max 側 (Gamma03)
        ==
    Mirror(C) の max 側 (Gamma03)

隣接ピッチ間の接続:
    Mirror(C) の min 側 (Gamma11 + Gamma12)
        ==
    次の C の min 側 (Gamma11 + Gamma12)

従って複数ピッチは

    C == Mirror(C) == C == Mirror(C) == C == Mirror(C) ...

となり、

    fin | gap | fin | gap | fin | ... | fin

のように左右端にも fin が残る。

--copies は symmetric モードでは
「GEO から生成された半コンポーネント数」
= 「GEDATSU domain 数」を表す。

1 domain = 1 GEO component とし、

    domain 0 : C
    domain 1 : Mirror(C)
    domain 2 : C
    domain 3 : Mirror(C)
    ...

のように 1 対 1 で割り当てる。

左右端を fin 側にするため、--copies は偶数とする。

例:
    --copies 8 --symmetric

なら

    C M C M C M C M

を作り、GEDATSU の domain 数も 8 とする。

これは旧仕様の

    --copies 8

（4個の C+M ピッチ、domain も4）

と同じ全体形状を作る。
新仕様では 1 GEO component = 1 domain なので 8 が必要になる。

============================================================
Surface connectivity
============================================================

surface connectivity は elem.dat と同様:

    Nface NEN
    n1 n2 ... nNEN
    ...

左右対称セルでは center interface を複数ファイルで指定できる。

thermal_fin_fig3.geo の場合:

    center:
        Gamma11_quad_connectivity.dat
        Gamma12_quad_connectivity.dat

    gap/cell boundary:
        Gamma03_quad_connectivity.dat

実行例:

python3 repeat_dat_mesh_with_bc.py \
    ./mesh_karman_vortex/node.dat \
    ./mesh_karman_vortex/elem.dat \
    --copies 8 \
    --axis x \
    --gap 0.0 \
    --symmetric \
    --center-interface ./mesh_tmp/Gamma11_quad_connectivity.dat \
    --center-interface ./mesh_tmp/Gamma12_quad_connectivity.dat \
    --gap-interface ./mesh_tmp/Gamma03_quad_connectivity.dat \
    --out-dir merged_mesh

scalar convdiff の全面 bounding-box Dirichlet BC も作る例:

python3 repeat_dat_mesh_with_bc.py \
    ./mesh_karman_vortex/node.dat \
    ./mesh_karman_vortex/elem.dat \
    --copies 8 \
    --axis x \
    --gap 0.0 \
    --symmetric \
    --center-interface ./mesh_tmp/Gamma11_quad_connectivity.dat \
    --center-interface ./mesh_tmp/Gamma12_quad_connectivity.dat \
    --gap-interface ./mesh_tmp/Gamma03_quad_connectivity.dat \
    --block-length 1 \
    --bc xmin:1 \
    --bc xmax:1 \
    --bc ymin:1 \
    --bc ymax:1 \
    --bc zmin:1 \
    --bc zmax:1 \
    --bc-value 0.0 \
    --out-dir merged_mesh

実形状の全外部境界に BC を付ける推奨例:

python3 repeat_dat_mesh_with_bc_symmetric.py \
    ./mesh_karman_vortex/node.dat \
    ./mesh_karman_vortex/elem.dat \
    --copies 8 \
    --axis x \
    --gap 0.0 \
    --symmetric \
    --center-interface ./mesh_tmp/Gamma11_quad_connectivity.dat \
    --center-interface ./mesh_tmp/Gamma12_quad_connectivity.dat \
    --gap-interface ./mesh_tmp/Gamma03_quad_connectivity.dat \
    --bc-surface ./mesh_tmp/All_external_boundaries_quad_connectivity.dat \
    --block-length 1 \
    --bc-surface-dofs 1 \
    --bc-value 0.0 \
    --out-dir merged_mesh

--bc-surface を使う場合、C/M の weld で内部化する Gamma03 および
Gamma11+Gamma12 の surface elements は自動的に BC 対象から除外する。
ただし最左端 C の center 面と最右端 Mirror(C) の center 面は外表面なので残す。

重要:
- symmetric モードでは --gap は 0 でなければならない。
- symmetric モードは現在 Gmsh 8-node hexahedron (NEN=8) を対象とする。
- 鏡像側 hexahedron は node 順を反転して orientation/Jacobian を維持する。
- center/gap interface は surface connectivity から unique node ID を取得する。
- center interface は元コンポーネントの repeat-axis min 面全体、
  gap interface は repeat-axis max 面全体であることを検査する。
- shared node は given_subdomain.dat で重複所有させない。
  セル間 shared node は前の domain が owner になる。
- D_bc.dat を生成しない実行では、古い D_bc.dat が残らないよう削除する。
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
# GEO component generation
# =============================================================================


def require_file(
    path: Path,
    label: str,
) -> Path:
    path = path.resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"{label} was not found: {path}"
        )

    return path


def require_executable(
    executable: str,
) -> str:
    """
    Resolve an executable name/path.
    """

    p = Path(executable)

    if p.parent != Path(".") or p.is_absolute():
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
    bc_group: str | None,
    clean_work: bool,
) -> Dict[str, object]:
    """
    Generate the half-component directly from a Gmsh .geo file.

    Workflow
    --------
    .geo
      -> gmsh -3 -format msh2 -0
      -> component.msh
      -> save_physical_groups.py
      -> mesh_tmp/*_connectivity.dat
      -> use these files as the component input

    The current project convention is assumed:
        Fluid_node_coordinates.dat
        Fluid_hexahedron_connectivity.dat
        GammaXX_quad_connectivity.dat
    """

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
    print(
        "Generating component from GEO..."
    )
    print(
        f"geo                 : {geo_path}"
    )
    print(
        f"component work dir  : {component_work_dir}"
    )
    print(
        f"gmsh                : {gmsh_executable}"
    )
    print()

    gmsh_cmd = [
        gmsh_executable,
        "-3",
        "-format",
        "msh2",
        "-0",
        str(geo_path),
        "-o",
        str(msh_out),
    ]

    subprocess.run(
        gmsh_cmd,
        cwd=component_work_dir,
        check=True,
    )

    if not msh_out.is_file():
        raise RuntimeError(
            "Gmsh completed but the requested MSH file "
            f"was not created: {msh_out}"
        )

    # Run the existing project extractor in the component work directory.
    #
    # The observed save_physical_groups.py convention writes:
    #     ./mesh_tmp/...
    subprocess.run(
        [
            sys.executable,
            str(physical_group_script),
            str(msh_out),
        ],
        cwd=component_work_dir,
        check=True,
    )

    node_dat = (
        mesh_tmp_dir
        / f"{volume_group}_node_coordinates.dat"
    )

    elem_dat = (
        mesh_tmp_dir
        / f"{volume_group}_hexahedron_connectivity.dat"
    )

    node_dat = require_file(
        node_dat,
        f"{volume_group} node coordinate file",
    )

    elem_dat = require_file(
        elem_dat,
        f"{volume_group} hexahedron connectivity file",
    )

    center_files = [
        require_file(
            physical_surface_file(
                mesh_tmp_dir,
                name,
            ),
            f"center physical surface '{name}'",
        )
        for name in center_groups
    ]

    gap_files = [
        require_file(
            physical_surface_file(
                mesh_tmp_dir,
                name,
            ),
            f"gap physical surface '{name}'",
        )
        for name in gap_groups
    ]

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
    print(
        "GEO component generation completed."
    )
    print(
        f"node input          : {node_dat}"
    )
    print(
        f"elem input          : {elem_dat}"
    )
    print(
        f"center surfaces     : {[str(p) for p in center_files]}"
    )
    print(
        f"gap surfaces        : {[str(p) for p in gap_files]}"
    )

    if bc_files:
        print(
            f"BC surface          : {bc_files[0]}"
        )

    return {
        "geo": geo_path,
        "msh": msh_out,
        "mesh_tmp": mesh_tmp_dir,
        "node_dat": node_dat,
        "elem_dat": elem_dat,
        "center_interface": center_files,
        "gap_interface": gap_files,
        "bc_surface": bc_files,
        "component_work_dir": component_work_dir,
    }


# =============================================================================
# Basic I/O
# =============================================================================


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

    if elems.size == 0:
        elems = np.empty((0, nen), dtype=np.int64)
    elif elems.ndim == 1:
        elems = elems.reshape(1, -1)

    if elems.shape != (nelem, nen):
        raise ValueError(
            f"{path}: header=({nelem}, {nen}), actual shape={elems.shape}"
        )

    return elems, nen


def detect_base(
    elems: np.ndarray,
    nnode: int,
    requested: str,
) -> int:
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
            raise ValueError(
                "Invalid 1-based connectivity"
            )
        return 1

    if 0 < mn and mx < nnode:
        raise ValueError(
            "Index base is ambiguous. "
            "Specify --input-base 0 or --input-base 1. "
            f"Connectivity range=[{mn}, {mx}], nnode={nnode}"
        )

    raise ValueError(
        f"Invalid connectivity range [{mn}, {mx}] for nnode={nnode}"
    )


def write_nodes(
    path: Path,
    nodes: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(nodes)}\n")

        for x, y, z in nodes:
            f.write(
                f"{x:.16g} {y:.16g} {z:.16g}\n"
            )


def write_elements(
    path: Path,
    elems0: np.ndarray,
    nen: int,
    base: int,
) -> None:
    elems_out = elems0 + base

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"{len(elems_out)} {nen}\n"
        )

        for row in elems_out:
            f.write(
                " ".join(
                    str(int(v))
                    for v in row
                )
                + "\n"
            )


def write_given_subdomain(
    path: Path,
    domain_nodes_zero: Sequence[Sequence[int]],
    base: int,
) -> None:
    """
    GEDATSU --given_subdomain graph format.

    domain_nodes_zero[k]:
        domain k が所有する zero-based global node IDs.

    shared node は複数 domain に重複して書かない。
    """

    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"{len(domain_nodes_zero)}\n"
        )

        for domain_id, node_ids in enumerate(domain_nodes_zero):
            ids0 = sorted(
                set(int(v) for v in node_ids)
            )

            if not ids0:
                raise ValueError(
                    f"Subdomain {domain_id} owns zero nodes."
                )

            point_id = (
                domain_id
                + base
            )

            ids_out = [
                v + base
                for v in ids0
            ]

            f.write(
                f"{point_id} {len(ids_out)}"
            )

            if ids_out:
                f.write(
                    " "
                    + " ".join(
                        str(v)
                        for v in ids_out
                    )
                )

            f.write("\n")


def write_bc(
    path: Path,
    bc_rows: Sequence[Tuple[int, int]],
    block_length: int,
    value: float,
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(
            f"{len(bc_rows)} {block_length}\n"
        )

        for node_id, block_id in bc_rows:
            f.write(
                f"{node_id} {block_id} {value:.16g}\n"
            )


# =============================================================================
# Boundary condition helpers
# =============================================================================


def parse_bc(
    spec: str,
    block_length: int,
) -> Tuple[str, List[int]]:
    if ":" not in spec:
        raise ValueError(
            f"Invalid --bc '{spec}'. "
            "Example: --bc ymin:1"
        )

    face, dof_text = spec.split(
        ":",
        1,
    )

    face = face.strip().lower()

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
            f"Unknown face '{face}'. "
            f"Choose from {sorted(allowed)}"
        )

    dofs = [
        int(s.strip())
        for s in dof_text.split(",")
        if s.strip()
    ]

    if not dofs:
        raise ValueError(
            f"No DOFs specified in '{spec}'"
        )

    for dof in dofs:
        if dof < 1 or dof > block_length:
            raise ValueError(
                f"DOF {dof} is outside "
                f"1..{block_length} in '{spec}'"
            )

    return (
        face,
        sorted(set(dofs)),
    )


def local_face_nodes(
    nodes: np.ndarray,
    face: str,
    tol: float,
) -> np.ndarray:
    axis = {
        "x": 0,
        "y": 1,
        "z": 2,
    }[face[0]]

    use_min = face.endswith(
        "min"
    )

    target = (
        nodes[:, axis].min()
        if use_min
        else nodes[:, axis].max()
    )

    return np.where(
        np.abs(
            nodes[:, axis]
            - target
        )
        <= tol
    )[0].astype(np.int64)


# =============================================================================
# Surface connectivity helpers
# =============================================================================


def extract_surface_node_ids(
    surface_elems: np.ndarray,
    base: int,
    nnode: int,
    name: str,
) -> np.ndarray:
    if surface_elems.size == 0:
        raise ValueError(
            f"{name}: surface connectivity is empty"
        )

    surface0 = (
        surface_elems
        - base
    )

    mn = int(
        surface0.min()
    )

    mx = int(
        surface0.max()
    )

    if mn < 0 or mx >= nnode:
        raise ValueError(
            f"{name}: node IDs are incompatible "
            "with volume node.dat/elem.dat. "
            f"zero-based range=[{mn}, {mx}], "
            f"nnode={nnode}, base={base}"
        )

    return np.unique(
        surface0.reshape(-1)
    ).astype(np.int64)


def read_surface_union(
    paths: Sequence[Path],
    base: int,
    nnode: int,
    label: str,
) -> np.ndarray:
    if not paths:
        raise ValueError(
            f"No surface file was specified for {label}"
        )

    union = set()

    for path in paths:
        surf, _ = read_elements(
            path
        )

        ids = extract_surface_node_ids(
            surf,
            base,
            nnode,
            f"{label}: {path}",
        )

        union.update(
            int(v)
            for v in ids
        )

    return np.array(
        sorted(union),
        dtype=np.int64,
    )


def compare_id_sets(
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

    msg = [
        f"{actual_name} does not match {expected_name}.",
        f"{actual_name} nodes={len(a)}, "
        f"{expected_name} nodes={len(e)}.",
    ]

    if missing:
        msg.append(
            "first missing local IDs="
            + str(missing[:10])
        )

    if extra:
        msg.append(
            "first extra local IDs="
            + str(extra[:10])
        )

    raise ValueError(
        " ".join(msg)
    )


# =============================================================================
# Standard repeated-copy welding
# =============================================================================


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


def build_pairing_from_two_surfaces(
    nodes: np.ndarray,
    left_elems: np.ndarray,
    right_elems: np.ndarray,
    base: int,
    axis_id: int,
    pitch: float,
    tol: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[int, int], float]:
    left_ids = extract_surface_node_ids(
        left_elems,
        base,
        len(nodes),
        "left interface",
    )

    right_ids = extract_surface_node_ids(
        right_elems,
        base,
        len(nodes),
        "right interface",
    )

    if len(left_ids) != len(right_ids):
        raise ValueError(
            "Cannot make conformal interface: "
            f"left nodes={len(left_ids)}, "
            f"right nodes={len(right_ids)}"
        )

    shift = np.zeros(
        3,
        dtype=float,
    )

    shift[axis_id] = pitch

    left_xyz = (
        nodes[left_ids]
        + shift
    )

    right_xyz = (
        nodes[right_ids]
    )

    lo = lexicographic_order(
        left_xyz
    )

    ro = lexicographic_order(
        right_xyz
    )

    left_sorted_ids = (
        left_ids[lo]
    )

    right_sorted_ids = (
        right_ids[ro]
    )

    err = np.max(
        np.abs(
            left_xyz[lo]
            - right_xyz[ro]
        ),
        axis=1,
    )

    max_err = (
        float(err.max())
        if len(err)
        else 0.0
    )

    if np.any(
        err > tol
    ):
        raise ValueError(
            "Interface coordinate pairing failed: "
            f"max error={max_err:.16g}, "
            f"weld_tol={tol:.16g}"
        )

    right_to_left = {
        int(r): int(l)
        for r, l in zip(
            right_sorted_ids,
            left_sorted_ids,
        )
    }

    return (
        left_ids,
        right_ids,
        right_to_left,
        max_err,
    )


def repeat_disconnected(
    nodes: np.ndarray,
    elems0: np.ndarray,
    copies: int,
    axis_id: int,
    pitch: float,
    start_offset: float,
):
    nnode = len(
        nodes
    )

    all_nodes = []
    all_elems = []
    copy_maps = []
    domain_nodes = []

    for k in range(copies):
        shift = np.zeros(
            3,
            dtype=float,
        )

        shift[axis_id] = (
            start_offset
            + k * pitch
        )

        nodes_k = (
            nodes
            + shift
        )

        node_offset = (
            k * nnode
        )

        local_to_global = (
            np.arange(
                nnode,
                dtype=np.int64,
            )
            + node_offset
        )

        all_nodes.append(
            nodes_k
        )

        all_elems.append(
            local_to_global[
                elems0
            ]
        )

        copy_maps.append(
            local_to_global
        )

        domain_nodes.append(
            local_to_global.tolist()
        )

    return (
        np.vstack(all_nodes),
        np.vstack(all_elems),
        copy_maps,
        domain_nodes,
    )


def repeat_continuous(
    nodes: np.ndarray,
    elems0: np.ndarray,
    copies: int,
    axis_id: int,
    pitch: float,
    start_offset: float,
    left_ids: np.ndarray,
    right_to_left: Dict[int, int],
):
    nnode = len(
        nodes
    )

    left_set = set(
        int(v)
        for v in left_ids
    )

    merged_nodes_list: List[np.ndarray] = []
    all_elems = []
    copy_maps = []
    domain_nodes = [
        []
        for _ in range(copies)
    ]

    previous_map = None

    for k in range(copies):
        shift = np.zeros(
            3,
            dtype=float,
        )

        shift[axis_id] = (
            start_offset
            + k * pitch
        )

        nodes_k = (
            nodes
            + shift
        )

        local_to_global = np.empty(
            nnode,
            dtype=np.int64,
        )

        if k == 0:
            for lid in range(nnode):
                gid = len(
                    merged_nodes_list
                )

                merged_nodes_list.append(
                    nodes_k[lid].copy()
                )

                local_to_global[lid] = gid

                domain_nodes[k].append(
                    gid
                )

        else:
            assert previous_map is not None

            for prev_right, this_left in right_to_left.items():
                local_to_global[this_left] = (
                    previous_map[
                        prev_right
                    ]
                )

            for lid in range(nnode):
                if lid in left_set:
                    continue

                gid = len(
                    merged_nodes_list
                )

                merged_nodes_list.append(
                    nodes_k[lid].copy()
                )

                local_to_global[lid] = gid

                domain_nodes[k].append(
                    gid
                )

        all_elems.append(
            local_to_global[
                elems0
            ]
        )

        copy_maps.append(
            local_to_global
        )

        previous_map = (
            local_to_global
        )

    return (
        np.asarray(
            merged_nodes_list,
            dtype=float,
        ),
        np.vstack(
            all_elems
        ),
        copy_maps,
        domain_nodes,
    )


# =============================================================================
# Symmetric-cell mode
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

# Orientation-reversing symmetry of the reference Hex8.
#
# 1 <-> 2
# 3 <-> 4
# 5 <-> 6
# 7 <-> 8
#
# zero-based:
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
    """
    Jacobian determinant at Hex8 reference center.
    Used only as an orientation sanity check.
    """

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


def mirrored_reference_component(
    nodes: np.ndarray,
    elems0: np.ndarray,
    axis_id: int,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Reflect original component around its repeat-axis minimum plane
    and fix Hex8 orientation.
    """

    if elems0.shape[1] != 8:
        raise ValueError(
            "--symmetric currently supports NEN=8 "
            "Gmsh hexahedra only."
        )

    axis_min = float(
        nodes[:, axis_id].min()
    )

    mirror_nodes = (
        nodes.copy()
    )

    mirror_nodes[:, axis_id] = (
        2.0 * axis_min
        - nodes[:, axis_id]
    )

    mirror_elems0 = (
        elems0[
            :,
            HEX8_MIRROR_PERM,
        ]
    )

    original_det = (
        hex8_center_jacobian_det(
            nodes,
            elems0,
        )
    )

    mirror_det = (
        hex8_center_jacobian_det(
            mirror_nodes,
            mirror_elems0,
        )
    )

    tiny = 1.0e-14

    if np.any(
        np.abs(original_det)
        <= tiny
    ):
        idx = int(
            np.argmin(
                np.abs(original_det)
            )
        )

        raise ValueError(
            "Original mesh contains a degenerate/near-zero "
            f"Hex8 center Jacobian. element={idx}, "
            f"det={original_det[idx]:.16g}"
        )

    if np.any(
        np.abs(mirror_det)
        <= tiny
    ):
        idx = int(
            np.argmin(
                np.abs(mirror_det)
            )
        )

        raise ValueError(
            "Mirrored mesh contains a degenerate/near-zero "
            f"Hex8 center Jacobian. element={idx}, "
            f"det={mirror_det[idx]:.16g}"
        )

    sign_mismatch = (
        np.sign(original_det)
        != np.sign(mirror_det)
    )

    if np.any(
        sign_mismatch
    ):
        idx = int(
            np.where(
                sign_mismatch
            )[0][0]
        )

        raise ValueError(
            "Mirrored Hex8 orientation check failed. "
            f"element={idx}, "
            f"original det={original_det[idx]:.16g}, "
            f"mirrored det={mirror_det[idx]:.16g}"
        )

    ratio = (
        mirror_det
        / original_det
    )

    orientation_error = float(
        np.max(
            np.abs(
                ratio - 1.0
            )
        )
    )

    return (
        mirror_nodes,
        mirror_elems0,
        orientation_error,
    )


def transform_original_half(
    nodes: np.ndarray,
    axis_id: int,
    axis_min: float,
    center_position: float,
) -> np.ndarray:
    """
    Original half:
        local axis_min -> center
        local axis_max -> right cell boundary
    """

    out = (
        nodes.copy()
    )

    u = (
        nodes[:, axis_id]
        - axis_min
    )

    out[:, axis_id] = (
        center_position
        + u
    )

    return out


def transform_mirror_half(
    nodes: np.ndarray,
    axis_id: int,
    axis_min: float,
    center_position: float,
) -> np.ndarray:
    """
    Mirror half:
        local axis_min -> center
        local axis_max -> left cell boundary
    """

    out = (
        nodes.copy()
    )

    u = (
        nodes[:, axis_id]
        - axis_min
    )

    out[:, axis_id] = (
        center_position
        - u
    )

    return out


def validate_symmetric_interfaces(
    nodes: np.ndarray,
    axis_id: int,
    center_ids: np.ndarray,
    gap_ids: np.ndarray,
    tol: float,
) -> Tuple[float, float]:
    """
    End-fin symmetric arrangement validation.

    Local component C:
        center_ids = repeat-axis MIN plane
                     thermal_fin: Gamma11 + Gamma12
        gap_ids    = repeat-axis MAX plane
                     thermal_fin: Gamma03

    Final arrangement:
        C == Mirror(C) == C == Mirror(C) ...

    Therefore:
      * inside one pitch:
            original GAP == mirror GAP
      * between adjacent pitches:
            previous mirror CENTER == next original CENTER

    The surface-connectivity node sets are also checked against the
    complete local min/max bounding-box planes.
    """

    axis_name = "xyz"[axis_id]

    bbox_center_ids = local_face_nodes(
        nodes,
        f"{axis_name}min",
        tol,
    )

    bbox_gap_ids = local_face_nodes(
        nodes,
        f"{axis_name}max",
        tol,
    )

    compare_id_sets(
        center_ids,
        bbox_center_ids,
        "center interface surface union",
        f"{axis_name}min bounding-box node set",
    )

    compare_id_sets(
        gap_ids,
        bbox_gap_ids,
        "gap interface surface union",
        f"{axis_name}max bounding-box node set",
    )

    axis_min = float(nodes[:, axis_id].min())
    axis_max = float(nodes[:, axis_id].max())
    half_length = axis_max - axis_min

    # Reference:
    #
    #   original C: center at x=0, gap at x=L
    #   mirror  M : gap at x=L, center at x=2L
    #   next C    : center at x=2L
    original0 = transform_original_half(
        nodes,
        axis_id,
        axis_min,
        0.0,
    )

    mirror0 = transform_mirror_half(
        nodes,
        axis_id,
        axis_min,
        2.0 * half_length,
    )

    original1 = transform_original_half(
        nodes,
        axis_id,
        axis_min,
        2.0 * half_length,
    )

    gap_err = float(
        np.max(
            np.abs(
                original0[gap_ids]
                - mirror0[gap_ids]
            )
        )
    )

    center_err = float(
        np.max(
            np.abs(
                mirror0[center_ids]
                - original1[center_ids]
            )
        )
    )

    if gap_err > tol:
        raise ValueError(
            "In-pitch gap interface coordinate pairing failed: "
            f"max error={gap_err:.16g}, "
            f"weld_tol={tol:.16g}"
        )

    if center_err > tol:
        raise ValueError(
            "Inter-pitch center interface coordinate pairing failed: "
            f"max error={center_err:.16g}, "
            f"weld_tol={tol:.16g}"
        )

    return (
        center_err,
        gap_err,
    )


def repeat_symmetric_cells(
    nodes: np.ndarray,
    elems0: np.ndarray,
    components: int,
    axis_id: int,
    start_offset: float,
    center_ids: np.ndarray,
    gap_ids: np.ndarray,
):
    """
    Build

        C == Mirror(C) == C == Mirror(C) == ...

    with ONE GEO half-component per GEDATSU domain.

    --copies / components therefore means:

        domain 0 : C
        domain 1 : Mirror(C)
        domain 2 : C
        domain 3 : Mirror(C)
        ...

    Adjacent interfaces alternate:

        C -> Mirror(C):
            GAP(max) == GAP(max)

        Mirror(C) -> next C:
            CENTER(min) == CENTER(min)

    Shared interface nodes are owned by the preceding domain.
    GEDATSU can then add them as overlap nodes for the following domain.

    For an even number of components, both global ends are CENTER/fin faces.

    cell_maps is returned as:
        [(C_map, Mirror_map), (C_map, Mirror_map), ...]
    only for BC generation; it does NOT define GEDATSU ownership.
    """

    if components < 2:
        raise ValueError(
            "--symmetric requires --copies >= 2"
        )

    if components % 2 != 0:
        raise ValueError(
            "--symmetric requires an even --copies value so that "
            "the array starts with C and ends with Mirror(C), leaving "
            "fin/center faces at both global ends."
        )

    nnode = len(nodes)

    axis_min = float(nodes[:, axis_id].min())
    axis_max = float(nodes[:, axis_id].max())
    half_length = axis_max - axis_min

    (
        _mirror_reference_nodes,
        mirror_elems0,
        orientation_error,
    ) = mirrored_reference_component(
        nodes,
        elems0,
        axis_id,
    )

    center_set = set(int(v) for v in center_ids)
    gap_set = set(int(v) for v in gap_ids)

    overlap = center_set & gap_set

    if overlap:
        raise ValueError(
            "center interface and gap interface share local node IDs. "
            "They must be opposite, disjoint planes. "
            f"Example shared IDs={sorted(overlap)[:10]}"
        )

    merged_nodes_list: List[np.ndarray] = []
    all_elems = []

    # Kept as C/M pairs for BC construction.
    cell_maps: List[Tuple[np.ndarray, np.ndarray]] = []

    # IMPORTANT:
    # one domain per GEO half-component.
    domain_nodes: List[List[int]] = [
        []
        for _ in range(components)
    ]

    previous_map = None
    pending_original_map = None

    for comp in range(components):
        is_original = (
            comp % 2 == 0
        )

        local_to_global = np.empty(
            nnode,
            dtype=np.int64,
        )

        if is_original:
            # C occupies:
            #   [start + comp*L, start + (comp+1)*L]
            #
            # its local MIN/center face is at the left end.
            center_position = (
                start_offset
                + comp * half_length
            )

            xyz = transform_original_half(
                nodes,
                axis_id,
                axis_min,
                center_position,
            )

            if comp == 0:
                # First component: everything is new and owned by domain 0.
                for lid in range(nnode):
                    gid = len(
                        merged_nodes_list
                    )

                    merged_nodes_list.append(
                        xyz[lid].copy()
                    )

                    local_to_global[lid] = gid
                    domain_nodes[comp].append(gid)

            else:
                assert previous_map is not None

                # Previous component is Mirror(C).
                #
                # previous Mirror center/min
                #        ==
                # current C center/min
                #
                # Shared nodes remain owned by previous domain.
                for lid in center_ids:
                    lid_i = int(lid)

                    local_to_global[lid_i] = (
                        previous_map[lid_i]
                    )

                for lid in range(nnode):
                    if lid in center_set:
                        continue

                    gid = len(
                        merged_nodes_list
                    )

                    merged_nodes_list.append(
                        xyz[lid].copy()
                    )

                    local_to_global[lid] = gid
                    domain_nodes[comp].append(gid)

            all_elems.append(
                local_to_global[elems0]
            )

            pending_original_map = (
                local_to_global
            )

        else:
            # Mirror(C) occupies:
            #   [start + comp*L, start + (comp+1)*L]
            #
            # its local MAX/gap face is at its LEFT end,
            # and local MIN/center face is at its RIGHT end.
            right_center_position = (
                start_offset
                + (comp + 1) * half_length
            )

            xyz = transform_mirror_half(
                nodes,
                axis_id,
                axis_min,
                right_center_position,
            )

            assert previous_map is not None
            assert pending_original_map is not None

            # Previous component is C:
            #
            # previous C gap/max
            #       ==
            # current Mirror gap/max
            #
            # Shared nodes stay owned by previous C domain.
            for lid in gap_ids:
                lid_i = int(lid)

                local_to_global[lid_i] = (
                    previous_map[lid_i]
                )

            for lid in range(nnode):
                if lid in gap_set:
                    continue

                gid = len(
                    merged_nodes_list
                )

                merged_nodes_list.append(
                    xyz[lid].copy()
                )

                local_to_global[lid] = gid
                domain_nodes[comp].append(gid)

            all_elems.append(
                local_to_global[mirror_elems0]
            )

            # Complete one physical C+Mirror pitch for BC processing.
            cell_maps.append(
                (
                    pending_original_map,
                    local_to_global,
                )
            )

            pending_original_map = None

        previous_map = local_to_global

    return (
        np.asarray(
            merged_nodes_list,
            dtype=float,
        ),
        np.vstack(all_elems),
        cell_maps,
        domain_nodes,
        orientation_error,
    )


# =============================================================================
# BC construction
# =============================================================================


def build_regular_bc(
    nodes: np.ndarray,
    copy_maps: Sequence[np.ndarray],
    bc_specs: Sequence[Tuple[str, Sequence[int]]],
    base: int,
    tol: float,
    continuous: bool,
    axis_id: int,
) -> List[Tuple[int, int]]:
    face_cache: Dict[str, np.ndarray] = {}

    for face, _ in bc_specs:
        if face not in face_cache:
            face_cache[
                face
            ] = local_face_nodes(
                nodes,
                face,
                tol,
            )

    axis_name = (
        "xyz"[axis_id]
    )

    min_face = (
        f"{axis_name}min"
    )

    max_face = (
        f"{axis_name}max"
    )

    result = set()

    for face, dofs in bc_specs:
        local_ids = (
            face_cache[
                face
            ]
        )

        face_axis = {
            "x": 0,
            "y": 1,
            "z": 2,
        }[face[0]]

        for k, local_to_global in enumerate(copy_maps):
            if continuous and face_axis == axis_id:
                if (
                    face == min_face
                    and k != 0
                ):
                    continue

                if (
                    face == max_face
                    and k != len(copy_maps) - 1
                ):
                    continue

            for lid in local_ids:
                gid = int(
                    local_to_global[
                        int(lid)
                    ]
                )

                for dof in dofs:
                    result.add(
                        (
                            gid + base,
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


def build_symmetric_bc(
    nodes: np.ndarray,
    cell_maps: Sequence[Tuple[np.ndarray, np.ndarray]],
    bc_specs: Sequence[Tuple[str, Sequence[int]]],
    base: int,
    tol: float,
    axis_id: int,
) -> List[Tuple[int, int]]:
    """
    Bounding-box BC for the end-fin arrangement

        C == Mirror(C) == C == Mirror(C) ...

    cell_maps contains:
        (original_map, mirror_map)

    Along the repeat axis:
        global MIN = first original C local MIN (center/fin face)
        global MAX = last mirror C local MIN (center/fin face)

    All internal center/gap cuts are excluded.

    Along transverse axes the same local face is applied to every C/M half.
    """

    axis_name = "xyz"[axis_id]
    min_axis_face = f"{axis_name}min"
    max_axis_face = f"{axis_name}max"

    local_center_ids = local_face_nodes(
        nodes,
        min_axis_face,
        tol,
    )

    face_cache: Dict[str, np.ndarray] = {}

    for face, _ in bc_specs:
        face_axis = {
            "x": 0,
            "y": 1,
            "z": 2,
        }[face[0]]

        if face_axis != axis_id:
            face_cache[face] = local_face_nodes(
                nodes,
                face,
                tol,
            )

    result = set()

    for face, dofs in bc_specs:
        face_axis = {
            "x": 0,
            "y": 1,
            "z": 2,
        }[face[0]]

        if face_axis == axis_id:
            if face == min_axis_face:
                original_map = cell_maps[0][0]

                for lid in local_center_ids:
                    gid = int(original_map[int(lid)])

                    for dof in dofs:
                        result.add(
                            (
                                gid + base,
                                int(dof),
                            )
                        )

            elif face == max_axis_face:
                mirror_map = cell_maps[-1][1]

                for lid in local_center_ids:
                    gid = int(mirror_map[int(lid)])

                    for dof in dofs:
                        result.add(
                            (
                                gid + base,
                                int(dof),
                            )
                        )

            continue

        local_ids = face_cache[face]

        for original_map, mirror_map in cell_maps:
            for local_to_global in (
                original_map,
                mirror_map,
            ):
                for lid in local_ids:
                    gid = int(
                        local_to_global[int(lid)]
                    )

                    for dof in dofs:
                        result.add(
                            (
                                gid + base,
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


def parse_dof_list(
    text: str,
    block_length: int,
    option_name: str,
) -> List[int]:
    dofs = [
        int(s.strip())
        for s in text.split(",")
        if s.strip()
    ]

    if not dofs:
        raise ValueError(
            f"{option_name}: no DOFs were specified"
        )

    dofs = sorted(set(dofs))

    for dof in dofs:
        if dof < 1 or dof > block_length:
            raise ValueError(
                f"{option_name}: DOF {dof} is outside "
                f"1..{block_length}"
            )

    return dofs


def read_surface_faces_union(
    paths: Sequence[Path],
    base: int,
    nnode: int,
    label: str,
) -> np.ndarray:
    """
    Read one or more surface-connectivity files and return unique
    zero-based surface elements.

    Duplicate faces are removed using a node-order-independent key.
    """

    if not paths:
        return np.empty(
            (0, 0),
            dtype=np.int64,
        )

    unique_rows = {}
    nen_ref = None

    for path in paths:
        surf, nen = read_elements(path)

        if nen_ref is None:
            nen_ref = nen
        elif nen != nen_ref:
            raise ValueError(
                f"{label}: mixed surface NEN values are not supported: "
                f"{nen_ref} and {nen} ({path})"
            )

        surf0 = surf - base

        if surf0.size:
            mn = int(surf0.min())
            mx = int(surf0.max())

            if mn < 0 or mx >= nnode:
                raise ValueError(
                    f"{label}: {path}: surface node IDs are outside "
                    f"0..{nnode - 1}: range=[{mn}, {mx}]"
                )

        for row in surf0:
            key = tuple(
                sorted(int(v) for v in row)
            )

            if key not in unique_rows:
                unique_rows[key] = np.asarray(
                    row,
                    dtype=np.int64,
                ).copy()

    if not unique_rows:
        return np.empty(
            (0, nen_ref or 0),
            dtype=np.int64,
        )

    return np.vstack(
        list(unique_rows.values())
    )


def surface_face_key_set(
    faces0: np.ndarray,
) -> set:
    return {
        tuple(sorted(int(v) for v in row))
        for row in faces0
    }


def build_symmetric_surface_bc(
    cell_maps: Sequence[Tuple[np.ndarray, np.ndarray]],
    boundary_faces0: np.ndarray,
    center_faces0: np.ndarray,
    gap_faces0: np.ndarray,
    base: int,
    dofs: Sequence[int],
) -> List[Tuple[int, int]]:
    """
    Build BC from actual boundary surface elements for

        C == Mirror(C) == C == Mirror(C) ...

    The candidate boundary_faces0 should normally be
    All_external_boundaries_quad_connectivity.dat.

    Artificial cuts that are welded are removed at the FACE level:
      * gap/Gamma03 is internal for every C and every Mirror(C)
      * center/Gamma11+Gamma12 is internal except:
            - first C center face  -> global left end
            - last Mirror center   -> global right end

    Face-level removal is important: edge/corner nodes shared with another
    true external surface remain in the BC node set.
    """

    center_keys = surface_face_key_set(
        center_faces0
    )

    gap_keys = surface_face_key_set(
        gap_faces0
    )

    bc_node_ids = set()
    last_cell = len(cell_maps) - 1

    for cell, (
        original_map,
        mirror_map,
    ) in enumerate(cell_maps):

        # ----------------------------------------------------------
        # Original C
        # ----------------------------------------------------------
        original_center_is_external = (
            cell == 0
        )

        for row in boundary_faces0:
            key = tuple(
                sorted(int(v) for v in row)
            )

            # C max / Gamma03 is welded to Mirror(C) max.
            if key in gap_keys:
                continue

            # C min is external only at the global left end.
            if (
                key in center_keys
                and not original_center_is_external
            ):
                continue

            for lid in row:
                bc_node_ids.add(
                    int(
                        original_map[
                            int(lid)
                        ]
                    )
                )

        # ----------------------------------------------------------
        # Mirror(C)
        # ----------------------------------------------------------
        mirror_center_is_external = (
            cell == last_cell
        )

        for row in boundary_faces0:
            key = tuple(
                sorted(int(v) for v in row)
            )

            # Mirror max / Gamma03 is welded to C max.
            if key in gap_keys:
                continue

            # Mirror min is external only at the global right end.
            if (
                key in center_keys
                and not mirror_center_is_external
            ):
                continue

            for lid in row:
                bc_node_ids.add(
                    int(
                        mirror_map[
                            int(lid)
                        ]
                    )
                )

    rows = []

    for gid in sorted(bc_node_ids):
        for dof in dofs:
            rows.append(
                (
                    gid + base,
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
    mn = int(
        elems0.min()
    )

    mx = int(
        elems0.max()
    )

    if mn < 0 or mx >= nnode:
        raise RuntimeError(
            "Generated element connectivity contains "
            "an invalid node ID: "
            f"range=[{mn}, {mx}], nnode={nnode}"
        )


def validate_domain_ownership(
    domain_nodes: Sequence[Sequence[int]],
    nnode: int,
) -> None:
    owner = np.full(
        nnode,
        -1,
        dtype=np.int64,
    )

    for domain_id, ids in enumerate(domain_nodes):
        for gid in ids:
            gid_i = int(
                gid
            )

            if owner[
                gid_i
            ] != -1:
                raise RuntimeError(
                    f"Global node {gid_i} is owned by more "
                    "than one domain: "
                    f"{int(owner[gid_i])} and {domain_id}"
                )

            owner[
                gid_i
            ] = domain_id

    missing = np.where(
        owner < 0
    )[0]

    if len(missing):
        raise RuntimeError(
            f"{len(missing)} generated nodes have no "
            "given_subdomain owner. "
            f"First missing={int(missing[0])}"
        )


# =============================================================================
# main
# =============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Repeat a component mesh, optionally as left-right symmetric "
            "continuous cells."
        )
    )

    parser.add_argument(
        "node_dat",
        nargs="?",
        type=Path,
        help=(
            "Existing component node.dat. "
            "Omit when --geo is used."
        ),
    )

    parser.add_argument(
        "elem_dat",
        nargs="?",
        type=Path,
        help=(
            "Existing component elem.dat. "
            "Omit when --geo is used."
        ),
    )

    # -----------------------------------------------------------------
    # Component source
    # -----------------------------------------------------------------

    parser.add_argument(
        "--geo",
        type=Path,
        default=None,
        help=(
            "Generate the component directly from this Gmsh .geo file. "
            "When specified, positional node_dat/elem_dat are not required."
        ),
    )

    parser.add_argument(
        "--gmsh",
        type=str,
        default="gmsh",
        help=(
            "Gmsh executable name/path. Default: gmsh"
        ),
    )

    parser.add_argument(
        "--physical-group-script",
        type=Path,
        default=Path("mesh_io/save_physical_groups.py"),
        help=(
            "Existing script that extracts Physical Groups from the generated "
            "MSH file. Default: mesh_io/save_physical_groups.py"
        ),
    )

    parser.add_argument(
        "--component-work-dir",
        type=Path,
        default=Path("generated_component"),
        help=(
            "Work directory for the Gmsh-generated half component. "
            "Default: generated_component"
        ),
    )

    parser.add_argument(
        "--volume-group",
        type=str,
        default="Fluid",
        help=(
            "Physical Volume used for node/hex connectivity. Default: Fluid"
        ),
    )

    parser.add_argument(
        "--center-group",
        action="append",
        default=[],
        metavar="PHYSICAL_SURFACE",
        help=(
            "Physical Surface forming the fin-center side of the half component. "
            "Can be repeated. GEO-mode default: Gamma11 and Gamma12."
        ),
    )

    parser.add_argument(
        "--gap-group",
        action="append",
        default=[],
        metavar="PHYSICAL_SURFACE",
        help=(
            "Physical Surface forming the half-pitch gap side. "
            "Can be repeated. GEO-mode default: Gamma03."
        ),
    )

    parser.add_argument(
        "--bc-group",
        type=str,
        default="All_external_boundaries",
        help=(
            "Physical Surface group used to generate D_bc.dat in GEO mode. "
            "Default: All_external_boundaries"
        ),
    )

    parser.add_argument(
        "--no-auto-bc",
        action="store_true",
        help=(
            "In GEO mode, do not automatically use --bc-group for D_bc.dat."
        ),
    )

    parser.add_argument(
        "--no-clean-component-work",
        action="store_true",
        help=(
            "Do not remove the generated mesh_tmp directory before GEO meshing."
        ),
    )

    parser.add_argument(
        "--copies",
        type=int,
        required=True,
        help=(
            "Normal mode: number of component copies. "
            "With --symmetric: number of GEO half-components, and also "
            "the number of GEDATSU domains. One domain = one component. "
            "Use an even number so both global ends are fin/center faces."
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
    )

    spacing = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    spacing.add_argument(
        "--gap",
        type=float,
        help=(
            "Normal mode: component bbox gap. "
            "With --symmetric it must be 0."
        ),
    )

    spacing.add_argument(
        "--pitch",
        type=float,
        help=(
            "Normal mode: translation pitch. "
            "Not supported with --symmetric."
        ),
    )

    parser.add_argument(
        "--start-offset",
        type=float,
        default=0.0,
        help=(
            "Normal mode: translation of first copy. "
            "Symmetric mode: global coordinate of the leftmost cell boundary "
            "along the repeat axis."
        ),
    )

    parser.add_argument(
        "--input-base",
        choices=(
            "auto",
            "0",
            "1",
        ),
        default="auto",
    )

    # -----------------------------------------------------------------
    # Normal continuous mode
    # -----------------------------------------------------------------

    parser.add_argument(
        "--continuous",
        action="store_true",
        help=(
            "Normal repeated-copy mode: weld explicit left/right interfaces."
        ),
    )

    parser.add_argument(
        "--interface-left",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--interface-right",
        type=Path,
        default=None,
    )

    # -----------------------------------------------------------------
    # Symmetric mode
    # -----------------------------------------------------------------

    parser.add_argument(
        "--symmetric",
        action="store_true",
        help=(
            "Build each domain/pitch as component+Mirror(component). "
            "The global left and right ends are fin/center-side faces."
        ),
    )

    parser.add_argument(
        "--center-interface",
        action="append",
        type=Path,
        default=[],
        metavar="SURFACE_DAT",
        help=(
            "Surface connectivity belonging to the symmetry-center interface. "
            "Can be repeated. thermal_fin: Gamma11 + Gamma12."
        ),
    )

    parser.add_argument(
        "--gap-interface",
        action="append",
        type=Path,
        default=[],
        metavar="SURFACE_DAT",
        help=(
            "Surface connectivity belonging to the outer/cell interface. "
            "Can be repeated. thermal_fin: Gamma03."
        ),
    )

    parser.add_argument(
        "--weld-tol",
        type=float,
        default=1.0e-10,
    )

    # -----------------------------------------------------------------
    # BC
    # -----------------------------------------------------------------

    parser.add_argument(
        "--block-length",
        type=int,
        default=1,
        help=(
            "D_bc.dat block length. Default=1, suitable for scalar convdiff. "
            "Use 3 explicitly for a 3-DOF vector problem."
        ),
    )

    parser.add_argument(
        "--bc",
        action="append",
        default=[],
        metavar="FACE:DOFS",
        help=(
            "Bounding-box Dirichlet BC. "
            "Example: --bc ymin:1. Can be repeated."
        ),
    )

    parser.add_argument(
        "--bc-surface",
        action="append",
        type=Path,
        default=[],
        metavar="SURFACE_DAT",
        help=(
            "Surface connectivity used to create Dirichlet BC nodes. "
            "Can be repeated. In --symmetric mode, welded center/gap "
            "surface elements are automatically removed. "
            "Recommended: All_external_boundaries_quad_connectivity.dat"
        ),
    )

    parser.add_argument(
        "--bc-surface-dofs",
        type=str,
        default="1",
        help=(
            "Comma-separated DOFs applied to --bc-surface nodes. "
            "Default: 1"
        ),
    )

    parser.add_argument(
        "--bc-tol",
        type=float,
        default=1.0e-10,
    )

    parser.add_argument(
        "--bc-value",
        type=float,
        default=0.0,
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(
            "merged_mesh"
        ),
    )

    args = parser.parse_args()

    # -----------------------------------------------------------------
    # Resolve component source
    # -----------------------------------------------------------------

    component_geo_info = None

    if args.geo is not None:
        if args.node_dat is not None or args.elem_dat is not None:
            raise ValueError(
                "Use either --geo OR positional node_dat/elem_dat, not both."
            )

        center_groups = (
            args.center_group
            if args.center_group
            else ["Gamma11", "Gamma12"]
        )

        gap_groups = (
            args.gap_group
            if args.gap_group
            else ["Gamma03"]
        )

        auto_bc_group = (
            None
            if args.no_auto_bc
            else args.bc_group
        )

        component_geo_info = generate_component_from_geo(
            geo_path=args.geo,
            gmsh_executable=args.gmsh,
            physical_group_script=args.physical_group_script,
            component_work_dir=args.component_work_dir,
            volume_group=args.volume_group,
            center_groups=center_groups,
            gap_groups=gap_groups,
            bc_group=auto_bc_group,
            clean_work=not args.no_clean_component_work,
        )

        args.node_dat = component_geo_info[
            "node_dat"
        ]

        args.elem_dat = component_geo_info[
            "elem_dat"
        ]

        # In GEO mode the Physical Groups define the intended interfaces.
        if not args.center_interface:
            args.center_interface = list(
                component_geo_info[
                    "center_interface"
                ]
            )

        if not args.gap_interface:
            args.gap_interface = list(
                component_geo_info[
                    "gap_interface"
                ]
            )

        # Preserve explicit --bc-surface if the user supplied one.
        # Otherwise use All_external_boundaries from the GEO model.
        if (
            not args.bc_surface
            and component_geo_info[
                "bc_surface"
            ]
        ):
            args.bc_surface = list(
                component_geo_info[
                    "bc_surface"
                ]
            )

    else:
        if args.node_dat is None or args.elem_dat is None:
            raise ValueError(
                "Specify either --geo GEO_FILE or both positional "
                "node_dat and elem_dat."
            )

        args.node_dat = require_file(
            args.node_dat,
            "component node.dat",
        )

        args.elem_dat = require_file(
            args.elem_dat,
            "component elem.dat",
        )

    # -----------------------------------------------------------------
    # Argument validation
    # -----------------------------------------------------------------

    if args.copies < 1:
        raise ValueError(
            "--copies must be >= 1"
        )

    if args.block_length < 1:
        raise ValueError(
            "--block-length must be >= 1"
        )

    if args.bc_tol < 0:
        raise ValueError(
            "--bc-tol must be >= 0"
        )

    if args.weld_tol <= 0:
        raise ValueError(
            "--weld-tol must be > 0"
        )

    if args.symmetric and args.continuous:
        raise ValueError(
            "Do not combine --symmetric and --continuous. "
            "--symmetric is already a continuous welded mode."
        )

    if args.bc_surface and not args.symmetric:
        raise ValueError(
            "--bc-surface is currently supported with --symmetric mode. "
            "For the thermal-fin array use --symmetric."
        )

    if args.symmetric:
        if args.copies < 2:
            raise ValueError(
                "--symmetric requires --copies >= 2"
            )

        if args.copies % 2 != 0:
            raise ValueError(
                "--symmetric requires an even --copies value. "
                "One GEO component = one GEDATSU domain, and an even count "
                "is needed to leave fin/center faces at both global ends."
            )

        if args.pitch is not None:
            raise ValueError(
                "--pitch is not used with --symmetric. "
                "Use --gap 0.0."
            )

        if args.gap is None:
            raise ValueError(
                "--symmetric requires --gap 0.0"
            )

        if abs(
            float(args.gap)
        ) > args.weld_tol:
            raise ValueError(
                "--symmetric creates a conformal continuous structure, "
                "so --gap must be 0.0. "
                f"Got gap={args.gap}"
            )

        if not args.center_interface:
            raise ValueError(
                "--symmetric requires one or more "
                "--center-interface files."
            )

        if not args.gap_interface:
            raise ValueError(
                "--symmetric requires one or more "
                "--gap-interface files."
            )

    else:
        if args.continuous:
            if (
                args.interface_left is None
                or args.interface_right is None
            ):
                raise ValueError(
                    "--continuous requires "
                    "--interface-left and --interface-right"
                )

    # -----------------------------------------------------------------
    # Input
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

    if nnode == 0:
        raise ValueError(
            "node.dat has zero nodes"
        )

    if nelem == 0:
        raise ValueError(
            "elem.dat has zero elements"
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

    axis_id = {
        "x": 0,
        "y": 1,
        "z": 2,
    }[args.axis]

    xyz_min = (
        nodes.min(
            axis=0
        )
    )

    xyz_max = (
        nodes.max(
            axis=0
        )
    )

    component_length = float(
        xyz_max[axis_id]
        - xyz_min[axis_id]
    )

    if component_length <= 0:
        raise ValueError(
            f"Component length along {args.axis} "
            f"is invalid: {component_length}"
        )

    print(
        f"Detected connectivity base: "
        f"{base}-based"
    )

    # -----------------------------------------------------------------
    # Build mesh
    # -----------------------------------------------------------------

    interface_count_center = 0
    interface_count_gap = 0
    pair_error_center = 0.0
    pair_error_gap = 0.0
    orientation_error = 0.0

    if args.symmetric:
        center_ids = read_surface_union(
            args.center_interface,
            base,
            nnode,
            "center interface",
        )

        gap_ids = read_surface_union(
            args.gap_interface,
            base,
            nnode,
            "gap interface",
        )

        (
            pair_error_center,
            pair_error_gap,
        ) = validate_symmetric_interfaces(
            nodes,
            axis_id,
            center_ids,
            gap_ids,
            args.weld_tol,
        )

        (
            merged_nodes,
            merged_elems0,
            cell_maps,
            domain_nodes,
            orientation_error,
        ) = repeat_symmetric_cells(
            nodes,
            elems0,
            args.copies,
            axis_id,
            args.start_offset,
            center_ids,
            gap_ids,
        )

        interface_count_center = (
            len(center_ids)
        )

        interface_count_gap = (
            len(gap_ids)
        )

        copy_maps = None

        physical_pitch = (
            2.0
            * component_length
        )

    else:
        if args.pitch is not None:
            physical_pitch = float(
                args.pitch
            )

            gap = (
                physical_pitch
                - component_length
            )

        else:
            gap = float(
                args.gap
            )

            physical_pitch = (
                component_length
                + gap
            )

        if physical_pitch <= 0:
            raise ValueError(
                f"pitch must be positive: {physical_pitch}"
            )

        if args.continuous:
            left_surf, _ = read_elements(
                args.interface_left
            )

            right_surf, _ = read_elements(
                args.interface_right
            )

            (
                left_ids,
                right_ids,
                right_to_left,
                max_pair_error,
            ) = build_pairing_from_two_surfaces(
                nodes,
                left_surf,
                right_surf,
                base,
                axis_id,
                physical_pitch,
                args.weld_tol,
            )

            (
                merged_nodes,
                merged_elems0,
                copy_maps,
                domain_nodes,
            ) = repeat_continuous(
                nodes,
                elems0,
                args.copies,
                axis_id,
                physical_pitch,
                args.start_offset,
                left_ids,
                right_to_left,
            )

            interface_count_center = (
                len(left_ids)
            )

            pair_error_center = (
                max_pair_error
            )

        else:
            (
                merged_nodes,
                merged_elems0,
                copy_maps,
                domain_nodes,
            ) = repeat_disconnected(
                nodes,
                elems0,
                args.copies,
                axis_id,
                physical_pitch,
                args.start_offset,
            )

            cell_maps = None

    # -----------------------------------------------------------------
    # Validate generated mesh
    # -----------------------------------------------------------------

    validate_connectivity(
        merged_elems0,
        len(merged_nodes),
    )

    validate_domain_ownership(
        domain_nodes,
        len(merged_nodes),
    )

    # -----------------------------------------------------------------
    # Boundary conditions
    # -----------------------------------------------------------------

    bc_specs = [
        parse_bc(
            spec,
            args.block_length,
        )
        for spec in args.bc
    ]

    bc_surface_dofs = parse_dof_list(
        args.bc_surface_dofs,
        args.block_length,
        "--bc-surface-dofs",
    )

    bc_row_set = set()

    # Optional bounding-box BC.
    if args.symmetric:
        if bc_specs:
            for row in build_symmetric_bc(
                nodes,
                cell_maps,
                bc_specs,
                base,
                args.bc_tol,
                axis_id,
            ):
                bc_row_set.add(row)
    else:
        if bc_specs:
            for row in build_regular_bc(
                nodes,
                copy_maps,
                bc_specs,
                base,
                args.bc_tol,
                args.continuous,
                axis_id,
            ):
                bc_row_set.add(row)

    # Actual external-surface BC.
    #
    # The interface surface files are read as FACE connectivity so that
    # internal faces can be removed without accidentally removing edge
    # nodes that also belong to a true external surface.
    if args.bc_surface:
        boundary_faces0 = read_surface_faces_union(
            args.bc_surface,
            base,
            nnode,
            "BC surface",
        )

        center_faces0 = read_surface_faces_union(
            args.center_interface,
            base,
            nnode,
            "center interface faces",
        )

        gap_faces0 = read_surface_faces_union(
            args.gap_interface,
            base,
            nnode,
            "gap interface faces",
        )

        for row in build_symmetric_surface_bc(
            cell_maps,
            boundary_faces0,
            center_faces0,
            gap_faces0,
            base,
            bc_surface_dofs,
        ):
            bc_row_set.add(row)

    bc_rows = sorted(
        bc_row_set,
        key=lambda x: (
            x[0],
            x[1],
        ),
    )

    generate_bc = bool(
        bc_specs
        or args.bc_surface
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
    else:
        # Do not leave a stale BC file from a previous run.
        if bc_out.exists():
            bc_out.unlink()

    # -----------------------------------------------------------------
    # Counts / report
    # -----------------------------------------------------------------

    if args.symmetric:
        n_full_pitches = (
            args.copies
            // 2
        )

        expected_nodes = (
            args.copies
            * nnode
            - n_full_pitches
            * interface_count_gap
            - (n_full_pitches - 1)
            * interface_count_center
        )

        expected_elements = (
            args.copies
            * nelem
        )

        if len(merged_nodes) != expected_nodes:
            raise RuntimeError(
                "Symmetric node count sanity check failed: "
                f"actual={len(merged_nodes)}, "
                f"expected={expected_nodes}"
            )

        if len(merged_elems0) != expected_elements:
            raise RuntimeError(
                "Symmetric element count sanity check failed: "
                f"actual={len(merged_elems0)}, "
                f"expected={expected_elements}"
            )

        disconnected_equivalent_nodes = (
            args.copies
            * nnode
        )

    else:
        expected_elements = (
            args.copies
            * nelem
        )

        disconnected_equivalent_nodes = (
            args.copies
            * nnode
        )

    welded_duplicates = (
        disconnected_equivalent_nodes
        - len(merged_nodes)
    )

    with report_out.open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write(
            "Repeated mesh report\n"
        )
        f.write(
            "====================\n"
        )

        if component_geo_info is not None:
            f.write(
                f"component source        : GEO\n"
            )
            f.write(
                f"component geo           : {component_geo_info['geo']}\n"
            )
            f.write(
                f"component msh           : {component_geo_info['msh']}\n"
            )
            f.write(
                f"component work dir      : {component_geo_info['component_work_dir']}\n"
            )
        else:
            f.write(
                f"component source        : node.dat / elem.dat\n"
            )

        f.write(
            f"input node file         : {args.node_dat}\n"
        )
        f.write(
            f"input elem file         : {args.elem_dat}\n"
        )
        f.write(
            f"input nodes             : {nnode}\n"
        )
        f.write(
            f"input elements          : {nelem}\n"
        )
        f.write(
            f"nodes/element           : {nen}\n"
        )
        f.write(
            f"index base              : {base}\n"
        )
        f.write(
            f"repeat axis             : {args.axis}\n"
        )
        f.write(
            f"component bbox min      : {xyz_min.tolist()}\n"
        )
        f.write(
            f"component bbox max      : {xyz_max.tolist()}\n"
        )
        f.write(
            f"half component length   : {component_length:.16g}\n"
        )
        f.write(
            f"components/domains      : {args.copies}\n"
        )
        if args.symmetric:
            f.write(
                f"full C+Mirror pitches  : {args.copies // 2}\n"
            )
        f.write(
            f"symmetric end-fin mode  : {args.symmetric}\n"
        )
        f.write(
            f"continuous regular mode : {args.continuous}\n"
        )
        f.write(
            f"effective pitch         : {physical_pitch:.16g}\n"
        )
        f.write(
            f"start offset            : {args.start_offset:.16g}\n"
        )

        if args.symmetric:
            f.write(
                f"center interface files  : "
                f"{[str(p) for p in args.center_interface]}\n"
            )
            f.write(
                f"gap interface files     : "
                f"{[str(p) for p in args.gap_interface]}\n"
            )
            f.write(
                f"center interface nodes  : {interface_count_center}\n"
            )
            f.write(
                f"gap interface nodes     : {interface_count_gap}\n"
            )
            f.write(
                f"center pair error       : {pair_error_center:.16g}\n"
            )
            f.write(
                f"gap pair error          : {pair_error_gap:.16g}\n"
            )
            f.write(
                f"mirror orientation error: {orientation_error:.16g}\n"
            )

        f.write(
            f"welded duplicate nodes  : {welded_duplicates}\n"
        )
        f.write(
            f"output nodes            : {len(merged_nodes)}\n"
        )
        f.write(
            f"output elements         : {len(merged_elems0)}\n"
        )
        f.write(
            f"given subdomain file    : {subdomain_out}\n"
        )

        f.write(
            "\nowned nodes/domain\n"
        )

        for domain_id, ids in enumerate(domain_nodes):
            f.write(
                f"  domain {domain_id}: "
                f"{len(ids)}\n"
            )

        if generate_bc:
            f.write(
                "\nBC\n"
            )
            f.write(
                f"block length            : {args.block_length}\n"
            )
            f.write(
                f"bbox specs              : {args.bc}\n"
            )
            f.write(
                f"surface files           : {[str(p) for p in args.bc_surface]}\n"
            )
            f.write(
                f"surface DOFs            : {bc_surface_dofs}\n"
            )
            f.write(
                f"value                   : {args.bc_value:.16g}\n"
            )
            f.write(
                f"D_bc rows               : {len(bc_rows)}\n"
            )
            f.write(
                f"D_bc file               : {bc_out}\n"
            )
        else:
            f.write(
                "\nD_bc                    : not generated\n"
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
            f"-n {args.copies} "
            "-in node.dat "
            "-ie elem.dat "
            "--given_subdomain given_subdomain.dat "
            "-d parted.0\n"
        )

    # -----------------------------------------------------------------
    # Console
    # -----------------------------------------------------------------

    print()
    print(
        "Done."
    )
    print()

    print(
        f"mode                : "
        f"{'symmetric GEO components (1 component/domain)' if args.symmetric else 'normal repetition'}"
    )
    print(
        f"component nodes     : {nnode}"
    )
    print(
        f"component elements  : {nelem}"
    )
    print(
        f"components/domains  : {args.copies}"
    )
    if args.symmetric:
        print(
            f"full C+Mirror pitches: {args.copies // 2}"
        )

    if args.symmetric:
        print(
            f"center interface    : {interface_count_center} nodes"
        )
        print(
            f"gap interface       : {interface_count_gap} nodes"
        )
        print(
            f"full fin pitch      : {physical_pitch}"
        )
        print(
            f"center pair error   : {pair_error_center}"
        )
        print(
            f"gap pair error      : {pair_error_gap}"
        )
        print(
            f"orientation error   : {orientation_error}"
        )

    print(
        f"welded duplicates   : {welded_duplicates}"
    )
    print(
        f"output nodes        : {len(merged_nodes)}"
    )
    print(
        f"output elements     : {len(merged_elems0)}"
    )

    print()
    print(
        f"node.dat            : {node_out}"
    )
    print(
        f"elem.dat            : {elem_out}"
    )
    print(
        f"given_subdomain.dat : {subdomain_out}"
    )

    if generate_bc:
        print(
            f"D_bc.dat            : {bc_out}"
        )
        print(
            f"D_bc rows           : {len(bc_rows)}"
        )
        print(
            f"BC block length     : {args.block_length}"
        )
    else:
        print(
            "D_bc.dat            : not generated"
        )

    print(
        f"report              : {report_out}"
    )

    print()
    print(
        "GEDATSU command:"
    )
    print()
    print(
        f"cd {args.out_dir}"
    )
    print(
        "mkdir -p parted.0"
    )
    print(
        "gedatsu_simple_mesh_partitioner "
        f"-n {args.copies} "
        "-in node.dat "
        "-ie elem.dat "
        "--given_subdomain given_subdomain.dat "
        "-d parted.0"
    )


if __name__ == "__main__":
    main()