#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert an ASCII Gmsh MSH2 mesh to the mixed-element BBFE/Monolis files.

Canonical outputs in --output-dir:
  node.dat          : shared compact 0-based node coordinates
  elem_tet.dat      : first-order tetra connectivity (NEN=4)
  elem_prism.dat    : first-order prism/wedge connectivity (NEN=6)
  surf.dat          : Cylinder_wall connectivity if it has one face type
  <group>_<type>_connectivity.dat
  <group>_node_ids.dat
  conversion_report.txt

The same compact node numbering is used by every output file.
"""
from __future__ import annotations

import argparse
import shlex
from collections import defaultdict
from pathlib import Path

ELEMENT_TYPES = {
    1: ("line", 1, 2),
    2: ("triangle", 2, 3),
    3: ("quad", 2, 4),
    4: ("tetra", 3, 4),
    5: ("hexahedron", 3, 8),
    6: ("wedge", 3, 6),  # Gmsh prism; meshio name is wedge
    7: ("pyramid", 3, 5),
    15: ("vertex", 0, 1),
}


def section(lines: list[str], name: str, required: bool = True) -> list[str]:
    start = f"${name}"
    end = f"$End{name}"
    try:
        i0 = next(i for i, line in enumerate(lines) if line.strip() == start)
        i1 = next(i for i, line in enumerate(lines[i0 + 1 :], i0 + 1) if line.strip() == end)
    except StopIteration:
        if required:
            raise ValueError(f"section {start} ... {end} not found")
        return []
    return lines[i0 + 1 : i1]


def parse_mesh_format(lines: list[str]) -> None:
    sec = section(lines, "MeshFormat")
    p = sec[0].split()
    if not p[0].startswith("2.") or int(p[1]) != 0:
        raise ValueError("converter requires ASCII MSH2; run gmsh with '-format msh2'")


def parse_physical_names(lines: list[str]) -> dict[tuple[int, int], str]:
    sec = section(lines, "PhysicalNames", required=False)
    if not sec:
        return {}
    n = int(sec[0])
    out: dict[tuple[int, int], str] = {}
    for raw in sec[1 : 1 + n]:
        p = shlex.split(raw)
        if len(p) < 3:
            raise ValueError(f"invalid PhysicalNames line: {raw}")
        out[(int(p[0]), int(p[1]))] = p[2]
    return out


def parse_nodes(lines: list[str]) -> dict[int, tuple[float, float, float]]:
    sec = section(lines, "Nodes")
    n = int(sec[0])
    nodes: dict[int, tuple[float, float, float]] = {}
    for raw in sec[1 : 1 + n]:
        p = raw.split()
        nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
    if len(nodes) != n:
        raise ValueError(f"node count mismatch: expected {n}, parsed {len(nodes)}")
    return nodes


def parse_elements(lines: list[str], physical_names: dict[tuple[int, int], str]):
    sec = section(lines, "Elements")
    n = int(sec[0])
    groups: dict[str, dict[str, list[list[int]]]] = defaultdict(lambda: defaultdict(list))
    unknown: dict[int, int] = defaultdict(int)
    no_physical = 0

    for raw in sec[1 : 1 + n]:
        p = raw.split()
        etype = int(p[1])
        ntags = int(p[2])
        tags = list(map(int, p[3 : 3 + ntags]))
        conn = list(map(int, p[3 + ntags :]))
        info = ELEMENT_TYPES.get(etype)
        if info is None:
            unknown[etype] += 1
            continue
        cell_name, dim, nen = info
        if len(conn) != nen:
            raise ValueError(f"Gmsh element type {etype} expects {nen} nodes, got {len(conn)}")
        physical_tag = tags[0] if tags else 0
        if physical_tag == 0:
            no_physical += 1
            continue
        name = physical_names.get((dim, physical_tag), f"Physical_{dim}_{physical_tag}")
        groups[name][cell_name].append(conn)
    return groups, dict(unknown), no_physical


def write_nodes(path: Path, tags: list[int], nodes: dict[int, tuple[float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(tags)}\n")
        for tag in tags:
            x, y, z = nodes[tag]
            f.write(f"{x:.16g} {y:.16g} {z:.16g}\n")


def write_connectivity(path: Path, elems: list[list[int]], node_map: dict[int, int]) -> None:
    if not elems:
        return
    nen = len(elems[0])
    if any(len(e) != nen for e in elems):
        raise ValueError(f"mixed NEN in {path}")
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(elems)} {nen}\n")
        for elem in elems:
            f.write(" ".join(str(node_map[n]) for n in elem) + "\n")


def write_node_ids(path: Path, ids: set[int]) -> None:
    vals = sorted(ids)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(vals)}\n")
        for n in vals:
            f.write(f"{n}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_mesh", type=Path)
    ap.add_argument("--output-dir", type=Path, default=Path("mesh_karman_vortex"))
    ap.add_argument("--fluid-name", default="Fluid")
    ap.add_argument("--wall-name", default="Cylinder_wall")
    ap.add_argument(
        "--project-cylinder-radius",
        type=float,
        default=None,
        help=(
            "Project nodes in the Cylinder_wall Physical Surface onto this exact "
            "xy radius before node.dat is written. Useful because Gmsh topological "
            "boundary-layer extrusion uses smoothed discrete normals."
        ),
    )
    args = ap.parse_args()

    lines = args.input_mesh.read_text(encoding="utf-8").splitlines()
    parse_mesh_format(lines)
    phys = parse_physical_names(lines)
    nodes = parse_nodes(lines)
    groups, unknown, no_physical = parse_elements(lines, phys)

    if args.fluid_name not in groups:
        raise SystemExit(
            f"Physical Volume '{args.fluid_name}' not found. Available groups: "
            + ", ".join(sorted(groups))
        )

    fluid = groups[args.fluid_name]
    tets = fluid.get("tetra", [])
    prisms = fluid.get("wedge", [])
    unsupported_vol = {
        k: len(v) for k, v in fluid.items() if k in {"hexahedron", "pyramid"} and v
    }
    if not tets:
        raise SystemExit("Fluid contains no first-order tetrahedra")
    if not prisms:
        raise SystemExit("Fluid contains no first-order prisms/wedges")
    if unsupported_vol:
        raise SystemExit(f"Fluid also contains unsupported volume types: {unsupported_vol}")

    fluid_tags = sorted({n for elems in (tets, prisms) for e in elems for n in e})
    node_map = {tag: i for i, tag in enumerate(fluid_tags)}

    # Gmsh topological boundary-layer extrusion follows a Gouraud-smoothed
    # discrete normal field. When an outer cylindrical interface is extruded
    # inward, the generated top surface is therefore only approximately r=d/2.
    # For the CFD geometry the wall position is more important than preserving
    # that small normal-field error, so optionally project only wall nodes onto
    # the exact cylinder. Connectivity and the shared compact node IDs are not
    # changed.
    nodes_out = dict(nodes)
    projection_stats = None
    if args.project_cylinder_radius is not None:
        target = args.project_cylinder_radius
        if target <= 0:
            raise SystemExit("--project-cylinder-radius must be positive")
        wall_group = groups.get(args.wall_name, {})
        wall_tags = {
            n
            for cell_name in ("triangle", "quad")
            for elem in wall_group.get(cell_name, [])
            for n in elem
            if n in node_map
        }
        if not wall_tags:
            raise SystemExit(
                f"Cannot project cylinder: Physical Surface '{args.wall_name}' "
                "has no triangle/quad nodes belonging to Fluid"
            )
        before = []
        max_shift = 0.0
        for tag in wall_tags:
            x, y, z = nodes_out[tag]
            r = (x*x + y*y) ** 0.5
            if r <= 1.0e-14:
                raise SystemExit(f"Cylinder wall node {tag} lies on the axis")
            before.append(r)
            scale = target / r
            xn, yn = x * scale, y * scale
            max_shift = max(max_shift, ((xn-x)**2 + (yn-y)**2) ** 0.5)
            nodes_out[tag] = (xn, yn, z)
        projection_stats = (len(wall_tags), min(before), max(before), max_shift, target)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    write_nodes(out / "node.dat", fluid_tags, nodes_out)
    write_connectivity(out / "elem_tet.dat", tets, node_map)
    write_connectivity(out / "elem_prism.dat", prisms, node_map)

    report: list[str] = []
    for group_name in sorted(groups):
        local_nodes: set[int] = set()
        for cell_name in sorted(groups[group_name]):
            elems = groups[group_name][cell_name]
            valid = [e for e in elems if all(n in node_map for n in e)]
            if not valid:
                continue
            path = out / f"{group_name}_{cell_name}_connectivity.dat"
            write_connectivity(path, valid, node_map)
            for e in valid:
                local_nodes.update(node_map[n] for n in e)
            report.append(f"{group_name:24s} {cell_name:12s} {len(valid):9d}  {path.name}")
        if local_nodes:
            write_node_ids(out / f"{group_name}_node_ids.dat", local_nodes)

    # Legacy surface reader supports one NEN only. For a wall-normal wedge layer,
    # Cylinder_wall is normally triangular. Keep mixed wall faces explicit instead.
    wall = groups.get(args.wall_name, {})
    wall_tri = [e for e in wall.get("triangle", []) if all(n in node_map for n in e)]
    wall_quad = [e for e in wall.get("quad", []) if all(n in node_map for n in e)]
    if wall_tri and not wall_quad:
        write_connectivity(out / "surf.dat", wall_tri, node_map)
    elif wall_quad and not wall_tri:
        write_connectivity(out / "surf.dat", wall_quad, node_map)
    elif wall_tri and wall_quad:
        write_connectivity(out / "surf_triangle.dat", wall_tri, node_map)
        write_connectivity(out / "surf_quad.dat", wall_quad, node_map)

    with (out / "conversion_report.txt").open("w", encoding="utf-8") as f:
        f.write("Karman mixed tetra/prism conversion\n")
        f.write("===================================\n")
        f.write(f"mesh       : {args.input_mesh}\n")
        f.write(f"nodes      : {len(fluid_tags)}\n")
        f.write(f"tetra      : {len(tets)}\n")
        f.write(f"prism      : {len(prisms)}\n")
        f.write(f"wall tri   : {len(wall_tri)}\n")
        f.write(f"wall quad  : {len(wall_quad)}\n")
        if projection_stats is not None:
            nw, rmin, rmax, max_shift, target = projection_stats
            f.write(
                f"wall projection: nodes={nw}, before_r=[{rmin:.16g}, {rmax:.16g}], "
                f"target_r={target:.16g}, max_shift={max_shift:.16g}\n"
            )
        f.write("\nPhysical-group connectivity\n")
        for line in report:
            f.write(line + "\n")
        if unknown:
            f.write(f"\nIgnored unsupported Gmsh element types: {unknown}\n")
        if no_physical:
            f.write(f"Elements without Physical tag ignored: {no_physical}\n")

    print(f"[convert] nodes={len(fluid_tags)} tetra={len(tets)} prism={len(prisms)}")
    if projection_stats is not None:
        nw, rmin, rmax, max_shift, target = projection_stats
        print(
            f"[convert] projected Cylinder_wall nodes={nw}: "
            f"r_before=[{rmin:.6g},{rmax:.6g}] -> {target:.6g}, "
            f"max_shift={max_shift:.6g}"
        )
    print(f"[convert] output={out}")


if __name__ == "__main__":
    main()
