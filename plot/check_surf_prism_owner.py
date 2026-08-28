#!/usr/bin/env python3

import sys


def read_elements(path):
    """
    対応する形式:

    Format A:
        Nelem NEN
        node0 node1 node2 ...
        node0 node1 node2 ...
        ...

    Format B:
        Nelem NEN
        elem_id node0 node1 ...
        ...

    Format C:
        Nelem
        elem_id NEN node0 node1 ...
        ...

    自動判別する。
    """

    with open(path, "r") as f:
        lines = [
            line.strip()
            for line in f
            if line.strip()
        ]

    if not lines:
        raise RuntimeError(
            f"empty file: {path}"
        )

    # ============================================================
    # Header
    # ============================================================

    header = lines[0].split()

    elems = []

    # ------------------------------------------------------------
    # Format:
    #
    #   Nelem NEN
    #
    # 例:
    #   33244 3
    # ------------------------------------------------------------

    if len(header) == 2:

        nelem = int(header[0])
        header_nen = int(header[1])

        data_lines = lines[1:]

        if len(data_lines) < nelem:
            raise RuntimeError(
                f"{path}: expected {nelem} elements, "
                f"but only {len(data_lines)} lines found"
            )

        for e in range(nelem):

            cols = list(
                map(
                    int,
                    data_lines[e].split()
                )
            )

            # --------------------------------------------
            # Case 1:
            #
            # node0 node1 ...
            # --------------------------------------------

            if len(cols) == header_nen:

                eid = e
                nodes = cols

            # --------------------------------------------
            # Case 2:
            #
            # elem_id node0 node1 ...
            # --------------------------------------------

            elif len(cols) == header_nen + 1:

                eid = cols[0]
                nodes = cols[1:]

            # --------------------------------------------
            # Case 3:
            #
            # elem_id NEN node0 node1 ...
            # --------------------------------------------

            elif (
                len(cols) == header_nen + 2
                and cols[1] == header_nen
            ):

                eid = cols[0]
                nodes = cols[2:]

            else:

                raise RuntimeError(
                    f"{path}: cannot understand element line "
                    f"{e + 2}:\n"
                    f"  {data_lines[e]}\n"
                    f"header NEN = {header_nen}, "
                    f"number of columns = {len(cols)}"
                )

            if len(nodes) != header_nen:

                raise RuntimeError(
                    f"{path}: invalid connectivity "
                    f"at element {e}: "
                    f"expected {header_nen}, "
                    f"got {len(nodes)}"
                )

            elems.append(
                (
                    eid,
                    header_nen,
                    nodes
                )
            )

        return elems


    # ------------------------------------------------------------
    # Format:
    #
    #   Nelem
    #   elem_id NEN node0 node1 ...
    # ------------------------------------------------------------

    elif len(header) == 1:

        nelem = int(header[0])

        data_lines = lines[1:]

        if len(data_lines) < nelem:

            raise RuntimeError(
                f"{path}: expected {nelem} elements, "
                f"but only {len(data_lines)} lines found"
            )

        for e in range(nelem):

            cols = list(
                map(
                    int,
                    data_lines[e].split()
                )
            )

            if len(cols) < 2:

                raise RuntimeError(
                    f"{path}: invalid line {e + 2}: "
                    f"{data_lines[e]}"
                )

            eid = cols[0]
            nen = cols[1]

            nodes = cols[2:]

            if len(nodes) != nen:

                raise RuntimeError(
                    f"{path}: invalid connectivity "
                    f"eid={eid}: "
                    f"NEN={nen}, "
                    f"nodes={len(nodes)}"
                )

            elems.append(
                (
                    eid,
                    nen,
                    nodes
                )
            )

        return elems


    else:

        raise RuntimeError(
            f"{path}: invalid header:\n"
            f"  {lines[0]}"
        )


# ================================================================
# Main
# ================================================================

if len(sys.argv) != 3:

    print(
        "Usage:\n"
        "  python3 check_surf_prism_owner.py "
        "surf.dat elem_prism.dat"
    )

    sys.exit(2)


surf_file = sys.argv[1]
prism_file = sys.argv[2]


print()
print("Reading:")
print(f"  surface : {surf_file}")
print(f"  prism   : {prism_file}")
print()


surf = read_elements(
    surf_file
)

prisms = read_elements(
    prism_file
)


print(
    f"Loaded surface elements : {len(surf)}"
)

print(
    f"Loaded prism elements   : {len(prisms)}"
)


# ================================================================
# Validate element types
# ================================================================

for eid, nen, nodes in surf:

    if nen != 3:

        raise RuntimeError(
            f"surface is not TRI3: "
            f"eid={eid}, NEN={nen}"
        )


for eid, nen, nodes in prisms:

    if nen != 6:

        raise RuntimeError(
            f"volume element is not PRI6: "
            f"eid={eid}, NEN={nen}"
        )


# ================================================================
# Build PRI6 triangular-face dictionary
#
#
# PRI6 reference:
#
#        3 -------- 4
#         \        /
#          \      /
#            5
#
#
#        0 -------- 1
#         \        /
#          \      /
#            2
#
#
# triangular faces:
#
#   face 0 = {0,1,2}
#   face 1 = {3,4,5}
#
# ================================================================

prism_faces = {}


for eid, nen, nodes in prisms:

    face0 = tuple(
        sorted(
            (
                nodes[0],
                nodes[1],
                nodes[2]
            )
        )
    )

    face1 = tuple(
        sorted(
            (
                nodes[3],
                nodes[4],
                nodes[5]
            )
        )
    )


    prism_faces.setdefault(
        face0,
        []
    ).append(
        (
            eid,
            0
        )
    )


    prism_faces.setdefault(
        face1,
        []
    ).append(
        (
            eid,
            1
        )
    )


# ================================================================
# Check TRI3 -> PRI6
# ================================================================

mapped = []

missing = []


for eid, nen, nodes in surf:

    key = tuple(
        sorted(nodes)
    )

    owners = prism_faces.get(
        key
    )


    if owners is None:

        missing.append(
            (
                eid,
                nodes
            )
        )

    else:

        mapped.append(
            (
                eid,
                nodes,
                owners
            )
        )


# ================================================================
# Report
# ================================================================

nsurf = len(surf)
nmapped = len(mapped)
nmissing = len(missing)


print()
print("========================================")
print("TRI3 -> PRI6 owner check")
print("========================================")

print(
    f"surface faces : {nsurf}"
)

print(
    f"mapped        : {nmapped}"
)

print(
    f"missing       : {nmissing}"
)


if nsurf > 0:

    ratio = (
        100.0
        * nmapped
        / nsurf
    )

    print(
        f"mapped ratio  : {ratio:.8f} %"
    )


# ================================================================
# Missing faces
# ================================================================

if nmissing > 0:

    print()
    print("----------------------------------------")
    print("First missing faces")
    print("----------------------------------------")

    for eid, nodes in missing[:50]:

        print(
            f"surf eid={eid:8d} "
            f"nodes={nodes}"
        )


    print()
    print("========================================")
    print("RESULT: FAILED")
    print("========================================")

    print()
    print(
        "Some Ahmed-body TRI3 faces are NOT "
        "triangular end faces of PRI6 elements."
    )

    print()
    print(
        "This means the problem already exists "
        "before MPI partitioning."
    )

    sys.exit(1)


# ================================================================
# Success
# ================================================================

print()
print("========================================")
print("RESULT: OK")
print("========================================")

print()
print(
    "Every surface TRI3 has a PRI6 owner "
    "before MPI partitioning."
)

print()
print(
    "Therefore, if MPI execution still reports "
    "'no local PRI6 owner',"
)

print(
    "the problem is in the partitioning of "
    "surf_graph.dat versus graph_elem_prism.dat."
)

sys.exit(0)
