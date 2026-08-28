#!/usr/bin/env python3

import sys


def iter_elements(path):
    with open(path, "r") as f:

        first = f.readline().strip()

        if not first:
            raise RuntimeError(f"empty file: {path}")

        header = first.split()

        # ========================================================
        # Format:
        #
        #   Nelem NEN
        # ========================================================

        if len(header) == 2:

            nelem = int(header[0])
            header_nen = int(header[1])

            for e in range(nelem):

                line = f.readline()

                if not line:
                    raise RuntimeError(
                        f"{path}: unexpected EOF at element {e}"
                    )

                cols = list(map(int, line.split()))

                # node0 node1 ...
                if len(cols) == header_nen:

                    eid = e
                    nodes = cols

                # elem_id node0 node1 ...
                elif len(cols) == header_nen + 1:

                    eid = cols[0]
                    nodes = cols[1:]

                # elem_id NEN node0 node1 ...
                elif (
                    len(cols) == header_nen + 2
                    and cols[1] == header_nen
                ):

                    eid = cols[0]
                    nodes = cols[2:]

                else:

                    raise RuntimeError(
                        f"{path}: invalid element line:\n"
                        f"{line}"
                    )

                yield eid, header_nen, nodes


        # ========================================================
        # Format:
        #
        #   Nelem
        #   elem_id NEN node0 ...
        # ========================================================

        elif len(header) == 1:

            nelem = int(header[0])

            for e in range(nelem):

                line = f.readline()

                if not line:
                    raise RuntimeError(
                        f"{path}: unexpected EOF at element {e}"
                    )

                cols = list(map(int, line.split()))

                if len(cols) < 2:

                    raise RuntimeError(
                        f"{path}: invalid element line:\n"
                        f"{line}"
                    )

                eid = cols[0]
                nen = cols[1]
                nodes = cols[2:]

                if len(nodes) != nen:

                    raise RuntimeError(
                        f"{path}: invalid connectivity "
                        f"eid={eid}, nen={nen}, nodes={len(nodes)}"
                    )

                yield eid, nen, nodes


        else:

            raise RuntimeError(
                f"{path}: invalid header: {first}"
            )


def key3(a, b, c):
    return tuple(sorted((a, b, c)))


def write_tri_file(path, faces):

    with open(path, "w") as f:

        f.write(f"{len(faces)} 3\n")

        for _, nodes in faces:

            f.write(
                f"{nodes[0]} "
                f"{nodes[1]} "
                f"{nodes[2]}\n"
            )


# ================================================================
# Arguments
# ================================================================

if len(sys.argv) != 4:

    print(
        "Usage:\n"
        "  python3 classify_surf_owner.py "
        "surf.dat elem_prism.dat elem_tet.dat"
    )

    sys.exit(2)


surf_file = sys.argv[1]
prism_file = sys.argv[2]
tet_file = sys.argv[3]


# ================================================================
# Read surface
# ================================================================

surface = []

surface_by_key = {}


for eid, nen, nodes in iter_elements(surf_file):

    if nen != 3:

        raise RuntimeError(
            f"surface is not TRI3: "
            f"eid={eid}, nen={nen}"
        )

    entry = (eid, nodes)

    surface.append(entry)

    k = tuple(sorted(nodes))

    surface_by_key.setdefault(
        k, []
    ).append(entry)


surface_keys = set(
    surface_by_key.keys()
)


print()
print("========================================")
print("Surface")
print("========================================")

print(
    f"surface faces = {len(surface)}"
)


# ================================================================
# Scan PRI6
#
# PRI6 triangular end faces:
#
#   face 0 = 0,1,2
#   face 1 = 3,4,5
# ================================================================

prism_owned = set()

nprism = 0


for eid, nen, nodes in iter_elements(prism_file):

    if nen != 6:

        raise RuntimeError(
            f"not PRI6: eid={eid}, nen={nen}"
        )

    nprism += 1

    f0 = key3(
        nodes[0],
        nodes[1],
        nodes[2]
    )

    f1 = key3(
        nodes[3],
        nodes[4],
        nodes[5]
    )

    if f0 in surface_keys:
        prism_owned.add(f0)

    if f1 in surface_keys:
        prism_owned.add(f1)


print()
print("========================================")
print("PRI6 scan")
print("========================================")

print(
    f"PRI6 elements = {nprism}"
)

print(
    f"surface keys owned by PRI6 = "
    f"{len(prism_owned)}"
)


# ================================================================
# Surface faces without PRI6 owner
# ================================================================

no_prism = (
    surface_keys
    - prism_owned
)


print(
    f"surface keys without PRI6 owner = "
    f"{len(no_prism)}"
)


# ================================================================
# Scan TET4
#
# TET4 faces:
#
#   0,1,2
#   0,1,3
#   0,2,3
#   1,2,3
#
# Only test faces which are already missing from PRI6.
# ================================================================

tet_owned = set()

ntet = 0


for eid, nen, nodes in iter_elements(tet_file):

    if nen != 4:

        raise RuntimeError(
            f"not TET4: eid={eid}, nen={nen}"
        )

    ntet += 1

    faces = (
        key3(nodes[0], nodes[1], nodes[2]),
        key3(nodes[0], nodes[1], nodes[3]),
        key3(nodes[0], nodes[2], nodes[3]),
        key3(nodes[1], nodes[2], nodes[3]),
    )

    for face in faces:

        if face in no_prism:

            tet_owned.add(face)


print()
print("========================================")
print("TET4 scan")
print("========================================")

print(
    f"TET4 elements = {ntet}"
)

print(
    f"missing PRI6 faces owned by TET4 = "
    f"{len(tet_owned)}"
)


# ================================================================
# Classification
# ================================================================

prism_faces = []

tet_faces = []

orphan_faces = []


for eid, nodes in surface:

    k = tuple(sorted(nodes))

    if k in prism_owned:

        prism_faces.append(
            (eid, nodes)
        )

    elif k in tet_owned:

        tet_faces.append(
            (eid, nodes)
        )

    else:

        orphan_faces.append(
            (eid, nodes)
        )


print()
print("========================================")
print("FINAL CLASSIFICATION")
print("========================================")

print(
    f"total      : {len(surface)}"
)

print(
    f"PRI6 owner : {len(prism_faces)}"
)

print(
    f"TET4 owner : {len(tet_faces)}"
)

print(
    f"orphan     : {len(orphan_faces)}"
)


# ================================================================
# Output
# ================================================================

write_tri_file(
    "surf_prism.dat",
    prism_faces
)

write_tri_file(
    "surf_tet.dat",
    tet_faces
)

write_tri_file(
    "surf_orphan.dat",
    orphan_faces
)


print()
print("Generated:")
print("  surf_prism.dat")
print("  surf_tet.dat")
print("  surf_orphan.dat")


# ================================================================
# First orphan faces
# ================================================================

if orphan_faces:

    print()
    print("----------------------------------------")
    print("First orphan faces")
    print("----------------------------------------")

    for eid, nodes in orphan_faces[:30]:

        print(
            f"surf eid={eid:8d} "
            f"nodes={nodes}"
        )


print()
