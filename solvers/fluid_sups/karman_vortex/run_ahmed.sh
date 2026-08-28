#!/usr/bin/env bash

run_ahmed_main() (
    # ============================================================
    # Safe shell settings
    # ============================================================

    set -euo pipefail


    # ============================================================
    # Basic settings
    # ============================================================

    ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


    # ------------------------------------------------------------
    # Existing Gmsh mesh
    # ------------------------------------------------------------

    MSH_FILE="${MSH_FILE:-}"


    # ------------------------------------------------------------
    # Output directory
    # ------------------------------------------------------------

    OUT_DIR="${OUT_DIR:-$ROOT/mesh_ahmed_groundprism}"


    # ------------------------------------------------------------
    # Optional solver executable
    # ------------------------------------------------------------

    SOLVER_BIN="${SOLVER_BIN:-}"


    # ------------------------------------------------------------
    # Inlet velocity
    #
    # U = (60, 0, 0)
    # ------------------------------------------------------------

    INLET_U="${INLET_U:-60.0}"


    # ------------------------------------------------------------
    # Existing validator compatibility
    #
    # Ahmed body does NOT use cylinder-radius checking.
    # ------------------------------------------------------------

    DIAMETER="${DIAMETER:-1.0}"


    # ============================================================
    # Tool paths
    # ============================================================

    SAVE_PHYSICAL_GROUPS="$ROOT/mesh_io/save_physical_groups.py"

    BUILD_GRAPH="$ROOT/mesh_io/build_graph.py"

    ELEM_TO_BC="$ROOT/mesh_io/elem_to_bc.py"

    FILTER_INDEX="$ROOT/mesh_io/fileter_index.py"

    MERGE_BC="$ROOT/mesh_io/merge_bc.py"

    VALIDATE_MESH="$ROOT/mesh_io/validate_mixed_mesh.py"

    CLASSIFY_SURF="$ROOT/plot/classify_surf_owner.py"


    # ============================================================
    # Usage
    # ============================================================

    usage() {

        cat <<EOF

Usage:
  bash run_ahmed.sh [options]

Options:

  --msh FILE
      Existing MSH2 mesh file.

  --out DIR
      Output directory.
      Default: $ROOT/mesh_ahmed_groundprism

  --solver EXE
      Solver executable.

  --inlet-u VALUE
      Inlet x velocity.
      Default: $INLET_U

  --skip-gmsh
      Accepted for compatibility.
      Gmsh generation is always skipped.

  -h, --help
      Show this help.


Boundary conditions:

  Inlet
      u = $INLET_U
      v = 0
      w = 0
      strong Dirichlet

  lowerWall
      u = 0
      v = 0
      w = 0
      strong no-slip

  upperWall
      uz = 0
      slip

  frontAndBack
      uy = 0
      slip

  Ahmed body / PRI6
      Nitsche no-slip

  Ahmed body / TET4
      Nitsche no-slip

  Outlet
      no velocity Dirichlet condition


Pipeline:

  existing .msh
        |
        v
  save_physical_groups.py
        |
        +--> node.dat
        +--> elem_tet.dat
        +--> elem_prism.dat
        +--> boundary TRI3 files
        |
        v
  ahmed_body_triangle_connectivity.dat
        |
        v
      surf.dat
        |
        v
  classify_surf_owner.py
        |
        +--> surf_prism.dat
        |
        +--> surf_tet.dat
        |
        +--> surf_orphan.dat

  TET4 + PRI6
        |
        v
  build_graph.py
        |
        v
      graph.dat

  Boundary files
        |
        v
      D_bc_v.dat


Environment variables:

  MSH_FILE   (optional; one *.msh under ROOT is auto-detected if unset)
  OUT_DIR    (optional; default: $ROOT/mesh_ahmed_groundprism)
  SOLVER_BIN
  INLET_U


Example:

  bash run_ahmed.sh \\
      --msh /path/to/your_mesh.msh \
      --out /path/to/output_mesh_dir

EOF
    }


    # ============================================================
    # Error trap
    # ============================================================

    trap '
        status=$?

        if [[ $status -ne 0 ]]; then

            echo >&2
            echo "==================================================" >&2
            echo "run_ahmed.sh FAILED" >&2
            echo "status  : $status" >&2
            echo "command : $BASH_COMMAND" >&2
            echo "==================================================" >&2

        fi
    ' EXIT


    # ============================================================
    # Utility functions
    # ============================================================

    require_file() {

        local file="$1"

        if [[ ! -f "$file" ]]; then

            echo "ERROR: required file not found:" >&2
            echo "  $file" >&2

            exit 1

        fi
    }


    require_command() {

        local cmd="$1"

        command -v "$cmd" >/dev/null 2>&1 || {

            echo "ERROR: command not found: $cmd" >&2

            exit 127
        }
    }


    # ------------------------------------------------------------
    # Read first integer in an element file.
    #
    # Supported headers:
    #
    #   N
    #
    # or:
    #
    #   N 3
    # ------------------------------------------------------------

    read_elem_count() {

        local file="$1"

        awk '
            NF > 0 {
                print $1
                exit
            }
        ' "$file"
    }


    # ============================================================
    # Parse command-line arguments
    # ============================================================

    while [[ $# -gt 0 ]]; do

        case "$1" in

            --msh)

                [[ $# -ge 2 ]] || {

                    echo "ERROR: --msh requires FILE" >&2

                    exit 2
                }

                MSH_FILE="$2"

                shift 2
                ;;


            --out)

                [[ $# -ge 2 ]] || {

                    echo "ERROR: --out requires DIR" >&2

                    exit 2
                }

                OUT_DIR="$2"

                shift 2
                ;;


            --solver)

                [[ $# -ge 2 ]] || {

                    echo "ERROR: --solver requires EXE" >&2

                    exit 2
                }

                SOLVER_BIN="$2"

                shift 2
                ;;


            --inlet-u)

                [[ $# -ge 2 ]] || {

                    echo "ERROR: --inlet-u requires VALUE" >&2

                    exit 2
                }

                INLET_U="$2"

                shift 2
                ;;


            --skip-gmsh)

                # Existing mesh only.
                # Kept for compatibility.

                shift
                ;;


            -h|--help)

                usage

                exit 0
                ;;


            *)

                echo "ERROR: unknown option: $1" >&2
                echo >&2

                usage

                exit 2
                ;;

        esac

    done


    # ============================================================
    # Basic checks
    # ============================================================

    require_command python3
    require_command awk

    # ------------------------------------------------------------
    # Backward-compatible defaults for sourced execution.
    #
    # exec_fom_fluid_sups_mixed.sh historically sources this script
    # without positional arguments. Keep that invocation valid.
    #
    # Input mesh priority:
    #   1. MSH_FILE environment variable / --msh
    #   2. exactly one *.msh file directly under ROOT
    #
    # OUT_DIR keeps the historical name because the existing execution
    # script consumes mesh files from mesh_ahmed_groundprism.
    # ------------------------------------------------------------

    if [[ -z "$MSH_FILE" ]]; then
        mapfile -t _msh_candidates < <(
            find "$ROOT" -maxdepth 1 -type f -name '*.msh' -print | sort
        )

        if [[ ${#_msh_candidates[@]} -eq 1 ]]; then
            MSH_FILE="${_msh_candidates[0]}"
            echo "Auto-detected MSH_FILE: $MSH_FILE"
        elif [[ ${#_msh_candidates[@]} -eq 0 ]]; then
            echo "ERROR: no .msh file found under: $ROOT" >&2
            echo "Set MSH_FILE=/path/to/mesh.msh before sourcing the execution script." >&2
            exit 2
        else
            echo "ERROR: multiple .msh files found under: $ROOT" >&2
            printf '  %s
' "${_msh_candidates[@]}" >&2
            echo "Set MSH_FILE=/path/to/the_mesh.msh explicitly." >&2
            exit 2
        fi
    fi

    require_file "$MSH_FILE"

    require_file "$SAVE_PHYSICAL_GROUPS"
    require_file "$BUILD_GRAPH"

    require_file "$ELEM_TO_BC"
    require_file "$FILTER_INDEX"
    require_file "$MERGE_BC"

    require_file "$CLASSIFY_SURF"


    mkdir -p "$OUT_DIR"


    # ============================================================
    # Header
    # ============================================================

    echo
    echo "============================================================"
    echo "Ahmed body mixed TET4 + PRI6 preprocessing"
    echo "============================================================"
    echo
    echo "MSH_FILE : $MSH_FILE"
    echo "OUT_DIR  : $OUT_DIR"
    echo "INLET_U  : $INLET_U"
    echo
    echo "Boundary conditions:"
    echo "  Inlet       : ($INLET_U, 0, 0) strong"
    echo "  lowerWall   : (0, 0, 0) strong no-slip"
    echo "  upperWall   : slip"
    echo "  front/back  : slip"
    echo "  Ahmed PRI6  : Nitsche no-slip"
    echo "  Ahmed TET4  : Nitsche no-slip"
    echo "  Outlet      : natural"
    echo
    echo "Gmsh generation:"
    echo "  skipped"
    echo


    # ============================================================
    # 1. Existing MSH
    # ============================================================

    echo
    echo "============================================================"
    echo "[1/6] Existing Gmsh mesh"
    echo "============================================================"
    echo

    echo "  MSH : $MSH_FILE"
    echo


    # ------------------------------------------------------------
    # Optional MSH version check
    # ------------------------------------------------------------

    if grep -q '^\$MeshFormat' "$MSH_FILE"; then

        mesh_version="$(
            awk '
                /^\$MeshFormat$/ {
                    getline
                    print $1
                    exit
                }
            ' "$MSH_FILE"
        )"

        echo "  Mesh format version : $mesh_version"

        case "$mesh_version" in

            2.*)
                ;;

            *)

                echo >&2
                echo "WARNING:" >&2
                echo "  This pipeline was designed for MSH2." >&2
                echo "  detected version = $mesh_version" >&2
                ;;

        esac

    fi


    # ============================================================
    # 2. Convert Physical Groups
    # ============================================================

    echo
    echo "============================================================"
    echo "[2/6] Convert Physical Groups"
    echo "============================================================"
    echo


    python3 \
        "$SAVE_PHYSICAL_GROUPS" \
        "$MSH_FILE" \
        --output-dir "$OUT_DIR"


    # ------------------------------------------------------------
    # Required volume files
    # ------------------------------------------------------------

    require_file "$OUT_DIR/node.dat"

    require_file "$OUT_DIR/elem_tet.dat"

    require_file "$OUT_DIR/elem_prism.dat"


    # ------------------------------------------------------------
    # Required boundary connectivity files
    # ------------------------------------------------------------

    require_file \
        "$OUT_DIR/Inlet_triangle_connectivity.dat"

    require_file \
        "$OUT_DIR/frontAndBack_triangle_connectivity.dat"

    require_file \
        "$OUT_DIR/lowerWall_triangle_connectivity.dat"

    require_file \
        "$OUT_DIR/upperWall_triangle_connectivity.dat"

    require_file \
        "$OUT_DIR/ahmed_body_triangle_connectivity.dat"


    echo
    echo "Volume mesh:"
    echo
    echo "  $OUT_DIR/node.dat"
    echo "  $OUT_DIR/elem_tet.dat"
    echo "  $OUT_DIR/elem_prism.dat"


    # ============================================================
    # 3. Split Ahmed wall into PRI6 / TET4 surfaces
    # ============================================================

    echo
    echo "============================================================"
    echo "[3/6] Split Ahmed-body surface"
    echo "============================================================"
    echo


    # ------------------------------------------------------------
    # Complete Ahmed surface
    # ------------------------------------------------------------

    cp -f \
        "$OUT_DIR/ahmed_body_triangle_connectivity.dat" \
        "$OUT_DIR/surf.dat"


    require_file "$OUT_DIR/surf.dat"


    # ------------------------------------------------------------
    # classify_surf_owner.py writes into the current directory.
    # Execute it inside OUT_DIR.
    # ------------------------------------------------------------

    (
        cd "$OUT_DIR"

        python3 \
            "$CLASSIFY_SURF" \
            "surf.dat" \
            "elem_prism.dat" \
            "elem_tet.dat"
    )


    require_file "$OUT_DIR/surf_prism.dat"

    require_file "$OUT_DIR/surf_tet.dat"

    require_file "$OUT_DIR/surf_orphan.dat"


    # ------------------------------------------------------------
    # Surface counts
    # ------------------------------------------------------------

    NSURF="$(
        read_elem_count "$OUT_DIR/surf.dat"
    )"

    NPRI_SURF="$(
        read_elem_count "$OUT_DIR/surf_prism.dat"
    )"

    NTET_SURF="$(
        read_elem_count "$OUT_DIR/surf_tet.dat"
    )"

    NORPHAN="$(
        read_elem_count "$OUT_DIR/surf_orphan.dat"
    )"


    echo
    echo "Ahmed wall classification:"
    echo
    echo "  total  : $NSURF"
    echo "  PRI6   : $NPRI_SURF"
    echo "  TET4   : $NTET_SURF"
    echo "  orphan : $NORPHAN"
    echo


    # ------------------------------------------------------------
    # Count consistency
    # ------------------------------------------------------------

    if [[ $((NPRI_SURF + NTET_SURF + NORPHAN)) -ne "$NSURF" ]]; then

        echo "ERROR:" >&2
        echo "  surface classification count mismatch." >&2
        echo >&2

        echo "  total  = $NSURF" >&2
        echo "  PRI6   = $NPRI_SURF" >&2
        echo "  TET4   = $NTET_SURF" >&2
        echo "  orphan = $NORPHAN" >&2

        exit 1

    fi


    # ------------------------------------------------------------
    # Orphan faces are forbidden
    # ------------------------------------------------------------

    if [[ "$NORPHAN" -ne 0 ]]; then

        echo "ERROR:" >&2
        echo "  Ahmed wall contains orphan TRI3 faces." >&2
        echo "  orphan = $NORPHAN" >&2
        echo >&2

        echo "Check:" >&2
        echo "  $OUT_DIR/surf_orphan.dat" >&2

        exit 1

    fi


    # Mesh-independent validation: only classification consistency and
    # orphan==0 are required. PRI6/TET4 face counts depend on the mesh.

    echo "Surface split: OK"
    echo

    echo "  PRI6 Nitsche surface:"
    echo "    $OUT_DIR/surf_prism.dat"

    echo
    echo "  TET4 Nitsche surface:"
    echo "    $OUT_DIR/surf_tet.dat"

    echo


    # ============================================================
    # 4. Unified TET4 + PRI6 nodal graph
    # ============================================================

    echo
    echo "============================================================"
    echo "[4/6] Build unified TET4 + PRI6 graph"
    echo "============================================================"
    echo


    python3 \
        "$BUILD_GRAPH" \
        --node "$OUT_DIR/node.dat" \
        --tet "$OUT_DIR/elem_tet.dat" \
        --prism "$OUT_DIR/elem_prism.dat" \
        --output "$OUT_DIR/graph.dat"


    require_file "$OUT_DIR/graph.dat"


    echo
    echo "Generated:"
    echo
    echo "  $OUT_DIR/graph.dat"


    # ============================================================
    # 5. Velocity Dirichlet BC
    # ============================================================

    echo
    echo "============================================================"
    echo "[5/6] Build velocity Dirichlet BC"
    echo "============================================================"
    echo

    echo "Strong velocity BC:"
    echo
    echo "  Inlet:"
    echo "    u = $INLET_U"
    echo "    v = 0"
    echo "    w = 0"
    echo
    echo "  lowerWall:"
    echo "    u = 0"
    echo "    v = 0"
    echo "    w = 0"
    echo "    no-slip"
    echo
    echo "  upperWall:"
    echo "    uz = 0"
    echo "    slip"
    echo
    echo "  frontAndBack:"
    echo "    uy = 0"
    echo "    slip"
    echo
    echo "Weak velocity BC:"
    echo
    echo "  Ahmed PRI6:"
    echo "    Nitsche no-slip"
    echo
    echo "  Ahmed TET4:"
    echo "    Nitsche no-slip"
    echo


    BC_TMP="$OUT_DIR/bc_tmp"


    # ------------------------------------------------------------
    # Remove stale BC files
    # ------------------------------------------------------------

    rm -rf "$BC_TMP"

    mkdir -p "$BC_TMP"


    # ------------------------------------------------------------
    # 5-1. Inlet
    #
    # u = INLET_U
    # v = 0
    # w = 0
    # ------------------------------------------------------------

    echo "  Building Inlet strong BC..."


    python3 \
        "$ELEM_TO_BC" \
        "$OUT_DIR/Inlet_triangle_connectivity.dat" \
        "$BC_TMP/Inlet_triangle_bc.dat" \
        3 "$INLET_U" 0 0


    # ------------------------------------------------------------
    # IMPORTANT:
    #
    # Ahmed body must NOT be added to D_bc_v.dat.
    #
    # Ahmed no-slip is imposed by:
    #
    #   PRI6/TRI3 Nitsche
    #   TET4/TRI3 Nitsche
    #
    # ------------------------------------------------------------


    # ------------------------------------------------------------
    # 5-2. frontAndBack
    #
    # Original full BC:
    #
    #   ux = 0
    #   uy = 0
    #   uz = 0
    #
    # Remove ux and uz:
    #
    #   uy = 0
    #
    # => slip wall
    # ------------------------------------------------------------

    echo "  Building frontAndBack slip BC..."


    python3 \
        "$ELEM_TO_BC" \
        "$OUT_DIR/frontAndBack_triangle_connectivity.dat" \
        "$BC_TMP/frontAndBack_triangle_bc.dat" \
        3 0 0 0


    python3 \
        "$FILTER_INDEX" \
        "$BC_TMP/frontAndBack_triangle_bc.dat" \
        "$BC_TMP/filtered_frontAndBack_triangle_bc.dat" \
        0 2


    # ------------------------------------------------------------
    # 5-3. lowerWall
    #
    # STRONG NO-SLIP
    #
    #   ux = 0
    #   uy = 0
    #   uz = 0
    #
    # IMPORTANT:
    # DO NOT call fileter_index.py here.
    # ------------------------------------------------------------

    echo "  Building lowerWall strong no-slip BC..."


    python3 \
        "$ELEM_TO_BC" \
        "$OUT_DIR/lowerWall_triangle_connectivity.dat" \
        "$BC_TMP/lowerWall_triangle_bc.dat" \
        3 0 0 0


    # ------------------------------------------------------------
    # 5-4. upperWall
    #
    # Original full BC:
    #
    #   ux = 0
    #   uy = 0
    #   uz = 0
    #
    # Remove ux and uy:
    #
    #   uz = 0
    #
    # => slip wall
    # ------------------------------------------------------------

    echo "  Building upperWall slip BC..."


    python3 \
        "$ELEM_TO_BC" \
        "$OUT_DIR/upperWall_triangle_connectivity.dat" \
        "$BC_TMP/upperWall_triangle_bc.dat" \
        3 0 0 0


    python3 \
        "$FILTER_INDEX" \
        "$BC_TMP/upperWall_triangle_bc.dat" \
        "$BC_TMP/filtered_upperWall_triangle_bc.dat" \
        0 1


    # ------------------------------------------------------------
    # 5-5. Merge
    #
    # merge_bc.py:
    # later files take priority for shared nodes.
    #
    # Priority:
    #
    #   front/back
    #       <
    #   upperWall
    #       <
    #   lowerWall
    #       <
    #   Inlet
    #
    # Inlet is last intentionally.
    # ------------------------------------------------------------

    echo
    echo "  Merging velocity BC..."


    python3 \
        "$MERGE_BC" \
        "$BC_TMP/filtered_frontAndBack_triangle_bc.dat" \
        "$BC_TMP/filtered_upperWall_triangle_bc.dat" \
        "$BC_TMP/lowerWall_triangle_bc.dat" \
        "$BC_TMP/Inlet_triangle_bc.dat" \
        -d ./ \
        -o "$BC_TMP/merged_bc.dat"


    require_file "$BC_TMP/merged_bc.dat"


    cp -f \
        "$BC_TMP/merged_bc.dat" \
        "$OUT_DIR/D_bc_v.dat"


    require_file "$OUT_DIR/D_bc_v.dat"


    echo
    echo "Generated:"
    echo
    echo "  $OUT_DIR/D_bc_v.dat"


    # ============================================================
    # 6. Validation
    # ============================================================

    echo
    echo "============================================================"
    echo "[6/6] Validate mixed mesh"
    echo "============================================================"
    echo


    if [[ -f "$VALIDATE_MESH" ]]; then

        # --------------------------------------------------------
        # Ahmed body:
        # DO NOT use --check-cylinder-radius
        # --------------------------------------------------------

        python3 \
            "$VALIDATE_MESH" \
            "$OUT_DIR" \
            --diameter "$DIAMETER"

    else

        echo "WARNING:"
        echo "  validator not found:"
        echo "  $VALIDATE_MESH"
        echo
        echo "  mixed mesh validator skipped."

    fi


    # ============================================================
    # Compatibility aliases
    # ============================================================

    echo
    echo "============================================================"
    echo "Compatibility aliases"
    echo "============================================================"
    echo


    cp -f \
        "$OUT_DIR/node.dat" \
        "$OUT_DIR/440_0_node.dat"


    cp -f \
        "$OUT_DIR/elem_tet.dat" \
        "$OUT_DIR/440_0_elem_10.dat"


    cp -f \
        "$OUT_DIR/elem_prism.dat" \
        "$OUT_DIR/440_0_elem_13.dat"


    cp -f \
        "$OUT_DIR/graph.dat" \
        "$OUT_DIR/440_0_graph_sort.dat"


    # ============================================================
    # Final summary
    # ============================================================

    cat <<EOF


============================================================
Ahmed mixed-mesh preprocessing finished
============================================================

Input mesh:

  $MSH_FILE


Volume mesh:

  $OUT_DIR/node.dat
  $OUT_DIR/elem_tet.dat
  $OUT_DIR/elem_prism.dat


Unified graph:

  $OUT_DIR/graph.dat


Ahmed complete surface:

  $OUT_DIR/surf.dat

  faces = $NSURF


PRI6-owned Ahmed surface:

  $OUT_DIR/surf_prism.dat

  faces = $NPRI_SURF

  BC:
      Nitsche no-slip


TET4-owned Ahmed surface:

  $OUT_DIR/surf_tet.dat

  faces = $NTET_SURF

  BC:
      Nitsche no-slip


Orphan surface:

  $OUT_DIR/surf_orphan.dat

  faces = $NORPHAN


Velocity Dirichlet BC:

  $OUT_DIR/D_bc_v.dat


============================================================
Boundary conditions
============================================================

Inlet:

  ux = $INLET_U
  uy = 0
  uz = 0

  strong Dirichlet


lowerWall:

  ux = 0
  uy = 0
  uz = 0

  strong no-slip


upperWall:

  uz = 0

  ux = free
  uy = free

  slip


frontAndBack:

  uy = 0

  ux = free
  uz = free

  slip


Ahmed body / PRI6:

  ux = 0
  uy = 0
  uz = 0

  weak Nitsche no-slip

  surface:
      surf_prism.dat


Ahmed body / TET4:

  ux = 0
  uy = 0
  uz = 0

  weak Nitsche no-slip

  surface:
      surf_tet.dat


Outlet:

  velocity Dirichlet BC:
      none


============================================================
Compatibility files
============================================================

  $OUT_DIR/440_0_node.dat

  $OUT_DIR/440_0_elem_10.dat

  $OUT_DIR/440_0_elem_13.dat

  $OUT_DIR/440_0_graph_sort.dat


============================================================

EOF


    # ============================================================
    # Optional solver
    # ============================================================

    if [[ -n "$SOLVER_BIN" ]]; then


        [[ -x "$SOLVER_BIN" ]] || {

            echo "ERROR: solver is not executable:" >&2
            echo "  $SOLVER_BIN" >&2

            exit 1
        }


        echo
        echo "============================================================"
        echo "[solver]"
        echo "============================================================"
        echo
        echo "Executable:"
        echo "  $SOLVER_BIN"
        echo
        echo "Directory:"
        echo "  $OUT_DIR"
        echo


        # --------------------------------------------------------
        # Temporarily disable set -e so the solver status can be
        # captured and propagated.
        # --------------------------------------------------------

        set +e


        "$SOLVER_BIN" "$OUT_DIR"


        solver_status=$?


        set -e


        echo
        echo "============================================================"
        echo "Solver finished"
        echo "exit status = $solver_status"
        echo "============================================================"


        exit "$solver_status"


    else


        cat <<EOF

Solver was not launched.


Example:

  bash run_ahmed.sh \\
      --msh "$MSH_FILE" \\
      --out "$OUT_DIR" \\
      --solver /path/to/solver


EOF

    fi
)


# ================================================================
# Run
# ================================================================

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then

    # ------------------------------------------------------------
    # Executed directly:
    #
    #   bash run_ahmed.sh ...
    # ------------------------------------------------------------

    run_ahmed_main "$@"


else

    # ------------------------------------------------------------
    # Sourced:
    #
    #   . ./run_ahmed.sh
    #
    # Do not inherit the caller shell's positional parameters.
    # ------------------------------------------------------------

    run_ahmed_main

fi

