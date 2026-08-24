#!/usr/bin/env bash

run_karman_main() (
    # ここからサブシェルなので、
    # set -e や exit が親の対話シェルを終了させない
    set -euo pipefail

    ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    GEO_FILE="${GEO_FILE:-$ROOT/gmsh/karman_vortex_prism_tet.geo}"
    OUT_DIR="${OUT_DIR:-$ROOT/mesh_karman_vortex}"
    MSH_FILE="${MSH_FILE:-$OUT_DIR/karman_vortex_prism_tet.msh}"
    SOLVER_BIN="${SOLVER_BIN:-}"
    INLET_U="${INLET_U:-1.0}"
    DIAMETER="${DIAMETER:-1.0}"
    CHECK_CYLINDER_RADIUS="${CHECK_CYLINDER_RADIUS:-1}"
    PROJECT_CYLINDER_WALL="${PROJECT_CYLINDER_WALL:-1}"
    GMSH_LOG="${GMSH_LOG:-$OUT_DIR/gmsh.log}"

    usage() {
        cat <<EOF
Usage: run_karman.sh [--geo FILE] [--out DIR] [--solver EXE] [--skip-gmsh]

Pipeline:
  Gmsh -> save_physical_groups.py -> node/elem_tet/elem_prism
       -> build_graph.py -> graph.dat -> BC -> validation -> optional solver

Environment alternatives:
  GEO_FILE
  OUT_DIR
  MSH_FILE
  SOLVER_BIN
  INLET_U
  DIAMETER
  CHECK_CYLINDER_RADIUS
  PROJECT_CYLINDER_WALL
  GMSH_LOG
EOF
    }

    trap '
        status=$?
        if [[ $status -ne 0 ]]; then
            echo >&2
            echo "========================================" >&2
            echo "run_karman.sh failed" >&2
            echo "status : $status" >&2
            echo "command: $BASH_COMMAND" >&2
            echo "========================================" >&2
        fi
    ' EXIT

    SKIP_GMSH=0

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --geo)
                GEO_FILE="$2"
                shift 2
                ;;

            --out)
                OUT_DIR="$2"
                MSH_FILE="$OUT_DIR/karman_vortex_prism_tet.msh"
                GMSH_LOG="$OUT_DIR/gmsh.log"
                shift 2
                ;;

            --solver)
                SOLVER_BIN="$2"
                shift 2
                ;;

            --skip-gmsh)
                SKIP_GMSH=1
                shift
                ;;

            -h|--help)
                usage
                exit 0
                ;;

            *)
                echo "Unknown option: $1" >&2
                usage
                exit 2
                ;;
        esac
    done


    mkdir -p "$OUT_DIR"


    # ============================================================
    # 1. Gmsh
    # ============================================================

    if [[ "$SKIP_GMSH" -eq 0 ]]; then

        command -v gmsh >/dev/null 2>&1 || {
            echo "ERROR: gmsh is not in PATH." >&2
            echo "Install Gmsh or use --skip-gmsh with an existing MSH2 file." >&2
            exit 127
        }

        echo
        echo "[1/5] Gmsh"
        echo "  GEO : $GEO_FILE"
        echo "  MSH : $MSH_FILE"
        echo "  LOG : $GMSH_LOG"
        echo

        gmsh \
            -3 \
            -format msh2 \
            "$GEO_FILE" \
            -o "$MSH_FILE" \
            2>&1 | tee "$GMSH_LOG"

    else

        echo
        echo "[1/5] Gmsh skipped"
        echo "  using: $MSH_FILE"

        [[ -f "$MSH_FILE" ]] || {
            echo "ERROR: mesh not found: $MSH_FILE" >&2
            exit 1
        }

    fi


    # ============================================================
    # 2. Convert Physical Groups
    # ============================================================

    echo
    echo "[2/5] Convert Physical Groups"

    CONVERT_ARGS=(
        "$MSH_FILE"
        --output-dir "$OUT_DIR"
    )

    if [[ "$PROJECT_CYLINDER_WALL" == "1" ]]; then

        CYL_RADIUS="$(
            awk -v d="$DIAMETER" \
                'BEGIN { printf "%.17g", 0.5*d }'
        )"

        echo "  project cylinder radius = $CYL_RADIUS"

        CONVERT_ARGS+=(
            --project-cylinder-radius "$CYL_RADIUS"
        )

    fi

    python3 \
        "$ROOT/mesh_io/save_physical_groups.py" \
        "${CONVERT_ARGS[@]}"


    # ============================================================
    # 3. Graph
    # ============================================================

    echo
    echo "[3/5] Build unified tetra+prism nodal graph"

    python3 \
        "$ROOT/mesh_io/build_graph.py" \
        --node "$OUT_DIR/node.dat" \
        --tet "$OUT_DIR/elem_tet.dat" \
        --prism "$OUT_DIR/elem_prism.dat" \
        --output "$OUT_DIR/graph.dat"


    # ============================================================
    # 4. Velocity BC
    # ============================================================

    echo
    echo "[4/5] Build velocity Dirichlet BC"

    python3 \
        "$ROOT/mesh_io/build_velocity_bc.py" \
        "$OUT_DIR" \
        --inlet-u "$INLET_U"


    # ============================================================
    # 5. Validation
    # ============================================================

    echo
    echo "[5/5] Validate mixed mesh"

    if [[ "$CHECK_CYLINDER_RADIUS" == "1" ]]; then

        python3 \
            "$ROOT/mesh_io/validate_mixed_mesh.py" \
            "$OUT_DIR" \
            --diameter "$DIAMETER" \
            --check-cylinder-radius

    else

        python3 \
            "$ROOT/mesh_io/validate_mixed_mesh.py" \
            "$OUT_DIR" \
            --diameter "$DIAMETER"

    fi


    # ============================================================
    # Compatibility aliases
    # ============================================================

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
    # Summary
    # ============================================================

    cat <<EOF

Generated mixed mesh:
  $OUT_DIR/node.dat
  $OUT_DIR/elem_tet.dat
  $OUT_DIR/elem_prism.dat
  $OUT_DIR/graph.dat
  $OUT_DIR/D_bc_v.dat
EOF

    if [[ -f "$OUT_DIR/surf.dat" ]]; then
        echo "  $OUT_DIR/surf.dat"
    fi


    # ============================================================
    # Optional solver
    # ============================================================

    if [[ -n "$SOLVER_BIN" ]]; then

        [[ -x "$SOLVER_BIN" ]] || {
            echo "ERROR: solver is not executable: $SOLVER_BIN" >&2
            exit 1
        }

        echo
        echo "========================================"
        echo "[solver]"
        echo "  executable : $SOLVER_BIN"
        echo "  directory  : $OUT_DIR"
        echo "========================================"
        echo

        # IMPORTANT:
        # Do NOT use exec here.
        "$SOLVER_BIN" "$OUT_DIR"

        solver_status=$?

        echo
        echo "========================================"
        echo "Solver finished"
        echo "exit status = $solver_status"
        echo "========================================"

        exit "$solver_status"

    else

        cat <<EOF

Solver was not launched.

To launch:
  SOLVER_BIN=/path/to/solver . ./run_karman.sh --skip-gmsh

or:
  . ./run_karman.sh --solver /path/to/solver
EOF

    fi
)


# ================================================================
# Run
# ================================================================

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    # bash run_karman.sh ... の場合
    run_karman_main "$@"
else
    # . run_karman.sh の場合
    # 呼び出し元シェルの $@ を引き継がない
    run_karman_main
fi