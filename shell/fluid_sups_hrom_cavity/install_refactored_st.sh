#!/usr/bin/env bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
mode="${2:-copy}"
package_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
solver_src="${package_root}/solvers/fluid_sups"
solver_dst="${repo_root}/solvers/fluid_sups"
shell_src="${package_root}/shell/fluid_sups_hrom_cavity"
shell_dst="${repo_root}/shell/fluid_sups_hrom_cavity"
[[ -d "$solver_dst" && -d "$shell_dst" ]] || { echo "ERROR: repository layout not found" >&2; exit 2; }
stamp="$(date +%Y%m%d_%H%M%S)"
backup="${repo_root}/refactor_backup_${stamp}"
mkdir -p "$backup/solvers/fluid_sups" "$backup/shell"
cp -a "$solver_dst/fluid_sups_core" "$backup/solvers/fluid_sups/"
[[ -d "$solver_dst/fluid_sups_main_cavity" ]] && cp -a "$solver_dst/fluid_sups_main_cavity" "$backup/solvers/fluid_sups/"
cp -a "$shell_dst" "$backup/shell/"
[[ -f "$solver_dst/Makefile_ST_cavity" ]] && cp "$solver_dst/Makefile_ST_cavity" "$backup/solvers/fluid_sups/"

# Copy only refactored ST/build files. Ordinary core_FOM/ROM/HROM/NR sources are untouched.
mkdir -p "$solver_dst/fluid_sups_core" "$solver_dst/fluid_sups_main_cavity" "$shell_dst" "$solver_dst/st_refactor_tools"
cp -a "$solver_src/fluid_sups_core/fluid_sups_st_model.c" "$solver_src/fluid_sups_core/fluid_sups_st_model.h" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/fluid_sups_st_model" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/st_fom_common.c" "$solver_src/fluid_sups_core/st_fom_common.h" "$solver_src/fluid_sups_core/st_fom_common" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/rom_std_stddrom.c" "$solver_src/fluid_sups_core/rom_std_stddrom.h" "$solver_src/fluid_sups_core/rom_std_stddrom" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/rom_std_stddrom_fluid.c" "$solver_src/fluid_sups_core/rom_std_stddrom_fluid.h" "$solver_src/fluid_sups_core/rom_std_stddrom_fluid" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/rom_std_stddrom_hlpod_bridge.c" "$solver_src/fluid_sups_core/rom_std_stddrom_hlpod_bridge.h" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_core/fluid_sups_jold.h" "$solver_src/fluid_sups_core/fluid_sups_nr_api_compat.h" "$solver_src/fluid_sups_core/generate_nr_api_from_core_fom.py" "$solver_dst/fluid_sups_core/"
cp -a "$solver_src/fluid_sups_main_cavity/." "$solver_dst/fluid_sups_main_cavity/"
cp "$solver_src/Makefile_ST_cavity" "$solver_dst/Makefile_ST_cavity"
cp -a "$solver_src/st_refactor_tools/." "$solver_dst/st_refactor_tools/"

# Preserve the three ordinary POD-Galerkin scripts; install the clean ST set.
for f in st_case.env st_pipeline.sh run_ST.sh run_ST_DDROM.sh run_validated_case.sh verify_validated_case.py validated_thresholds.json inspect_ST_DDROM_snapshots.py meshgen.sh; do cp "$shell_src/$f" "$shell_dst/$f"; done
chmod +x "$shell_dst"/*.sh "$shell_dst"/*.py || true

if [[ "$mode" == "prune" ]]; then
  find "$solver_dst/fluid_sups_core" -maxdepth 1 -type f \( -name '*.o' -o -name '*.d' -o -name '*.before_*' -o -name '*.bak' -o -name '*.cpy' \) -delete
  rm -f "$solver_dst/fluid_sups_core/core_FOM_ST.c" "$solver_dst/fluid_sups_core/core_FOM_ST.h" "$solver_dst/fluid_sups_core/fluid_sups_nr_refactored.c"
  keep='^(execution\.sh|execution_offline\.sh|execution_online\.sh|merge_graph\.sh|meshgen\.sh|st_case\.env|st_pipeline\.sh|run_ST\.sh|run_ST_DDROM\.sh|run_validated_case\.sh|verify_validated_case\.py|validated_thresholds\.json|inspect_ST_DDROM_snapshots\.py)$'
  find "$shell_dst" -maxdepth 1 -type f | while read -r f; do b="$(basename "$f")"; [[ "$b" =~ $keep ]] || rm -f "$f"; done
fi

echo "refactor installed; backup: $backup"
echo "build: cd $solver_dst && make -f Makefile_ST_cavity dev_build"
echo "static check: cd $solver_dst && make -f Makefile_ST_cavity verify_refactor"
