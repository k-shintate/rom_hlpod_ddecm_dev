#!/usr/bin/env bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2; exit "$rc"' ERR

if [[ "$#" -ne 6 ]]; then
    echo "usage: $0 <e> <domain-size> <nm> <nd> <np> <pa>" >&2
    exit 2
fi

e="$1"; ep="$2"; nm="$3"; nd="$4"; np="$5"; pa="$6"
: "$ep" "$pa"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
default_result_dir="${project_root}/result_fluid_sups_cavity/${nm}-${np}-${nd}"
result_dir="${FLUID_ST_RESULT_DIR:-$default_result_dir}"
test_thermal="${TEST_THERMAL_ROOT:-${project_root}/../test_thermal}"

rm -rf "$result_dir"
mkdir -p "$result_dir"
cd "$result_dir"

"${test_thermal}/bin/cmd2cond" "#density_array" double 1 100.0
mv cond.dat density.dat
"${test_thermal}/bin/cmd2cond" "#viscosity_array" double 1 1
mv cond.dat viscosity.dat
"${test_thermal}/bin/cmd2cond" "#target_density" double 1 100.0
mv cond.dat target_density.dat
"${test_thermal}/bin/cmd2cond" "#target_viscosity" double 1 1
mv cond.dat target_viscosity.dat
"${test_thermal}/bin/cmd2cond" \
    "#snapshot_interval" int 1 1 \
    "#rom_finish_time" double 1 0.5 \
    "#rom_output_interval" int 1 1
mv cond.dat rom_cond.dat
"${test_thermal}/bin/cmd2cond" "#inc_svd_interval" int 1 50
mv cond.dat hrom_cond.dat
"${test_thermal}/bin/cmd2cond" \
    "#time_spacing" double 1 "${FLUID_ST_DT:-0.001}" \
    "#output_interval" int 1 "${FLUID_ST_OUTPUT_INTERVAL:-10}" \
    "#finish_time" double 1 "${FLUID_ST_FINISH_TIME:-0.02}"

"${test_thermal}/bin/meshgen_hex" "$e" "$e" "$e" 1.0 1.0 1.0
"${test_thermal}/submodule/monolis/bin/monolis_extract_all_surf_hex"
"${test_thermal}/bin/mesh_surf_extract" 0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_z.dat -ov surf_z.vtk
"${test_thermal}/bin/mesh_surf_remove"  0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_wall.dat -ov surf_wall.vtk
"${test_thermal}/bin/surf_dbc" 3 1.0 0.0 0.0 -ie surf_z.dat -o D_bc_z.dat
"${test_thermal}/bin/surf_dbc" 3 0.0 0.0 0.0 -ie surf_wall.dat -o D_bc_wall.dat
"${test_thermal}/bin/surf_bc_merge" D_bc_z.dat D_bc_wall.dat -o D_bc_v.dat
"${test_thermal}/bin/mesh_surf_extract" 0.48 0.48 0.0 0.52 0.52 0.0 -oe surf_p.dat -ov surf_p.vtk
"${test_thermal}/bin/surf_dbc" 1 0.0 -ie surf_p.dat -o D_bc_p.dat
