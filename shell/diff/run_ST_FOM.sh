#!/bin/bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2' ERR

# mesh
# 一方向分割数
e=30
# 解析領域の大きさ
ep=5

# 既存ディレクトリ命名との互換用。FOM-onlyではPODを実行しない。
num_modes=(10)
pa=0
st=1

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd "${script_dir}/../.." && pwd -P)"
cd "${project_root}"

for nm in "${num_modes[@]}"; do
    echo "============================================================"
    echo "ST-FOM only: e=${e} ep=${ep} nm(label only)=${nm}"
    echo "============================================================"

    # executable scriptはsourceしない。
    bash shell/diff/meshgen_ST.sh "$e" "$ep" "$nm" 4 4 "$pa"
    bash shell/diff/execution_ST_FOM.sh "$e" "$ep" "$nm" 4 4 "$pa" "$st"
done

echo
echo "ST-FOM driver completed."
