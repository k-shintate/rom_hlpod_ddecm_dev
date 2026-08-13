#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
target="${repo_root}/shell/fluid_sups_hrom_cavity/execution_ST.sh"

if [ ! -f "${target}" ]; then
    echo "ERROR: not found: ${target}" >&2
    exit 2
fi

python3 - "${target}" <<'PY'
from pathlib import Path
import shutil
import sys

path = Path(sys.argv[1])
text = path.read_text()

old = """if [ "${nr_api_linkage}" = "missing" ] || [ "${nr_api_linkage}" = "mixed" ]; then
    # The current fluid ST-FOM links the verified NR kernels from
    # fluid_sups_core/core_FOM_st_provider.o.  The old detector
    # requires all six kernels to be exported by libBBFE_elemmat,
    # which is not the current link architecture.
    if [[ "${FLUID_ST_REQUIRE_BBFE_NR_API:-0}" == "1" ]]; then
        bash "${detector}" ../../../test_thermal
    else
        echo "Skipping legacy libBBFE_elemmat NR API detector; using local core_FOM ST provider."
    fi
    exit 1
fi
"""

new = """if [ "${nr_api_linkage}" = "missing" ] || [ "${nr_api_linkage}" = "mixed" ]; then
    # The current fluid ST-FOM links the verified NR solve-path kernels from
    # fluid_sups_core/core_FOM_st_provider.o, not from libBBFE_elemmat.
    if [[ "${FLUID_ST_REQUIRE_BBFE_NR_API:-0}" == "1" ]]; then
        bash "${detector}" ../../../test_thermal
        exit 1
    else
        echo "Skipping legacy libBBFE_elemmat NR API detector; using local core_FOM ST provider."
        nr_api_linkage="local_core_FOM_cpp"
    fi
fi
"""

if old in text:
    backup = path.with_suffix(path.suffix + ".before_stale_exit_fix")
    if not backup.exists():
        shutil.copy2(path, backup)
    text = text.replace(old, new, 1)
    path.write_text(text)
    print(f"patched: {path}")
    print(f"backup : {backup}")
elif 'nr_api_linkage="local_core_FOM_cpp"' in text:
    print("execution_ST.sh is already corrected.")
else:
    raise SystemExit(
        "ERROR: expected preflight block was not found exactly. "
        "Please inspect execution_ST.sh lines 30-50."
    )
PY

echo
echo "Corrected preflight block:"
sed -n '28,58p' "${target}"
