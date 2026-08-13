#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
file="${repo_root}/solvers/fluid_sups/fluid_sups_core/core_HROM.c"

if [ ! -f "${file}" ]; then
    echo "ERROR: core_HROM.c not found: ${file}" >&2
    exit 2
fi

python3 - "${file}" <<'PY'
from pathlib import Path
import re
import shutil
import sys

path = Path(sys.argv[1])
text = path.read_text()
original = text

backup = path.with_suffix(path.suffix + ".before_st_compile_fix")
if not backup.exists():
    shutil.copy2(path, backup)

# ----------------------------------------------------------------------
# 1. hlpod_mat is used in solver_hrom_NR() but is not declared.
#    The existing fluid ROM/HROM code stores the matrix at sys->rom.hlpod_mat.
# ----------------------------------------------------------------------
m = re.search(
    r'(void\s+solver_hrom_NR\s*\([^)]*\)\s*\{)',
    text,
    flags=re.S,
)

if not m:
    raise SystemExit(
        "ERROR: solver_hrom_NR() was not found; core_HROM.c layout differs."
    )

# Limit the alias check to the function vicinity so declarations elsewhere
# do not suppress the required local alias.
near = text[m.end():m.end()+1200]
if not re.search(r'\bHLPOD_MAT\s*\*\s*hlpod_mat\b', near):
    insertion = (
        m.group(1)
        + "\n"
        + "    /* Canonical HLPOD storage owned by FE_SYSTEM. */\n"
        + "    HLPOD_MAT* hlpod_mat = &(sys->rom.hlpod_mat);"
    )
    text = text[:m.start()] + insertion + text[m.end():]

# ----------------------------------------------------------------------
# 2. nrm_v/nrm_dv already exist earlier in solver_hrom_NR().
#    Later lines assign state/update norms, so remove only the second type
#    specifier instead of creating a second variable with the same name.
# ----------------------------------------------------------------------
text, count_v = re.subn(
    r'\bdouble\s+nrm_v\s*=\s*sqrt\s*\(\s*norm_v\s*\)\s*;',
    'nrm_v = sqrt(norm_v);',
    text,
    count=1,
)

text, count_dv = re.subn(
    r'\bdouble\s+nrm_dv\s*=\s*sqrt\s*\(\s*norm_delta_v\s*\)\s*;',
    'nrm_dv = sqrt(norm_delta_v);',
    text,
    count=1,
)

if count_v != 1:
    raise SystemExit(
        "ERROR: expected exactly one 'double nrm_v = sqrt(norm_v);' declaration."
    )
if count_dv != 1:
    raise SystemExit(
        "ERROR: expected exactly one 'double nrm_dv = sqrt(norm_delta_v);' declaration."
    )

if text == original:
    print("core_HROM.c already appears patched.")
else:
    path.write_text(text)
    print(f"patched: {path}")
    print(f"backup : {backup}")
PY

echo
echo "Relevant patched locations:"
grep -n -E \
    'HLPOD_MAT\* hlpod_mat|nrm_v[[:space:]]*=[[:space:]]*sqrt\(norm_v\)|nrm_dv[[:space:]]*=[[:space:]]*sqrt\(norm_delta_v\)' \
    "${file}" || true
