#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/core_NR.c"
header="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_linalg_compat.h"

if [ ! -f "${src}" ]; then
    echo "ERROR: not found: ${src}" >&2
    exit 2
fi

if [ ! -f "${header}" ]; then
    echo "ERROR: compatibility header not found: ${header}" >&2
    exit 2
fi

python3 - "${src}" <<'PY'
from pathlib import Path
import re
import shutil
import sys

path = Path(sys.argv[1])
text = path.read_text()

backup = path.with_suffix(path.suffix + ".before_mkl_portability")
if not backup.exists():
    shutil.copy2(path, backup)

mkl_include = re.compile(
    r'^[ \t]*#[ \t]*include[ \t]*[<"]mkl\.h[>"][ \t]*$',
    flags=re.M,
)

compat_include = '#include "fluid_linalg_compat.h"'

if compat_include in text:
    print("core_NR.c already uses fluid_linalg_compat.h")
elif mkl_include.search(text):
    text, count = mkl_include.subn(compat_include, text, count=1)
    if count != 1:
        raise SystemExit("ERROR: failed to replace exactly one mkl.h include.")
    path.write_text(text)
    print(f"patched: {path}")
    print(f"backup : {backup}")
else:
    raise SystemExit(
        "ERROR: neither <mkl.h> nor fluid_linalg_compat.h was found. "
        "core_NR.c layout differs from the diagnosed source."
    )
PY

echo
echo "Linear algebra includes:"
grep -n -E 'mkl\.h|fluid_linalg_compat\.h|Accelerate|cblas\.h' "${src}" || true

echo
echo "Remaining MKL call sites are compatibility macros, not a hard MKL link:"
grep -n -E '\bmkl_(malloc|free)\b' "${src}" | head -n 20 || true
