#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/core_NR.c"

python3 - "${src}" <<'PY'
from pathlib import Path
import re, shutil, sys

path = Path(sys.argv[1])
text = path.read_text()

# Strip C/C++ comments approximately for dependency classification.
scan = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
scan = re.sub(r'//.*', '', scan)

# Ignore the include itself while looking for actual MKL/LAPACKE/CBLAS usage.
scan_without_include = re.sub(
    r'^\s*#\s*include\s*[<"]mkl\.h[>"]\s*$',
    '',
    scan,
    flags=re.M
)

patterns = {
    "MKL-specific": r'\b(?:MKL_[A-Za-z0-9_]*|mkl_[A-Za-z0-9_]*)\b',
    "LAPACKE": r'\bLAPACKE_[A-Za-z0-9_]+\b',
    "CBLAS": r'\bcblas_[A-Za-z0-9_]+\b',
}

found = {
    name: sorted(set(re.findall(pattern, scan_without_include)))
    for name, pattern in patterns.items()
}

active = {k:v for k,v in found.items() if v}
if active:
    print("core_NR.c contains active linear-algebra API dependencies:")
    for kind, names in active.items():
        print(f"  {kind}:")
        for name in names:
            print(f"    {name}")
    print()
    print("No automatic source edit was made.")
    print("Run diagnose_core_NR_linalg.sh and convert the actual calls explicitly.")
    raise SystemExit(3)

matches = list(re.finditer(
    r'^\s*#\s*include\s*[<"]mkl\.h[>"]\s*$',
    text,
    flags=re.M
))
if not matches:
    print("No #include <mkl.h> remains. Nothing to patch.")
    raise SystemExit(0)

backup = path.with_suffix(path.suffix + ".before_mkl_include_removal")
if not backup.exists():
    shutil.copy2(path, backup)

text = re.sub(
    r'^\s*#\s*include\s*[<"]mkl\.h[>"]\s*\n?',
    '',
    text,
    flags=re.M
)
path.write_text(text)

print(f"Removed stale mkl.h include from: {path}")
print(f"Backup: {backup}")
PY
