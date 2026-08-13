#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
makefile="${repo_root}/solvers/fluid_sups/Makefile_ST_cavity"

if [ ! -f "${makefile}" ]; then
    echo "ERROR: Makefile_ST_cavity not found: ${makefile}" >&2
    exit 2
fi

python3 - "${makefile}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()

if "ACCELERATE_CXXFLAGS := -flax-vector-conversions" not in text:
    old = """CFLAGS += -DACCELERATE_NEW_LAPACK=1
FLUID_CBLAS_LIBS := -framework Accelerate"""
    new = """CFLAGS += -DACCELERATE_NEW_LAPACK=1
ACCELERATE_CXXFLAGS := -flax-vector-conversions
FLUID_CBLAS_LIBS := -framework Accelerate"""
    if old not in text:
        raise SystemExit("ERROR: Darwin Accelerate block not found.")
    text = text.replace(old, new, 1)

rule = """
./fluid_sups_core/core_NR.o: ./fluid_sups_core/core_NR.c
\t$(CC) $(CFLAGS) $(ACCELERATE_CXXFLAGS) \\
\t\t$(INCLUDES) $(INCLUDES_BB) $(INCLUDES_ROM) \\
\t\t-c $< -o $@

"""

if "./fluid_sups_core/core_NR.o:" not in text:
    generic = "%.o: %.c\n"
    pos = text.find(generic)
    if pos < 0:
        raise SystemExit("ERROR: generic %.o rule not found.")
    text = text[:pos] + rule + text[pos:]

path.write_text(text)
print(f"patched: {path}")
PY

echo
echo "Relevant Makefile settings:"
grep -n -E 'ACCELERATE_NEW_LAPACK|ACCELERATE_CXXFLAGS|core_NR\.o:' "${makefile}" || true
