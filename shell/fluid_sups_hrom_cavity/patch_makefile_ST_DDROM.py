#!/usr/bin/env python3
from pathlib import Path
import shutil
import sys

if len(sys.argv) != 2:
    raise SystemExit(
        "usage: patch_makefile_st_ddrom.py <repository-root>"
    )

root = Path(sys.argv[1]).resolve()
makefile = root / "solvers/fluid_sups/Makefile_ST_cavity"
fragment = root / "solvers/fluid_sups/Makefile_ST_DDROM.fragment"

if not makefile.is_file():
    raise SystemExit(f"ERROR: not found: {makefile}")
if not fragment.is_file():
    raise SystemExit(f"ERROR: not found: {fragment}")

text = makefile.read_text()
marker = "# Fluid nonlinear all-at-once ST-DDROM"

if marker in text:
    print("Makefile_ST_cavity already contains the ST-DDROM target.")
    raise SystemExit(0)

backup = makefile.with_suffix(
    makefile.suffix + ".before_st_ddrom")
if not backup.exists():
    shutil.copy2(makefile, backup)

with makefile.open("a") as fp:
    fp.write("\n")
    fp.write(fragment.read_text())
    fp.write("\n")

print(f"patched: {makefile}")
print(f"backup : {backup}")
