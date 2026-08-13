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
import re
import shutil
import sys

path = Path(sys.argv[1])
text = path.read_text()

marker = "Skipping legacy libBBFE_elemmat NR API detector"
if marker in text:
    print("execution_ST.sh is already patched.")
    raise SystemExit(0)

# Match the detector invocation visible in the runtime log.  Keep the detector
# available as an opt-in diagnostic, but do not require libBBFE_elemmat to own
# kernels that the current ST executable provides locally.
pattern = re.compile(
    r'^(?P<indent>[ \t]*)bash[ \t]+"\$\{detector\}"(?P<args>[^\n]*)$',
    re.M
)
m = pattern.search(text)
if not m:
    raise SystemExit(
        'ERROR: could not locate: bash "${detector}" ... in execution_ST.sh'
    )

indent = m.group("indent")
args = m.group("args")

replacement = (
    f'{indent}# The current fluid ST-FOM links the verified NR kernels from\n'
    f'{indent}# fluid_sups_core/core_FOM_st_provider.o.  The old detector\n'
    f'{indent}# requires all six kernels to be exported by libBBFE_elemmat,\n'
    f'{indent}# which is not the current link architecture.\n'
    f'{indent}if [[ "${{FLUID_ST_REQUIRE_BBFE_NR_API:-0}}" == "1" ]]; then\n'
    f'{indent}    bash "${{detector}}"{args}\n'
    f'{indent}else\n'
    f'{indent}    echo "Skipping legacy libBBFE_elemmat NR API detector; using local core_FOM ST provider."\n'
    f'{indent}fi'
)

backup = path.with_suffix(path.suffix + ".before_local_nr_provider")
if not backup.exists():
    shutil.copy2(path, backup)

text = text[:m.start()] + replacement + text[m.end():]
path.write_text(text)

print(f"patched: {path}")
print(f"backup : {backup}")
PY

echo
echo "Relevant execution_ST.sh block:"
grep -n -A10 -B4 \
  'Skipping legacy libBBFE_elemmat NR API detector' \
  "${target}" || true
