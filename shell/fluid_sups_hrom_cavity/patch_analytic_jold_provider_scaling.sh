#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
target="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_jold.h"

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

old = """    const double previous_transient_trial =
        -density * N_j / dt;

    const double momentum_test =
        N_i + tau * density * advective_test;
"""
new = """    /*
     * Provider-consistent transient scaling: the linked residual receives
     * du_time = v - v_old as the raw time increment.
     */
    const double previous_transient_trial_momentum =
        -density * N_j;

    const double previous_transient_trial_pspg =
        -N_j;

    const double momentum_test =
        N_i + tau * density * advective_test;

    (void)dt;
"""

old2 = """    for(b = 0; b < 3; b++){
        J_old[b][b] =
            previous_transient_trial * momentum_test;

        J_old[3][b] =
            previous_transient_trial
            * tau
            * grad_N_i[b];
    }
"""
new2 = """    for(b = 0; b < 3; b++){
        J_old[b][b] =
            previous_transient_trial_momentum
            * momentum_test;

        J_old[3][b] =
            previous_transient_trial_pspg
            * tau
            * grad_N_i[b];
    }
"""

if old not in text:
    if "previous_transient_trial_momentum" in text:
        print("analytic J_old provider-scaling fix is already applied.")
        raise SystemExit(0)
    raise SystemExit("ERROR: expected old J_old coefficient block not found.")

if old2 not in text:
    raise SystemExit("ERROR: expected old J_old assignment block not found.")

backup = path.with_suffix(path.suffix + ".before_provider_scaling_fix")
if not backup.exists():
    shutil.copy2(path, backup)

text = text.replace(old, new, 1)
text = text.replace(old2, new2, 1)
path.write_text(text)

print(f"patched: {path}")
print(f"backup : {backup}")
PY

echo
echo "Rebuild:"
echo "  cd ${repo_root}/solvers/fluid_sups"
echo "  rm -f fluid_sups_core/fluid_sups_st_model.o fluid_sups_st_fom"
echo "  make -f Makefile_ST_cavity dev_build"
