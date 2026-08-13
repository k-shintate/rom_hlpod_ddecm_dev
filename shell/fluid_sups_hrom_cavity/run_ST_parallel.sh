#!/usr/bin/env bash
set -Eeuo pipefail

# Compatibility wrapper.
# The spatial domain decomposition is now explicitly prepared with the same
# GEDATSU workflow as the convection-diffusion ST-DDROM launcher.
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_ST_GEDATSU_MPI.sh" "$@"
