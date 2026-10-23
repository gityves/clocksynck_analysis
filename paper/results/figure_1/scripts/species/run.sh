#!/usr/bin/env bash
# Figure S12 (species and known-time validations). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_colon_injury.py --jobs 4 "$@"            # A21 (CellBender counts)
"$PY" compute_colon_injury.py --jobs 4 --counts raw "$@"   # A21 sensitivity (Cell Ranger raw counts)
"$PY" compute_sheep.py --jobs 4 "$@"                   # A34
"$PY" compute_baboon.py --jobs 4 "$@"                  # A36, A40, A41
"$PY" compute_baboon_stability.py "$@"             # A36 stability, gate, leave-one-animal-out
"$PY" compute_sheep_study_gates.py "$@"            # A34 per-study gates
"$PY" compute_mouse_skin_timed.py "$@"             # A72
