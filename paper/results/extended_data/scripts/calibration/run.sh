#!/usr/bin/env bash
# Figure S01 (calibration). Usage: run.sh [--subsample N]   (full run by default; small_n and the coverage part take hours)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_small_n.py --jobs 4 "$@"
"$PY" compute_anchored_simulations.py --jobs 4 "$@"
"$PY" compute_vonmises_calibration.py "$@"
"$PY" compute_orpheus_bias.py "$@"
"$PY" compute_orpheus_cells.py --jobs 3 "$@"
"$PY" compute_age_acceptance.py --jobs 4 "$@"
