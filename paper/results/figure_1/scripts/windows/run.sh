#!/usr/bin/env bash
# Figure S02 (sampling windows, internal time). Usage: run.sh [--subsample N]   (full run by default; the cohort
# bootstraps and the V25.1 grid take hours)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_window_check.py "$@"
"$PY" compute_internal_time.py --jobs 4 "$@"
"$PY" compute_baboon_windows.py --jobs 4 "$@"
"$PY" compute_chen_windows.py "$@"
"$PY" compute_breast_timed.py "$@"
