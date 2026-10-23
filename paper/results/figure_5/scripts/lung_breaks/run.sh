#!/usr/bin/env bash
# Figure 7 (breaks). Usage: run.sh [--subsample N]   (full run by default; the permutation steps take about an hour)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_modules.py --jobs 4 "$@"
"$PY" compute_decomposition.py --jobs 3 "$@"
"$PY" compute_a18_followups.py --jobs 4 "$@"
