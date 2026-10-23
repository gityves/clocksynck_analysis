#!/usr/bin/env bash
# Figure 4 (cell types). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_oscillator_gtex.py --jobs 4 "$@"
"$PY" compute_family_summary.py --jobs 4 "$@"
"$PY" compute_per_cell.py --jobs 2 "$@"
"$PY" compute_cancer_split.py --jobs 4 "$@"
"$PY" compute_positive_control.py "$@"
"$PY" compute_followups_a18.py --jobs 2 "$@"
