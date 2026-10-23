#!/usr/bin/env bash
# Figure S13 (simulations and validations). Usage: run.sh [--subsample N]   (full runs: V16 / V18 about 1 h, V13 and R26 hours)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_hierarchy_validation.py --jobs 4 "$@"
"$PY" compute_celltype_truth.py --jobs 4 "$@"
"$PY" compute_a22_colon.py "$@"
"$PY" compute_r26_cohorts.py --jobs 4 "$@"
