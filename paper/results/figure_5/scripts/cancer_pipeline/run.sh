#!/usr/bin/env bash
# Figure S10 (cancer details: the cohort pipeline). Usage: run.sh [--subsample N]   (full run by default: hours)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_cohort_pipeline.py "$@"
