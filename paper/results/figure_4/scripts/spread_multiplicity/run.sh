#!/usr/bin/env bash
# Figure S06 (spread, hierarchy, multiplicity). Usage: run.sh [--subsample N]   (full run by default)
# compute_multiplicity.py reads resilience/a17_contrast.csv: run that folder's
# compute_contrasts_a17.py first (with the same --subsample).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_spread_inputs.py --jobs 4 "$@"
"$PY" compute_hierarchy.py --jobs 4 "$@"
"$PY" compute_multiplicity.py "$@"
"$PY" compute_phase_cancer.py --jobs 4 "$@"
