#!/usr/bin/env bash
# Figure S05 (atlas fits). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_atlas_fits.py --jobs 3 "$@"
"$PY" compute_step7_confirmations.py "$@"   # after figure_S04 compute_coherence.py and organs
