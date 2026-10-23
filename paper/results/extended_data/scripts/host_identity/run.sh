#!/usr/bin/env bash
# Figure S04 (host identity, anchors). Usage: run.sh [--subsample N]   (full run by default)
# compute_coherence.py first: compute_anchors.py and compute_host_identity.py read its tables (same output folder).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_coherence.py --jobs 4 "$@"
"$PY" compute_anchors.py "$@"
"$PY" compute_host_identity.py --jobs 4 "$@"
"$PY" compute_cortex_modules.py "$@"
"$PY" compute_group_anchors.py --jobs 4 "$@"
