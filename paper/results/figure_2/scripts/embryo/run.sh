#!/usr/bin/env bash
# Figure 2 (development). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_embryo_timed.py "$@"
"$PY" compute_embryo_emergence.py --jobs 4 "$@"
"$PY" compute_embryo_lineages.py "$@"
"$PY" compute_fetal_atlases.py "$@"
"$PY" compute_placenta.py "$@"
