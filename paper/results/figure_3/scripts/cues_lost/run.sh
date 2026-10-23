#!/usr/bin/env bash
# Figure 4 (cues lost: ventilated and slowly dying donors). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_walltime_hosts.py "$@"
"$PY" compute_ventilator_tissues.py --jobs 4 "$@"
"$PY" compute_ventilator_gc.py --jobs 4 "$@"
# WP1-WP3 (A109-A112, docs/WORKPLAN.md) add their scripts here.
