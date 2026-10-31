#!/usr/bin/env bash
# Extended Data, V5 spike-ins in real biopsies (SRP445031). Usage: run.sh [--subsample N]   (full run: 280 fits, about 1 h on 4 workers)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_v5_spikeins.py --jobs 4 "$@"
