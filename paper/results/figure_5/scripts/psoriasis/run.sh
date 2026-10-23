#!/usr/bin/env bash
# Figure 6 (psoriasis). Usage: run.sh [--subsample N]   (full run by default)
# compute_replication.py first: the later scripts read its per-series tables (gated series, weights) from the same
# output folder (compute_cues_gtex.py too; it needs the S07 cue reference file, gtex_cues).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_replication.py --jobs 4 "$@"
for s in compute_dose compute_treatment compute_response compute_cues compute_cues_gtex compute_sc_status; do
  "$PY" $s.py --jobs 4 "$@"
done
