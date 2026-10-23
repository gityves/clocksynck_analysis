#!/usr/bin/env bash
# Figure 6 (cancer). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_tcga_paired.py "$@"
"$PY" compute_tumour_gate.py "$@"
"$PY" compute_decomposition.py "$@"
"$PY" compute_malignant_gate.py "$@"
"$PY" compute_delay.py "$@"
"$PY" compute_proliferation.py "$@"
"$PY" compute_anchored_cancer.py --jobs 4 "$@"
# Cancer follow-ups (work package of 2026-10-22). Needs host_identity/compute_coherence.py and
# gtex_cues/compute_cue_reference.py run first (same --subsample); --build-missing builds absent cell tables.
"$PY" compute_anchors_cancer.py --jobs 4 "$@"
"$PY" compute_cycling_a24.py --jobs 4 "$@"
"$PY" compute_cycling_a24.py --jobs 4 --definition quartile "$@"
"$PY" compute_cue_following.py "$@"
"$PY" compute_matched_cues_tcga.py "$@"
