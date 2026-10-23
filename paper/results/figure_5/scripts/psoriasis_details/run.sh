#!/usr/bin/env bash
# Figure S09 (psoriasis details). Usage: run.sh [--subsample N]   (full run by default)
# compute_other_diseases.py before compute_cues_companion.py (the companion reads its gated series and weights).
# compute_sc_attribution.py reads the A54 per-gene shifts and compute_cues_matched.py / compute_cue_context.py the gated
# A54 series of psoriasis (run psoriasis/run.sh first, same --subsample). The TCGA arms of
# compute_cues_diseases / compute_cues_matched / compute_cue_context read cancer/compute_tcga_paired.py outputs
# (a52_per_cancer.csv; skipped with a message when missing); the cue reference and selection tables are outputs of
# gtex_cues/compute_cue_reference.py.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
for s in compute_other_diseases compute_cues_companion compute_layers compute_timed_layers compute_sc_attribution \
  compute_differentiation compute_bin_mix compute_taurus_amplitude compute_gtex_skin compute_cues_diseases compute_cues_matched \
  compute_cue_context; do
  "$PY" $s.py --jobs 4 "$@"
done
