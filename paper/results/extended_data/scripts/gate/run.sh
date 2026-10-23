#!/usr/bin/env bash
# Figure S03 (gate power and failures, waveform sensitivity, extended clock). Usage: run.sh [--subsample N]
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_gate_power.py --jobs 4 "$@"
"$PY" compute_extended_clock.py --jobs 4 "$@"
"$PY" compute_baboon_reference.py --jobs 4 "$@"
"$PY" compute_waveform_sensitivity.py --jobs 4 "$@"
"$PY" compute_gene_screen.py "$@"   # after figure_S07 compute_internal_external.py and compute_dip_arm.py
"$PY" compute_extended_e4.py "$@"
"$PY" compute_sheep_reference_primary.py "$@"   # after figure_S12 compute_sheep.py
