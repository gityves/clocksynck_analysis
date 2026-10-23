#!/usr/bin/env bash
# Figure 5 (ageing and resilience). Usage: run.sh [--subsample N]   (full run by default)
# Not in the default run (hours): `compute_contrasts_a17.py --mde`, `compute_procurement.py --mde`,
# `compute_masld_gut.py --mde` (minimum detectable effects).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_contrasts_a17.py --jobs 4 "$@"
"$PY" compute_procurement.py "$@"
"$PY" compute_age.py --jobs 4 "$@"
"$PY" compute_masld_gut.py --jobs 4 "$@"
"$PY" compute_mouse_ageing.py "$@"
"$PY" compute_gtex_age.py "$@"
"$PY" compute_gtex_chiral_age.py --jobs 4 "$@"
"$PY" compute_gtex_sex_bmi.py "$@"
"$PY" compute_a45_learned.py "$@"
"$PY" compute_gtex_age_tiebreak.py "$@"
"$PY" compute_gtex_age_continuous.py "$@"
