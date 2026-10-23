#!/usr/bin/env bash
# Figure S11 (natural experiments). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_covid_gates.py "$@"
"$PY" compute_covid_couplings.py "$@"
"$PY" compute_covid_ventilation.py "$@"
"$PY" compute_fetal_secondary.py "$@"
"$PY" compute_placenta_sensitivity.py "$@"
"$PY" compute_gtex_reference.py --jobs 4 "$@"
