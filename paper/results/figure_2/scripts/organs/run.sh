#!/usr/bin/env bash
# Figure 3 (organs). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_organs_timed.py "$@"
"$PY" compute_brain_artefact.py "$@"
"$PY" compute_brain_artefact_d4.py "$@"
"$PY" compute_brain_periphery_species.py "$@"
"$PY" compute_chen_age.py "$@"
"$PY" compute_chen_age.py --adjust ph rin pmi "$@"
# Not in the default run (hours): compute_species_rstar_full.py --jobs 4  (A42 readout 1 / A103.3 with the registered bootstraps, in parallel)
