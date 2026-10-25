#!/usr/bin/env bash
# Figure S08 (gene rhythms of cell types). Usage: run.sh [--subsample N]   (full run by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_gtex_truth.py "$@"                       # answer key first (the others read its tables)
"$PY" compute_atlas_gene_rhythms.py "$@"               # A27 part 1, A31b
"$PY" compute_donor_curve.py --jobs 3 "$@"             # A39
"$PY" compute_gene_map.py --jobs 4 "$@"                # A50 (--n-boot 200 by default in a full run)
"$PY" compute_blood_clock.py "$@"                      # A32, A35
"$PY" compute_atlas_gene_rhythms_modes.py --jobs 3 "$@"   # A27 / A31b covariate modes, composition, z > 2
