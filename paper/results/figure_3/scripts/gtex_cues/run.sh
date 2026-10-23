#!/usr/bin/env bash
# Figure S07 (GTEx covariates, cues, internal vs external time). Usage: run.sh [--subsample N]   (full run by default)
# compute_cue_reference.py first: it writes the cue reference (a81b_cue_reference.csv) and the 26-tissue selection tables.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_cue_reference.py "$@"
"$PY" compute_cue_hierarchy.py "$@"
"$PY" compute_internal_external.py "$@"
"$PY" compute_drug_targets.py "$@"
"$PY" compute_chiral_age.py "$@"
"$PY" compute_covariate_sensitivity.py "$@"
"$PY" compute_ventilation_dose.py "$@"
"$PY" compute_metabolic.py "$@"
"$PY" compute_coupling_sources.py "$@"
"$PY" compute_glucocorticoid.py "$@"
"$PY" compute_gonogo_diagnostics.py --jobs 3 "$@"
"$PY" compute_dip_arm.py --jobs 2 "$@"
"$PY" compute_gc_atlas.py "$@"
"$PY" compute_gc_li_breast.py "$@"
"$PY" compute_chiral_prep_counts.py --jobs 2 "$@"
