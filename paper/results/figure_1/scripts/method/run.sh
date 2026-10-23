#!/usr/bin/env bash
# Figure 1 (method and truth validations). Usage: run.sh [--subsample N]   (full run by default; the GTEx benchmark and CYCLOPS take hours)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$(conda info --base)/envs/sccc/bin/python}
"$PY" compute_gtex_benchmark.py --ncores 4 "$@"
"$PY" compute_celltype_truth.py "$@"
"$PY" compute_celltype_truth.py --arms C-DC,R-PC --tag sens_RPC "$@"
"$PY" compute_celltype_truth.py --window "$@"
"$PY" compute_baboon_truth.py "$@"
"$PY" compute_skin_truth.py "$@"
"$PY" compute_tcga_cofe.py --ncores 4 "$@"
"$PY" compute_cyclops_anchoring.py --ncores 4 "$@"
