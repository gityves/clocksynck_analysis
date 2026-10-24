#!/usr/bin/env bash
# A133 V133.1 (simulation gate). Usage: run.sh [--subsample 1] [--seed S]. Capped at 30 GB through common/run.sh.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A133_axial_decomposition/scripts/compute_v133_sims.py "$@"
