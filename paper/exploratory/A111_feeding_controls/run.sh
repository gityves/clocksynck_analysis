#!/usr/bin/env bash
# A111 / A112 / A129 (WP3 feeding and genetic controls). Usage: run.sh [--subsample N] [--only a111 a112 clockloss]
# --subsample N: quick end-to-end check (fewer permutations), outputs in output/subsample/. Capped at 30 GB (common/run.sh).
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A111_feeding_controls/scripts/compute_feeding_controls.py "$@"
