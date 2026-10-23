#!/usr/bin/env bash
# A109 (coupling nulls). Usage: run.sh [--scout | --subsample N]   (full scoring by default; do not run before the
# data-dependent part of A109 is committed). Heavy: capped at 30 GB through common/run.sh.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A109_coupling_nulls/scripts/compute_coupling_nulls.py "$@"
