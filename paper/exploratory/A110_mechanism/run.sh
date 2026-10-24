#!/usr/bin/env bash
# A110 (what keeps organs coupled). Usage: run.sh [--scout | --subsample N] [--only h1 h2 h3 h4 h5]
# (full scoring by default; do not run before the scout is committed). Heavy: capped at 30 GB through common/run.sh.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A110_mechanism/scripts/compute_mechanism.py "$@"
