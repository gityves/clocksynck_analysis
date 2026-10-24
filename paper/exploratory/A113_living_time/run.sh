#!/usr/bin/env bash
# A113 and A130 (WP4). Usage: run.sh [--subsample N] [--seed S]   (full scoring by default). Capped at 30 GB through common/run.sh.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A113_living_time/scripts/compute_living_time.py "$@"
