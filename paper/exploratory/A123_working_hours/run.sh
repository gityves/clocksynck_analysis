#!/usr/bin/env bash
# A123 (GTEx organ map from working-hours deaths only). Usage: run.sh [--subsample N] [--boot B] [--perm P]. Capped at 30 GB.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A123_working_hours/scripts/compute_working_hours.py "$@"
