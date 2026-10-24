#!/usr/bin/env bash
# A135 part 1 (control-only brain reference choice). Usage: run.sh [--subsample N] [--jobs J]   (full run by default).
# Capped at 30 GB through common/run.sh; at most 4 jobs.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A135_brain_reference/scripts/compute_brain_reference.py --jobs 4 "$@"
