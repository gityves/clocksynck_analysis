#!/usr/bin/env bash
# A132 (brain death or ventilation). Usage: run.sh [--scout | --subsample N]   (full scoring by default; do not run
# before the data-dependent part of A132 is committed). Heavy: capped at 30 GB through common/run.sh.
set -euo pipefail
cd "$(dirname "$0")/../../.."
RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh paper/exploratory/A132_brain_death/scripts/compute_brain_death.py "$@"
