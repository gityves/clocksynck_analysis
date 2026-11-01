#!/usr/bin/env bash
# sequential (one GPU job at a time) CoPhaser env; usage: run_all.sh STUDY...
cd "$(dirname "$0")/../../../../.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1 RUN_CAP_GB=30
for s in "$@"; do paper/common/run.sh paper/exploratory/C4_single_cell/pilot/scripts/fit_cellcycle.py "$s"; tr "\r" "\n" < paper/results/logs/fit_cellcycle.log | grep -v "^Training:" > paper/exploratory/C4_single_cell/pilot/output/log_$s.txt; done
