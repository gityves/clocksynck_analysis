#!/usr/bin/env bash
# sequential (one GPU job at a time), CoPhaser env. Reuses the pilot's fit_cellcycle.py.
# usage: run_all.sh "STUDY [extra args]"...   e.g. run_all.sh "GSE131907 --extra --seeds 5,6,7,8,9" GSE161529
cd "$(dirname "$0")/../../../../.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1 RUN_CAP_GB=30
OUT=paper/exploratory/C4_single_cell/fits/output
for spec in "$@"; do
  set -- $spec; s=$1; shift
  tag=$s; [[ "$*" == *--extra* ]] && tag=${s}_extra
  paper/common/run.sh paper/exploratory/C4_single_cell/pilot/scripts/fit_cellcycle.py "$s" --out $OUT "$@"
  tr "\r" "\n" < paper/results/logs/fit_cellcycle.log | grep -v "^Training:" > $OUT/log_$tag.txt
done
