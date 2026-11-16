#!/usr/bin/env bash
# Genome-wide (all genes) Census group pseudobulks of PsychAD RADC and MSSM (cell-type program C3: third and fourth
# cortex cohorts for H274.3 and the cortex offset sign). Same route and minimum cells as psychad_hbcc in run.sh.
#   RUN_CAP_GB=24 paper/common/run.sh paper/exploratory/celltype_map_builds/run_psychad.sh
set -u
PY=/home/maxine/miniconda3/envs/census_env/bin/python
BUILD=/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/common/build_pseudobulk.py
WG=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg
build() {  # name dataset_id min_cells
  if [ -f "$WG/$1/pseudobulk.h5ad" ]; then echo "$(date '+%F %T') skip $1 (exists)"; return; fi
  mkdir -p "$WG/$1"
  echo "$(date '+%F %T') start $1 ($2, min cells $3)"
  "$PY" -u "$BUILD" "$2" "$WG/$1" --threads 6 --min-cells "$3" 2>&1 | awk '{ print strftime("%H:%M:%S"), $0; fflush() }'
  echo "$(date '+%F %T') end $1 (exit ${PIPESTATUS[0]}), $(du -sh "$WG/$1" | cut -f1)"
}
build psychad_radc  5e57cd50-8e42-42d6-940d-5c1660d06864 10
build psychad_mssm  37a17b78-4864-4a42-b67b-31c00962795a 10
echo "$(date '+%F %T') all done"
