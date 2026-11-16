#!/usr/bin/env bash
# Genome-wide (all genes) Census group pseudobulks of the atlases that only have a clock-gene panel table, for the
# cell-type clock map (A272 follow-up). Same grouping and minimum cells as the panel build of each atlas (see
# $SCCC_DATA/atlas/<name>/PROVENANCE.md), no gene subset. Output: $SCCC_DATA/atlas_wg/<name>/ (pseudobulk.h5ad,
# donors.csv, PROVENANCE.md written by the builder). Processed counts streamed from the CELLxGENE Census 2025-11-08.
#   RUN_CAP_GB=24 paper/common/run.sh paper/exploratory/celltype_map_builds/run.sh
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
build psychad_hbcc  5c97eeeb-7e52-44b3-b010-b832b1f5424c 10
build optic_nerve   3c364ffe-d599-44d5-a77f-6c5fd7038cf6 1
build eye_anterior  978a566a-dc27-478d-b306-26daba116c1f 1
build kpmp_sc       dea717d4-7bc0-4e46-950f-fd7e1cc8df7d 1
build islets_t2d    e51bae9a-c747-4b64-904a-4da7cda218ab 1
build seaad_dlpfc   6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3 10
build seaad_mtg     c2876b1b-06d8-4d96-a56b-5304f815b99a 10
echo "$(date '+%F %T') all done"
