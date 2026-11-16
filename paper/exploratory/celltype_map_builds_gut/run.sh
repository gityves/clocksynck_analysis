#!/usr/bin/env bash
# Second lane of the genome-wide Census builds (see ../celltype_map_builds/run.sh): the Kong 2023 gut atlases (ileum and
# colon: epithelial, immune and stromal datasets; groups also split by disease, as the panel builds) and the Natri 2024
# lung (groups also split by lineage, annotation_level_1 and Processing_site). Output: $SCCC_DATA/atlas_wg/<name>/.
#   RUN_CAP_GB=16 paper/common/run.sh paper/exploratory/celltype_map_builds_gut/run.sh
set -u
PY=/home/maxine/miniconda3/envs/census_env/bin/python
BUILD=/home/maxine/Documents/paychere/scCircadianCancer/legacy/desync/common/build_pseudobulk.py
WG=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg
build() {  # name dataset_id min_cells source-group columns...
  local name=$1 id=$2 mc=$3; shift 3
  if [ -f "$WG/$name/pseudobulk.h5ad" ]; then echo "$(date '+%F %T') skip $name (exists)"; return; fi
  mkdir -p "$WG/$name"
  echo "$(date '+%F %T') start $name ($id, min cells $mc, source groups: $*)"
  "$PY" -u "$BUILD" "$id" "$WG/$name" --threads 3 --min-cells "$mc" --source-group "$@" 2>&1 | awk '{ print strftime("%H:%M:%S"), $0; fflush() }'
  echo "$(date '+%F %T') end $name (exit ${PIPESTATUS[0]}), $(du -sh "$WG/$name" | cut -f1)"
}
build kong_co_str  6cf3634d-e911-44ad-bf52-c747a9af3c01 1 disease
build kong_co_epi  5ce42b38-d867-487f-9b40-e8bb00b21d0b 1 disease
build kong_co_imm  518d9049-2a76-44f8-8abc-1e2b59ab5ba1 1 disease
build kong_ti_str  0f4865d5-8000-4f68-8ac7-f5efea9e5e70 1 disease
build kong_ti_epi  fe4b89d5-461e-440c-a5a8-621b37b122c0 1 disease
build ileum_kong   a37f857c-779f-464e-9310-3db43a1811e7 1 disease
build natri2024_ild f14bc322-1322-4184-8d16-409557525ea5 10 lineage annotation_level_1 Processing_site
echo "$(date '+%F %T') all done"
