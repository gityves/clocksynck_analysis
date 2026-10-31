#!/usr/bin/env bash
# A172: download recount3 TCGA gene sums and metadata of the projects not among the 11 A97 projects (same mirror and layout as
# $SCCC_DATA/tcga/recount3/, see its PROVENANCE.md). Writes MANIFEST_A172.tsv (path, url, size, sha256). Nothing downloaded is executed.
set -euo pipefail
R=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga/recount3
BASE=https://duffel.rail.bio/recount3/human/data_sources/tcga
M=$R/MANIFEST_A172.tsv
[ -f "$M" ] || printf "path\turl\tsize_bytes\tsha256\n" > "$M"
for P in ACC CESC CHOL DLBC ESCA GBM HNSC KICH LAML LGG LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THYM UCS UVM; do
  S=${P: -2}
  mkdir -p "$R/$P"
  for f in "gene_sums/$S/$P/tcga.gene_sums.$P.G029.gz" "metadata/$S/$P/tcga.tcga.$P.MD.gz" "metadata/$S/$P/tcga.recount_qc.$P.MD.gz" \
           "metadata/$S/$P/tcga.recount_project.$P.MD.gz"; do
    out="$R/$P/$(basename "$f")"
    [ -s "$out" ] || curl -fsSL --retry 3 -o "$out" "$BASE/$f"
    grep -q "^$P/$(basename "$f")	" "$M" || printf "%s\t%s\t%s\t%s\n" "$P/$(basename "$f")" "$BASE/$f" "$(stat -c %s "$out")" "$(sha256sum "$out" | cut -d' ' -f1)" >> "$M"
  done
  echo "$(date '+%F %T') $P done"
done
