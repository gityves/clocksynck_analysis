#!/usr/bin/env bash
# A243 downloads (2026-11-15). One new directory per source under $SCCC_DATA/natexp/<source>/raw/, each with a
# PROVENANCE.md (URL, date, sha256). Files are only parsed afterwards (python -I), never executed.
set -u
D=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp
GEO=https://ftp.ncbi.nlm.nih.gov/geo
get() {  # source dir, file name, url, [description]
  local d=$D/$1/raw; mkdir -p "$d"; local f=$d/$2
  [ -s "$f" ] || curl -sSL --retry 3 --fail -o "$f" "$3" || { echo "FAILED $3"; return 1; }
  local p=$D/$1/PROVENANCE.md
  [ -f "$p" ] || printf '# %s (A243 download, %s; parsed with python -I only)\n\n| file | bytes | sha256 | URL |\n|---|---|---|---|\n' "$1" "$(date +%F)" > "$p"
  grep -q "| raw/$2 |" "$p" || printf '| raw/%s | %s | %s | %s |\n' "$2" "$(stat -c %s "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$3" >> "$p"
}
for x in "GSM4008686 Fetal-Heart1" "GSM4008687 Fetal-Heart2" "GSM4008693 Fetal-Kidney3" "GSM4008694 Fetal-Kidney4" \
         "GSM4008695 Fetal-Kidney5" "GSM4008696 Fetal-Kidney6" "GSM4008699 Fetal-Lung1" "GSM4008700 Fetal-Lung2"; do
  set -- $x; get han2020_organs "${1}_${2}_dge.txt.gz" $GEO/samples/GSM4008nnn/$1/suppl/${1}_${2}_dge.txt.gz
done
get han2020_organs GSE134355-GPL20795_series_matrix.txt.gz $GEO/series/GSE134nnn/GSE134355/matrix/GSE134355-GPL20795_series_matrix.txt.gz
for f in ENCFF646WRS ENCFF292LKT ENCFF564JYO ENCFF193UXG ENCFF012IBV ENCFF598STD ENCFF589BBG ENCFF181TFL ENCFF675DVW \
         ENCFF678NWU ENCFF780YUS ENCFF846SOK; do
  get encode_fetal_organs $f.tsv https://www.encodeproject.org/files/$f/@@download/$f.tsv
done
ser() {  # source dir, GSE, files...
  local s=$1 g=$2; shift 2; local n=${g:0:${#g}-3}nnn
  get $s ${g}_series_matrix.txt.gz $GEO/series/$n/$g/matrix/${g}_series_matrix.txt.gz
  for f in "$@"; do get $s $f $GEO/series/$n/$g/suppl/$f; done
}
ser cui2019 GSE106118 GSE106118_UMI_count_merge.txt.gz GSE106118_barcode_information.txt.gz
ser knight_schrijver2023 GSE216019 GSE216019_RAW.tar
ser hochane2019 GSE114530 GSE114530_RAW.tar
ser lindstrom_gse124472 GSE124472 GSE124472_RAW.tar
ser lindstrom_gse112570 GSE112570 GSE112570_RAW.tar
ser lindstrom_gse139280 GSE139280 GSE139280_RAW.tar
ser lindstrom_gse102596 GSE102596 GSE102596_RAW.tar
ser sountoulidis2023 GSE215895 GSE215895_RAW.tar
ser renesme2025 GSE280880 GSE280880_UMI_count_matrix.csv.gz GSE280880_meta_data_GEO.csv.gz
ser roost2015_organs GSE66302   # series matrix only (ages of the on-disk read table)
echo DONE
