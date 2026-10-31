#!/usr/bin/env bash
# Download GEO series matrices and platform annotation (ftp.ncbi.nlm.nih.gov only) and write PROVENANCE.md with sha256. Files are only parsed later (python -I).
set -euo pipefail
D=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_gastric_lung; OLD=/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_luad_lihc
declare -A WHY=([GSE62254]="ACRG gastric, 300 tumours" [GSE15459]="Singapore gastric, 200 tumours" [GSE84437]="Samsung gastric, 483 tumours" [GSE26253]="Yonsei FFPE gastric, 432 tumours"
  [GSE30219]="lung histology (ADC 85, squamous 82)" [GSE37745]="lung histology (adeno 106, squamous 66)" [GSE81089]="lung RNA-seq (Uppsala, 199 NSCLC)")
mkdir -p $D/platforms
for g in "${!WHY[@]}"; do
  mkdir -p $D/$g; f=${g}_series_matrix.txt.gz
  [ -s $D/$g/$f ] || curl -sf --max-time 900 -o $D/$g/$f https://ftp.ncbi.nlm.nih.gov/geo/series/${g:0:${#g}-3}nnn/$g/matrix/$f
  printf '# %s\n\nSource: GEO series matrix (processed), downloaded %s from https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/\nScouting for the gastric / lung-histology tests (paper/exploratory/scouts_2026-10-31/gastric_lung_histology; %s). Parsed only (python -I), never executed.\n\n```\n%s\n```\n' \
    $g $(date +%F) ${g:0:${#g}-3} $g "${WHY[$g]}" "$(cd $D/$g && sha256sum $f)" > $D/$g/PROVENANCE.md
done
# GSE41271: copy from the LUAD/LIHC scout (not re-downloaded)
mkdir -p $D/GSE41271; cp -n $OLD/GSE41271/GSE41271_series_matrix.txt.gz $D/GSE41271/
printf '# GSE41271\n\nSource: GEO series matrix (processed), originally downloaded 2026-10-31 from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41271/ into scout_luad_lihc and COPIED here (not re-downloaded). Lung histology scout (275 NSCLC: 183 ADC, 80 squamous per the previous scout).\nParsed only (python -I), never executed.\n\n```\n%s\n```\n' "$(cd $D/GSE41271 && sha256sum GSE41271_series_matrix.txt.gz)" > $D/GSE41271/PROVENANCE.md
# platforms: copy what exists, download the rest
for p in GPL570.annot.gz GPL96.txt GPL10558.annot.gz GPL6884.annot.gz; do cp -n $OLD/platforms/$p $D/platforms/; done
[ -s $D/platforms/GPL6947.annot.gz ] || curl -sf -o $D/platforms/GPL6947.annot.gz https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz
[ -s $D/platforms/GPL8432_family.soft.gz ] || curl -sf -o $D/platforms/GPL8432_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8432/soft/GPL8432_family.soft.gz
{ printf '# platforms\n\nGEO platform annotation tables from https://ftp.ncbi.nlm.nih.gov/geo/platforms/ (GPL570, GPL6884 copied from scout_luad_lihc/platforms; GPL6947 annot, GPL8432 family soft downloaded %s). GPL16791 (HiSeq 2500, RNA-seq) has no probe table.\n```\n' $(date +%F); (cd $D/platforms && sha256sum GPL*); printf '```\n'; } > $D/platforms/PROVENANCE.md
du -sh $D
