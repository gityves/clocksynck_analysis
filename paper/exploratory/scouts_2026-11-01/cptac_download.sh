#!/usr/bin/env bash
# Usage: cptac_download.sh LUAD [LSCC CCRCC ...]   (open-access GDC STAR gene-count files, no token needed)
# Run by the user. Destination $SCCC_DATA/cptac_<type>/ ; resumable (skips existing files), md5 checked.
set -euo pipefail
SCCC_DATA=${SCCC_DATA:-/mnt/max8tb/data/CoPhaser_data/cancer_circadian}
HERE=$(cd "$(dirname "$0")" && pwd)
for T in "$@"; do
  M=$HERE/cptac_manifests/${T}_manifest.tsv; OUT=$SCCC_DATA/cptac_$T; mkdir -p "$OUT"
  cp "$M" "$HERE/cptac_manifests/${T}_samples.tsv" "$OUT"/
  tail -n +2 "$M" | while IFS=$'\t' read -r id fn md5 size state; do
    f=$OUT/$id.$fn
    if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then continue; fi
    curl -sS --fail -o "$f" "https://api.gdc.cancer.gov/data/$id"
    [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ] || { echo "md5 mismatch $id" >&2; exit 1; }
  done
  echo "$T done: $(du -sh "$OUT" | cut -f1)"
done
# Alternative: gdc-client download -m cptac_manifests/<T>_manifest.tsv -d $SCCC_DATA/cptac_<T>/
