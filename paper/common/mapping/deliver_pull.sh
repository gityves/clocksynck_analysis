#!/usr/bin/env bash
# Run on YOUR machine (not on Jed): pull finished datasets from Jed over ssh into the NAS folder you have mounted,
# verify every file's sha256 on the NAS copy, write DELIVERED_<date>.ok, append DELIVERIES.tsv, copy REPORT.md,
# and send a receipt back to Jed (so the working files there can be deleted). No password is stored anywhere.
# Usage: deliver_pull.sh <jed ssh host> <mounted NAS folder, e.g. /Volumes/upnae/sccc_mapping> [dataset ...]
set -euo pipefail
H=$1; NAS=$2; shift 2
W=/work/upnae/paychere_sccc
if command -v sha256sum >/dev/null; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
mkdir -p "$NAS"; D=$(date +%F)
DS=("$@"); [ ${#DS[@]} -gt 0 ] || DS=($(ssh "$H" "ls $W/deliverables_meta | sed 's/.tsv$//'"))
for d in "${DS[@]}"; do
  if ls "$NAS/$d"/DELIVERED_*.ok >/dev/null 2>&1; then echo "$d: already delivered, skipped"; continue; fi
  echo "== $d"
  mkdir -p "$NAS/$d"
  rsync -rt --checksum --progress "$H:$W/deliverables/$d/" "$NAS/$d/"
  chk=$(cd "$NAS/$d" && $SHA -c SHA256SUMS 2>&1) || { echo "$chk"; echo "$d: VERIFY FAILED"; continue; }
  n=$(grep -c . "$NAS/$d/SHA256SUMS"); ok=$(echo "$chk" | grep -c ': OK$' || true)
  [ "$n" = "$ok" ] || { echo "$chk"; echo "$d: VERIFY FAILED $ok/$n"; continue; }
  echo "verified $ok/$n files (sha256) on the NAS copy, $(date)" > "$NAS/$d/DELIVERED_$D.ok"
  meta=$(ssh "$H" "cat $W/deliverables_meta/$d.tsv")
  size=$(du -sk "$NAS/$d" | cut -f1)
  [ -s "$NAS/DELIVERIES.tsv" ] || printf "dataset\tdate\tn_files\ttotal_bytes\tCHF\tpipeline_commit\n" > "$NAS/DELIVERIES.tsv"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$d" "$D" "$(ls "$NAS/$d" | wc -l | tr -d ' ')" "$((size*1024))" "$(echo "$meta" | cut -f2)" "$(echo "$meta" | cut -f3)" >> "$NAS/DELIVERIES.tsv"
  scp -q "$NAS/$d/DELIVERED_$D.ok" "$H:$W/delivery_receipts/$d.ok"
  echo "$d: delivered and verified"
done
scp -q "$H:/home/paychere/mapping/REPORT.md" "$NAS/REPORT.md" && echo "REPORT.md copied"
