#!/usr/bin/env bash
# Hourly (cron): pull finished Jed mapping deliveries to the lab NAS, then copy them into $SCCC_DATA/sccc_mapping.
# Both copies are sha256-verified against the dataset's SHA256SUMS. Uses the reviewed deliver_pull.sh pinned next to
# this file, never a fresh copy from Jed. Needs the Jed key loaded in ~/.ssh/agent_jed.sock (ssh alias jed-batch)
# and the NAS mounted at /mnt/upnae; otherwise it logs why and exits. Log: paper/results/logs/mapping_sync.log.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
H=jed-batch
W=/work/upnae/paychere_sccc
NAS=/mnt/upnae/sccc_mapping
DATA=/mnt/max8tb/data/CoPhaser_data/cancer_circadian
DEST=$DATA/sccc_mapping
CAP_KB=$((500 * 1024 * 1024))
LOG=$HERE/../../results/logs/mapping_sync.log
NAME_RE='^[a-z0-9_]+$'

mkdir -p "$(dirname "$LOG")"
exec >>"$LOG" 2>&1 9>"$HOME/.cache/sccc_mapping_sync.lock"
flock -n 9 || exit 0
log() { echo "$(date '+%F %T') $*"; }
banner='SCITAS EPFL|Authorized use|scitas-doc'

mountpoint -q /mnt/upnae || { log "NAS not mounted at /mnt/upnae, skipped"; exit 0; }
meta=$(ssh "$H" "ls $W/deliverables_meta" 2>/dev/null) || { log "Jed unreachable (key not loaded?), skipped"; exit 0; }

# 1. Jed -> NAS for datasets without a receipt on the NAS
new=()
for f in $meta; do
  d=${f%.tsv}
  [[ $d =~ $NAME_RE ]] || { log "unexpected dataset name '$d', ignored"; continue; }
  ls "$NAS/$d"/DELIVERED_*.ok >/dev/null 2>&1 || new+=("$d")
done
if [ ${#new[@]} -gt 0 ]; then
  log "pulling to NAS: ${new[*]}"
  bash "$HERE/deliver_pull.sh" "$H" "$NAS" "${new[@]}" 2>&1 | grep -v -E "$banner" | grep -E 'delivered|FAILED|REPORT|: FAILED'
fi

# 2. NAS -> $SCCC_DATA for delivered datasets not yet copied here
mkdir -p "$DEST"
for ok in "$NAS"/*/DELIVERED_*.ok; do
  [ -e "$ok" ] || continue
  src=$(dirname "$ok"); d=$(basename "$src")
  [[ $d =~ $NAME_RE ]] || continue
  [ -d "$DEST/$d" ] && continue
  need=$(du -sk "$src" | cut -f1); used=$(du -sk "$DATA" | cut -f1)
  if [ $((used + need)) -gt "$CAP_KB" ]; then log "$d: would exceed the 500 GB cap, not copied"; continue; fi
  rm -rf "$DEST/.incoming/$d"; mkdir -p "$DEST/.incoming"
  rsync -rt "$src" "$DEST/.incoming/" \
    && (cd "$DEST/.incoming/$d" && sha256sum -c --quiet SHA256SUMS) \
    || { log "$d: local copy FAILED verification, left in .incoming"; continue; }
  mv "$DEST/.incoming/$d" "$DEST/$d"
  row=$(awk -F'\t' -v d="$d" '$1 == d {print $2" (pipeline "$6", "$5" CHF)"}' "$NAS/DELIVERIES.tsv" | tail -1)
  echo "| $d | $row | $(date +%F) |" >> "$DEST/PROVENANCE.md"
  log "$d: copied to \$SCCC_DATA/sccc_mapping, $(wc -l < "$DEST/$d/SHA256SUMS") files sha256 OK"
done
cp "$NAS/DELIVERIES.tsv" "$NAS/REPORT.md" "$DEST/" 2>/dev/null
log "checked: $(echo $meta | wc -w) on Jed, ${#new[@]} pulled"
exit 0
