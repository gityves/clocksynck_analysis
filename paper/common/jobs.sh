#!/usr/bin/env bash
# Status of the jobs launched through paper/common/run.sh.
#   paper/common/jobs.sh [n_finished=15] [silent_min=60]
# Running: one line per live job (elapsed, current and peak memory, minutes since its log was last written; SILENT
# when longer than silent_min). LOST: a state file whose runner is gone without a ledger line (killed from outside,
# machine reboot). Finished: the last n lines of logs/jobs.tsv, any status other than ok marked with "!!".
set -u
repo=$(cd "$(dirname "$0")/../.." && pwd); logs="$repo/paper/results/logs"
n=${1:-15}; silent=${2:-60}; now=$(date +%s)
echo "== running"
for f in "$logs"/jobs/*.state; do
  [ -e "$f" ] || continue
  rp=$(basename "$f" .state); step=$(sed -n 's/^step=//p' "$f"); log=$(sed -n 's/^log=//p' "$f")
  if ! kill -0 "$rp" 2>/dev/null; then echo "!! LOST runner $rp  $step  (log $log)"; continue; fi
  age=$(( (now - $(stat -c %Y "$log" 2>/dev/null || echo "$now")) / 60 ))
  flag=""; [ "$age" -gt "$silent" ] && flag="  SILENT ${age} min"
  printf '%s  pid %s  since %s  mem %s/%s GB (cap %s)  log %s min ago%s\n' "$step" "$rp" "$(sed -n 's/^start=//p' "$f")" \
    "$(sed -n 's/^pss_gb=//p' "$f")" "$(sed -n 's/^peak_gb=//p' "$f")" "$(sed -n 's/^cap_gb=//p' "$f")" "$age" "$flag"
done
echo "== finished (last $n)"
if [ -f "$logs/jobs.tsv" ]; then tail -n "$n" "$logs/jobs.tsv" | grep -v '^start' | \
  awk -F'\t' '{m = ($3 == "ok") ? "  " : "!!"; printf "%s %s  %-13s exit %-4s peak %s GB  %ss  %s %s\n", m, $2, $3, $4, $5, $6, $8, $9}'
else echo "(no runs recorded yet)"; fi
