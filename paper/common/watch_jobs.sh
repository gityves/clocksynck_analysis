#!/usr/bin/env bash
# Wait quietly until a job launched through paper/common/run.sh needs attention, then print why and exit. Meant to run
# in the background so that its exit wakes the agent: exit 1 = a run ended failed / killed_memory / interrupted, a
# runner was LOST, or a log stayed silent longer than the limit; exit 0 = nothing is running any more (all done).
#   paper/common/watch_jobs.sh [silent_min=120] [extra PIDs to wait for (jobs started without run.sh)...]
# WATCH_IGNORE="pid ..." skips those runners (e.g. another session's LOST runner already reported).
# Polls every 60 s; reads only logs/jobs.tsv, logs/jobs/*.state and log modification times.
set -u
repo=$(cd "$(dirname "$0")/../.." && pwd); logs="$repo/paper/results/logs"
silent=${1:-120}; shift || true; extra=("$@")
ledger="$logs/jobs.tsv"; seen=$( [ -f "$ledger" ] && wc -l < "$ledger" || echo 0 )
while true; do
  sleep 60
  if [ -f "$ledger" ] && [ "$(wc -l < "$ledger")" -gt "$seen" ]; then
    new=$(tail -n +"$((seen + 1))" "$ledger"); seen=$(wc -l < "$ledger")
    bad=$(awk -F'\t' '$1 != "start" && $3 != "ok" {print $3 " exit " $4 " peak " $5 " GB: " $8 " " $9 " (log " $10 ")"}' <<< "$new")
    [ -n "$bad" ] && { echo "ATTENTION run ended badly:"; echo "$bad"; exit 1; }
  fi
  running=0; now=$(date +%s)
  for f in "$logs"/jobs/*.state; do
    [ -e "$f" ] || continue
    rp=$(basename "$f" .state); [[ " ${WATCH_IGNORE:-} " == *" $rp "* ]] && continue; step=$(sed -n 's/^step=//p' "$f"); log=$(sed -n 's/^log=//p' "$f")
    if ! kill -0 "$rp" 2>/dev/null; then
      sleep 5; [ -e "$f" ] || continue  # it may just have finished
      echo "ATTENTION runner LOST (no ledger line): $step (log $log)"; exit 1
    fi
    running=$((running + 1))
    age=$(( (now - $(stat -c %Y "$log" 2>/dev/null || echo "$now")) / 60 ))
    [ "$age" -gt "$silent" ] && { echo "ATTENTION log silent for $age min: $step (log $log)"; exit 1; }
  done
  for p in "${extra[@]}"; do
    if kill -0 "$p" 2>/dev/null; then running=$((running + 1)); else
      echo "NOTE extra PID $p ended (started without run.sh: check its log)"; exit 1; fi
  done
  [ "$running" -eq 0 ] && { echo "DONE nothing running"; "$repo/paper/common/jobs.sh" 30; exit 0; }
done
