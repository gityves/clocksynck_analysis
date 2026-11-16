#!/usr/bin/env bash
# Run a paper step with a RAM cap (the machine is shared): a script (`python`, env `sccc`), a notebook (executed in
# place with papermill) or a group's `run.sh` (log `logs/m3_<group>.log`).
#   [RUN_CAP_GB=30] paper/common/run.sh paper/results/figure_1/scripts/method/compute_gtex_benchmark.py [script args...]
# Log: paper/results/logs/<name>.log (or <RUN_LOG_NAME>.log, for concurrent runs of the same script). Memory (PSS, so pages shared between joblib workers count once) is summed
# over the whole process tree (papermill starts the kernel in its own session, and joblib workers under it); on
# breach every process of the tree is killed.
# Every run is recorded: a live state file logs/jobs/<pid>.state while it runs, and one line in logs/jobs.tsv when it
# ends (status ok | failed | killed_memory | interrupted, exit code, peak GB, runtime). `paper/common/jobs.sh` lists
# running jobs (with log silence) and the last finished ones. Exit codes: the step's own; 137 if killed for memory;
# 130 if interrupted.
set -u
step=$1; shift; cap_gb=${RUN_CAP_GB:-30}; args="$*"
repo=$(cd "$(dirname "$0")/../.." && pwd)
name=$(basename "$step"); log="$repo/paper/results/logs/${name%.*}.log"
[[ $step == *.sh ]] && log="$repo/paper/results/logs/m3_$(basename "$(dirname "$step")").log"
[ -n "${RUN_LOG_NAME:-}" ] && log="$repo/paper/results/logs/${RUN_LOG_NAME}.log"
jobs_dir="$repo/paper/results/logs/jobs"; ledger="$repo/paper/results/logs/jobs.tsv"
mkdir -p "$jobs_dir"
[ -f "$ledger" ] || printf 'start\tend\tstatus\texit\tpeak_gb\truntime_s\trunner_pid\tstep\targs\tlog\n' > "$ledger"
cd "$repo"

tree() {  # pid and all its descendants
  ps -eo pid=,ppid= | awk -v root="$1" '{kids[$2] = kids[$2] " " $1} END {
    n = split(root, q, " "); for (i = 1; i <= n; i++) print q[i];
    while (n > 0) { p = q[n--]; m = split(kids[p], c, " "); for (j = 1; j <= m; j++) { print c[j]; q[++n] = c[j] } } }'
}

t0=$(date +%s); start=$(date '+%F %T'); peak_kb=0; state="$jobs_dir/$$.state"
finish() {  # status exit_code: append to the ledger, drop the state file, exit
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$start" "$(date '+%F %T')" "$1" "$2" \
    "$(LC_ALL=C awk -v k="$peak_kb" 'BEGIN {printf "%.1f", k / 1048576}')" "$(( $(date +%s) - t0 ))" "$$" "$step" "$args" "$log" >> "$ledger"
  echo "run: $1 (exit $2, peak $((peak_kb / 1048576)) GB)" >> "$log"
  rm -f "$state"; exit "$2"
}
on_signal() { kill $(tree "$pid" | sort -u) 2>/dev/null; finish interrupted 130; }

if [[ $step == *.sh ]]; then  # a group's run.sh (it picks the sccc python itself)
  setsid bash "$step" "$@" > "$log" 2>&1 < /dev/null &
else
  if [[ $step == *.ipynb ]]; then cmd=(papermill -k sccc --cwd "$(dirname "$step")" "$step" "$step"); else cmd=(python -u "$step" "$@"); fi
  setsid conda run --no-capture-output -n sccc "${cmd[@]}" > "$log" 2>&1 < /dev/null &
fi
pid=$!
trap on_signal INT TERM HUP
while kill -0 "$pid" 2>/dev/null; do
  pids=$(tree "$pid" | sort -u | tr '\n' ' ')
  pss_kb=$(for p in $pids; do cat "/proc/$p/smaps_rollup" 2>/dev/null; done | awk '/^Pss:/ {s += $2} END {print s + 0}')
  [ "$pss_kb" -gt "$peak_kb" ] && peak_kb=$pss_kb
  printf 'start=%s\nstep=%s\nlog=%s\nchild=%s\npss_gb=%s\npeak_gb=%s\ncap_gb=%s\n' "$start" "$step" "$log" "$pid" \
    "$((pss_kb / 1048576))" "$((peak_kb / 1048576))" "$cap_gb" > "$state"
  if [ "$pss_kb" -gt $((cap_gb * 1024 * 1024)) ]; then
    echo "run: PSS $((pss_kb / 1048576)) GB > $cap_gb GB, killing" >> "$log"
    kill $pids 2>/dev/null; sleep 2; kill -9 $pids 2>/dev/null
    finish killed_memory 137
  fi
  sleep 5
done
wait "$pid"; status=$?
[ "$status" -eq 0 ] && finish ok 0 || finish failed "$status"
