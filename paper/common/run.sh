#!/usr/bin/env bash
# Run a paper step with a RAM cap (the machine is shared): a script (`python`, env `sccc`) or a notebook (executed in
# place with papermill).
#   [RUN_CAP_GB=30] paper/common/run.sh paper/results/figure_1/scripts/method/compute_gtex_benchmark.py [script args...]
# Log: paper/results/logs/<name>.log. Memory (PSS, so pages shared between joblib workers count once) is summed
# over the whole process tree (papermill starts the kernel in its own session, and joblib workers under it); on
# breach every process of the tree is killed.
set -u
step=$1; shift; cap_gb=${RUN_CAP_GB:-30}
repo=$(cd "$(dirname "$0")/../.." && pwd)
name=$(basename "$step"); log="$repo/paper/results/logs/${name%.*}.log"; mkdir -p "$(dirname "$log")"
cd "$repo"

tree() {  # pid and all its descendants
  ps -eo pid=,ppid= | awk -v root="$1" '{kids[$2] = kids[$2] " " $1} END {
    n = split(root, q, " "); for (i = 1; i <= n; i++) print q[i];
    while (n > 0) { p = q[n--]; m = split(kids[p], c, " "); for (j = 1; j <= m; j++) { print c[j]; q[++n] = c[j] } } }'
}

if [[ $step == *.ipynb ]]; then cmd=(papermill -k sccc --cwd "$(dirname "$step")" "$step" "$step"); else cmd=(python -u "$step" "$@"); fi
setsid conda run --no-capture-output -n sccc "${cmd[@]}" > "$log" 2>&1 < /dev/null &
pid=$!
while kill -0 "$pid" 2>/dev/null; do
  pids=$(tree "$pid" | sort -u | tr '\n' ',' | sed 's/,$//')
  pss_kb=$(for p in ${pids//,/ }; do cat "/proc/$p/smaps_rollup" 2>/dev/null; done | awk '/^Pss:/ {s += $2} END {print s + 0}')
  if [ "$pss_kb" -gt $((cap_gb * 1024 * 1024)) ]; then
    echo "run: PSS $((pss_kb / 1048576)) GB > $cap_gb GB, killing" >> "$log"
    kill $(echo "$pids" | tr ',' ' ') 2>/dev/null; exit 1
  fi
  sleep 5
done
wait "$pid"; status=$?
echo "run: exit $status" >> "$log"; exit $status
