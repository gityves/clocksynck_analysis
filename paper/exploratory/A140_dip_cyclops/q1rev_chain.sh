#!/usr/bin/env bash
# A140 Q1 revision, autonomous chain (no agent needed): wait for the running gate, score it, and run the primaries
# (H140.1-3) only if the registered simulator gate passes at factor scale 1 (core clock genes excluded), then score them.
#
#   setsid nohup paper/exploratory/A140_dip_cyclops/q1rev_chain.sh <gate run.sh PID> > /dev/null 2>&1 &
# Log: paper/results/logs/a140_chain.log. Outputs: output/q1rev/{simulator_gate,primary_scores,primary_arm_B_descriptive}.csv
set -uo pipefail
cd "$(dirname "$0")/../../.."
GATE_PID=${1:?gate run.sh PID}
LOG=paper/results/logs/a140_chain.log
PY=$(conda info --base)/envs/sccc/bin/python
export A140_OUT=q1rev A140_CALIB=1 A140_N_MODELS=80
D=paper/exploratory/A140_dip_cyclops
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

say "waiting for gate run $GATE_PID"
while kill -0 "$GATE_PID" 2>/dev/null; do sleep 300; done
status=$(awk -F'\t' -v p="$GATE_PID" '$7 == p {s = $3} END {print s}' paper/results/logs/jobs.tsv)
say "gate run ended: status '$status'"
if [ "$status" != "ok" ]; then say "STOP: gate run not ok; nothing launched"; exit 1; fi

"$PY" "$D/q1_score.py" gate >> "$LOG" 2>&1 || { say "STOP: gate scoring failed"; exit 1; }
pass=$("$PY" -c "
import pandas as pd
t = pd.read_csv('$D/output/q1rev/simulator_gate.csv')
r = t[(t.factor_scale == 1.0) & (~t.core.astype(bool))]
print('yes' if len(r) == 1 and bool(r.gate_passes.iloc[0]) else 'no')")
say "simulator gate at scale 1 passes: $pass"
if [ "$pass" != "yes" ]; then
  say "registered revision: H140.1-3 not scored; report the largest passing scale (simulator_gate.csv). Done."
  exit 0
fi

say "launching primaries (H140.1-3), 50 replicates, 10 workers"
RUN_CAP_GB=30 paper/common/run.sh "$D/q1_run.py" --step primary --jobs 10 --reps 50 &
PRIM_PID=$!
say "primary run.sh PID $PRIM_PID"
wait "$PRIM_PID"
status=$(awk -F'\t' -v p="$PRIM_PID" '$7 == p {s = $3} END {print s}' paper/results/logs/jobs.tsv)
say "primary run ended: status '$status'"
if [ "$status" != "ok" ]; then say "STOP: primary run not ok (memory kill: rerun with --jobs 1 by hand)"; exit 1; fi
"$PY" "$D/q1_score.py" primary >> "$LOG" 2>&1 || { say "STOP: primary scoring failed"; exit 1; }
say "done: output/q1rev/primary_scores.csv"
