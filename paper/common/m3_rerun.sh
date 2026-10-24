#!/usr/bin/env bash
# M3 full reruns (provisory_plan_paper.md, migration step M3): every figure group's run.sh without --subsample, one
# group at a time (each group uses up to 4 workers), RAM-capped by common/run.sh. The current outputs are copied first
# to $SCCC_DATA/m3_snapshot_backup/<date>/ so the archive-snapshot comparison stays possible.
#   nohup paper/common/m3_rerun.sh [groups...] > paper/results/logs/m3_queue.log 2>&1 &
# M3_JOBS=J passes --jobs J to every script (overrides the run.sh defaults); M3_AFTER=PID waits for that process to end
# first (two lanes can then share the machine without overlapping a running group).
# Per-group logs: paper/results/logs/m3_<group>.log; queue log lists start, end and exit status per group (137 =
# killed for memory); every group is also in logs/jobs.tsv (`paper/common/jobs.sh`).
set -u
repo=$(cd "$(dirname "$0")/../.." && pwd)
cd "$repo"
backup="${SCCC_DATA:-/mnt/max8tb/data/CoPhaser_data/cancer_circadian}/m3_snapshot_backup/$(date +%F)"
if [ -d "$backup" ]; then echo "backup exists, kept: $backup"; else  # never overwrite the pre-rerun outputs
  mkdir -p "$backup"
  cp -a paper/results/*/ "$backup/" 2>/dev/null; rm -rf "$backup"/logs "$backup"/*/scripts "$backup"/*/data
  echo "backup: $backup"
fi
groups_arg=("$@")  # optional: rerun only these groups (e.g. figure_5/scripts/psoriasis)

# Order: the psoriasis positive and cancer first (Figure 5), then the method (Figure 1), cues lost (Figure 3),
# backbone and organs (Figure 2), resilience (Figure 4), the rest, Extended Data last.
groups=(
  figure_5/scripts/psoriasis figure_5/scripts/psoriasis_details figure_5/scripts/cancer
  figure_1/scripts/method figure_3/scripts/cues_lost figure_3/scripts/gtex_cues
  figure_2/scripts/backbone figure_2/scripts/organs figure_2/scripts/atlas_fits
  figure_4/scripts/resilience figure_4/scripts/natural_experiments figure_4/scripts/spread_multiplicity
  figure_1/scripts/species figure_1/scripts/windows figure_2/scripts/gene_rhythms figure_2/scripts/embryo
  figure_5/scripts/lung_breaks figure_5/scripts/cancer_pipeline
  extended_data/scripts/calibration extended_data/scripts/gate extended_data/scripts/host_identity
  extended_data/scripts/simulations
)
[ ${#groups_arg[@]} -gt 0 ] && groups=("${groups_arg[@]}")
if [ -n "${M3_AFTER:-}" ]; then while kill -0 "$M3_AFTER" 2>/dev/null; do sleep 60; done; fi
extra=(); [ -n "${M3_JOBS:-}" ] && extra=(--jobs "$M3_JOBS")
for g in "${groups[@]}"; do
  echo "$(date '+%F %T') start $g"
  RUN_CAP_GB=${RUN_CAP_GB:-30} paper/common/run.sh "paper/results/$g/run.sh" "${extra[@]}"; rc=$?
  echo "$(date '+%F %T') end $g exit $rc"
done
echo "$(date '+%F %T') queue done"
paper/common/jobs.sh 40
