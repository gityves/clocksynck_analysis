#!/usr/bin/env bash
# Reruns of 2026-11-01 (HANDOFF_2026-10-31: "Reruns of possibly stale results" and T4), on clocksync c266e20, where
# anchor_readout="loo" became the default (Q13). Two lanes of 4 workers, each step through common/run.sh.
#   lane a: T4, true-time groups whose anchor amplitudes/offsets are quoted (baboon A36/A40/A41, skin A37, V25.5 + A25)
#   lane b: results not rerun after a default-changing change (STATE 2026-10-31, section 2.2)
# The previous CSV outputs are copied first to $SCCC_DATA/m3_snapshot_backup/2026-11-01_pre_loo/ (old vs new compare).
#   setsid nohup paper/exploratory/WP0_integrity/rerun_2026-11-01/queue.sh a > /dev/null 2>&1 &
set -u
repo=$(cd "$(dirname "$0")/../../../.." && pwd); cd "$repo"
backup="${SCCC_DATA:-/mnt/max8tb/data/CoPhaser_data/cancer_circadian}/m3_snapshot_backup/2026-11-01_pre_loo"
log="paper/results/logs/rerun_2026-11-01_$1.log"
run() { echo "$(date '+%F %T') start $*" >> "$log"; RUN_CAP_GB=${RUN_CAP_GB:-20} paper/common/run.sh "$@";
        echo "$(date '+%F %T') end $1 exit $?" >> "$log"; }
backup() { for d in "$@"; do (cd "$repo/paper" && find "$d" -name '*.csv' -print0 | cpio -0pdm --quiet "$backup"); done; }
R=paper/results; E=paper/exploratory
case "$1" in
a)
  backup results/figure_1/output/species results/figure_1/output/method results/figure_1/output/windows \
         exploratory/A25_prior_baseline/output
  run $R/figure_1/scripts/species/compute_baboon.py --jobs 4
  run $R/figure_1/scripts/species/compute_baboon_stability.py
  run $R/figure_1/scripts/method/compute_baboon_truth.py
  run $R/figure_1/scripts/method/compute_skin_truth.py
  run $R/figure_1/scripts/windows/compute_baboon_windows.py --jobs 4
  run $R/figure_1/scripts/windows/compute_breast_timed.py
  run $E/A25_prior_baseline/scripts/compute_prior_baseline.py ;;
b)
  backup results/figure_3/output results/figure_4/output results/figure_2/output/backbone \
         exploratory/A150_atlas_specificity/output exploratory/A144_strong_clocks/output exploratory/WP0_integrity/output
  run $R/figure_3/scripts/gtex_cues/compute_ventilation_dose.py
  run $R/figure_3/scripts/gtex_cues/compute_metabolic.py
  run $R/figure_4/scripts/natural_experiments/compute_covid_ventilation.py
  run $R/figure_4/scripts/resilience/compute_contrasts_a17.py --jobs 4
  run $R/figure_4/scripts/spread_multiplicity/run.sh --jobs 4
  run $R/figure_2/scripts/backbone/compute_oscillator_gtex.py --jobs 4
  run $R/figure_2/scripts/backbone/compute_cancer_split.py --jobs 4
  run $R/figure_2/scripts/backbone/compute_positive_control.py
  run $R/figure_2/scripts/backbone/compute_followups_a18.py --jobs 2
  run $E/A150_atlas_specificity/scripts/compute_a150.py --jobs 4
  run $E/A144_strong_clocks/scripts/compute_split_half.py --jobs 4
  run $E/WP0_integrity/scripts/a21_shift_without_lymphatic.py --jobs 4 ;;
esac
echo "$(date '+%F %T') lane $1 done" >> "$log"
