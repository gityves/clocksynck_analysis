# Roadmap 4.4 (descriptive): prior-only baseline for V25.5 (breast per-patient internal time)

Script `scripts/compute_prior_baseline.py` (seconds; `python` of env `sccc`). Input: `paper/results/figure_1/output/windows/v25_5_patients.csv`
(42 patients, primary run). Output `output/prior_baseline.csv` (aggregates only).

Choices beyond the roadmap item:
- Same offset logic as A123c (`compute_truth_iii.py`, prior baseline): method scored with (a) the offset minimising its median error
  (0.05 h grid; the A123c convention, 1.07 h here) and (b) the circular mean of its hours set to 13:00 (no recorded time).
  The published V25.5 1.12 h (`v25_5_summary.csv`, `mae_day`) uses the mean-error offset (-2.2 h) and is not recomputed here.
- Constants: always 13:00 (the hospital-prior window centre used by the estimator); circular mean of the recorded hours and the median-error
  best constant (both use the truth, so are favourable to the baseline).
