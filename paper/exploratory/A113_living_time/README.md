# A113 (living-biopsy internal time) and A130 (internal minus wall time against age)

Registered in `paper/docs/registered_experiments.md` (A113, A130, WP4 note). Script `scripts/compute_living_time.py`, run with
`run.sh [--subsample N]` (log: `paper/results/logs/compute_living_time.log`; outputs in `output/`, git-ignored).

Choices beyond the note (see the script docstring): muscle projector refit on GTEx with 13 genes (A28 recipe, all fast-death donors,
15 tissues; sensitivity: ARNTL row dropped from the published waveform); permutation = time labels shuffled within subject;
estimators JS and Fisher-Lee (JS is ill-defined for the evenly spread muscle times, Fisher-Lee reported as the valid one, post hoc);
muscle window subsets by the A25 estimator on sample states (single-tissue host, shift 0); V25.5 ages from SI Table S1 added to
`common/datasets/breast_timed_gse233242/resection_times.csv` by `add_age.py`.

Outputs: `a113_summary.csv`, `a113_subjects.csv`, `a113_windows.csv`, `a113_phases.csv`, `waveform_refit.csv`, `a130_subjects.csv`,
`a130_summary.csv`.
