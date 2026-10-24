# A123: GTEx organ map from working-hours deaths only

Registered in `paper/docs/registered_experiments.md` (A113/A123/A130/A131 section and the WP4 note). Script:
`scripts/compute_working_hours.py` (`run.sh [--subsample N]`, 30 GB cap, log `paper/results/logs/compute_working_hours.log`).

Outputs (git-ignored, `output/`; per-tissue aggregates and summary statistics only, no per-donor time, phase or coordinate):
`organ_map.csv` (per tissue: estimate with 95 % donor-bootstrap interval, truth (i), truth (ii), CHIRAL), `tests.csv`
(H123.1 / H123.2 per method and truth), `donor_time.csv` (H123.3 and the CHIRAL donor ordering), `truth_i.csv`, `truth_ii.csv`.

Choices beyond the registration:
- Truth (ii) is re-derived from CHIRAL's per-sample CT (`sample_metadata.txt`): the same cosinor + rank-1 complex SVD as truth (i),
  on all samples with a CT in the labelled GTEx table (all Hardy classes, one per donor and tissue); whole blood has no CT.
- The waveform of the untimed fit is fitted on the fast deaths outside 08-18 h (the benchmark's cross-fit); a run with the
  all-donor `cs.gtex()` waveform is reported as a sensitivity (`clocksync_gtex_waveform`).
- All tissues are anchors (free host), as in the benchmark and the A25 GTEx windows; bootstrap 200 donor resamples.
- Offsets and log2 relative amplitudes are centred over the tissues common to estimate and truth (reference free); H123.2 uses the
  log2 scale (the linear-scale Pearson is also in `tests.csv`).
- Tissue-permutation p: 10 000 permutations of the truth across tissues, one-sided, (1 + #>=)/(1 + n).
- H123.3: tangent coordinate of `host_estimate` with the hospital-prior flip, Spearman against recorded time among window donors;
  donor-bootstrap interval resamples the pairs (host not refitted).
- CHIRAL: pychiral per tissue (seed 0) on log2-normalised counts of the window donors, per-donor phases merged by the benchmark
  consensus, tissue estimands by cosinor on the consensus times + SVD.
