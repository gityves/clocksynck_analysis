# A157: JS correlation on an evenly spread 4-h grid

Registered in `paper/docs/registered_experiments.md` (A157). `scripts/compute_a157.py` (4 s, single thread). Outputs: `output/a157_calibration.csv`
(+ `a157_calibration_{muscle,adipose}.csv`, bisection trace), `a157_replicates.csv`, `a157_summary.csv`.

Choices beyond the registration:
- True internal time = recorded time + subject offset N(0,1 h) + sample noise N(0,1 h) (noise SD not registered; = `clocksync.simulate` jitter_h).
- Generator: local (the times are given by the design, which `cs.simulate` cannot take): clock-plane state A exp(i W t), gene log expression =
  X beta + N(0, tau^2) (+ mesor, per-gene offset sd 0.5, common shift sd 0.2, Poisson at each sample's observed depth for adipose; muscle
  is an array: log2 values, unweighted projector). A = observed clock-plane amplitude (complex regression of the observed sample
  states on exp(iW t_recorded): muscle 1.56, adipose 2.33). "Observed coupling" is thus carried by A and tau.
- Muscle waveform: published GTEx waveform with the ARNTL row dropped (the A113 sensitivity arm), not the 13-gene refit (needs controlled-access tod).
- Calibration: tau by bisection on a pilot of 100 replicates (seed 1157), criterion median |error| of the estimated phase vs the true internal
  time after the common offset (A113 `score` convention), tolerance 0.01 h; no correlation was evaluated in the pilot.
- Correlations are of estimated phase vs TRUE INTERNAL time (registered); vs recorded time also in the replicates (`*_recorded`).
- The adipose design is run through the same estimator and criterion (error 1.11 h).
