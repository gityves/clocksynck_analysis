# A124: simulated night shift in human muscle (GSE106800)

Registration: `paper/docs/registered_experiments.md` (A124 and the dated data rule). Blind inventory: `inventory.md`.

    scripts/parse_inputs.py                       # stage 1, python -I: samples.csv, probe_map.csv.gz (no expression)
    scripts/compute_a124.py                       # stage 2: output/
    scripts/compute_a124.py --alias-sensitivity   # post hoc, output/alias_sensitivity/

Outputs: `tests.csv` (H124.1-2), `mde.csv`, `null_distribution.csv.gz`, `descriptive_deltas.csv`,
`descriptive_time_state_difference.csv`, `sensitivity_12_subjects.csv`, `gene_coverage.csv`, `gene_mapping.csv`.

## Choices beyond the registration
- Primary set: the 11 subjects with 4 arrays (44 arrays); one projector per gene set over these arrays (log2 RMA, arrays
  convention: no depth, equal gene weights, `A115 _lib.states_for`). States centred within subject over the 4 arrays, then
  `S.ccorr` per condition (which recentres over that condition's arrays, no Poisson term).
- Cue and locked sets from `a81b_sets.csv` (amp, acro_h); core = `cs.gtex()`.
- Gene names: `clocksync.genes.canonical` only (the project table holds 4 aliases). Coverage obtained: core 14/14,
  cue 35/40, locked 36/40; the data-rule expectation (40/40, 39/40) is not met because the table lacks the platform
  names FAM46B, MTERFD3, HIST2H2BE, C12orf23, SEPT4, C10orf116, CENPC1, FAM134B. These runs are the primary result.
- Post hoc sensitivity (not registered): the 8 aliases above added locally (`EXTRA` in the script) gives core 14, cue 40,
  locked 39. Reported separately; it does not replace the primary.
- Probe rule as in the data rule (main category; probe sets listing the symbol alone if any; mean of log2). NR1D1 uses its
  two THRA;NR1D1 read-through probe sets.
- Permutation: all 2^11 = 2,048 within-subject swaps of the two condition labels (a swap exchanges the subject's control and
  misaligned arrays); projector and states fixed. p = fraction of assignments with statistic >= observed (includes the
  identity); H124.1 statistic is |Δ| so two-sided; H124.2 one-sided (greater). Holm over the two p values.
- Delta = offset(misaligned) - offset(control), offset = arg of the complex correlation (hours, wrapped to [-12, 12)).
- MDE (not specified numerically in the registration): circular SD of the signed null Δ(core-cue) and
  (1.96 + 0.8416) x SD (two-sided 0.05, 80 % power, Gaussian approximation); the 95th percentile of null |Δ| is reported too.
- Bootstrap: 2,000 resamples of the 11 subjects (all their arrays), projector fixed, seed 124; interval = 2.5 and 97.5
  percentiles of the wrapped deviations around the circular mean.
- Time-of-day state difference: per subject and condition, state(19:00) - state(07:00) on the uncentred states of the same
  projectors; mean over subjects (complex). Reported: amplitude, direction (h, arg x 24 / 2 pi), direction of misaligned
  minus control (12 h would be a full reversal), a bootstrap interval of that angle, the cosine between the two mean
  difference vectors, and amplitude over standard error. Direction is in the projector's clock-plane convention, not clock time.
- Coherence = |ccorr| per condition (in `descriptive_deltas.csv`).
- RAS9 (3 arrays) enters only `sensitivity_12_subjects.csv` (descriptive, projector over all 47 arrays).
