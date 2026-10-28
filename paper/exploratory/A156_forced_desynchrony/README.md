# A156: forced desynchrony in human blood (GSE48113)

Registration and dated data rule: `paper/docs/registered_experiments.md` (A156). Blind inventory: `inventory.md`.

    scripts/parse_inputs.py    # stage 1 (python -I): output/samples.csv, probe_map.csv.gz, probe_missing.csv.gz
    scripts/compute_a156.py    # stage 2: everything else in output/ (about 3 min, one thread)

Outputs: `tests.csv`, `mde.csv`, `primary_summary.csv`, `sensitivities.csv`, `descriptive_deltas.csv`, `coherence_per_day.csv`,
`descriptive_timed.csv`, `cosinor_per_subject_day_anonymised.csv`, `per_subject_delta_anonymised.csv`, `gene_coverage.csv`,
`gene_probe_counts_alias.csv`, `imputation_by_gene_primary.csv`.

## Choices beyond the registration and data rule
- Estimator imported from A115 (`_lib.states_for`, arrays convention, no depth). One projector per set over the analysis arrays
  (primary: 259 arrays of the 19 subjects). Per subject-day offset: `S.ccorr` of the core state with the cue (locked) state
  (centred within subject x day, no Poisson term). Delta_i = wrap(offset D4 - offset D1).
- Null: a day-label swap within a subject negates Delta_i exactly (the per-subject-day offsets do not depend on labels), so the
  null is implemented as sign flips of Delta_i. n = 19: 10,000 random flips, p = (1 + #{null >= obs}) / (1 + 10,000); n = 15
  (7+7): all 2^15 flips, p = fraction. Statistic of H156.2 = |circ mean Delta cue| - |circ mean Delta locked|, one-sided.
- MDE as in A124: circular SD of the signed null mean Delta(core-cue), (1.96 + 0.8416) x SD; the 95th percentile of the null
  |mean Delta| is also reported.
- Probe mapping: non-control probes (CONTROL_TYPE FALSE), `canonical` plus the 8 A124 aliases (`scripts/_common.py`, copied
  because compute_a124.py runs on import); gene = nanmean of probes' log2 values. Imputation of a gene missing in an array: mean
  within the subject-day; where the subject-day has the gene missing in every array (frequency in
  `imputation_by_gene_primary.csv`), fallback = mean over arrays of the same day (not in the data rule; counted in
  `imputed_fallback`).
- Cosinor (descriptive): per subject-day, state c_t ~ A exp(i (w t + theta)), w = 2 pi / 24 h on the recorded clock time;
  theta (h) = arg of sum (c_t - mean) exp(-i w t), projector and states as above; samples without a clock time dropped.
  Phases are in the clock-plane convention (state angle = w t + theta), compared D4 minus D1 within subject at equal recorded
  clock time. Timed core-minus-cue change = (theta_core - theta_cue)(D4) - (theta_core - theta_cue)(D1).
- Circular correlation (Jammalamadaka) between untimed Delta core-cue and the timed change; bootstrap interval over subjects
  (2,000) and a two-sided permutation p (10,000, subject pairing shuffled). Seed 156 (+1 for descriptives).
- Descriptive bootstrap: subjects resampled (2,000), percentile interval of wrapped deviations around the circular mean.
- Coherence per day: mean over subjects of |z| (per subject-day complex correlation) and |mean z|.
- Subject identifiers in per-subject outputs are anonymised ranks (S01...).
- Sensitivities (same code): 15 subjects with 7+7 arrays (exact), all 21 subjects, primary subjects with canonical names only.
