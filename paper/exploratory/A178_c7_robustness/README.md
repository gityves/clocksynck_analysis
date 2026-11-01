# A178: robustness of the coherence losses to single genes and pseudo-pool size (C7)

Registered in `paper/docs/registered_experiments.md` (A178). Level L1.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A178_c7_robustness/scripts/run_a178.py   # 486 s, 4 workers

`scripts/a176_defs.py` reuses (executes, without copying) the definitions of `A176_c3_plane_variance/scripts/run_a176.py`
(registered contrasts, strata, label permutations, loaders). Before anything else the script asserts that the full-gene statistic of every
contrast equals the registered table value (1e-6; H166.2 from `A166_metabric/output/a166_tests.csv`).

- `output/logo.csv`: V178.1, contrast, dropped gene, statistic, one-sided p (add-one, 2,000 registered permutations, waveform subset of the remaining genes).
- `output/pseudo1000.csv`: V178.2 and descriptive, percentile of the observed statistic among 1,000 new pseudo-clock sets (A168 construction, pool
  `A166.pool_genes`, seed crc32("A178"); one set collection per cohort, shared by the contrasts of that cohort), exact (Clopper-Pearson) 95 % interval.
  Percentile = share of sets with statistic <= observed (as in A171/A173).
- `output/verdict.csv`: V178.1 robust = all 14 statistics < 0 and >= 12 with p < 0.05; V178.2 clock-specific = percentile < 5, `at_the_line` = interval contains 5.
  `pseudo_sets_<cohort>.csv` are the sets.

Choices where the registration was silent: V178.1 passes overall only if all four contrasts are robust; V178.2 passes only if both cohorts are
clock-specific; the statistic of each contrast is the registered coherence difference (dR).
