# Psoriasis clock-delay robustness to gene sets (post hoc, descriptive, 2026-10-28)

Not pre-registered: no verdict, no p-values. Imports `figure_5/scripts/psoriasis/_shared.py`, `_cues.py` and `clocksync.paired`
unchanged. Pooled offset = weighted resultant of per-series complex correlations (lesional vs uninvolved), weights 1/se^2 of the
A54 series, held fixed. 15 gated series (A78 / A81b), 6 held-out series (A54). Seed 2026.

- `scripts/compute_robustness.py` (peak 1 GB, 3 min): reproduction (asserted to 1e-6 against the stored tables), leave-one-gene-out,
  random 80% subsets, pseudo sets, platform coverage. Reads only the needed parquet columns (all-gene tables at once exceed 10 GB).
- `scripts/compute_alias_coverage.py`: effect of mapping platform names with `clocksync.genes.canonical`.
- Outputs in `output/`: `reproduction.csv`, `leave_one_gene_out{,_summary}.csv`, `random_subsets{,_summary}.csv`,
  `pseudo_sets{,_summary}.csv`, `platform_coverage.csv`, `a124_alias_check.csv`, `alias_mapping_*.csv`, `coverage_after_canonical_mapping.csv`.

Notes: responsive-gene sets with the cue waveform were not run (the responsive pool has at most 27 genes present, fewer than the 37
cue genes, so sets would repeat genes, which the projector cannot take); the responsive control uses the 14-gene core waveform.
Pseudo sets are drawn per series from that series' uninvolved-sample means and pooled by set index, as in `clocksync.paired`.
Matched pseudo sets are not screened against the cue genes. Numbers: see the csv files.
