# A211: liver regeneration (GSE180012) and the clock

Registration: `paper/docs/registered_experiments.md`, section A211. Data: `$SCCC_DATA/gse180012_phx/` (PROVENANCE.md).

- `scripts/parse_gse180012.py <data dir>` (run with `python -I`): sample table and clock-gene / full expression tables from the GEO matrices and the processed xls files.
- `scripts/compute_a211.py --stage gate|power|tests`: manipulation gate, null / MDE of H211.1 (no observed coherence), then the tests
  (observed T, pseudo-clock companion, 4-16 h / 20-32 h splits, cosinor, day/night series). Reuses `a168_lib` (coherence, pseudo sets, Hallmark) and `clocksync.paired.reference_gate`.
- Data notes: file 1 holds log2(RPKM + 1)-like values (14,620 genes, 8 time points x 3 replicates per arm); Excel turned some gene symbols into dates (dropped);
  mouse symbols are upper-cased to match the human waveform, Hallmark and pool definitions. Surgery 08:00-09:00 (GEO extract protocol).
- Output: `output/` (samples, gate, power, null_draws; tests files after the tests stage).
