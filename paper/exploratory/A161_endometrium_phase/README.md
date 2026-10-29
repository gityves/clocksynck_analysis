# A161: endometrial clock coherence, secretory vs proliferative (GSE51981)

Registered: `paper/docs/registered_experiments.md` (A161). Data: `$SCCC_DATA/endometrium_gse51981/` (PROVENANCE.md). Scripts (python -I, in order):
`s00_convert_xlsx.py` (GSE141549 supplementary xlsx to parquet), `s01_metadata.py` (phase, endometriosis status, id overlap; before any readout),
`s02_map.py` (14-gene probe map, gene-level matrices), `s02b_overlap_expr.py` (all-gene near-duplicate check), `s03_run.py` (H161.1, pseudo set, MDE, sensitivities; about 2 min).
Results: `output/a161_results.csv` (one row per test), logs in `output/logs/`.

Choices and deviations from the registration
- GSE51981 id overlap: titles are the sequential indices 1-148, so the id check against GSE6364 / GSE4888 is uninformative; an all-gene correlation check added (not registered).
- GSE6364 series matrix holds raw intensities: log2 with a floor of 1. GSE141549 series matrices hold no values: supplementary normalised xlsx used (HT-12 = GPL10558, WG-6 = GPL13376).
- Clock probes: any probe whose symbol list contains the gene, highest mean (PER1 and NR1D1 probes are shared with MIR6883 and THRA).
- Pseudo pool: all single-symbol genes except clock-like, clock and the 8 composition markers; 20 nearest in standardised (mean, sd), 200 sets, one draw reused in (i)-(iii).
- Composition scores and the residualisation use all 148 samples; the replication series have no pseudo companion (not registered).
- GSE141549: technical replicates ("CE Replicate", 4 secretory) dropped (24 vs 40); the registered counts (24 vs 43, replicates kept, one stage-less secretory dropped) are reported alongside. z-scored within platform.
- Sensitivity rows also carry a permutation p and pseudo percentile, descriptive only.
