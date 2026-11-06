# A219: clock coherence between molecular subtypes at matched proliferation

Registration: `paper/docs/registered_experiments.md`, section A219 (power note after it). Script: `scripts/compute_a219.py`
(`python -I scripts/compute_a219.py --stage power|tests [--jobs 4]`; power takes about 30 s). Mirrors A217 (TCGA A210 cache and A217 subtype table
`A217_subtype_matched/output/subtypes_tcga.csv`; coherence `a168_lib` `L.coh` / `L.WF`).

- `--stage power`: `eligibility.csv`, `counts.csv`, `power.csv`, `null_draws.csv` (null of S / D from permuted labels only; 2,000 permutations, 20 chunks, seeds crc32("A219")).
  No coherence of any real group is evaluated.
- `--stage tests` (not yet run on real data): `tests.csv` (S / D, p one-sided greater, Holm over H219.1-3, pseudo companion), `pseudo_stat.csv`, `descriptive_*.csv`
  (subtype coherences per cell and group, per cohort, without proliferation matching with its own null, proliferation imbalance, ER within LumB / Her2). Tested only with
  `--smoke --out <scratch>` (random noise expression, real labels and cells).

Cohorts: H219.1 the 11 A217 projects (BRCA excluded); H219.2 METABRIC (`L.A166.load("metabric")`, `pam50`, A166 five) and TCGA BRCA (A217 subtypes, `BRCA_` prefix dropped);
H219.3 SCAN-B (`a127_lib.load()`: replicates dropped, PAM50 and ER recorded, prolif = G2M score).

Choices left open by the registration:
- Proliferation tertile (rank "first", A210 / A166 function) within project / cohort over all tumours with a subtype label, before the 15-tumour rule; one cell per tertile.
- A subtype enters a cell with >= 15 tumours; H219.1 cell needs >= 2 entering subtypes (all labels as given, small levels drop out); no minimum on a pooled group beyond it.
- D: LumA, LumB, Her2, Basal each enter by their own >= 15; a cell needs LumA and >= 1 other, the entering others are pooled; Normal-like left out of the contrast (but kept in the tertile).
- S: tumour-weighted SD around the weighted mean (population form, denominator N), unweighted mean over cells. D and S pooled over cohorts as the mean over all cells.
- Permutation of the subtype labels within cell, over the tumours of entering subtypes only.
- Pseudo companion (tests): TCGA 200 sets per project (A166 construction, seeds crc32("A219")); METABRIC the stored A166 sets (`pseudo_sets_metabric.csv`); SCAN-B the stored A127 sets matched on
  mean and variance (`pseudo_sets_variance_matched.csv`). Pseudo S / D is the same statistic on the same cells; clock-specific = observed S (D) above the 95th percentile of the pseudo values.
- "Without proliferation matching": one cell per project / cohort, same entry rule, own 2,000-permutation null. Imbalance: SD of the subtypes' mean proliferation within cell. ER within LumB / Her2:
  METABRIC and SCAN-B only (TCGA BRCA has no ER column in the cache), coherence of ER+ and ER- when >= 15, no test.
- Holm over the scored hypotheses among H219.1-3.
