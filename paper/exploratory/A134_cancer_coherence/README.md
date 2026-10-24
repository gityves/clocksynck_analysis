# A134 scout (before scoring; nothing observed)

Scripts (`scripts/`): `scout_study.py STUDY` (whole-transcriptome per-patient covariates, one cohort), `scout_null.py`
(feasibility tables, null spread), `scout_misc.py` (loci, gene-set versions, D134.3 counts). All outputs in `output/scout_*.csv`.
The observed H134.1 / H134.2 statistic was not computed (the null permutes the copy-number rows).

- **Copy number.** The stored `cnv/<study>.csv.gz` calls are per-cell malignant / normal / ambiguous labels (clonality of an
  inferCNV-like profile), not locus states, and the h5ad obs hold no author copy-number annotation. A locus state is
  therefore derived here from the same profile (`cancer_common.cnv`): patient mean of the profile at the locus over the
  malignant cells (own gene removed from its 101-gene window; unexpressed gene -> nearest expressed gene on the chromosome,
  <= 392 kb). Patient-level score SD, malignant / normal: 0.033-0.059 / 0.010-0.026.
- Files: `scout_cn_feasibility.csv`, `scout_cn_distribution.csv`, `scout_clock_loci.csv`, `scout_patient_covariates.csv`
  (per patient and group: cells, `cn_<gene>_mean/_median`, cycling fraction, S, G2M, MYC, hypoxia), `scout_geneset_*.csv`,
  `scout_geneset_versions.csv`, `scout_d134_3.csv`, `scout_null_spread.csv`, `scout_patients79.csv`.
