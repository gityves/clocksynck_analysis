# A170: is the ER− breast-tumour clock-coherence loss carried by the malignant cells?

Registered in `paper/docs/registered_experiments.md` (A170). Data: Wu 2021 (GSE176078) and Pal 2021 (GSE161529, figshare
17058077), `$SCCC_DATA/sc_breast_er/`; scout in `exploratory/scouts_2026-10-31/sc_breast_er/`.

```
scripts/prep_a170.py                 per-patient pseudobulks (malignant, whole tumour; Wu CAF, endothelial) -> output/units.csv, counts_*.parquet
scripts/run_a170.py --stage power    permutation null and MDE (power.csv, null_draws.csv)
scripts/run_a170.py --stage tests    H170.1, pseudo-clock companion (tests.csv, pseudo_sets.csv), descriptive.csv
scripts/rds_parse.py                 pure-python .rds reader (scout copy), meta.data only
```
Coherence is the A97 / A127 statistic (`a127_lib.coh` through `A168_endocrine/scripts/a168_lib.py`) on log fractions
centred within cohort. Run through `paper/common/run.sh`; prep about 2 min, each stage under a few minutes.
