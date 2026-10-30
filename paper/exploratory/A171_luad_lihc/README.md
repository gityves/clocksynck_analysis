# A171: proliferation clock-coherence loss in lung adenocarcinoma and hepatocellular carcinoma

Registered in `paper/docs/registered_experiments.md` (A171). Data: GEO series matrices under `$SCCC_DATA/scout_luad_lihc/`;
scout in `exploratory/scouts_2026-10-31/luad_lihc/`. Reuses the A168 parser and library (`A168_endocrine/scripts/`).

```
scripts/prep_a171.py                 series matrices -> output/{series}_meta.parquet, {series}_expr.parquet, prep_log.csv
scripts/run_a171.py --stage power    permutation nulls, MDE (power.csv, null_draws.csv); TCGA LUAD / LIHC (tcga.csv)
scripts/run_a171.py --stage tests    H171.1-4, Holm, pseudo-clock companion (tests.csv), strata.csv, descriptive.csv
```
Run through `paper/common/run.sh`; each step takes a few minutes.
