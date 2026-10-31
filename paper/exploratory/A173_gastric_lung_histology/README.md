# A173: proliferation clock-coherence loss in gastric cancer and by lung histology

Registered in `paper/docs/registered_experiments.md` (A173). H173.1-3: the H171 statistic in three gastric cohorts
(GSE62254, GSE84437, GSE15459; Holm). H173.4: adenocarcinoma minus squamous difference of the statistic in three lung
cohorts (GSE37745, GSE41271, GSE30219), on purity-proxy halves within histology.

| file | content |
|---|---|
| `scripts/prep_a173.py` | series matrices → `output/{key}_meta.parquet`, `{key}_expr.parquet` (A171 parser; GSE62254 log10 → log2) |
| `scripts/run_a173.py --stage power` | permutation nulls, MDE (`null_draws.csv`, `power.csv`) |
| `scripts/run_a173.py --stage tests` | `tests.csv`, `descriptive.csv`, `pseudo_sets_*.csv`, `pseudo_stats_*.csv` |

Data: `$SCCC_DATA/scout_gastric_lung/` (scout `paper/exploratory/scouts_2026-10-31/gastric_lung_histology/`).
