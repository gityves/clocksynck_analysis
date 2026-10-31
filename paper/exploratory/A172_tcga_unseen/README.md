# A172: proliferation clock-coherence loss in the TCGA projects not yet read

Registered in `paper/docs/registered_experiments.md` (A172). Pan-cancer test of the H166.3 / H171 statistic (coherence of
the high- minus the low-proliferation tertile, within purity-proxy strata) in the 22 TCGA projects outside the 11 of A97.

| file | content |
|---|---|
| `scripts/fetch_a172.sh` | recount3 download (gene sums + metadata) into `$SCCC_DATA/tcga/recount3/<P>/`; `MANIFEST_A172.tsv` |
| `scripts/run_a172.py --stage power` | eligibility (≥ 150 primary tumours), per-project permutation nulls, MDE |
| `scripts/run_a172.py --stage tests` | per-project T, H172.1 (mean), H172.2 (Holm), pseudo-clock companion, without PER1/NR1D1 |
| `output/` | `eligibility.csv`, `power.csv`, `null_draws.csv`, `per_project.csv`, `pseudo_T.csv`, `pooled.csv` |

Data are loaded through `paper/common/datasets/tcga_paired/build.py` in the parent process (the builder module, loaded
from file, does not pickle to joblib workers).
