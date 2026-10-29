# A160: psoriatic core-cue delay within keratinocytes (GSE173706)

Spec: `paper/docs/registered_experiments.md` (A160). Scripts, run in order (thread env vars = 1):

| script | content | outputs |
|---|---|---|
| `scripts/a160_01_build_pseudobulk.py` | all-gene pseudobulks per library x group (whole, keratinocyte, A63 bins) from the processed count matrices and the archived cell annotations; writes `$SCCC_DATA/psoriasis_sc/a160_pseudobulk.h5ad` | `output/a160_cells_per_library.csv` |
| `scripts/a160_02_blind_power.py` | feasibility, null SD and MDE from PN + healthy libraries only | `a160_feasibility.csv`, `a160_blind_power.csv`, `a160_null_draws.csv`, `a160_plan.json` |
| `scripts/a160_03_readout.py` | exact within-pair label swaps, Holm over qualifying variants, descriptives | `a160_tests.csv`, `a160_per_bin.csv`, `a160_descriptives.csv` |

`scripts/_lib.py`: estimator on `clocksync.paired` (`projector`, `states`, `ccorr`) and `clocksync.stats`. Neither variant
reaches MDE <= 1.7 h, so no confirmatory test is scored; all estimates in `a160_tests.csv` are descriptive.
