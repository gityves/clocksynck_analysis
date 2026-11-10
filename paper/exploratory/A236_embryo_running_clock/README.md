# A236 is the time-of-day state of the late mouse embryo a running clock? (Qiu 2024)

Registered in `paper/docs/registered_experiments.md` (A236, d51cc24; power note with three deviations before scoring).
Outputs in `output/` (git-ignored); results in `results_preliminary.md`.

| script | content | output |
|---|---|---|
| `scripts/score_a236.py power` | gene-level plasmode on the 77 complete cell types: size and power of H236.1 (S1) and H236.2 (rho), simulated critical values | `power.csv`, `power_info.csv`, `null_t.npy` |
| `scripts/score_a236.py score` | H236.1 / H236.1b gene content of the time-locked harmonics (S1, plasmode null), per-gene harmonics, leave-one-gene-out; H236.2 polarisation rho (jackknife t, simulated critical value), per trajectory | `h236_1_content.csv`, `h236_1_genes.csv`, `h236_1_arms.csv`, `h236_1b_timed.csv`, `h236_2_rotation.csv` |
| `scripts/build_genome.py` | genome-wide cell type x embryo pseudobulks of the late window, streamed from GEO (parse only); clock-gene columns checked against the stored panel | `$SCCC_DATA/natexp/qiu2024/genome_late/pseudobulk_celltype.h5ad` |
| `scripts/score_a236.py companions` | H236.3: H108.1's R against 200 matched and 200 responsive pseudo-clock sets, raw and residualised on embryonic day | `h236_3_sets.csv`, `h236_3_companions.csv` |

Run with the sccc Python; `power` takes about 160 s (2,000 data sets per scenario, 4 threads).
