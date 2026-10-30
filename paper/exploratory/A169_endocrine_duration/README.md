# A169: progressive loss of tumour clock coherence with the duration of letrozole

Registration: `paper/docs/registered_experiments.md`, A169. Library and definitions: `A168_endocrine/scripts/a168_lib.py`.

| file | role |
|---|---|
| `scripts/prep_a169.py` | GSE111563, GSE59515, GSE55374 series matrices → `output/<series>_{meta,expr}.parquet` (A168 parser; time 0 pre, 1 early, 2 late) |
| `scripts/run_a169.py` | `--stage power` (null SDs, MDE; no observed statistic), then `--stage tests` (H169.1-2, pseudo-clock companion, progression, sensitivities) |

Outputs (git-ignored, `output/`): `power.csv`, `null_draws.csv`, `tests.csv`, `coherence_tp.csv`, `progression.csv`,
`treatment_effect_scores.csv`, `sensitivity.csv`, `prep_log.csv`, `sample_counts.csv`.

## Choices not in the registration
1. GSE59515 is deposited linear (log2-transformed); GSE111563 and GSE55374 are log2. GSE111563 keeps 12 clock genes.
2. The progression slope is the least-squares slope of coherence on the time rank (0, 1, 2); with equal spacing it
   equals half of late − pre on the complete triplets, so the 2-week point is shown in `coherence_tp.csv`.
3. Pseudo sets are built on all samples of the series; seeds crc32("A169").
4. The "NR1D1" sensitivity is identical to "without PER1" in GSE111563 (no NR1D1).
