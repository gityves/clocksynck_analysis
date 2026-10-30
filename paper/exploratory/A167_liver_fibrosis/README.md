# A167: advanced liver fibrosis and the stellate / endothelial / hepatocyte clocks

Registration: `paper/docs/registered_experiments.md`, A167, amendment 1 and the power step (2026-10-30).

| file | role |
|---|---|
| `scripts/compute_a167.py` | `--stage power`: MDE on d* of the primaries (`power.mde`, A153 conversion); `--stage contrast`: gates and contrasts (D2 primary; D1 and leave-one-study-out secondaries) with the liver-wide companion |
| `scripts/score_a167.py` | H167.1 (Holm over 3), H167.2 (Holm over 9), secondaries |

Outputs (git-ignored, `output/`): `power.csv`, `analyses.csv`, `gates.csv`, `contrast_<analysis>.csv`, `h1671.csv`,
`h1672.csv`, `secondary.csv`, `summary.csv`.

## Choices not in the registration
1. The liver-wide companion uses every compartment with permutation draws (leave-self-out), via
   `A165_disease_readouts/scripts/compute_a165.py::companion_all`.
2. The 95 % interval of d* is d* ± 1.96 × the bootstrap SE (clocksync's `lo`/`hi` are for d).
3. The specificity label qualifies effects only; for other rows the companion p is reported without a label.
4. The two power runs (stellate/endothelial, then hepatocyte after amendment 1) share one run.sh log name; the CSV
   output was merged by compartment and checked.
