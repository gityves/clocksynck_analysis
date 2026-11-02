# C8 / A191: tumour clock state and outcome

Draft and decisions: `DRAFT_C8.md`. Registered as A191 in `paper/docs/registered_experiments.md` (ff1c1b4). Outputs in
`output/` (git-ignored); results in `results_preliminary.md`.

| script | content | output |
|---|---|---|
| `scripts/score_a191.py` | power step; H191.1 METABRIC per-tumour in-plane log-odds vs RFS, H191.2 GSE14520 tumour-liver phase distance vs recurrence (Cox, adjusted), 200 matched pseudo-clock sets each; secondary OS, all METABRIC, unadjusted, internal time at surgery | `power.csv`, `tests.csv`, `secondary.csv` |

Run through `paper/common/run.sh` (66 s). Uses A166 (METABRIC tables) and A177 (GSE14520 pairs) outputs.
