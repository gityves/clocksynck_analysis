# A190 within-adrenal clock coupling in the fetus (Cao 2020)

Registered in `paper/docs/registered_experiments.md` (A190, 0e88acc); follows A189 (`../A189_fetal_within/`), whose
genome-wide class pseudobulks it reads. Power: `../D1b_fetal_feasibility/scripts/feasibility.py --adrenal`
(`passrates_adrenal.csv`). Outputs in `output/` (git-ignored); results in `results_preliminary.md`.

| script | content | output |
|---|---|---|
| `scripts/score_a190.py` | H190.1 gate over adrenal parenchyma, stromal, endothelial (10,000 permutations; 200 matched / responsive pseudo sets); L1 coherence per class (descriptive) | `h190_1_gate.csv`, `l1_classes.csv` |

Run through `paper/common/run.sh` (168 s, 8 workers).
