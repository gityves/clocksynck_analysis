# D1b fetal feasibility (Cao 2020)

Feasibility of three fetal clock readouts that do not need fetuses to share organs (after A181). Not registered; no
fetal clock statistic computed. Summary and numbers: `FEASIBILITY.md`.

| script | content | output (git-ignored) |
|---|---|---|
| `scripts/feasibility.py` | adult calibrations at the fetal design (L1 per compartment, L2 within organ), fetal marker detection | `design.csv`, `passrates.csv`, `subsamples.csv`, `markers.csv` |
| `scripts/feasibility.py --markers` | counts expected at fetal depth under adult expression | `markers_expected.csv` |

Run: `RUN_CAP_GB=16 paper/common/run.sh paper/exploratory/D1b_fetal_feasibility/scripts/feasibility.py --jobs 8` (297 s).
