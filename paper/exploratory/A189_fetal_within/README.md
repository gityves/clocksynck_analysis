# A189 fetal clock organisation within compartments (Cao 2020)

Registered in `paper/docs/registered_experiments.md` (A189, ded6066); feasibility in `../D1b_fetal_feasibility/`.
Outputs in `output/` (git-ignored); results in `results_preliminary.md`.

| script | content | output |
|---|---|---|
| `scripts/build_genome.py` | genome-wide organ and class pseudobulks from the GEO loom (temporary 21 GB decompression, deleted); clock-gene columns checked against `fetal_cao` | `$SCCC_DATA/natexp/cao2020/genome/pseudobulk_{organ,class}.h5ad` |
| `scripts/score_a189.py` | H189.1 coherence per compartment (10,000 per-gene permutations, 200 matched / responsive pseudo sets, leave-one-gene-out); H189.2-3 within-organ gates with pseudo-set companions | `h189_1_coherence.csv`, `h189_2_3_gates.csv` |

Run through `paper/common/run.sh` (build 1,305 s, peak 24 GB with 4 workers; scoring 1,128 s).
