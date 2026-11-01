# A181 fetal depth control and responsive-gene shared state

Registered in `paper/docs/registered_experiments.md` (A181 and the A181 design note). Outputs in `output/` (git-ignored).

| script | content | output |
|---|---|---|
| `scripts/step0_repro.py` | reproduces the A107b gates H107b.1 and H107b.2 | `repro_a107b.csv` |
| `scripts/v181_1_depth_control.py` | V181.1: organ-presence-matched GTEx draws (23 donors, 08-18 h, any Hardy, BA9), binomial thinning to Cao depth, gate, 200 subsamples, unthinned reference | `v181_1_subsamples.csv`, `v181_1_passrates.csv`, `v181_1_depths.csv`, `v181_1.log` |
| `scripts/h181_2_responsive.py` | H181.2: pairing null on 15 responsive genes (identity projector, in-script) | `h181_2_responsive.csv` |

Controlled access: time of death is used in memory to select donors; only aggregates are written.
Thinning: p = min(1, Cao depth / library size); depth scaled by p. Cerebrum has p = 1 (GTEx library < Cao depth).
Run: `python scripts/<script>.py` (V181.1: 449 s with 4 workers; the others about 15 s).
