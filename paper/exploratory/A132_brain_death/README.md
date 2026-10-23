# A132: brain death or ventilation?

Registered in `paper/docs/registered_experiments.md` (A132). `scripts/compute_brain_death.py` splits Hardy 0 by
`COHORT` (organ donors, brain dead, vs ventilated post-mortem cases), compares their pair-offset relations with fast
deaths (Fisher-Lee rho, folded distance to 0 / 12 h, p_closer), tests the folded-distance difference at equal n
(H132.1), the ventilation-duration dose (H132.2, `DTHVNTD` hours, top vs bottom tertile) and long-term steroids
(D132.3, `MHSTRDLT`). Reuses `A109_coupling_nulls/scripts` (`_helpers`, `compute_coupling_nulls`). Pairs: >= 20 shared
donors in both groups compared. MDE = 2.49 x SD of the donor-label permutation null.

```
run.sh --scout           group sizes, pairs, tertiles, steroids, MDE (no observed contrast); needs no registration data
run.sh --subsample 500   end-to-end check, outputs in output/subsample/
run.sh                   full scoring (10,000 label permutations, 1,000 donor-label permutations)
```
Outputs (`output/`, aggregates only): `a132_groups.csv`, `a132_contrast.csv`, `a132_steroids.csv`, `scout_*.csv`.
Controlled access: phenotypes and time of death stay in memory.
