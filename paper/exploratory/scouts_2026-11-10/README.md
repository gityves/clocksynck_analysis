# Scout 2026-11-10: host-spread bias of group contrasts (post hoc)

Post hoc feasibility look made during the review `paper/docs/REVIEW_2026-11-10.md` (§1, M1). Not a registered analysis;
no result of it is entered in `results_preliminary.md`. The registered version is proposal D1 of the review (PI question Q27).

- `posthoc_spread_bias.py <n_seeds> <out.csv>`: two arms of 80 simulated donors with identical loadings and noise, reference
  window 6-18 h, case window 6-18 / 8-16 / 10-14 h, anchor noise 0.002 / 0.05 / 0.15; `Dataset.contrast` without
  permutations. Truth: d* = 0, log z2 ratio = 0.
- `posthoc_spread_bias.csv`: the 12-seed run quoted in the review (one row per seed, design and compartment).
