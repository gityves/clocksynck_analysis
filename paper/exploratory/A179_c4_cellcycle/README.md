# A179: is the proliferation coherence loss carried by cell-cycle-associated clock-gene variation? (C4, bulk)

Registered in `paper/docs/registered_experiments.md` (A179). Level L1.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A179_c4_cellcycle/scripts/run_a179.py   # 322 s, 4 workers

Definitions of the five contrasts (H127.1, H171.1, H171.2, H173.1, H173.3) come from A176 via `A178_c7_robustness/scripts/a176_defs.py`; the registered
dR is asserted to 1e-6 first. Within each arm (high / low proliferation tertile) x stratum every clock gene is regressed on [1, G1/S, G2M]
(G1/S: HALLMARK_E2F_TARGETS without clock genes, mean gene z-score as the A97 G2M score; G2M: the registered score; GSE68465 scores z-scored within site) and
coherence is recomputed on the residuals. S = dR_conditioned - dR. Null: 200 random score pairs, each built from A168 expression-matched
pseudo sets of the G1/S and of the G2M gene lists (pool without clock-like, G2M, purity and E2F genes), S_random computed with the same arms and strata.

- `output/tests.csv`: dR, dR_conditioned, S, random-S mean / SD / q95, p (add-one, share of random S >= observed), Holm over five,
  MDE_S, shrinkage_share = S / |dR|, `MDE_below_half_abs_dR`.
- **MDE of S = 2.49 x random-S SD** (2.49 = z_0.95 + z_0.80). This is the operational reading of the registration's "S's MDE".
- `output/descriptive.csv`: per-gene R^2 of the two scores (mean over strata, high and low arm), dlog_plane, dlog_off and their difference before and after conditioning.
- `output/verdict.csv`: H179.1 supported = >= 3 of 5 pass (Holm p < 0.05) and median shrinkage share >= 0.5. "Not explained by linear cell-cycle composition"
  is called when H179.1 is not supported and MDE_S < 0.5 |dR| in at least 3 of 5 contrasts (majority rule chosen here; the registration did not say how to aggregate).
