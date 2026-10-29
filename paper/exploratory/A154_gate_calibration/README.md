# A154: gate calibration under shared non-clock structure

Registered in `paper/docs/registered_experiments.md` (A154). Code: `scripts/compute_a154.py` (simulation; `--kinds N0,GNULL,COVER,N1,P`),
`scripts/score_a154.py` (V154.1-3 and descriptive tables). Raw per-replicate files in `output/raw/<kind>/`, tables in `output/`.
Depth values: `../A155_anchor_self_inclusion/output/depth_values.csv` (see that README).

## Choices not in the registration
- **Free gate, anchors for fit/contrast.** The gate uses all compartments (free host, as the A140/A145 gates). For the coverage and the global null, the first K/2
  compartments (2 for K = 4) are anchors and the remaining ones the targets; p and intervals are those of the targets (K/2 per replicate, pooled;
  the first target alone is also reported for the null, to remove within-replicate dependence).
- **Clock.** All compartments z = 1, c by psi = V (1 - c) / (2 c) with the window 8-18 h, jitter 1 h. N0 / N1: amplitude 0 (no clock).
- **Factor.** Applied as in A140 Q1 (`thin_clock`, binomial thinning, per-gene rescale to a maximum of 1) to the 14 clock-gene counts. Loadings
  are drawn over a 2,000-gene panel with unit-norm normal loadings on 10 or 30 % of the genes (clock genes are the first 14 and eligible);
  only the 14 clock-gene rows are needed, since the Dataset reads only those columns. Scores N(0, lambda_k) x scale per donor, shared by all
  compartments, lambda = 411.4, 243.6, 195.2 (`A140_dip_cyclops/output/q1/tumour_immune_genome/calibration_lambda.csv`); new loadings per replicate.
  The depth column is not changed by the thinning (as A140).
- **Global null.** One cohort of 2n donors (same compartment offsets in both arms; a first version simulated two independent cohorts, whose compartment offsets differed between arms, and is kept as `output/raw/GNULL_v1_two_cohorts_flawed`), the first n donors are arm 1, c = 0.6, `Dataset.contrast` with n_perm 1,000 (default), `n_boot = 0`
  (the bootstrap does not enter the permutation p; compute saving). p = the table's two-sided `p`.
- **Coverage.** `Dataset.fit(gate=False, n_boot=200)` percentile interval of c against the nominal c; (iv) `Dataset.lower_bound(n_sim=200)` for
  fits at the boundary, bound holds if nominal c >= c_lo. The gate is not run before the fit (conditioning on the gate is not part of V154.3).
- Seeds: SeedSequence([154, design, cell index, replicate, arm]); the COVER arm uses the same data as the P arm without factor.
