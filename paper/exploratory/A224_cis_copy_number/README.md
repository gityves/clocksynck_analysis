# A224 cis copy number of the clock genes (code and power stage)

Spec: `paper/docs/registered_experiments.md`, A224. `scripts/compute_a224.py --stage power|tests [--jobs 8] [--nperm 2000] [--npseudo 200]`,
`--smoke --out DIR` runs both stages on noise expression. Power stage run 2026-11-06 (45 s, peak 3 GB, `RUN_CAP_GB=30 paper/common/run.sh`).
The tests stage is implemented and was smoke-tested (noise expression; pseudo-gene path also checked on real counts with 2 sets); it has NOT been run on real data.
Expected tests runtime: loading 23 TCGA count tables (~1 min), pseudo sets, pooled beta of the pseudo genes (some 10^4 genes), 200 sets x 23 projects: of order 20-40 min, so run it through run.sh.

Outputs: `dosage.csv` (per gene x cohort TCGA / METABRIC, plus TCGA_BRCA descriptive), `power.csv`, `null_draws.csv`, `counts.csv`, `eligibility.csv`, `eligibility_h224_3.csv`.
Not computed before scoring: D_adj, T_adj, H224.3 mutant - wild-type difference, attenuation.

## Choices not fixed by the registration
- Seeds: the permutation nulls reuse the original ones (A219 H219.2: crc32("A219"), tag 1; A215 H215.1: crc32("A215"), tag 0), so the adjusted null shares draws with the unadjusted one.
  Pseudo sets reuse the A219 (BRCA, index in TP + BRCA) and A215 (index in the 23 projects) generators; METABRIC uses the stored A166 sets. crc32("A224") is used for the H224.3 label permutations only.
- TCGA beta pooled over the 23 projects of A215 (contains the A219 TCGA projects); cells and tumours of the fit: all tumours with an expression row and a CN call. METABRIC: all 1,980 rows.
- Phase: leave-gene-out projection y_c[others] @ (amp exp(-i psi)) of the project-centred natural-log levels, no standardisation of genes; one pooled cos / sin slope per cohort.
  Pseudo genes use the phase of the 14 real clock genes (nothing left out).
- Fit on natural-log Y (A219 / A215 expression definition); beta, se, phase coefficients in dosage.csv are in log2 units (beta_ln, se_ln also given).
- cn_var_share = var(beta CN, project-centred) / var(Y, project-centred); also cn_partial_r2 (added to project + phase) and cn_alone_r2 (CN alone with project).
- Interaction: CN x cos and CN x sin added to the model (int_cos, int_sin, SE, joint F p, int_amp); descriptive.
- CN names: CIART = C1orf51 in Xena. METABRIC data_cna.txt has no CIART row: CIART not adjusted in METABRIC (keeps its level); duplicated Hugo rows: first kept. Missing CN: no adjustment (CN treated as 0 in the offset).
- H224.3: eligible tumours = TP53 mutation call (A215 MC3 FILTER-passing non-silent) and GISTIC >= 0 at TP53 and PER1 and == 0 at PER2; projects need >= 10 mutant and >= 10 wild type (MIN_ARM as A215).
  Outcome = mean over PER1, PER2 of (adjusted log2 level - pooled phase cos / sin effect of that gene), difference with project fixed effects (within-project OLS slope of the label); labels permuted within project, 2,000.
  Null mean / SD from these draws; two-sided permutation p in the tests stage.
- H224.2 one-sided (T_adj < 0, as A215), H224.1 one-sided (D_adj > 0, as A219); Holm over the scored ones.
- Attenuation (D - D_adj) / D and (T - T_adj) / T: ratio of the mean over cells, leave-one-cell-out jackknife SE (95 % interval = +/- 1.96 SE).
- Descriptives (tests stage): naive mutant - wild-type (phase-adjusted, no dosage adjustment, project FE, all tumours with a mutation call), 2 x 2 means (mutation x TP53 GISTIC < 0).
