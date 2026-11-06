# A225 copy number at the clock loci: readout artefact or clock biology? (code, simulation gate and power stage)

Spec: `paper/docs/registered_experiments.md`, A225. `scripts/compute_a225.py --stage sim|power|tests [--jobs 8] [--nrep 200] [--nperm 2000] [--npseudo 200]`;
`--smoke --out DIR` runs every stage on noise expression (real labels, cells, copy numbers). `scripts/a225_lib.py`: coherence from a pairwise-complete covariance (batched over
arms / permutations), the simulator and the replicate worker. A224 objects (cells, GISTIC, dosage fits, adjusted expression) are obtained by executing `compute_a224.py` up to its stage
branch in a scratch directory (as `posthoc_a224.py`).
Run 2026-11-06: sim (200 replicates, 4 min, peak 6 GB), power (10 s, peak 3 GB), both through `run.sh`. The tests stage is implemented and smoke-tested only (incl. the masked H225.1 branch on a
hand-edited smoke power file); it has NOT been run on real data. Not computed: observed H225.1 D, observed D_g, pseudo companions.

Outputs: `sim_calibration.csv`, `sim_replicates.csv`, `sim_summary.csv` (gate, V225.1), `power.csv`, `null_draws.csv`, `eligibility.csv`.

## Coherence from covariance (a225_lib)
The A97 coherence (`clocksync.decompose.components`) is a function of the gene covariance only (state regression coefficients (M'S M)^-1 M'S_g). `arm_cov` builds covariances for many
row subsets at once with missing values (each pair uses the rows where both genes are observed, means over those rows); `coherence_from_cov` evaluates the coherence. Without missing
values it equals `a168_lib.coh` (max abs difference 3e-16 on the real A219 arms; also on a 13-gene sub-waveform). Kept in the A225 folder (no clocksync change).

## Choices not fixed by the registration
- Simulation scale: natural-log expression (A219: ln((count + 0.5)/depth) TCGA, ln2 x log2 intensity METABRIC), mean 0 (coherence is mean-invariant), library-depth noise ignored. Simulated tumours:
  all rows of the cohort table (METABRIC 1,980; TCGA BRCA 1,093), real GISTIC (NaN = no call, treated as 0), real subtype labels (unlabelled / Normal-like tumours get the "others" lam).
  Clock term lam_i amp_g cos(phi_i - psi_g) with the GTEx waveform of `L.WF`; phi uniform. Dosage: `beta_inter` and the CN x (cos, sin) coefficients of the cohort-specific row of A224
  `dosage.csv` (METABRIC; TCGA_BRCA = BRCA alone, not the pooled 23 projects), converted to ln; genes without a fit (METABRIC CIART) have no dosage term.
- Residual SD per gene: refit of A224's model (y ~ cohort + beta CN + cos + sin of the leave-gene-out phase) with our own OLS (`fit_cn`; betas asserted equal to `dosage.csv`), per cohort; the SD
  is of that residual (it contains un-modelled real clock variation). Gaussian noise independent across genes.
- Calibration of lam per cohort on LumA tumours of the cells only: bisection (log scale, 20 common-random-number simulations per evaluation, tolerance 0.003) so that the weighted mean over the 3 cells of the
  raw LumA coherence equals the A219 `descriptive_per_group.csv` value (METABRIC 0.6474, BRCA 0.4885). Achieved values in `sim_calibration.csv`. Non-LumA tumours: same lam (N) or 0.7 x lam (D).
- **Planted true D (not fixed by the registration).** The population coherence of a rank-2 clock is ~1 at any lam (the phase is consistent), so a population ("noise-free") covariance gives true D ~ 0.
  True D is therefore the mean, over 200 copy-number-free simulations (clock + Gaussian noise, no dosage terms, same cells and arm sizes, raw estimator), of the per-cell coh(LumA) - coh(others); under N it is 0.
  Pooled truth = mean over the cells the estimator can use. Gate "recovery" = mean D / true D within 0.75-1.25.
- Estimators per replicate: raw; adjusted = Y - beta CN with betas re-fitted in the replicate (our `fit_cn`, one cohort, leave-gene-out phase of the simulated data; NOT the pooled 23-project beta of A224);
  masked = GISTIC != 0 set missing, pairwise-complete covariance (no CN call = observed).
- Gate rule on the pooled level (mean over the usable cells of both cohorts): N bias |mean D| <= 0.03, rejection <= 7 %, D recovery within +-25 %. The permutation test is A219's (labels permuted within cell,
  2,000, one-sided greater, `perm_p`), recomputed in every replicate for every estimator from the same label permutations. Same clock/noise draws for N and D within a replicate.
- Pseudo genes in the simulation (descriptive, not part of the gate): one A219 pseudo set per replicate (replicate r uses set r mod 200; stored A166 sets for METABRIC, A219-seed sets for BRCA), no clock term, own CN,
  own beta and residual SD (cohort fit, phase of the 14 real clock genes). `pseudo_mean_D` / `pseudo_sd_D` in `sim_summary.csv`.
- **Masked estimator estimability (forced by the data).** A masked arm needs every gene pair observed together in >= 5 tumours of each arm (`MIN_PAIR`, `min_pair_n`). TCGA BRCA fails in all 3 A219 cells (the
  "others" arm of cell 1 has 15 tumours and 1 complete pair; LumA arm of cell 3 has 42 and 1), so the masked estimator uses METABRIC only (3 cells); its truth and pooled level are on those cells.
  In permutations a cell whose permuted arm has a pair with < 3 complete rows is left out of that draw's mean (nanmean).
- Power stage, H225.1: raw, adjusted and masked nulls are all computed (seed crc32("A225"), tag 1; labels permuted, no observed value); the gate-passing estimator (masked if it passes, else adjusted) is `primary`.
  The adjusted null uses A224's adjusted expression (pooled betas), so H225.1 with `adjusted` is identical to A224 H224.1; the tests stage then does not recompute it (prints a pointer to A224 `tests.csv`).
- Trans test cells (`trans_cells`): tumours with GISTIC(g) <= -1 (label 1) or == 0 (label 0) and FGA (A213 `fga_tcga.csv`) and proliferation; FGA and proliferation tertiles are taken within project over those tumours;
  cells = project x FGA tertile x proliferation tertile; eligible if both arms >= 10 tumours and the masked pair rule above holds. Outcome: coherence of the other 13 genes (13-gene sub-waveform), each masked where its own
  GISTIC != 0 (no call = observed). All 23 projects of A215 (A224 `PROJECTS`).
- H225.2 scored iff MDE <= 0.15 (not computable if there is no eligible cell); H225.3 `scored` = >= 5 eligible cells (registered); an MDE <= 0.15 flag (`mde_le_0.15`) is given for information only (no registered
  MDE rule for H225.3). Tests stage runs only the `scored` genes: pseudo sets (13 genes, 200, A215-style construction with its own seeds SeedSequence([crc32("A225"), 3, gene, project])), masked by their own CN, same arms.
  Holm over H225.1-2 and separately within H225.3; verdict pass = Holm p < 0.05 and beyond the pseudo q95 (H225.1) / q05 (H225.2-3).
- Expected tests runtime: dominated by loading the count tables of the projects with eligible cells and generating pseudo sets (minutes; A224 tests took 4 min); the number of eligible cells is small (see `power.csv`).
