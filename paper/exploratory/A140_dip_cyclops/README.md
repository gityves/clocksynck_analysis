# A140 (CYCLOPS 2.0 donor internal phase vs clocksync): code and results

## Q2 (i), (iii)

Registered in `paper/docs/registered_experiments.md` (A140 Q2 readouts, dated note 2026-10-26). Part (ii) (CYCLOPS
pipeline) is not run here (waits for Q1).

**Run.** `build_q2.py` (`--jobs 4`, log `paper/results/logs/build_q2.log`, peak 1 GB) streams the five GSE266919 objects
with `rds_parse.py` and sums raw UMIs per sample x compartment (`output/q2/pb_*.csv.gz`, `cells_per_sample.csv`,
`samples.csv`; git-ignored). `q2_readouts.py` (`--jobs 4`, log `paper/results/logs/q2_readouts.log`) runs (a) the orientation
check, (b) part (i) on the primary set (44 Pre-treatment biopsies; 43 enter, one has no compartment with >= 20 cells;
strata breast vs other) and on the all-78 set (77 enter), (c) part (iii) only if gate passes and stability >= 0.8.

**Orientation (a).** Simulated timed cohort, 5 compartments, 60 donors, 24 h window, 5 seeds: the resultant length of
(phase - t) is 0.984-0.989 and of (phase + t) 0.025-0.188; a larger `Dataset.host()` phase_h is later (no sign flip).
<!-- src: paper/exploratory/A140_dip_cyclops/output/q2/orientation.csv -->

**Part (i).**

| set | samples | gate p (shared signal) | gate p (clock shape) | clock share (null) | gate | anchors | stability | free agreement |
|---|---|---|---|---|---|---|---|---|
| primary | 43 | 0.000999 | 0.4585 | 0.0246 (0.0283) | fails | B, CD4_T, CD8_T, NK, myeloid (learned, status ok) | 1.000 (median similarity 0.980) | 1.000 |
| all 78 | 77 | 0.000999 | 0.5405 | 0.0249 (0.0311) | fails | same | 1.000 (0.989) | 1.000 |

<!-- src: paper/exploratory/A140_dip_cyclops/output/q2/part_i_summary.csv -->

The pairing p is at the permutation floor (1,000 permutations) but the clock-shape p is not: the shared covariance is
not clock-shaped (clock share below its null mean). The gate fails in both sets, so `Dataset.fit` couplings are withheld
(uninformative; not "decoupled"). The fixed stromal-vascular rule gives no anchor (no such compartment deposited).
Learned families (tissue scope, one breast study in the learning set): myeloid pooled 0.199 (lower 0.113), lymphocyte
0.123 (lower 0.040); both anchor families.
<!-- src: output/q2/primary_families.csv -->
These differ from the note's values (0.202 / 0.136 and 0.091 / 0.027), which came from `a17_anchor_families.csv` (target
organ_breast, all-tissue scope: the breast study itself was the target there); the anchor outcome is the same.
Stability and free agreement are computed on the anchored host regardless of the gate; the free host equals the anchored
host here because all five compartments are anchors, so the agreement of 1.0 is by construction. Primary pair
coherence excess is 0.46 (CD4_T-CD8_T) down to negative (B with CD4_T, NK, myeloid): `output/q2/primary_pairs.csv`.

**Part (iii).** Not estimable (gate failed, both sets): P140.1-3 rows read "not estimable: gate failed".
<!-- src: output/q2/part_iii_readouts.csv -->
The code path (host refit bootstrap, within-stratum permutations, Holm) was exercised on a check run only (40 bootstraps,
100 permutations on the primary set, forced past the gate; not a result and not kept).

### Choices not in the registration
- Compartments with >= 20 cells enter the clock table; the clock genes absent from an object would be set to 0 (none were).
- Depth = sum of all genes of the object (gene sets differ between objects).
- Stability uses the package default threshold 0.7 (share of bootstrap hosts above it); "identified" = share >= 0.8.
- Learning tables rebuilt as `compute_anchors.learn` (gate-passing healthy fits, tissue not "multi"); target pairs from
  `Dataset.coherence(n_perm=200, n_boot=200)` of each set (the all-78 set uses its own pairs).
- Gate seed 0, 1,000 permutations; fit would use n_boot 200, seed 0 (not run).
- Orientation: 60 donors, uniform 24 h window, no jitter, psi 0.05, z = 0.8-1.2 real amplitudes, anchors = all five; sign
  decided by the larger resultant length.
- Part (iii) specifics not exercised (fixed in code): fraction logit = log((n+0.5)/(N-n+0.5)) with N all deposited immune
  cells of the sample; scores z-scored once over the readout's samples (not per bootstrap replicate); bootstrap hosts
  unwrapped about their own circular mean; percentile CIs; interferon-response gene set (ISG15, IFIT1, IFIT2, IFIT3, IFI6,
  IFI44L, MX1, OAS1, STAT1, IRF7) chosen by me, T-cell scores in CD4_T and CD8_T separately, predicted sign > 0; the
  all-78 permutations run within Group x site strata without patient blocks.

## Q1 (simulation of the design)

Registered in `paper/docs/registered_experiments.md` (A140, "Q1"). Code: `q1_lib.py` (simulator, published pipeline,
clocksync gate), `q1_run.py` (steps `timing`, `gate`, `primary`, `grid`, `diag`), `q1_score.py` (gate, H140.1-3 with Holm,
grid summary). Outputs (git-ignored): `output/q1/`. Runs through `paper/common/run.sh`, 8 workers (the lead's
`compute_r26_cohorts.py` was running), 1 Julia worker per CYCLOPS fit, logs `paper/results/logs/q1_run.log`.

**Outcome: the simulator gate failed, Q1 stopped. Step 3 (H140.1-3, Holm) and the descriptive grid were not run.**
In arm D (clock coupled at 0.8, sampled over 24 h, K = 1, 1,500 seed genes, 50 replicates) the CYCLOPS phases have a
median circular correlation with the true time of 0.067 (maximum over replicates 0.259, 0 of 50 at or above 0.7).
<!-- src: output/q1/simulator_gate.csv; output/q1/gate/D_K1_r*.csv (truth_r_median) -->
The cross-compartment circular correlation of the same phases is |r| >= 0.89 in every replicate (median absolute value 0.986 over
replicates; the sign is arbitrary across replicates): CYCLOPS orders the cohort along the shared non-clock factor, and the
DIPs cover a median 0.29 of the circle (Rayleigh R 0.93). The clocksync gate (free host, package default) passes in 43 of 50
of these arm-D replicates, so the clock is there and observable; the CYCLOPS pipeline does not recover it.
<!-- src: output/q1/gate/D_K1_r*.csv (gate_pass, cc_median, dip_cover, dip_R) -->

**Cost (timing step, 2 replicates x arms A-D, K = 1, seed 500, 8 fits in parallel, 1 Julia worker each).** One CYCLOPS fit
(`train_n_models` = 10): 602 s on average (473-761 s); one clocksync gate (1,000 permutations): 2.6 s; simulation 0.3 s.
With the template's 80 models one fit took 1,227 s on 4 workers (about 80 worker-minutes), which made the registered 50
replicates unaffordable (about 100 fits in the budget), hence 10 models for all runs (choice below).
<!-- src: output/q1/timing/*.csv (fit_s, gate_s) -->
Gate step: 50 fits, 4,327 s wall on 8 workers, peak 12.4 GB (`paper/results/logs/jobs.tsv`).

**Post hoc diagnostic (not registered; descriptive).** Same arm D, seed 1,500, 8 replicates, with the factor's standard
deviation multiplied by 0.1, 0.25, 0.5 (calibrated value: 1). Median truth correlation 0.30 (1 of 8 at or above 0.7),
0.16 and 0.16 (0 of 8); at scale 0 CYCLOPS crashed (Julia error in its own code, 8 of 8), so a factor-free reference is
not available at this seed size. In a single check (arm D, K = 0 equivalent, seeds = clock genes plus followers, 8
models, one replicate) the truth correlation was 0.81, so the pipeline does recover a clock when the factor is absent.
<!-- src: output/q1/diag_summary.csv; the K = 0 check was interactive and is not kept -->
So under the registered calibration the factor, not the clock, is what CYCLOPS finds, and the gate cannot pass at the
registered factor size; the registration's own rule (fix the pipeline first) applies. The question for the PI is
whether the calibration (below) overstates the shared factor.

### Simulator as built (choices in the registration, made concrete)
- 76 patients x 3 compartments, 2,000 genes: 14 GTEx clock genes, 103 followers, 1,883 further genes (top by mean
  fraction in GSE266919 CD4_T/CD8_T/myeloid). Per-gene mesor = log mean fraction in those pseudobulks; depth 3.97e6 UMI
  (their median total). `clocksync.simulate` with default `tau`, `shift_sd`, `mesor_sd`, jitter 1 h, amplitude 1, z = 1 for the three
  compartments; coupling 0.8 by psi = z^2 V (1 - c) / (2 c) (V of the 08-18 h sampling spread, psi = 0.061; arm D keeps that psi).
  Arm B: waveform amplitudes of the clock genes and followers set to 0. Arm C: `Compartment(coupled=False)`.
- Non-clock factor: K axes, each with unit-norm normal loadings on a random 10-30 % of the 2,000 genes, patient scores
  N(0, lambda_k) shared by the three compartments, applied as binomial thinning (`thin_clock`; per-gene rescale to a maximum of 1).
  lambda_k = the k-th eigenvalue of the covariance of the log fractions of the 1,986 non-clock panel genes in the 44
  Pre-treatment samples (3 compartments centred and stacked, 128 rows): 78.7, 50.6, 47.8 (total 523.7).
<!-- src: output/q1/panel.npz (lam, lam_total, n_rows, n_follow, depth); output/q1/panel_genes.csv -->
- Followers: genes with GTEx (timed fast deaths, `figure_2/output/gene_rhythms/gtex_cosinor.csv.gz`) p < 0.01 in at least
  3 of 13 tissues (141, 103 in the panel); amplitude = mean over those tissues, phase = circular mean of the cosinor coefficients.
- Seed genes: ranked by max(multiple R with cos/sin of the true time, |r| with each factor score) in a simulated timed
  reference (arm D, 76 patients, same loadings, other seed), compartments centred; with core = the 14 clock genes plus the
  top n - 14 others, without = the top n among non-clock genes.

### Choices not in the registration
- **Calibration data.** The A17 cancer pseudobulk tables hold the 14 clock genes only, so non-clock PCs cannot be computed from
  the gate-passing cohorts; fallback as allowed: GSE266919 pseudobulks (output/q2), 44 Pre-treatment samples. These
  are not gate-passing (A140 Q2 part (i): gate fails), and the log-fraction PCs contain Poisson noise (the simulator adds
  its own), so lambda is if anything too large.
- CYCLOPS: the port with `covariates` = cell type (the CYCLOPS 2.0 eigengene regression on the discrete covariate, template
  options), CPM input, template hyperparameters except `train_n_models` = 10 (template 80; cost) and `seed_mth_Gene` = 1000
  (half of the genes; the template's 10,000 exceeds the 2,000 genes). `phase_raw`, seed of the sample order 0.
- Gate statistic: median over replicates of the median over the three compartments of |circular correlation| of the phase with
  the true time (orientation free); gate passes at >= 0.7. Evaluated on the primary variant (without core clock genes).
  Cross-compartment r is the sign-sensitive circular correlation (single joint model).
- DIP: circular mean, iteratively dropping phases > pi/2 from the mean (up to 5 iterations, a drop that would empty the set is
  not applied). Six DIP bins = equal-count groups along the arc obtained by cutting the circle at the largest DIP gap
  (equal-width bins of a DIP confined to a quarter of the circle left 2 occupied bins and no cosinor; changed after the timing
  step and before the gate step was scored). Bin time = circular mean DIP; cosinor on bin means (log2 CPM + 1), F test (2, 3 df), p < 0.05, all 2,000 genes.
- Planned primary conventions (not reached): seed 500 without core genes; H140.2 = the smallest of the three
  compartments' rhythmic-gene counts >= 500; H140.3 p = max of the two lower-tail binomial p-values (arms B, C) and arm A
  required at >= 35 as well; Holm over (p1, p2, p3).
- The data seeds are arbitrary (base seed 1500 for the gate step); loadings are shared by the arms of one replicate.
- The `timing` outputs were produced with equal-width DIP bins (their `n_rhy_*` are not meaningful); timings are unaffected.
