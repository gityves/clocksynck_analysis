# A141: are the brain disease contrasts stable enough to read?

Registered in `paper/docs/registered_experiments.md` (A141, 2026-10-27, before any readout). Descriptive, no disease
hypothesis. Brain contrast sets of WP0 item 8 (12 contrasts, 5 cohorts, A17 learned anchors); observed flips recomputed
here equal the WP0 table for all 12 (max |difference| of d 3e-11).

```
RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A141_brain_stability/scripts/step3_simulate.py --jobs 6   # 30 min, 200 sims
python scripts/step12_coverage_amplitude.py      # steps 1, 2 (needs step 3 output)
python scripts/step4_map.py                      # step 4 (only contrasts not noise-driven)
```
Scripts: `_a141.py` (contrast sets, fits, host statistics, simulator), `step3_simulate.py`,
`step12_coverage_amplitude.py`, `step4_map.py`. Outputs (`output/`): `step3/` (`step3_summary.csv`, `sims.csv`,
`observed.csv`, `observed_d_*.csv`, `params_*.csv`), `step1_coverage.csv`, `step1_gtex.csv`, `step2_amplitude.csv`,
`step2_summary.csv`, `step4_map.csv`, `step4_summary.csv`.

## Step 3 (decision): 200 simulations per contrast, 95th-percentile rule
Statistic: target compartments whose case-minus-control d changes sign between the A17 host and the all-LOO host.
`p_upper` = share of simulations with a count >= the observed one.

| contrast | targets | observed flips | sim mean | sim median | sim p95 | sim max | p_upper | noise-driven |
|---|---|---|---|---|---|---|---|---|
| AD_mssm | 9 | 3 | 1.43 | 1 | 3 | 4 | 0.170 | yes |
| DLBD_mssm | 9 | 1 | 1.525 | 1.5 | 3 | 4 | 0.835 | yes |
| PD_mssm | 9 | 7 | 1.58 | 1 | 3 | 5 | 0.000 | **no** |
| SCZ_mssm | 9 | 3 | 1.375 | 1 | 3 | 6 | 0.150 | yes |
| VaD_mssm | 9 | 5 | 1.47 | 1 | 4 | 5 | 0.010 | **no** |
| ASCVD_hbcc | 9 | 6 | 1.32 | 1 | 3 | 6 | 0.005 | **no** |
| BD_hbcc | 9 | 1 | 1.38 | 1 | 3 | 4 | 0.790 | yes |
| SCZ_hbcc | 9 | 0 | 1.835 | 2 | 4 | 6 | 1.000 | yes |
| AD_radc | 9 | 1 | 1.69 | 2 | 4 | 5 | 0.820 | yes |
| DM_radc | 9 | 1 | 1.14 | 1 | 3 | 4 | 0.690 | yes |
| AD_seaad_dlpfc | 8 | 2 | 1.005 | 1 | 3 | 4 | 0.290 | yes |
| pd | 34 | 15 | 2.195 | 2 | 5 | 6 | 0.000 | **no** |

Eight contrasts are noise-driven (not estimable at this design: nothing can be said about their compartments). Four are
not: PD_mssm, VaD_mssm, ASCVD_hbcc, pd. Wall time 30 min, so the registered 200 simulations were kept.
<!-- src: output/step3/step3_summary.csv; output/step3/observed.csv -->

## Step 1: coverage of the cycle
Rayleigh R of the host phases (A17-anchor host; the all-compartment host is in the csv). Reference: 2.5 / 97.5 % of R in
the uniform-sampling simulations of step 3 at the same donors, fitted amplitudes and noise. "All donors" = the labelled
donors of the contrast with the most controls (control set of A135). Full cycle = R inside the simulated range.

| cohort | donors | n | R | Rayleigh p | sim R 2.5 % / median / 97.5 % | full cycle |
|---|---|---|---|---|---|---|
| psychad_mssm | control | 257 | 0.027 | 0.83 | 0.006 / 0.021 / 0.047 | yes |
| psychad_mssm | all | 370 | 0.032 | 0.69 | 0.004 / 0.019 / 0.042 | yes |
| psychad_hbcc | control | 93 | 0.052 | 0.78 | 0.008 / 0.042 / 0.095 | yes |
| psychad_hbcc | all | 132 | 0.051 | 0.71 | 0.007 / 0.035 / 0.076 | yes |
| psychad_radc | control | 26 | 0.195 | 0.37 | 0.014 / 0.079 / 0.173 | **no** (above) |
| psychad_radc | all | 67 | 0.094 | 0.56 | 0.012 / 0.046 / 0.103 | yes |
| seaad_dlpfc | control | 41 | 0.063 | 0.85 | 0.011 / 0.056 / 0.137 | yes |
| seaad_dlpfc | all | 80 | 0.027 | 0.94 | 0.009 / 0.042 / 0.101 | yes |
| parkinson | control | 23 | 0.071 | 0.89 | 0.024 / 0.084 / 0.178 | yes |
| parkinson | all | 96 | 0.149 | 0.12 | 0.006 / 0.046 / 0.111 | **no** (above) |

No Rayleigh test rejects uniformity (smallest p 0.12). With the all-compartment host: psychad_mssm all donors (R 0.046,
range up to 0.042), psychad_radc all (0.157 vs 0.097) and parkinson all (0.132 vs 0.105) fall above the range;
parkinson all and radc control are above on the A17 host too.

GTEx fast-death brain (host on all 13 brain tissues with >= 40 donors, 212 donors): R 0.048 (Rayleigh p 0.62); on the
three cortical tissues (197 donors) R 0.105 (p 0.12). The cohort Rs (0.027-0.195) are of the same order.

Host spread V is frame-dependent: the factor model fixes the scale only through its reference compartment and a weakly
coupled host trades V against |z|^2 (V in the simulations varied 3-fold between replicates at the same truth), so V is
reported (`V`, `sim_V_median`, `ecc` columns of `step1_coverage.csv`) but not used for a verdict; the invariant is
|z|^2 V (step 2). Host eccentricity (smallest / largest eigenvalue of the host covariance; 1 = round): real 0.22-0.77
(A17 host), simulation medians 0.56-0.84; GTEx brain 0.21 (all tissues), 0.14 (cortical). Real hosts are less round than the
uniform simulations in every cohort, as GTEx timed brain also is.
<!-- src: output/step1_coverage.csv; output/step1_gtex.csv -->

## Step 2: absolute amplitude against GTEx brain
Absolute coupled amplitude sqrt(z2 V) (host-referenced amplitude times host amplitude, clock-plane units of the timed
GTEx fit) of each compartment on the control donors, A17 host. GTEx fast brain regions (13, `abs_amp`): min 0.438,
median 0.576, max 0.668. Same scale = inside that range. Reported, not tested.

| cohort | compartments | median amplitude | min-max | within GTEx range | above | below | median / GTEx median |
|---|---|---|---|---|---|---|---|
| psychad_mssm | 9 | 0.444 | 0.338-1.147 | 1 | 4 | 4 | 0.77 |
| psychad_hbcc | 9 | 0.760 | 0.470-1.624 | 3 | 6 | 0 | 1.32 |
| psychad_radc | 9 | 0.325 | 0-0.716 | 1 | 2 | 6 | 0.56 |
| seaad_dlpfc | 8 | 0.605 | 0.477-1.030 | 5 | 3 | 0 | 1.05 |
| parkinson | 34 | 0.443 | 0-1.327 | 3 | 14 | 17 | 0.77 |
| all | 69 | 0.558 | 0-1.624 | 13 | 29 | 27 | 0.97 |

The median amplitude over all 69 compartments is close to the GTEx median (0.97 x), but only 13 of 69 compartments lie
in the GTEx range; compartments scatter on both sides (vascular and microglia high, oligodendrocytes and OPC low).
Ten compartments have amplitude 0 (debiased |z|^2 floored at 0).
<!-- src: output/step2_amplitude.csv; output/step2_summary.csv; paper/results/figure_2/output/organs/tissues.csv -->

## Step 4: descriptive map (PD_mssm, VaD_mssm, ASCVD_hbcc, pd)
Agreement = same sign and |d_A17 - d_LOO| <= pooled bootstrap SE. Share agreeing (`step4_summary.csv`):

| contrast | all targets | expected agree (neurons, astro/oligo) | neurons | astro/oligo | expected noise | expected agree, not noise | expected noise, not agree |
|---|---|---|---|---|---|---|---|
| PD_mssm | 0/9 | 0/6 | 0/4 | 0/2 | 0/5 | 0/2 | 0/1 |
| VaD_mssm | 3/9 | 2/6 | 2/4 | 0/2 | 2/3 | 0/3 | 0/0 |
| ASCVD_hbcc | 1/9 | 1/6 | 0/4 | 1/2 | 0/3 | 1/3 | 0/0 |
| pd | 13/34 | 7/18 | 2/8 | 5/10 | 5/22 | 4/7 | 2/11 |

The expectations are only partly met: in pd, expected-agree compartments agree in 7 of 18 (neurons 2 of 8, astrocytes
and oligodendrocytes 5 of 10) against 5 of 22 for the expected-noise set (4 of 7 vs 2 of 11 once the two sets are
separated); in the three PsychAD contrasts at most 3 of 9 compartments agree, neurons included (0-2 of 4), and
expected-noise compartments are not less in agreement than the others (VaD_mssm 2 of 3). Reading: the neuron classes do not
agree among themselves; astrocyte/oligodendrocyte agree only in pd and ASCVD_hbcc.
<!-- src: output/step4_summary.csv; output/step4_map.csv -->

The lead should draft the time-of-death field request (outward-facing: PI decides). Not drafted or sent here.

## Choices not in the registration
- **Simulator.** `clocksync.simulate.simulate` draws Poisson counts at one scalar depth, which cannot keep the measurement covariance V of the real rows
  (depth varies by orders of magnitude between donors and compartments). The simulation is at the level of the clock-plane state (the level of
  `clocksync.simulate.Like`): beta = a_k exp(i theta_k) exp(i phi_s) + N(0, V_row + psi_k I), host phase phi_s uniform on 24 h
  (shared by a donor's compartments), a_k = sqrt(z2 V), theta_k = `offset_h` and psi_k from `Dataset.fit` on the contrast's
  donors (control + case pooled, so no case effect, A17 host; the anchor targets on their leave-one-out host), real rows, V,
  strata, intercept groups and donors per group. Fits on simulated data replace `Dataset.data` in a shallow copy (uses
  `clocksync._core.host` internals for the host statistics only). Compartments without an estimate (a handful of tiny `pd`
  compartments, not targets) get amplitude 0 and the median psi; per-compartment intercepts are re-estimated.
  Recovery was checked on SCZ_hbcc: two refits of one simulated set recovered most amplitudes to within about 0.15 of the input (microglia less well: 0.47 and 0.84 against 0.99) and psi to within about 0.2 (checked once, not tabulated).
- **Cohort = dataset; simulations per contrast** (12 contrasts x 200), not per dataset: each contrast has its own donors and groups.
  Observed counts are those of WP0 item 8 (equal to the recomputed ones).
- **Targets** = the compartments of the WP0 table of the contrast (>= 5 donors per group, non-"other" family).
- **d** = `Dataset.contrast` with no permutations or bootstrap (d does not depend on them).
- **95th percentile** = the smallest count with cumulative frequency >= 0.95 (`inverted_cdf`), so the threshold is an attained
  count; "noise-driven" is observed <= it.
- **Step 1.** Cohort controls = controls of the contrast with the most controls (A135); "all donors" = the labelled donors of that
  contrast (not every donor of the dataset). Host = A17 anchors of that contrast, fitted on those donors; the all-compartment
  host (no leave-one-out) as sensitivity. Reference simulations are the step 3 simulations of that contrast, hosts refitted on
  the control donors and on all donors. R = |mean exp(i phase)|, phase = arg(x + iy) of the host mean; "full cycle" = R inside
  the two-sided 2.5-97.5 % range (no R fell below it). Spread V is reported with eccentricity and not used for the verdict (frame-dependent, see above).
- **GTEx reference** = clocksync host (all brain tissues as anchors) on the fast-death cohort (`organs/_shared.cohort("fast")`,
  tissues with >= 40 donors); only n, R, V, eccentricity written. The registered range of `tissues.csv` is used for step 2.
- **Step 2** uses the control donors of the cohort's control-set contrast with the A17 host, on the contrast's targets; the amplitude is
  the coupled part sqrt(z2 V) (the own noise excluded).
- **Step 4** definitions (agreement, lowest third of control coupling, neuron and glia classes) are in `scripts/step4_map.py`, set before scoring.
