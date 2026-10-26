# windows: sampling windows and internal time

Does the method's variance scale with the sampling window, and can the host place a cohort on the internal day? Public `clocksync`
(`Dataset.fit / timed`, `timed_agreement`, `simulate`, `waveform.fit_waveform`) plus the local A25 estimator `_internal_time.py`. No plotting.
Outputs: `paper/results/figure_1/output/windows/` (`.../subsample/`). GTEx rows are windows (aggregates), never donors.

| script | amendment | content | legacy |
|---|---|---|---|
| `compute_window_check.py` | A18 part 1 (1a GTEx, 1b HLCA) | ratio of shared covariance between donor sets vs 1 - R_W^2 | `common/window_check.py` |
| `compute_internal_time.py` | A25 (V25.1 sim, V25.2 GTEx, V25.4 colon, R25 cohorts, T1-T3) | window centre / width from host spread | `common/a25_internal_time.py` |
| `compute_nontemporal.py` | A26 | g of the non-temporal state, score R^2 on time, V26.2, V26.3 (counts corrected by thinning) | `common/a26_nontemporal.py` |
| `compute_baboon_windows.py` | A40 | 3-6 animal windows: tissue offsets, ordering, centre | `validation/baboon/scripts/a40_window.py` |
| `compute_chen_windows.py` | A25 V25.3 | frame transfer in Chen 2016 prefrontal cortex (microarray, ZT known): axis / curvature centres, ratio r, internal time per donor | `common/a25_internal_time.py` (chen_step) |
| `compute_breast_timed.py` | A25 note V25.5 (post hoc) | Li 2024 breast (normal + tumour, resection hour): per-patient internal time vs recorded time, hospital-prior and curvature flips, 2,000 shuffles, patient bootstrap | `validation/scripts/v25_5_breast_timed.py` |

## Checks
| readout | here | legacy | note |
|---|---|---|---|
| window check GTEx 08-18 h: ratio, observed / predicted | 0.764, 1.37 [1.05, 1.80] (M3, 500 boots) | 0.764, 1.37 [1.05, 1.80] | other windows and HLCA (0.521, 1.15) identical point estimates |
| A25 V25.1 acceptance (M3, 50 replicates per cell) | PASS (n 40, coupling 0.8, w 6/9/12 h: axis error 0.34/0.27/0.40 h, flip 0.91/1.00/1.00) | FAIL at w 6 h (flip 0.69) | different generator (`clocksync.simulate` counts); see results_preliminary 2026-10-25 |
| A25 V25.2 (M3) | median axis error 1.26 h, flip 0.78, 18 windows: PASS | PASS | |
| A25 R25 T1 / T2 / T3 (M3, cancer main tables, rerun 2026-10-26) | 17 living cohorts: p 0.0014 PASS / 14 of 17, p 0.006 PASS / p 0.176 not supported | 17: 0.010 / 0.001 / 0.214 | all 50 cohorts fit; anchor sets as legacy |
| A25 GTEx all donors: c_axis12, c_curv, r, w, lambda, kappa, shift | 11.30, 23.30, 0.67, 20.16, 1.90, 0.389, -0.98 | 11.30, 23.30, 0.67, 20.20, 1.91, 0.389, -0.98 | `host_cov` reproduces the model's Sigma_f |
| A40 full-day frame F, c_curv, r_axis | 6.56, 23.68, 0.52 | 6.56, 23.68, 0.52 | |
| A40 m = 5 (12 windows): centre errors, flips | identical; flips 11/12; median Fisher-Lee r 0.552, median abs Spearman 1.0 | r 0.516; 11/12; 1.0 | offsets read on one full host (legacy: leave-one-out per tissue); primary (1) passes here (> 0.5), narrowly as legacy |
| A26 score R^2 on time: Chen (g_GTEx), colon (g_GTEx) | 0.19, 0.62 (fail) | 0.18, 0.64 | rejection reproduced |
| V25.3 Chen (33 rows): c_axis12, c_curv, r, flip_ok, err_curv, window `all` | 9.98, 9.98, 0.17, False, -6.66 h | 10.03, 10.03, 0.16, False, -6.62 h | all window rows agree to the first decimal; microarray as pseudo-counts (below); acceptance: median abs error of c_curv 4.15 h (> 3 h), flips 21 of 32: V25.3 fails, as the run note |
| V25.5 primary (30 M reads): n, Spearman, circular corr., MAE (hospital prior), offset, flip agreement | 42, 0.558, 0.565, 1.124 h, -2.209 h, False | 42, 0.558, 0.565, 1.124 h, -2.213 h, False | 10 M / 60 M / breast waveform / joined patient: 0.556, 0.559, 0.548, 0.536 vs identical; intervals Spearman [0.30, 0.77] vs [0.28, 0.77] (200 patient bootstraps) |
| A28-related E2 etc. | see figure_S03 | | |

## Deviations
- **Cohort n (R25).** `n` is every donor of the cohort table, and the donor bootstrap draws from all of them; legacy counted and
  resampled the donors with at least one anchor compartment (25 of 50 cohorts differ, e.g. GSE178341 60 vs 34). Donors without an anchor
  do not enter the host; in the bootstrap they only vary the number of host donors per replicate.
- Sigma_f is rebuilt (`host_cov`) from posterior host means plus mean posterior covariance; the 12 h flip, frame and bootstrap as A25 (bootstrap rebuilds a `Dataset` per draw).
- A26 correction acts on counts (binomial thinning), not on states.
- V25.1 uses `clocksync.simulate` counts (legacy simulated states). Defaults: 50 replicates; windows' tissues >= 8 donors, >= 20 donors.
- **V25.3** (`compute_chen_windows.py`): `Dataset` takes counts, so the 14 clock-gene log2 intensities enter as pseudo-counts 1e6 exp(ln2 (y - region median)) at a constant depth: equal gene weights, negligible Poisson variance (the legacy states had V = 0). The two-region host (< 3 anchors) uses the local `compute_breast_timed.two_host_estimate`.
- **V25.5** (`compute_breast_timed.py`): `Dataset.fit` needs >= 3 anchors, so the 2-compartment (normal, tumour) frame shift is read from their cross-covariance (local: not in clocksync); counts are approximated from TPM x length (builder `breast_timed_gse233242`).
- Cohorts with < 3 anchors are not estimable (`Dataset.fit`).

## Not ported / not rerun
- V26.1 and R26 (corrected cohorts: `spread_multiplicity` has `compute_r26_cohorts.py`).
- Before M3 (2026-10-25, 5.6 h, all parts full) run only on subsamples: V25.1 grid, R25 cohorts with bootstraps, A40 sizes 3, 4, 6.

## API gaps
`Dataset.fit` for two compartments (frame shift and states; V25.3 / V25.5 / A28 E4 use a local estimator); host covariance Sigma_f (`Fit.host` has posterior means only); `Dataset.fit` with 2 anchors; array / state-level `Dataset` input; a per-sample state accessor.
