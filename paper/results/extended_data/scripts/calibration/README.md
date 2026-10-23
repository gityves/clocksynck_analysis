# calibration: small n, coverage, boundary, estimator checks

Calibration of the estimators at small donor numbers and by simulation: the donor-subsampling rule of A10, the acceptance simulations of the
anchored host (A17 step 5, V12), the von Mises per-cell moment model (V11) and the ORPHEUS phase-dispersion bias (V14, amendment A19 part 3a).
Estimates come from `clocksync` (`Dataset.fit`, `Dataset.gate`, `Dataset.resample`, `simulate`, `power.mde`, `cells.vonmises`); data through
`common.data.load` (heart, kidney, lung atlases, cancer single-cell cohorts). No plotting code. Outputs: `paper/results/extended_data/output/calibration/`
(`.../subsample/` for checks).

```
run.sh [--subsample N]               all scripts (full run by default: small_n and the coverage part take hours)
compute_small_n.py
compute_anchored_simulations.py
compute_vonmises_calibration.py
compute_orpheus_bias.py
compute_orpheus_cells.py
compute_age_acceptance.py
```

## Analyses
| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_small_n.py` | A10 step 1 (rule: smallest n with gate pass >= 0.9, coverage >= 0.85, median error <= 0.15 in heart and kidney, else 15) | healthy heart (45 donors), kidney (26), lung (198): 40 random donor subsets of n = 8, 10, 12, 15; gate (200 permutations) and coupling of every compartment with a 200-draw donor bootstrap, against the full-data coupling: gate pass rate, median absolute error, coverage, boundary share | `desync/common/small_n.py` |
| `compute_anchored_simulations.py` | A17 step 5 (V12): identity, coverage, boundary, power | simulated cohorts with known couplings: (1) free host follows its reference compartment, anchored host follows host time (similarity >= 0.9, stability); (2) coverage >= 0.9 and |bias| <= 0.05 of the 95 % interval, c in {0.3, 0.6, 0.9}, precise or noisy anchors; (3) boundary share at c = 0.95, 1; (4) power curve under damping, level at s = 1 <= 0.07 | `scripts/v12_anchored.py` (identity, coverage, bound, power; bias: not ported) |
| `compute_vonmises_calibration.py` | post hoc (V11, 2026-10-08) | recovery of per-cell amplitude lam, synchrony R, dispersion D = lam^2 (1 - R^2) from simulated single cells at human depths (3000 UMI per cell, 2.4 or 8 clock UMI, 600 cells, 100 groups, metacells of 1 or 10) | `scripts/v11_vm_calibration.py` |
| `compute_orpheus_cells.py` | A19 part 3b | within-sample dispersion D and synchrony R of excitatory neurons, AD minus control, `AD_mssm`, `AD_radc`, `AD_seaad_dlpfc` (A18 per-cell estimator of `backbone`); inverse-variance pool of d_D / D_control | `cell_split.py` (AD rows) |
| `compute_age_acceptance.py` | A19 part 2, validations 1-2 (null, planted -0.1 / decade locked amplitude and +0.15 gamma) | `Dataset.covariate` at the real designs of `psychad_mssm` and `lung_hlca` (highest, median, lowest coupling compartment): null by permuting ages within stratum, planted effects by thinning; rejection <= 0.07, median `_c` within 25 % | `20_age/scripts/acceptance.py` |
| `compute_orpheus_bias.py` | A19 part 3a (V14) | ORPHEUS moment estimator of cell-to-cell phase dispersion when subject phases carry error: simulation and closed form | `scripts/v14_orpheus_bias.py` |

## Outputs
`small_n_{draws,summary,full}.csv`, `anchored_sim_{identity,coverage,boundary,power,summary}.csv`, `vonmises_{runs,summary}.csv`,
`orpheus_bias_{simulation,summary,analytic}.csv`, `orpheus_cells_{split,contrasts,pooled}.csv`, `age_acceptance_{runs,summary}.csv`, `vonmises_D_recovery.csv`. Column names follow the legacy files (`coverage`, `median_abs_error`, `gate_pass`,
`boundary`; `check, value, criterion` in the simulation summary). `small_n_draws.csv` additionally holds the package's `boundary` flag
(ML psi at its floor) next to the legacy `coupling >= 0.999` share.

## Checks
Every script runs with `--subsample N` (`compute_small_n.py`: N donors per set, sizes 8 and 12, 2 draws, 5 bootstraps; `compute_anchored_simulations.py`: N cohorts per
part, 10 bootstraps, 20 permutations; `compute_vonmises_calibration.py`: N groups of 100 cells, 1 replicate; `compute_orpheus_bias.py`: N subjects, 40 cells, 10
genes). Point estimates against `legacy/` outputs:

| readout | here | legacy | note |
|---|---|---|---|
| V12 identity, 4 cohorts: similarity of free host (ref in anchor module) to host time; of free host (ref in competing module) to the competing state; of anchored host to host time | max difference 2e-16, 0, 0 | `outputs/v12_anchored/identity.csv` | free host with the reference in the competing module is closer to the competing state in 4/4 cohorts (legacy 100/100), with the reference in the anchor module in 0/4 (0/100) |
| V12 coverage, 24 cohort fits (c x anchor noise x 4 reps), point estimate of the target | max difference 0.0 | `coverage.csv` | intervals from 10 bootstraps in the check, not compared |
| V12 power: damped coupling under s = 0.6, 0 and the implied drop | c_case 0.9355, 0.8393, 0; drop 0.0962 | identical | rejection at s = 1: 0.00 with 4 replicates (legacy 0.07 with 100); full 100-replicate power run not repeated |
| V14 ORPHEUS, full run (400 subjects, 300 cells, 60 genes, 5 replicates) | median sigma_hat 1.13 h at (sigma_c, sigma_s) = (1, 0), 2.47 h at (1, 3), 3.03 h at (3, 0) (12 cells of the grid) | identical (max difference 0) | same generator and seed |
| V11 von Mises, full run | see Deviations: the legacy summary cannot be reproduced by the current estimator (legacy `sccc` gives the same unstable values) | `summary.csv` | median D_hat / D 0.86 (IQR -0.56 to 1.49), 88 % of the 96 replicates degenerate (R_hat at 0 or 1, lam_hat > 20) |
| A19 3b, full run (30 bootstraps; excitatory D control / case, d_D) | MSSM 0.904 / 0.859, -0.045; RADC 0.835 / 0.774, -0.061; SEA-AD 1.168 / 1.162, -0.007 | identical point estimates (`a18_cell_split_contrasts.csv`) | intervals from 30 draws are narrower (MSSM -0.26 to 0.26 against -0.41 to 0.40); pooled relative change -2.1 % (-18 to +14) against -2 % (-23 to +19) from 100 draws |
| A19 acceptance, `--subsample 80`, 3 simulations | runs (6 compartments, null and planted); numbers not interpretable at 3 simulations | `desync/20_age/output/acceptance_summary.csv` | full run (200 simulations x 200 null draws x 6 compartments x 2 scenarios) not made |
| V11 re-derivation (the RESULTS statements, from the stored 96 runs, `--summarise-only`) | D_hat - D, SD over runs: 2.4 UMI, m = 1: -1.18, 1.67; 8 UMI, m = 1: -0.07, 0.78; metacells m = 10: SD 7.8 and 2.5 (3-5 times larger) | RESULTS: D unbiased, SD 2.5 and 0.7; metacells inflate the variance 2-5x | lam and R not identified: 75-96 % degenerate fits (R_hat at 0 or 1, lam_hat > 20), as RESULTS (split rests on Assumption W). The stored legacy summary shows D_hat / D 1.2-5.6 (upward bias): the superseded free-phase-per-group estimator that RESULTS says was dropped |
| A10, full-data couplings (heart 7, kidney 5, lung 18 compartments) | heart 0.51-0.79, kidney 0.34-0.79, lung 0.07-0.66 | 0.58-0.90, 0.44-1.00, 0.60-1.00 | different estimator, see Deviations |
| A10 reduced run (kidney, n = 15, 8 draws) | not completed within the time limit (> 10 min); only the `--subsample 30` run finished | kidney n = 15: gate pass 1.00, median error 0.054, coverage 0.875, boundary 0.36 | |

Not rerun in full: the A19 acceptance simulations, the A10 draws (40 subsets x 4 sizes x 3 sets x 200 bootstraps: hours), the V12 coverage (1,200 cohorts x 200 bootstraps), identity
stability on 20 cohorts and the 100-replicate power curve (hours).

## Deviations
- **A10 estimator.** The legacy small-n calibration used the pre-A17 free-host fit (`Study.fit`, one host from all compartments, no small-sample
  correction of |z|^2 and psi). `Dataset.fit` with every compartment an anchor reads each compartment on the host of the others with the A17
  corrections ("free_corrected" in V13): couplings are lower (full-data heart 0.51-0.79 against 0.58-0.90) and the legacy ranking is not
  expected to repeat. The uncorrected estimator is not in the public API. The A10 rule must be re-derived with the full run (not made here).
- **Lung set.** `lung_hlca` (healthy donors, the A17 lung builder, strata `study|assay`) replaces the legacy A4 healthy set; the subsets draw from
  `sorted(donors)` (the legacy drew from the table order), so no draw is identical to a legacy draw.
- **V12 coverage target.** `Dataset.fit(targets=["k"], gate=False)` gives the same anchored coupling and the same interval design (donor bootstrap,
  host refitted, 200 draws) as `anchored.estimates`; the bootstrap random stream differs.
- **V12 stability.** Local similarity of bootstrap hosts to the full host (share of 100 `Dataset.resample` hosts with similarity > 0.7).
- **V11.** The legacy baseline came from the 21 cohorts passing the main gate (`outputs/studies/summary_main_gate.csv`); here every `cancer_sc`
  cohort with a `tumour:T_NK` compartment enters the baseline (normalised clock-gene fractions 0.001-0.51, near the legacy ones). The legacy summary was
  made with an earlier estimator; the current rotation-invariant one is unstable at these depths (note in the script docstring).
- **V14.** None (no clocksync call).

## Not ported
- **V12 part 3 (lower bound at the boundary)** and **part 5 (per-fit bias prediction `design_bias`)**: they need the simulator at a fit's real design (host
  means and posterior covariances, fitted intercepts, debiased loading) and the estimator on simulated rows; both live in `clocksync._core.host` and the host
  posterior covariance is not public (`Fit.host` has means only). Not rewritten locally (it would duplicate the package estimator).
- **V13**: ported in `simulations` (builders `colon_prepub`, `skin_dermis_zt`); only the A1-A16 free estimator is missing.
- **V15 (variance decomposition of a contrast at real designs):** same obstacle as V12 parts 3 and 5 (simulator at the real design, `on_host` on given hosts).
- **V16 / V18**: ported in `simulations` (import of the local hierarchical model of `spread_multiplicity`).
- **A19 part 3b host-uncertainty proxy** (tr S_u of the host, AD against control): `Fit.host` has posterior means only.
- **A17 simulations, `small_n` step 2** (small-n secondary analyses of A8): they are the contrasts of figure 5.

## Deviations of the A19 scripts
- **Acceptance (`compute_age_acceptance.py`).** Count-level simulation at the real design (permuted ages, thinning) instead of the legacy state-level simulation from the fitted model; the null
  is therefore a true null of the real clock data. Legacy numbers are not reproduced draw by draw. z and psi of the planted compartment come from `species/_shared.known_host` on the fit's host.
- **3b.** The excitatory pool is the inverse-variance (fixed-effect) pool of d_D / D_control with SE = interval width / 3.92, which reproduces the RESULTS value from the legacy contrasts (-2 %, -23 to +19 %). The "+22 %" implied by ORPHEUS is a hand calculation, not in code.

## Local helpers (not in clocksync)
- `similarity` (|complex correlation| of two host states) and the bootstrap stability loop (`compute_anchored_simulations.py`).
- `split_projectors` (gene-half projectors of the von Mises estimator; the package's `cells.vonmises._split_projectors` is private).

## API gaps
- `Dataset.fit` with the pre-A17 free estimator (or a flag to switch the small-sample corrections off), for A10 and V13.
- Public `lower_bound` and `design_bias` of an anchored coupling (V12 parts 3 and 5), a simulator at a fitted design (`simulate.like(dataset)`) and host posterior covariances for V15.
- `cells.vonmises.split_projectors` as a public function.
