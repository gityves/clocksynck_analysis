# simulations: validations of the pooling model, of the estimators against known time, and of the A22 / A26 rules

Ports of pre-registered validations that did not fit the other folders. Estimates come from `clocksync` (`Dataset.fit`, `Dataset.resample`,
`Dataset.coherence` through `host_identity/compute_group_anchors.confirm`); the hierarchical model is the local copy
`spread_multiplicity/_hierarchy.py` (imported, not copied). No plotting code. Outputs: `paper/results/extended_data/output/simulations/`
(`.../subsample/` for checks).

```
run.sh [--subsample N]            all scripts
compute_hierarchy_validation.py   V16, V18
compute_celltype_truth.py         V13 (mouse colon, mouse skin)
compute_a22_colon.py              A22 V17.2
compute_r26_cohorts.py            A26 R26
```

## Analyses
| script | amendment, hypothesis | content | legacy script replaced |
|---|---|---|---|
| `compute_hierarchy_validation.py` | A20 V16, A23 V18 | null (theta = 0; calls lfsr < 0.05 <= 5 %, cell intervals cover 0 in >= 93 %) and planted data sets at the real structure of the A20 / A23 pooled tables; 2 chains x 1000 draws | `scripts/v16_hierarchical.py`, `a23_phase.py v18` |
| `compute_celltype_truth.py` | A19 part 4 (V13; method check, not a test) | coupling of every cell type against the coupling to known ZT (truth, same estimator with the host = time), free_corrected and anchored estimators; animal bootstrap | `scripts/v13_celltype_truth.py` |
| `compute_a22_colon.py` | A22 V17.2 | coherence group-confirmation rule on C-DC vs R-DC colon (must drop lymphatic EC) | `a22_group_anchors.py colon` |
| `compute_r26_cohorts.py` | A26 R26 (T1-T4) | the A25 cohorts corrected with g_GTEx (g_Chen sensitivity), then the A25 window estimator | `a26_nontemporal.py apply` |

New builders: `data/colon_prepub` (pre-publication mouse colon h5ad, 8 animals x 2 regions) and `data/skin_dermis_zt` (GSE223109 skin, 18 animals), listed in
`data/INDEX_natexp_cancer.md`. Both agree exactly with the legacy V13 pseudobulks (counts, depth, n_cells; 195 colon rows after the >= 30 cells and 75 % presence
rules, 53 skin rows).

## Checks (subsample runs; point estimates against `legacy/` outputs)
| readout | here | legacy | note |
|---|---|---|---|
| V16 / V18, data sets 0-2 (y1, y2, y3 null; y1 planted; V18 null and planted): `called`, `cells_cover0`, planted medians | identical to 3 decimals in all 15 runs | `outputs/v16_hierarchical/runs.csv`, `common_output/a23_v18.csv` | same simulated data and sampler seeds. `scenario = null` is read back as NaN by pandas (as legacy) |
| V13 colon (full point estimates, 3 bootstraps): mean error against truth, free_corrected / anchored; Spearman | +0.121 / +0.108; 0.928 / 0.958 | +0.12 / +0.11; 0.93 / 0.96 (RESULTS part 4) | intervals not compared (3 bootstraps) |
| V13 skin: immune coupling, truth / anchored | 0.80 / 0.77 | 0.80 / 0.75 | 3 fibroblast parts as anchors (see Deviations) |
| A22 V17.2 | lymphatic EC not dropped (delta -0.40, SE 0.16, e_case 0.048): fails, as legacy | not dropped (delta -0.11, SE 0.07, e_case 0.044) | e_control / delta differ, see Deviations |
| R26 (40 donors per cohort, 8 cohorts, 20 bootstraps) | runs; T1-T3 uninterpretable (2 interpretable cohorts) | `a26_tests_*.csv`: T1 pass, T3 p 0.058 | full run not made |

## Deviations
- **V13 estimators.** The A1-A16 free estimator (`Study.fit`) is not public: only free_corrected and anchored are compared (the legacy column "free" is not ported; RESULTS quotes +0.29 for it).
  Skin: `Dataset.fit` needs >= 3 anchors, so fibroblasts are split into 3 random parts (legacy 2 halves).
- **A22 V17.2 numbers.** The legacy `a22_colon_confirm.csv` was computed on a pseudobulk csv whose empty `stratum` was read as NaN; with that column the legacy pair coherence
  is not centred within stratum (e_control 0.16-0.32). Here the colon is one stratum (centred): e_control 0.45-0.55, delta -0.17 to -0.40. Reproduced with the legacy table and its NaN
  stratum: the legacy values exactly. The verdict (lymphatic EC not dropped, V17.2 fails) is the same. Other legacy coherence runs on tables with an empty stratum may carry the same artefact.
  Recomputing the A21 contrasts on a confirmed host is only needed when an anchor is dropped; not made.
- **V16 / V18 inputs.** The A20 / A23 tables of `spread_multiplicity` when a full run exists there, else the frozen tables of `legacy/desync/common_output` (same structure).
- **R26.** The correction acts on counts by thinning (as `windows/compute_nontemporal.py`); the cohorts and the estimator are those of `windows/compute_internal_time.py`.

## Not ported
- **V26.1** (Chen 2016 corrected with g_GTEx): the arrays cannot enter `Dataset` (no count-level input); it is the same obstacle as V25.3.
- **V15** (variance decomposition at real designs) and **V12 parts 3 and 5**: see `calibration/README.md`.
- Not rerun in full: V16 / V18 (100 null x 3 estimands + 50 planted; about an hour at 4 jobs), V13 with 200 bootstraps (one bootstrap of the colon takes about 4 min for the free_corrected fit),
  R26 with 200 bootstraps per cohort, A22 V17.2 at 200 / 200 (done: seconds).

## Local helpers (not in clocksync)
`species/_shared.known_host` (coupling to a known host), imported; the hierarchical model (`_hierarchy.py`, imported).

## API gaps
Public `Dataset.fit(flag)` for the pre-A17 free estimator; array-level input to `Dataset`.
