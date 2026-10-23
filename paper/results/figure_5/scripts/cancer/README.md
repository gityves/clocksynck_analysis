# cancer: the cancer decomposition

Tumours keep (or lose) the clock-shaped structure of their clock genes; the decomposition says which way. Everything
runs on `clocksync` (paired designs, `decompose`, `Dataset`, `stats`) plus `common/` and the builders of
`paper/common/datasets/{tcga_paired,cancer_sc}`. No plotting code. Outputs: `paper/results/figure_5/output/cancer/`
(`.../subsample/` for checks; `cache/` holds the reduced TCGA tumour tables). Per-cell inputs: `data/cells_cancer`, `data/cells_cancer_cycle`, `data/cue_units`.

```
run.sh [--subsample N]      all scripts in order (full run by default; --subsample N: N random patients per cohort, few resamples)
_shared.py                  loaders (TCGA, single cell, all primary tumours) and the local helpers listed below
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_tcga_paired.py` | A52: H52.1, H52.2, secondaries | TCGA tumour vs adjacent normal of the same patient (10 cancers): r = log clock-plane variance, tumour minus normal; pooled over cancers whose normals pass the gate, swap null, pseudo-clock specificity; without LIHC / KIRC, by stage, pooled offset | `paired_damping/scripts/a52_tcga.py`, `a52_secondary.py` |
| `compute_tumour_gate.py` | A83: H83 | is a clock observable in tumours? clock-plane share, tumour vs normal of the same patients, patient-level label swaps; all primary tumours of 11 projects (descriptive) | `a83_tumour_gate.py` |
| `compute_decomposition.py` | A87 (descriptive, reading rule) | why is the tumour clock less observable: plane variance, off-plane variance, gene coherence, alignment; TCGA (vs paired normal) and single-cell malignant vs normal epithelium (alignment vs the stromal host); patient bootstrap | `a87_decomposition.py` |
| `compute_malignant_gate.py` | A86: H86 | gate in the malignant cells: pooled Poisson-corrected clock-plane share per compartment; malignant minus normal epithelium share in the same patients | `a86_malignant_gate.py` |
| `compute_delay.py` | A97: H97.1 | phase of malignant vs normal epithelium relative to the stromal host (delta < 0 = delayed); swap null, bootstrap, host-permutation informativeness | `a97_delay.py` |
| `compute_proliferation.py` | A97: H97.2 (+ pre-specified exploratory survival and breast ER) | gene coherence in high vs low proliferation tertile of TCGA tumours, with a purity companion on the same permutations | `a97_tcga_load.py`, `a97_proliferation.py`, `a97_survival.py` |
| `compute_anchored_cancer.py` | A17 (cancer targets), A18 part 3 (copy number) | anchored coupling of tumour vs adjacent-normal compartments per cell class (21 gate-passing cohorts, learned and fixed A17 anchors), pooled; malignant (copy-number) vs normal epithelium | `desync/common/run_anchored.py` (`cancer`), `followups.py` (`cnv_one`) |
| `compute_anchors_cancer.py` | A17 step 2 (the 21 cancer targets) | anchors of every cancer cohort learned leave one study out (coherence of the cohort's own pairs can only drop an anchor); checked against the frozen table | `desync/common/anchors.py` (cancer targets) |
| `compute_cycling_a24.py` | A24 cancer arm (+ note of 2026-10-13, `--definition quartile`) | y5 = phase of S, G2M and cycling epithelial cells minus G1 of the same compartment on the host of learned anchors minus epithelium, tumour and adjacent normal, pooled by DerSimonian-Laird; tumour minus normal | `desync/common/a24_cycling.py` |
| `compute_cycling_scores.py` | cycling follow-ups T1, T2, T3, C1 (`docs/cancer/ANALYSIS_PLAN.md`, amendment 2026-10-04), exploratory, BH 10 % | G2M-hi and S-hi vs non-cycling tumour epithelium (offset, log amplitude; with / without ERBB2-high patients), clock-gene classes (RORE vs E-box), G2/M fraction vs host-predicted NR1D1, leave-one-study-out of T1 | `scripts/s1_cycling_checks.py`, `s1_cycling_genes.py`, `s1_states.py` |
| `compute_cycling_cophaser.py` | C2, C2b (post hoc), C3, C4 and the decision rule of the same amendment | CoPhaser-defined G2/M vs non-cycling offset of tumour epithelium, half-life gradient, T/NK and myeloid controls, landmark-window G2/M | `scripts/s4_g2m_checks.py`, `s5_cc_diagnostics.py`, `s6_c2b.py` |
| `compute_cue_following.py` | A85: H85 | cue-following index |b| / (|a| + |b|) of each class in 9 genome-wide healthy atlases; immune minus structural classes | `a85_cue_following.py` |
| `compute_matched_cues_tcga.py` | A84, TCGA arm: H84 | tumour minus normal offset of the core clock, the tissue-matched GTEx cue set and the locked set (10 cancers, label swaps) | `a84_matched_cues.py` (TCGA part) |

The single-cell cohort pipeline (gate, E1-E3, H1-H3, own rhythm E5, gene loadings E4, microenvironment) is in
`cancer_pipeline`.

## Outputs
`a52_{per_cancer,pooled,pseudo,secondary,stage}.csv`, `a83_{paired,all_tumours,verdict}.csv`,
`a87_{tcga,tcga_mean,single_cell}.csv`, `a86_{compartments,verdict}.csv`, `a97_{delay,delay_per_study}.csv`,
`a97_{proliferation,proliferation_per_cancer,tertiles,survival_per_cancer,survival_pooled,breast_er}.csv`,
`a17_cancer{,_tumour_minus_normal,_pooled}.csv`, `a18_cnv{,_pooled}.csv`, `a24_cancer{,_pooled}[_quartile].csv`, `a84_*`, `a85_{classes,verdict}.csv`, `anchors_cancer_{relearned,vs_frozen}.csv`, `cycling_gates.csv`, `t1_{by_study,summary}.csv`, `t2_{by_gene,genes,summary}.csv`, `t3_{by_study,pooled}.csv`, `c1_loso.csv`, `c2_c4_{by_study,pooled}.csv`, `c3_genes{,_by_study}.csv`, `c2b_{landmarks,by_study,pooled}.csv`, `decision.csv`. Column names follow the legacy files, except:
`alignment` (A87 TCGA, `align_corr` in the legacy file), `align_case` / `align_ref` (`align_malignant` / `align_normal`),
`coherence_case` / `coherence_ref`. TCGA are open-access recount3 counts; nothing here is controlled access.

## Checks
Every script runs with `--subsample` (patients per cohort; minutes). Full point estimates against the legacy output files
(`legacy/desync/validation/paired_damping/output/`, `legacy/desync/common_output/`) with the full resample sizes except the
last row:

| readout | result here | legacy | note |
|---|---|---|---|
| A52 per-cancer r, offset, coherence, n pairs (10 cancers) | identical to the printed digits | | |
| A52 per-cancer bootstrap SE | 0.134-0.387 | 0.136-0.370 | bootstrap draws differ (one seed per script) |
| A52 pooled r (fixed) / random | -0.0580 / -0.0518 | -0.0566 / -0.0512 | weights from the bootstrap SEs |
| A52 swap p (lower), H52.1 | 0.191, FAIL | 0.198, FAIL | 10,000 swaps, Monte Carlo |
| A52 pseudo-clock: pooled p, H52.2 | 0.005, PASS | 0.005, PASS | 200 pseudo sets; per-cancer p within Monte Carlo (e.g. BRCA 0.050 vs 0.020) |
| A52 secondaries (w/o LIHC, KIRC; stage; weighted offset) | -0.089; 0.007 / 0.139; -1.066 h | -0.086; 0.011 / 0.141; -1.066 h | bootstrap SE 500 vs 500, other seed |
| A83 all (10 cancers: shares, p, diff, swap p; 11 projects) | identical (1,000 swaps, same seed order) | | mean diff -0.316, p 0.002, tumours gated 6/10, normals 10/10 |
| A86 all compartments (16 rows) and H86 | identical (max share difference 2e-16) | | -0.164, p 0.004, 110 patients, 15 studies |
| A87 single cell (79 patients, 12 studies) | identical (estimates, intervals) | | d_coherence -0.554 [-0.726, -0.139] |
| A87 TCGA (10 cancers x 6 estimates) | estimates identical (1e-16); interval ends differ by < 0.06 | | bootstrap draws differ |
| A97 H97.1 delay (79 patients) | identical: delta -0.430 h, p 0.311, \|c_m\| 0.273 vs null q95 0.246 (uninformative) | | |
| A97 H97.2 | -0.1022, p 0.002 (PASS) | -0.0951, p 0.007 (PASS) | the Hallmark G2M list: here MSigDB v2024.1 from `$SCCC_DATA/reference`, the legacy run used v2026.1; with v2026.1 (`SCCC_HALLMARK_GMT=...`) every row is identical |
| A97 purity companion, purity score, survival pooled HR, breast ER | identical | | HR 0.976 [0.904, 1.055] per SD (9 cancers) |
| A17 cancer coupling, psi, V, n (SRP273482, learned and fixed anchors) | identical to the printed digits | | point estimates only; the bootstrap (500, host refitted) was not run in full |
| A18.3 copy-number coupling, ERP149257 (malignant 0.2556, normal epithelium 0.1922, non-malignant 0.2128) | identical | | point estimates only |
| A17 anchors of the 21 cancer cohorts, learned and fixed rule | identical to the frozen table (21 of 21, both rules) | | frozen table written before any anchored coupling |
| A24 cancer, pooled y5 of cycling tumour epithelium (h) | -0.49 [-1.02, 0.03], 7 studies | -0.35 [-0.83, 0.13], 10 | verdict (interval includes 0) the same; studies lost to the 3-anchor rule (see deviations) |
| A24 tumour minus normal, cycling | -0.17 [-4.0, 3.7], 2 studies | -1.39 [-4.6, 1.8], 2 | uninformative in both |
| A85 H85 (9 atlases, 70 classes) | immune 0.182, structural 0.095, diff +0.087, one-sided p 0.002, **PASS** | 0.496, 0.421, +0.074, p 0.160, FAIL | **legacy bug, not a deviation**: the legacy table read the Ensembl-indexed `var_names` with gene symbols, so every gene resolved to column -1 (the last gene) and every core and cue gene held the same counts (checked: legacy ARNTL equals CLOCK equals CAMKK1 in the lung table). The legacy A85 numbers (RESULTS A85) are invalid. Here the genes are matched on `feature_name`; 3 of the 54 genes are absent from the atlases (SLC9A3R2, C20ORF27, FAM126B) |
| A84 TCGA: core minus matched cue, p two-sided | -0.83 h, 0.996 (no offset) | -0.76 h, 0.994 | set-wise offsets core -1.05 / cue -0.22 / locked -0.31 h vs -1.04 / -0.28 / -0.35; matched-set cue \|corr\| 0.151 (legacy 0.151); global A81b set 0.201 |
| C2 CoPhaser G2/M vs non-cycling, pooled (9 studies; 200 -> 40 bootstraps) | -0.02 h [-0.55, 0.51] | -0.16 [-0.91, 0.60] | same verdict (C2 interval includes 0; positive finding: no) |
| C2b (10 studies) | -0.01 h [-0.61, 0.59] | -0.10 [-0.95, 0.75] | |
| T1-T3, C1, T2 | run with `--subsample` (outputs in `subsample/`); full run not made | -0.28 h (G2M_hi), T2 RORE minus E-box +0.06 (p 0.025), T3 slope -0.02, leave-one-out | 200 bootstraps x 16 studies x 2 ERBB2 variants did not finish in 40 minutes on the loaded machine; with 40 bootstraps the CoPhaser script (C2-C4) finished in 17 minutes, the score-based script did not in 30 |

## Deviations
- **Pseudo-clock sets (A52).** `paired` takes the pseudo-gene counts from the sample table, which holds the 14 clock genes
  only. `compute_tcga_paired.py` therefore reads the genome-wide counts of the paired patients (builder `_pairs`) and calls
  the public building blocks (`pseudo_waveforms`, `cohort_states`, `log_var_ratio`) per cancer; the match is on the normals'
  mean log fraction as in the legacy run. The per-cancer pseudo p is the lower tail (H52.2 tests below the 5th percentile).
- **Hallmark list (A97).** See the check table; set `SCCC_HALLMARK_GMT` to the v2026.1 file to reproduce the legacy numbers.
- **A87 / A97 single-cell cohort** is the table `cancer_sc.from_legacy()` (the pooled legacy pseudobulks), as in the legacy
  scripts; `build()` from the counts equals it on slices (`data/INDEX_natexp_cancer.md`).
- **Anchored cancer (A17, A18.3).** Only the coupling and its replicate differences are used; `Dataset.fit` reports amplitudes
  and offsets relative to the mean over compartments, not to the study reference as the legacy table, so those two columns
  differ by construction. Gate: recomputed on the study table (1,000 permutations); the 21 cohorts are the targets of the
  frozen A17 anchor table. The A17 cancer step used one seed per target, here `--seed`.
- **Seeds.** Each script has one fixed seed; sequences differ from the legacy scripts wherever a legacy script drew several
  cancers from one generator.

- **A24 cancer.** `Dataset.fit` needs 3 anchors after removing the epithelial compartments (the registration: 2), so the rows of
  studies with fewer than 3 anchors are not estimable here (7 rows), 7 rows have fewer than 6 patients, and 41 of the 60 study x part rows are
  estimable here against 64 of 69 in the legacy run. The host is refitted in every bootstrap replicate
  (legacy: fixed). The quartile definition (`--definition quartile`) runs on the subsample only; its legacy file is `a24_cancer_pooled_quartile.csv`.
- **Cycling follow-ups.** The states and the joint model are those of `cancer_pipeline/_shared.py` (`Dataset.study`); the
  cell-cycle fits (`data/cells_cancer_cycle`) are the CoPhaser outputs and are not recomputed. Bootstrap sizes 40 here (legacy 200). With `--subsample` the cycling scripts use 3 small studies (cophaser: `--subsample 40`; at 12 patients no contrast is estimable) and the cophaser script pools them without the gate table.
- **A85.** The legacy gene matching bug (see check table); the pooled estimate is therefore a new result and is not the registered
  verdict of the legacy run. Tell the user before anything quotes A85.
- **A84 TCGA.** The Holm over IBD and TCGA of the registration needs the IBD arm (figure_5); the TCGA row alone is reported.

## Not in this folder
- T1-T3 and C1 at full size (see check table); the IBD arm of A84 (belongs to figure_5 / figure_S07).

## Local helpers (not in clocksync)
`_shared.py`: `plane` and `share_cov` (clock-plane projector and share, private in the package), `pooled_share` and
`centred_units` (Poisson-corrected pooled gate with strata, A86), `load_tumours` (all primary tumours with proliferation and
purity scores), `tcga_full_counts`; `compute_proliferation.py`: `off_scores` (A97 survival score), the Cox fit (statsmodels).
