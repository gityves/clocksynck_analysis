# gene_rhythms: gene rhythms of cell types and their truth tests

Untimed gene-level rhythms per cell type (the loading of every gene on the host state, `Dataset.gene_rhythms`) against timed truths:
GTEx bulk (per-gene cosinor on the time of death), the mouse colon chrono-atlas (known ZT) and living blood (GSE39445). Also the
pooled gene-rhythm maps (A50), the donor-count curve (A39), the antiphase of non-clock genes (A30), the bulk clock gain
(A29) and the blood composition and blood-clock tests (A32, A35). Clock quantities come from `clocksync` (`Dataset`
`fit` / `gene_rhythms` / `timed` / `contrast`, `paired.projector`, `waveform.gtex`, `stats`); genome-wide tables from the builders
`paper/common/datasets/*_genome`, `gene_map_{lung,dlpfc}`, `colon_chrono_genome`, `blood_timed_gse39445`; GTEx through `common/gtex.py` and
the labelled h5ad. No plotting code. Outputs: `paper/results/figure_2/output/gene_rhythms/` (`.../subsample/` for checks).
GTEx time of death is controlled access: only per-gene and per-tissue aggregates are written, never per-donor rows.

```
run.sh [--subsample N]      all scripts in order (full run by default; --subsample N: N donors, few permutations)
_shared.py                  genome tables, GTEx truth, recombination, scores (local helpers listed below)
_classes.py                 A30 gene classes (fixed before scoring; copy of the legacy classes.py)
anchors_a31b.csv            A17 learned / fixed anchors re-run on the A31b atlases (frozen copy of the legacy anchors_new.csv)
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_gtex_truth.py` | A27 part 1 (truths v1-v3), A31b, A50 | answer key: per gene and tissue, cosinor of log(CPM + 0.5) on the time of death of fast deaths (v1), with RIN > 6 and covariates (v2), on the CHIRAL DIP of all deaths (v3); 2-df F p, BH q; the 11 matched tissues plus BA9 and Cortex | `gonogo/scripts/gtex_truth{,_v2,_v3}.py`, `gonogo_new/scripts/gtex_truth_new.py`, `gene_map/scripts/gtex_truth_gm.py` |
| `compute_atlas_gene_rhythms.py` | A27 part 1 (post hoc, no acceptance criterion); A31b (primary, pre-registered) | per atlas (heart nuclei, muscle, breast, colon; HLCA lung, HBCA breast, Xu liver, Ganier skin, Emont adipose): untimed class loadings recombined with RNA shares into a bulk rhythm, scored against GTEx per gene set (32 frozen followers with p < 0.05, GTEx p < 0.01 non-clock-plane genes, clock plane): R = abs(mean exp(i(pred - true))), offset, rise / fall at 13 h, 50 host permutations; A31b primary pooled over the five atlases (follower R above the pooled null 95 % and rise / fall > 0.5 with p < 0.05) | `gonogo/scripts/{sc_predict,compare}.py`, `gonogo_new/scripts/{atlases,predict,score}.py` |
| `compute_antiphase.py` | A30 (pre-registered); replication on A31b; `--source pb` exploratory | share of fasting + feeding genes in antiphase (12 +- 4 h) vs matched random genes (same atlas x tissue block, same GTEx p decile, 2000 draws), clock outputs as control, pooled over 4 atlases | `antiphase/scripts/{antiphase,classes,run,replicate_a31b,pb_predict}.py` |
| `compute_atlas_gene_rhythms_modes.py` | A27 part 1 and A31b grids (covariate modes `none` / `basic` / `pcs`), composition term, sc z > 2 sets | `compute_atlas_gene_rhythms.py` with the A27 covariate fits of the per-class regression (basic: age, sex, log depth, log cells; pcs: + 3 PCs) and `pred_total` (class composition loading added to every gene) | `gonogo/scripts/{sc_predict,compare}.py`, `gonogo_new/scripts/{predict,compare_new}.py` |
| `compute_antiphase_modes.py` | A30 on the covariate modes (+ A31b replication) | the A30 antiphase test on the predictions of each mode | `antiphase/scripts/run.py`, `replicate_a31b.py` |
| `compute_blood_gtex_arms.py` | A32 full design (run note arms) | GTEx whole blood: frame from the clock genes of other tissues, oracle and cross-fitted frames, blood's own clock, ventilated donors, donor bootstrap | `blood_composition/scripts/blood_comp.py` |
| `compute_donor_curve.py` | A39 (pre-registered) | HLCA lung donors subsampled to 25, 45, 100 (10 draws, seed 39) and 200: median R vs the median null 95 % at n = 45 | `gonogo_new/scripts/a39_subsample.py` |
| `compute_gene_map.py` | A50: H50.1-H50.4 (pre-registered) | pooled normal adult lung (HLCA + 21 datasets, 372 donors) and DLPFC (452 donors): R vs GTEx Lung / BA9 / Cortex (50 permutations), pooled vs HLCA alone with a donor bootstrap (`--n-boot`), split-half replication of the gene phases per class (HLCA vs others, PsychAD vs others), the map (class x gene: amplitude, phase, se, p, q) | `gene_map/scripts/{pool,fit,score_gm,boot,h503}.py` |
| `compute_bulk_clock.py` | A29 (pre-registered) | B = lam_host x sum s a exp(i delta) of 12 matched units (19 with sensitivities) vs the timed GTEx tissue gain: Spearman of absolute values (>= 0.5, p < 0.05), phase agreement, loss to misalignment, family shares | `bulk_clock/scripts/{gtex_gain,predict,score}.py` |
| `compute_celltype_truth.py` | A31 (pre-registered; E1 windows post hoc) | mouse colon chrono-atlas, 12 classes: untimed gene phases vs the timed cosinor on ZT, pooled follower R vs its host-permutation null; sensitivity arms C-DC + R-PC; host refitted on 3 or 2 consecutive ZTs | `celltype_truth/scripts/{build,score,window}.py` |
| `compute_blood_clock.py` | A32 (pre-registered), A35 (pre-registered) | marker-score rhythms in GTEx whole blood (timed vs untimed on the host of 14 tissues) and in living blood (GSE39445); blood's own 14-gene clock vs clock time (arms a, b, c, c+b, c+b_w, subject bootstrap), per-gene cosinor raw and composition-adjusted | `blood_composition/scripts/{blood_comp,external}.py`, `blood_clock/scripts/blood_clock.py` |

A27 part 2 (ventilator cases, glucocorticoid axis, D4-D9) is in `gtex_cues`; the A27 ventilator-as-atlas diagnostic D4 is not ported (below).

## Outputs
`gtex_cosinor{,_v2,_v3}.csv.gz`, `gtex_cosinor_summary.csv`; `atlas_{peratlas,pooled,primary,classes,fits}.csv`, `pred_<atlas>.csv.gz` (all genes,
complex loading as `c_re`, `c_im`), `null_<atlas>.csv.gz` (scored genes); `antiphase_{classes,class_sizes,replication,genes}[_pb]`;
`a39_{draws_scores,summary}.csv`; `h50_{scores,sets,boot,boot_summary,h503_lung,h503_dlpfc,h503_*_summary}.csv`,
`map_<fit>.csv.gz`, `map_counts.csv`; `a29_{gtex_gains,units,compartments,families,readouts}.csv`;
`a31_{scores,hosts,genes,window_E1}_<tag>.csv`; `a32_{gtex,external,external_vs_gtex}.csv`, `a35_{arms,per_block,genes,phases}.csv`, `a35_verdict.txt`.
Column names follow the legacy files, except `R` (legacy `R_err`), `offset_h` (legacy `R_rot_h`) in the atlas tables of `h50_scores`.

## Checks
Every script runs with `--subsample` (minutes). Full point estimates against `legacy/desync/validation/*/output/` (all full data; 50
host permutations, 200 for A31; bootstrap of H50.2 and the A35 / A32 bootstraps reduced):

| readout | result here | legacy | note |
|---|---|---|---|
| GTEx truths v1 / v2 / v3, the 7 A27 tissues (79,951 / 80,570 / 80,225 gene x tissue rows) | identical | | max diff cos coefficient 8e-16, p 2e-13 |
| A31b lung_hlca, GTEx p < 0.01 non-clock genes: n, R, offset, null 95 % | 98, 0.560, 1.37 h, 0.331 (p 0.020) | 98, 0.560, 1.37 h, 0.272 (p 0.020) | donors 200, 18 classes, 28 strata equal |
| A31b breast_hbca / liver_xu / skin_ganier (sun-exposed) / adipose: R | 0.223 / 0.051 / 0.029 / 0.169 | 0.223 / 0.045 / 0.025 / 0.170 | n 94 / 155 / 415 / 418 vs 94 / 156 / 418 / 421 |
| A27 heart LV / AA, muscle, breast, colon transverse / sigmoid: R | 0.252 / 0.169 / 0.076 / 0.408 / 0.226 / 0.026 | 0.240 / 0.174 / 0.074 / 0.408 / 0.220 / 0.024 | offsets within 0.3 h except sigmoid; donors 45 / 32 / 126 / 36 equal |
| A31b primary, pooled followers (n 33): R, null 95 %, p; rise / fall | 0.436, 0.430, 0.059; 0.455 (p 0.75) | 0.436, 0.376, 0.039; 0.455 (p 0.69) | fails on rise / fall in both; null quantiles differ by Monte Carlo (50 permutations, 33 genes) |
| A30 pooled, fasting + feeding (n 37): frac anti, null 95 %, p | 0.216, 0.432, 0.956 | 0.216, 0.432, 0.949 | all rhythmic genes 1398 vs 1409 |
| A30 clock outputs control (n 15); replication on A31b (34 genes) | 0.267, p 0.78; 0.500 vs null 0.382, p 0.12 | 0.267, p 0.80; 0.50 vs 0.38, p 0.12 | |
| A39 n = 45: median R, median null 95 %, draws passing | 0.199, 0.297, 3 / 10 | 0.20, 0.32, 2 / 10 | primary fails in both; n = 100: 0.282, 5 / 10 (0.28, 5 / 10); n = 200: 0.560; n = 25 infeasible |
| A50 H50.1: lung without HLCA, R (null 95 %, p) | 0.425 (0.299, 0.020) | 0.419 (0.343, 0.020) | pooled 0.462 vs 0.447; HLCA alone 0.560 vs 0.560; non-disease controls 0.285 / 0.486 vs 0.301 / 0.468 |
| A50 H50.4: DLPFC pooled vs BA9, R (p) | 0.538 (0.039) | 0.538 (0.020) | PsychAD 0.551 vs 0.553; other collections 0.482 vs 0.482; Cortex 0.131 vs 0.129 |
| A50 H50.3: classes replicating (lung; DLPFC) | 15 / 17; 7 / 8 | 15 / 17; 7 / 8 | all pass; same donors per set (372, 172, 452) |
| A29 primary Spearman(abs B, abs g), 12 units | 0.119 (p 0.36) | 0.168 (p 0.30) | fails in both; donors per unit identical, abs g within 3 %, lam_host 0-12 % lower; phase R 0.877 vs 0.898 |
| A31 pooled followers (n 33): R, offset, A; all rhythmic (n 2,673) | 0.982, 0.2 h, 0.98; 0.987 | 0.98, 0.2 h, 0.98; 0.99 | null 95 % 0.64 vs 0.82 (deviation: within-arm null); p 0.005 in both; 12 classes |
| A31 E1 windows, median R followers: 3 ZTs, 2 ZTs | 0.936, 0.259 | 0.937, 0.271 | |
| A35 arms b, c, c+b (rho with clock time) and (b)-(a) R, rho | 0.220, 0.309, 0.331; +0.034, +0.137 | identical | bootstrap 20 replicates in the check, intervals not compared |
| A32 external (GSE39445): L-N amplitude, peak, p; neutrophil peak | 0.454, 2.667 h, 2e-22; 15.08 h | identical | peak interval [1.44, 3.77] identical (500 bootstraps) |
| A32 GTEx whole blood (212 donors): L-N timed p, untimed p and peak | 0.459, 0.001, 8.9 h | 0.46, 0.001, 9.1 h | neutrophil timed p 0.085 in both |

| A27 / A31b covariate modes: R_err of the 42 tissue x set x prediction rows (4 A27 atlases, truth v1), within 0.02 of the legacy | none 32, basic 29, pcs 21 of 42 | `gonogo_<atlas>[_basic\|_pcs].csv` | max difference 0.086 / 0.14 / 0.23: the gene panel (below) and, for pcs, the PCs of the changed panel; breast basic and none identical to the digit (R 0.392, offset 23.39 h, ccorr -0.649, within 3 h 0.413); composition loadings (amp, acrophase) of breast none / basic identical (1e-15), pcs 0.002 / 0.17 h, heart basic 0.007 / 0.81 h |
| A30 on the modes, 1,188 class x scope x truth x mode rows | fasting + feeding, pooled: p 0.40-1.00 in none, basic, pcs | p 0.42-1.00 (`classes.csv`) | median abs difference of p 0.0015; gene counts n differ in 216 rows (panel); no excess in any mode |
| A32 GTEx whole blood arms (primary L-N): timed p, peaks (timed, untimed), d [CI] | 0.46, 8.68 h, 9.12 h, 0.44 h [-5.52, 6.59] | 0.459, 8.70, 9.13, 0.43 [-5.62, 6.51] | FAIL (uninformative: no timed rhythm), as legacy; oracle, cross-fit, own-clock and ventilated arms in `a32_arms.csv` (16 rows, same columns as the legacy `results.csv`) |

**M3 full rerun (2026-10-25, `run.sh --jobs 4`, 54 min, peak 13.3 GB).** Against the archive copies: 33 tables identical or
within 1e-9, A31 sensitivity (`sens_RPC`) and the H50.2 bootstrap new. The H50.2 donor bootstrap (200 replicates):
pooled lung R 0.462, bootstrap median 0.465 [0.310, 0.564], 3.0 % of replicates reach HLCA alone (0.560): H50.2 fails, as
legacy. The atlas tables (`atlas_{classes,peratlas,pooled,primary}`; `modes_peratlas` 168 of 1,080 rows, `modes_composition` 3 of 213) moved in the skin_ganier rows (plus two muscle v3 covariate-mode rows by < 0.007 h, not isolated): A31b
primary pooled followers R 0.435, null 95 % 0.430, p 0.078 (archive copy 0.436, 0.430, 0.059; fails on rise / fall as
before); 3 genes change rhythmic class in `atlas_classes`. Cause (checked, full `compute_atlas_gene_rhythms.py` with the
clocksync before the EM speed-up 09d0142 reproduces the archive copies to 2e-15): the leave-one-anchor-out host of
skin endothelial (`gene_rhythms` refits the host without each anchor, `max_iter=3000`) does not converge in 3,000
iterations (also breast_hbca endothelial); a non-converged EM path depends on rounding, so the two codes stop at
different points (loading differs by 0.19, log-likelihood -167.350 vs -167.339). Every converged host agrees to 3e-8.
Not isolated: the A39 rows n = 45 / 25 (median R 0.199 / 0.237, third decimal) and `a31_window_E1` (0.006); the
A39 verdict is unchanged. Open: `gene_rhythms` should report non-converged held-out hosts (API gap below).

## Deviations
- **Gene panel.** The legacy tables kept the genes expressed in the 7 A27 tissues; here the genome tables are cut to the genes expressed (median CPM >= 5) in any matched
  GTEx tissue of the three truths (a few hundred more genes), which changes the gene counts by 1-3 per set.
- **Atlas builders.** The genome-wide tables are rebuilt from the Census group pseudobulks (`data/*_genome`, `gene_map_*`); donors,
  classes and strata equal the legacy ones. The lung rule of `lung_hlca_genome` (A4 classes, stratum assay | study) is that of A31b.
  The builder `data/_atlas/census.py` raises a pandas 3 read-only-array error; `_shared._patch_census` replaces `sum_groups_sparse`
  at run time with a one-line corrected copy (needs a fix in `data/_atlas/census.py`, line 248: `.to_numpy().copy()`).
- **Null.** `gene_rhythms` permutes host states within stratum. The legacy A27 / A31b / A39 / A50 nulls were the same; the A31 null permuted animals across arms (here within arm).
- **Covariates and composition.** `compute_atlas_gene_rhythms.py` runs `mode none`; the `basic` / `pcs` fits and the composition term
  are done by the local covariate regression of `compute_atlas_gene_rhythms_modes.py` (`covariate_table`, `predict_all`), because `gene_rhythms` has no covariates.
  Null 50 permutations; the legacy GTEx split-half ceiling (`compare.py`) is not computed.
- **A30 `--source pb`** (exploratory): the atlas pseudobulk without the age, sex and log-depth covariates of the legacy D3 regression; min_mean 10.
- **A29.** lam_host is the standard deviation of the host state of `Dataset.fit` (3-12 % below the legacy trace of Sigma_f); the learned A17
  anchors come from the frozen `common/anchors_a17.csv`; the GTEx gain uses one sample per donor (`common/gtex.load`).
- **A32 untimed GTEx arm.** The frame shift from the host to clock time is read from a `gene_rhythms` call (PER1 acrophase against its raw
  loading), not from the waveform fit of the 14 other tissues' clock genes (legacy `frame_rot`); no bootstrap, oracle frame, cross-fit,
  blood own-clock or ventilator arms.
- **A50 H50.2 bootstrap.** Donors resampled within stratum with the legacy seed stream ([50, b]); the replicates refit the host and genes.
- **Exact `p` of the A31b primary** differ from the legacy by Monte Carlo (50 permutations).

## Not ported
- A27 D1-D4 diagnostics (the untimed step on GTEx bulk itself, degraded host, atlas pseudobulk, ventilator cases as an atlas); D3 is
  covered by `compute_antiphase.py --source pb`.
- GTEx split-half ceiling of `compare.py` (a descriptive reliability bound).
- A50: feasibility / selection scripts (`feasibility.py`, `select_datasets.py`), per-class composition tables, the Hollis / Anafi and post-mortem comparisons (not run in the legacy).
- A32, A35: nothing (A32 arms in `compute_blood_gtex_arms.py`).
- A29: the learned-rule re-fit per unit (the frozen anchors are used).

## Local helpers (not in clocksync), `_shared.py` unless noted
- `gtex_truth`: per-gene timed cosinor with covariates, 2-df F and BH (genome-wide; `timed.cosinor` has no test and one series).
- `shares`, `predict`: RNA-share recombination of class loadings into a bulk rhythm; `score`, `score_vs_null`, `st`: R, offset, rise / fall, null.
- `fit_atlas`: glue of `Dataset`, `gene_rhythms` and the recombination; `genome_table`, `class_rule`, `a31b_anchors`, `_patch_census`.
- `compute_antiphase.py`: `score_atlas`, `matched_null_test`; `compute_bulk_clock.py`: timed tissue gain, shares, family shares;
  `compute_celltype_truth.py`: timed cosinor with arm intercepts, rotation by the clock genes; `compute_blood_clock.py`: marker scores, block cosinor, arms.

- `compute_atlas_gene_rhythms_modes.py`: `covariate_table`, `predict_all` (covariate modes, composition, single-cell z); `compute_blood_gtex_arms.py`: `frame_rot` and the cross-fit.

## API gaps
- `Dataset.gene_rhythms` does not report the convergence of its leave-one-anchor-out hosts; a class read on a
  non-converged host is rounding-dependent (M3, skin_ganier endothelial).
- `Dataset.gene_rhythms(covariates=..., composition=True)`: the covariate (age, sex, log depth, log cells, PCs) and composition fits of the A27 grid.
- `Dataset.gene_rhythms` / `Fit`: expose the frame shift (host frame to the mean-anchor frame) and a way to fit an arbitrary score (not a count gene) on the host, and a null with permutation across strata.
- `Dataset.timed_loadings(hours)` and a known-host `on_host` (coupling to the known time, the A21 / A36 truths).
- `Dataset.contrast`: expose the bootstrap draws (R21.3 family intervals).
