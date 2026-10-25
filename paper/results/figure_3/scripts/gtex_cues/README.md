# gtex_cues: GTEx covariates, the cue reference, internal vs external time

GTEx v10 tissues as compartments of one donor, in fast deaths, ventilated donors (Hardy 0) and the slower death classes: the
dose-response and covariate tests of the host-locked clock (ventilation days, BMI, T2D), what keeps organs coupled in ventilated
donors, the GTEx **cue reference** (genes whose rhythm follows the body's timing cues rather than the clock), its hierarchy over
organs, clock-following vs time-following genes, drug targets, the CHIRAL age comparison and covariate sensitivity, and the
glucocorticoid gene classes. Hosts and couplings come from `clocksync` (`Dataset.fit(...).host`, `covariate`, `contrast`, `timed`,
`paired.projector / states / ccorr`); genome-wide per-gene regressions on the host are local helpers. No plotting code.
Outputs: `paper/results/figure_3/output/gtex_cues/` (`.../subsample/` for checks). GTEx time of death and protected
phenotypes are controlled access: aggregates only.

```
run.sh [--subsample N]      all scripts in order (full run by default); compute_cue_reference.py first
_shared.py                  tissue lists, per-tissue CPM, D2 covariates, partial F, host states (Dataset.fit(...).host)
_genes.py                   genome-wide cue pipeline: selection table (26 tissues), lock, sets, frame check   (cue family)
_covariates.py              cohorts, protected phenotypes in memory, module and permutation helpers          (covariate family)
```
Export `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` for joblib runs (`_genes.py` sets them; `limit_threads` does not reach loky workers).

## The cue reference (output for psoriasis_details and cancer)
`paper/results/figure_3/output/gtex_cues/a81b_cue_reference.csv` (columns gene, set = cue / locked, nt, fisher, Rw, lock, amp,
acro_h in the GTEx time-of-death frame) with `a81b_cue_reference.README.txt` (columns, frame, how to build a `clocksync.Waveform`
from `amp` and `acro_h`). These are the A81b sets (40 cue and 40 locked genes; the sets are identical to the legacy `a81b_sets.csv`) used by the
A81 H81.2 / H81.3 psoriasis arms, A82, A84 and the A88 TCGA arms. The frame check of A81b fails its pre-specified guard (0.481 < 0.8, shift -1.08 h), as in the legacy,
and the override (post hoc, labelled in A81b) applies. The 26-tissue selection table `selection_genes_26.csv.gz` also feeds the tissue-matched sets of A84.

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_cue_reference.py` | A81 (H81.1, E1); A81b (frame check, E1; H81.2/3 are scored in figure_S09 / psoriasis with the reference file) | per gene and tissue: timed cosinor in fast deaths and lock = host loading in ventilated donors / timed loading, relative to the tissue's clock-plane genes; cue (low lock) and locked (high lock) sets; core-vs-set offsets and correlations per tissue (fast vs ventilated); H81.1 drop of the core-cue correlation in ventilated donors minus the core-locked drop, Hardy labels permuted within tissue (2,000) | `paired_damping/scripts/a81_gtex_cues.py`, `a81b_cues.py` (selection part of `a84_matched_cues.py` for the added tissues) |
| `compute_cue_hierarchy.py` | A89 (H89); A90 | body-wide / group- / organ-specific cue-like genes over 26 tissues, within-group phase coherence vs 1,000 random tissue sets; p53 targets locked or cue-like (lock intervals, ventilated donors resampled 200) | `a89_cue_hierarchy.py`, `a90_p53_lock.py` |
| `compute_internal_external.py` | A91 (H91); A27 D5 (post hoc) | per tissue and gene: clock-following / time-following / both (partial F of host vs cos/sin of death time), noise floors, H91 lock difference, age ratio; D5: the same in ventilated donors, host-vs-time alignment | `a91_internal_external.py`, `glucocorticoid/scripts/diag5_int_ext.py` |
| `compute_drug_targets.py` | A100 (H100) | share of cue-like rhythmic pairs, drug targets minus other genes (gene labels permuted, 2,000); wall-time retention in ventilated / Hardy 3-4 donors | `a100_drug_targets.py` |
| `compute_chiral_age.py` | A101 (H101.1, H101.2); A101b (H101b) | age dampening of output genes vs core on the host, pooled / Hardy-matched / fast; CHIRAL's rhythmic-gene counts on CHIRAL's phases, pooled vs Hardy-matched | `a101_chiral_age.py`, `a101b_chiral_measure.py` |
| `compute_covariate_sensitivity.py` | A102 (A102.1-A102.6); A102b (A91, A101, A100 on CHIRAL preprocessing) | leakage of covariates into the host, A91 with extended covariates, interactions with ischaemic time and RIN, debiased age amplitudes, Hardy x RIN-matched counts, debiased retention; CHIRAL-preprocessed reruns | `a102_covariates.py`, `a102b_chiral_prep.py` |
| `compute_ventilation_dose.py` | A46: H46.1, H46.2 (Holm over both); variants cap7, log2, excl_steroid_sepsis | does the host-locked amplitude of {pancreas, muscle, liver} minus {tibial nerve, skin x2} fall with days of ventilation (Hardy 0)? `Dataset.covariate` per tissue on the leave-one-out host, bias-corrected `locked_amp_c`; and does the host alignment with the time of death fall with days (cos slope, permutation within cause of death)? tertile secondary | `gtex_covariates/scripts/a46_ventilation.py`, `common.py` |
| `compute_metabolic.py` | A47: H47.1; T2D and liver disease secondary | BMI (per SD) trend of the host-locked amplitude of {liver, adipose x2, pancreas, muscle} minus {nerve, skin x2}, Hardy 0 and Hardy 1-2 fitted separately, pooled by inverse null variance; T2D coupling loss (control minus case, d at own spread and d* at common spread) per tissue and module with the companion (`Dataset.contrast`) | `a47_metabolic.py` |
| `compute_coupling_sources.py` | A92: H92.1, H92.2 (+ GC adjustment, wall time); A93: H93.1, H93.2, H93.3 | inter-organ coherence C and wall-time alignment R of per-tissue clock states (`paired.projector/states/ccorr`): by ventilation tertile, brain-lesion proxy, GC adjustment; fast / H3 / H4 / ventilated at equal n; season of death vs host phase | `a92_coupling_sources.py`, `a93_disturbance.py` (+ `gtex_states.py`) |
| `compute_glucocorticoid.py` | A27 part 2 (post hoc): D6, D6b, D6c | a priori gene classes among fast-death rhythmic genes: retention, phase kept, amplitude, time following in ventilated donors vs matched random sets (host arm); GC score vs Hardy class, host and time of death; GC effect on the 14 clock genes (in-plane fraction, pull phase) | `glucocorticoid/scripts/diag6_classes.py`, `diag6b_gc.py`, `diag6c_gc_zeitgeber.py` (+ the host arm of `diag5_int_ext.py` for the gene table) |
| `compute_gonogo_diagnostics.py` | A27 part 1 diagnostics D1, D2, D4 (post hoc, exploratory) | the untimed per-gene step on GTEx bulk (every gene on the host of the other 14 tissues, no times) scored against the timed cosinor: D1 same donors, all day vs the 08-18 h window; D2 degraded (fewer donors, 3-tissue host, window, disjoint truth, timed oracle); D4 ventilated donors as the "atlas" scored against fast-death truth, with a fast-death control of the same n | `diag1_gtex.py`, `diag2_degrade.py`, `diag4_ventilator.py` |
| `compute_dip_arm.py` | A27 part 2: D5b, D6 (DIP arm), D6b (DIP), D6c part 1 (post hoc) | CHIRAL's donor internal phase replaces the host in ventilated donors: genes following it vs time of death, host-DIP agreement, a priori gene classes (retention, phase kept, amplitude ratio, comparable on the DIP), GC score on host / DIP / time of death, GC profile along the DIP | `diag5_dip.py`, `diag6_classes.py` (DIP arm), `diag6b_gc.py`, `diag6c_gc_zeitgeber.py` (part 1) |
| `compute_gc_atlas.py` | A27 part 2: D7 (post hoc) | GC score in the heart, breast, muscle and colon atlases: R^2 on the host, effect on the 14 clock genes (in-plane fraction, pull phase), the atlas-pseudobulk regression with the GC score as covariate vs GTEx | `diag7_gc_atlas.py` (+ `diag3_atlas_pb.py`) |
| `compute_gc_li_breast.py` | A27 part 2: D8 part B (post hoc) | GC score of Li 2024 breast (tumour, normal) against recorded resection time, V25.5 internal time and its error | `diag8_gc_cancer.py` (part B) |
| `compute_chiral_prep_counts.py` | A102b item 3 (A101b on CHIRAL's preprocessing, NaN rule) | rhythmic-gene counts (q < 0.2, log2 peak-to-trough > 1) on CHIRAL's phases with entries > 5 log2 from the gene mean excluded: R1, pooled, Hardy-matched, Hardy x RIN-matched | `a102b_chiral_prep.py` (`step_a101b`) |

## Outputs
**Cue reference (for figure_S09 and cancer, A81/A82/A84/A88 TCGA arms): `paper/results/figure_3/output/gtex_cues/a81b_cue_reference.csv`**
(+ `a81b_cue_reference.README.txt`: columns, frame, how to build the `Waveform`). Also `a81b_selection_genes.csv.gz`
(11 tissues, legacy columns tissue, gene, p_fast, amp_fast, acro_h, lock), `selection_genes_26.csv.gz` (26 tissues; the A84 and A89
tables of the legacy are subsets), `a81_sets.csv`, `a81_e1.csv`, `a81_h1_verdict.csv`, `a81b_sets.csv`, `a81b_e1.csv`,
`a81b_frame_check.csv`, `a89_*.csv`, `a90_*.csv`, `a91_*.csv(.gz)`, `diag5_{summary,host_vs_time}.csv`, `diag5_genes.csv.gz`,
`a100_*`, `a101_*`, `a101b_*`, `a102_*`, `a102b_*` (names as the legacy files). Second pass: `diag1_{summary,genes}`, `diag2_reps`, `diag4_reps`, `diag5_dip_{agreement,summary,overlap,genes}`, `diag6_classes_dip`, `diag6b_gc_dip`, `diag6c_gc_dip_profile`, `diag7_{<atlas>,<atlas>_gc,summary}`, `diag8_li_gc`, `a102b_a101b_{tissues,summary}`.

`a46_{tissues,modules,h462,primary,tertiles,permutation,fast_context}.csv`; `a47_{design_counts,bmi_tissues,bmi_cohort_modules,bmi_pooled,t2d_tissues,t2d_cohort_modules,t2d_pooled,t2d_companion_{vent,fast},primary}.csv`; `a92_{summary,walltime}.csv`; `a93_{groups,summary}.csv`; `diag6_classes.csv`, `diag6_tissues.csv`, `diag6b_gc.csv`, `diag6c_gc_summary.csv`, `diag6c_gc_clock_effects.csv`. Column names follow the legacy files except: A46/A47 tissue tables have `est_c` (bias-corrected) and `p_two` (package p, two-sided around the null median); A46 modules carry `perm_p_lower`; A47 T2D tables report the loss (control minus case) as `loss`, `loss_star`. Aggregates only.

## Checks
Every script runs with `--subsample`. Full point estimates against `legacy/desync/validation/{paired_damping,gtex_covariates,glucocorticoid}/output/`:

| readout | result here | legacy | note |
|---|---|---|---|
| A81 sets (selection halves), H81.1 | 3 cue + 3 locked genes; D -0.0566, p 0.975, FAIL | identical (D -0.0566, p 0.975) | 2,000 permutations; E1 table identical |
| A81b frame check | 14 genes, circular correlation 0.481, shift -1.08 h | identical | guard (> 0.8) fails, as legacy |
| A81b selection table (126,440 gene x tissue rows), sets | identical (max difference 1e-14); cue 40 / locked 40, same genes; median lock 0.527 / 0.992 | identical | host (EM) identical to machine precision |
| A89 levels (1,960 genes); H89 | body-wide 6, group 58, organ-specific 1,394; within-group R 0.915 vs random 0.870, p 0.269, FAIL | identical counts; 0.92 vs 0.87, p 0.27 | |
| A90 p53 (918 gene x tissue rows, 17 genes) | lock and 95 % limits identical; no target cue-like, CDKN1A median lock 1.10 | identical | whole blood is in the A90 tissue list (legacy quirk) |
| A91 (26 tissues): class counts, noise floors, H91 | identical; 12 tissues, median difference +1.232, p 0.0005, PASS | identical | 312,018 gene calls identical |
| A27 D5 (14 tissues) | counts, f-vs-t R, phase of internal-only genes identical; host vs time R 0.63 (fast) / 0.14 (ventilated) | identical | permutation-null medians differ (the legacy seeded with `hash(tissue)`) |
| A100 | 5,433 rhythmic pairs, 186 target pairs (142 genes), cue share 0.285 vs 0.346, diff -0.0615, p 0.087, FAIL | identical | |
| A101 point estimates | D_pooled 0.1301, D_core 0.0803, D_rel 0.0693, D_matched 0.1426, D_fast 0.2277 | identical | `--boot 0` (no bootstrap), so H101.1/2 not re-evaluated |
| A101b point estimates | 37 tissues; R1 counts old 6,200 / young 13,792; pooled log ratio -0.702, matched -0.200, matched - pooled +0.502 | identical | `--boot 0` |
| A102.1 leakage, A102.5, A102.2 | median excess R^2 0.0326 (not material); clock excess +42.7 %, H91 p 0.0005 (not robust); interaction medians identical | identical | |
| A102.3, A102.4, A102.6, A102b A91 / A101 / A100 | run on `--subsample 150` only (A102b A91: clock excess +42.5 % in the check, as legacy +42 %) | not compared | not rerun in full, see below |
| A27 D1 (full): R, median error (h) of the non-clock genes with timed p < 0.01, all donors / window, six tissues | e.g. heart 0.950, 0.76 h / 0.722, 2.33 h; thyroid window 0.122, 5.03 h | heart 0.950, 0.76 / 0.709, 2.37; thyroid window 0.128, 5.12 | donors and rhythmic gene counts identical in all 12 rows; R within 0.013 (window host fit, other EM start) |
| A27 D2 (5 replicates, six configurations), mean R and oracle R | 0.53 / 0.51 (reference), 0.56 / 0.38 (atlas-like) | 0.56 / 0.49, 0.57 / 0.40 | same pattern (untimed at least as good as the timed oracle under atlas-like conditions); other random streams (legacy used a process-dependent `hash`) |
| A102b A101b (NaN rule): log(total old / young) R1, pooled, Hardy-matched, Hardy x RIN-matched; tissues | -0.795, -0.727, -0.128, -0.023; 27, 26, 25, 23 | -0.795, -0.687, +0.210, -0.306; 27, 28, 23, 21 | R1 identical (same donors, phases, NaN rule); the 10 equal-n draws differ (seeds follow the tissue index of a different tissue list), so pooled / matched totals differ by 0.04 / 0.34; median tissue log ratios +0.26, +0.19 vs +0.31, +0.03: the matched counts are noisy, matched minus pooled stays positive (+0.60 here, +0.90 legacy) |
| A27 D4 (five replicates per tissue, six tissues): untimed R / offset (h) of the non-clock genes, ventilator vs fast-death control; clock-plane genes R | 0.34 / 23.1 h vs 0.48 / 0.6 h; frame 0.74 vs 0.86; timed 0.24 vs 0.46 | identical to 2 decimals (same seeds) | |
| A27 D5b: host vs DIP R (fast, ventilated), host vs time of death | 0.680, 0.577; 0.631, 0.142 | identical | |
| A27 D5b per tissue: genes following the DIP, internal-only, phase R | 1,440 / 1,137 / 0.350 (muscle), 2,433 / 2,057 / 0.453 (thyroid); all 14 tissues | identical | counts null-calibrated with 50 permutations |
| A27 D6 DIP arm: retention (clock outputs, SREBP, all rhythmic), p | 0.571 (0.001), 0.286 (0.025), 0.135 | identical (0.001, 0.023, 0.135) | all eight classes identical to 3 decimals; p within Monte Carlo error (2,000 draws) |
| A27 D6b DIP: R^2 of GC on host / DIP / time of death (ventilated); host-DIP R after removing GC; by GC tertile | 0.225 / 0.206 / 0.010; 0.529; 0.57-0.60 (time of death 0.13-0.19) | RESULTS: 23 % / 21 % / 1 %; 0.53; 0.56-0.60 (0.13-0.19) | |
| A27 D7 (heart, breast, muscle, colon): GC R^2 on host, pull phase, non-clock R with GC covariate | heart 0.029, 15.5 h, LV 0.22 (0.29 before); breast, muscle, colon likewise | identical GC summaries; R 0.22 vs 0.23 (heart LV), 0.49 / 0.49, 0.12 / 0.12, 0.25 / 0.24 | non-clock R within 0.03, offsets within 0.5 h (GTEx truth keeps one sample per donor; the legacy truth used every library, 844 vs 853 genes) |
| A27 D8 part B (Li breast, 42 patients): GC vs recorded time, internal time, error | -0.006, +0.128, +0.160 (p 0.31) | RESULTS: -0.01, +0.13, +0.16 (0.31) | |
| A46 H46.1 module D (per day) | +0.0167 | +0.0168 | identical to 3 decimals; hosts refitted by clocksync |
| A46 H46.1 p | 0.81 (within-cause permutation, 100) | 0.70 (independent simulated nulls); permutation 0.80 | different null, see deviations; same verdict |
| A46 H46.2 slope, p, fast-death R | -0.0012, 0.4895, 0.6305 | -0.0012, 0.4895, 0.6305 | identical; Holm 0.979 both |
| A46 tissue estimates (Pancreas, muscle, nerve) | 0.0135, -0.0364, -0.0102 | 0.0135, -0.0364, -0.0102 | liver 0.129 vs 0.126 (host EM) |
| A46 variants cap7, log2, excl. steroids/sepsis (D) | 0.014, 0.045, 0.047 | 0.023, 0.047, 0.043 | descriptive; none significant |
| A46 tertile host spread | 0.41, 0.31, 0.32 | 0.90, 0.68, 0.71 | different host scale (clocksync host units); mean cos 0.153, 0.140, 0.135 identical |
| A47 H47.1 pooled BMI module, p | +0.0253, 0.79 | +0.0207, 0.72 | verdict FAIL both; null here is a BMI permutation |
| A47 T2D module d* minus companion (pooled) | -0.0615 (p 0.88) | -0.0470 (p 0.87) | loss sign; own-spread d -0.0542 vs -0.0481 |
| A92 C_fast, C_vent, tertiles | 0.599, 0.451, 0.482 / 0.422 / 0.452 | identical | same random stream |
| A92 H92.1, H92.2 | +0.029 (p 0.296), -0.032 (p 0.831) | identical | both FAIL; GC drop 0.052, wall-time R 0.133 vs 0.528 identical |
| A93 R fast, H3, H4, vent | 0.526, 0.225, 0.096, 0.150 | identical | H93.1 diff -0.054 (p 0.905), H93.2 -0.050 (p 0.231), position index -0.14 identical |
| A93 H93.3 summer - winter | +0.137 h (p 0.838) | identical | |
| D6 clock outputs: n, retained, phase kept, amp_rel | 35, 0.657, 0.571, 0.882 (p 0.001) | 35, 0.657, 0.571, 0.882 | host arm identical; heat shock 9, 0.556, 0.444; all rhythmic 4,568, 0.173 |
| D6c in-plane fraction (fast, vent), pull h | 0.074, 0.382; 15.18, 17.69 h | 0.073, 0.381; 15.08, 17.67 h | null 0.766 / 0.639 vs 0.809 / 0.590 (200 random scores) |

## M3 full run (2026-10-26)
`run.sh --jobs 4`, 1.0 h, peak 3.3 GB. Every table of the archive outputs (`outputs/paper/figure_S07_gtex_covariates_cues`)
is reproduced to 1e-6 (A46, A47, A81, A81b, A90, A93, D1-D7; D4 to 0.006 h). First full runs, against
`legacy/desync/validation/paired_damping/output/`: A101 bootstraps (500) and A101b (200) equal the legacy intervals (H101b lower
bound -0.0117, FAIL, as legacy); `a102_age_summary` (A102.3), `a102_targets_summary`, A102b A101 and A100 summaries identical. A102.4 (counts)
differs as A102b A101b above (equal-n draws seeded by a different tissue list): pooled / matched total log ratio -0.737 /
-0.136 (legacy -0.706 / +0.029), 26 / 23 tissues (legacy 27 / 22). A102b A91 clock excess +41.0 % (legacy +42.5 %).
<!-- src: figure_3/output/gtex_cues/{a101_summary,a101b_summary,a102_*_summary,a102b_*_summary}.csv; legacy desync/validation/paired_damping/output/ -->

## Deviations
- The selection table of the legacy is split over three files (A81b 11 tissues, A84 7 added, A89 8 added); here one 26-tissue table (same code, same halves) and the A81b file is its 11-tissue subset.
- A101: as in the legacy, the 14 core genes are appended to the output genes, so core clock genes that A91 called clock-following count twice (kept for equivalence).
- A102b nan-OLS: genes with the same pattern of excluded entries are solved together (same numbers, faster).
- Hallmark enrichment (A89) uses the gmt present under `$SCCC_DATA/reference` (v2024.1; the legacy used v2026.1 from another directory); not compared.
- `--subsample` picks N fast and N ventilated donors; the 26-tissue table, A91 classes and A100 pairs of the subsample folder are therefore different objects from the full ones.
- **A46/A47 module null.** `Dataset.covariate` returns the bias-corrected estimates and a per-tissue p but not the simulated null draws, so the legacy module p (independent per-tissue null draws combined by Monte Carlo) cannot be rebuilt from the public API. The module p is the within-cause permutation of the ventilation hours (the legacy sensitivity check; 100 permutations, the host refitted in every call) and, for BMI and T2D, a permutation of BMI or of the T2D labels (100); pooling over the cohorts uses the permutation variances as weights. The module estimate is the mean of the package's `locked_amp_c`.
- **A46 design.** The covariate is centred per tissue (own donors), as the legacy run; permutations use one global design (centred over all donors), so observed and permuted statistics are computed identically.
- **A46/A47 bootstrap.** The donor bootstrap refits the host in every replicate (legacy: host fixed), 50 replicates (legacy 200). A47 gives null-based intervals (D +/- 1.96 null sd) instead of a bootstrap.
- **A47 T2D.** d and d* come from `Dataset.contrast` (common spread V0 of the two groups of a tissue), the legacy d* used the median own spread over tissues. Pooled over cohorts as the legacy; liver disease is descriptive counts only (below 5 cases per stratum, not estimable).
- **A92 / A93** run as written (A92 keeps the legacy "bootstrap" that resamples donors as a set, a 63 % subsample).
- **D6** gene table recomputed for the host arm (partial F of host and time of death, fast-death cosinor) with the 12 body tissues; the Hallmark file is the 2026.1 release (`CoPhaser/test/hazel/data`), fallback 2024.1 in `$SCCC_DATA/reference`.

- **A27 D1 / D2 / D4.** The legacy D2 replicate seed used Python's process-dependent `hash`; here a CRC of the configuration name, so replicate draws differ (summaries agree). The host of an arm is fitted on the tissues that have rows among the arm's donors (ventilated donors have no pituitary or brain rows). The timed truth keeps one sample per donor (`_shared.expression`).
- **A27 D5b / D6 DIP arm.** The DIP table is `data/chiral_DIP.csv` of the repository (per-donor, controlled; in memory only); the DIP enters as cos / sin of its phase in the same partial-F design as the host.
- **A27 D7.** The atlas genome tables, class rule and anchors come from `gene_rhythms/_shared.py` (read-only); the truth is one sample per donor.
- **A102b A101b.** Same preparation as A102.4 (`compute_covariate_sensitivity.prep` rules) with the NaN rule of the legacy `counts_nan`; the old / young / matched draws use `compute_chiral_age.draw_b`.

## Local helpers (not in clocksync)
New in the not-ported pass: `compute_gonogo_diagnostics` (`loadings`, `score`, `d1_one`, `d2_one`, `d4_one`; per-gene regressions on the host), `compute_dip_arm` (`read_dip`, DIP as a two-column host in `compute_internal_external.d5_tissue`), `compute_gc_atlas` (genome tables through `gene_rhythms/_shared.py`, read-only), `compute_chiral_prep_counts` (`counts_nan` on `ols_nan`).

Per-gene regressions on the host: `_shared.partial_f`, `timed_cosinor`, `lstsq_fit`; `_genes.select_tissue` (lock), `per_gene`, `cue_sets`, `frame_check`;
`compute_covariate_sensitivity` (`add_pcs`, `ols_nan`, debiased amplitudes); CHIRAL harmonic regression and BH (`compute_chiral_age`);
genome-wide CPM loader `_genes.load_rows` (gtex.load has no RIN / centre / all-gene CPM per tissue).

`_covariates.py`: `cause_group`, `protected`, `module`, `perm_within`, `p_lower/p_upper`; A46 `make_x`, `boot_D` (own donor bootstrap: `Dataset.resample` replicates cannot be passed to `covariate`), `cos_stat`/`slope`; A93 `pair_corr` (legacy Poisson-corrected |complex correlation| of the A93 form); D6 `gene_table`, `class_tests`, `d6b`, `d6c`.

## Not ported / not rerun in full
- A89 Hallmark enrichment compared with legacy; the A84 matched cue sets and their disease scoring (figure_S09 / cancer).
- Not rerun in full (written, run on `--subsample 150`): A102.3 (age, debiased), A102.4 (counts), A102.6 (retention), A102b A91 / A101 / A100 (A102b A101b was run in full), and the bootstraps of A101 / A101b (500 / 200).
- A27 D3 (atlas pseudobulk against the cell-type recombination, timed oracle at the atlas n) is `gene_rhythms/compute_antiphase.py --source pb` (no timed oracle there); the D7 rerun of it with the GC covariate is `compute_gc_atlas.py`.
- A27 D8 part A (the 21 single-cell cancer cohorts: GC against the A25 per-patient internal time, effect on PER1 / ARNTL / NR1D1 / DBP) and D9 (A25 per-patient times refitted without PER1, cohorts and V25.5) need the per-patient A25 times and the cohort panel pseudobulks of `tme_time.py` / `panel_pb.py`; windows has the cohort-level A25 estimator and the V25.5 breast fit only (no per-patient time for the 21 cohorts, no 13-gene waveform option): not ported (D8 part B is `compute_gc_li_breast.py`).
- `feasibility.py` (design power), A93.3 spring/fall descriptive beyond the season table, A46 simulated-null module p, the D6c in-plane fraction on atlas data beyond D7.
- Ported in the second pass: A27 D1, D2, D4 (`compute_gonogo_diagnostics.py`), the CHIRAL-DIP arms of D5 / D6 / D6b / D6c part 1 (`compute_dip_arm.py`), D7 (`compute_gc_atlas.py`), D8 part B, A102b A101b (`compute_chiral_prep_counts.py`).

## API gaps
`Dataset.covariate` should return the simulated null draws (or a module/contrast combination with companion) so the A46/A47 module test can be written without permutations; `Dataset.covariate` cannot be called on `Dataset.resample` replicates (the covariate index is matched against `<id>#j` labels while `on_host` looks up the base id); a covariate with several host-adjusted columns would cover the GC-adjusted couplings of A33.
Genome-wide per-gene regression on the host (partial F, timed cosinor, lock) and a genome-wide CPM loader with RIN / centre would remove most local helpers of the cue family.
