# Inventory of the cue line: what the project knows about what internal time controls (2026-11-15)

Roadmap item T0 (a). The existing analyses that ask which daily processes follow the body's internal time (host state),
the tissue's own clock, behaviour or the external day. No new statistic, apart from the A91 rerun (item A5). Numbers
were re-read in the tables named (paths relative to `paper/`; `gtex_cues/` = `results/figure_3/output/gtex_cues/`).
Archive registrations: `scCircadianCancer/docs/PLAN_internal_desync.md`; archive results: `.../RESULTS_preliminary.md`.
"Current code": migrated to `results/figure_3/scripts/gtex_cues/` and rerun on `clocksync` (M3, 2026-10-25/26), unless
stated otherwise.

## 1. One row per analysis

Times: **B** body (host state from the donor's other tissues), **L** the tissue's own clock, **X** external (recorded
hour), **H** behaviour. Companions: **M** matched pseudo-clock sets, **R** stress-responsive sets.

| analysis | question | design | verdict and number (table) | companions | code | contributes to the partition | why it failed |
|---|---|---|---|---|---|---|---|
| **A91** (reg.) | per gene: B or X? | GTEx fast deaths, 26 tissues; partial F of host vs cos/sin of hour of death; noise floors by permutation | **PASS**: H91 median lock difference 1.232, p 0.0005, 12 tissues. Muscle 1,480 B / 113 X / 676 both (floors 37 / 31.5 / 8.5); liver 2,717 / 5 / 326 (floors 26 / 19.5 / 10) (`gtex_cues/a91_verdict.csv`, `a91_tissues.csv`) | none (permutation floors only) | **rerun 2026-11-15 on current code: identical** (312,018 gene calls, lock to 1e-9) | the B-vs-X axis, two-way; L and H not separated ("internal time" includes the person's behaviour, as RP states) | - |
| A102 / A102b (reg.) | is A91 robust to covariates / CHIRAL preprocessing? | A91 with RIN, centre, 10 PCs; CHIRAL pipeline | robustness rule FAIL in the direction of more B: clock excess +42.7 % (`gtex_cues/a102_fast_summary.csv`); H91 holds (p 0.0005). X class fragile: 1,250.5 -> 380 excess with covariates; robust only in non-sun-exposed skin | none | current | B class robust; X class is mostly preprocessing-dependent | rule set on stability, B grew |
| A27 D5 (post hoc) | B vs X in ventilated donors | GTEx Hardy 0, 14 tissues | host vs hour of death R 0.142 (ventilated) vs 0.631 (fast) (`gtex_cues/diag5_host_vs_time.csv`); X-following genes at null in most tissues (3 of 12 above the migrated null, 2 in the archive: permutation seed) | none | current | the body keeps one state after losing the day (B without X) | - |
| A46 (reg.) | does alignment with the day decay with days of ventilation? | GTEx Hardy 0, ventilation hours (protected) | FAIL: H46.1 +0.0199, p 0.81; H46.2 slope -0.0012, p 0.49 (`gtex_cues/a46_primary.csv`) | nerve + skin module | current (2026-11-01; README check still quotes the legacy +0.0167) | loss of X alignment is a step, not a drift | no dose-response |
| A92 / A93 (reg.) | what keeps organs coupled without the day? | GTEx ventilated, Hardy 3 / 4 | H92.1 p 0.296, H92.2 p 0.831, H93.1 p 0.905 FAIL; coherence 0.599 fast vs 0.451 ventilated; wall-time R 0.528 vs 0.133 (`gtex_cues/a92_summary.csv`, `a93_summary.csv`); season +0.137 h, p 0.838 | none | current | B persists without X; behavioural cue loss (Hardy 4) suffices to lose X | no decay; season underpowered (SE about 1 h) |
| A27 D6 / D6b / D6c (post hoc) | which gene classes stay locked to B when X is lost? glucocorticoids? | GTEx ventilated vs fast | clock outputs retained 0.657 vs matched 0.257, p 0.001; all rhythmic 0.173 (`gtex_cues/diag6_classes.csv`); GC score R^2 on host 0.234 ventilated, 0.014 fast; on hour 0.011 / 0.015 (`diag6b_gc.csv`); GC in-plane share below its random-score null (erratum) | random matched sets (D6); random scores (D6c) | current | clock outputs follow B (or L: not separable); most other rhythms need cues; GC is a candidate systemic reader in ventilated donors, not shown clock-specific | - |
| A81 / A81b (reg.; override post hoc) | a gene class that follows cues rather than the clock | GTEx: lock = host loading in ventilated / timed loading | A81 FAIL, uninformative (3-gene sets; D -0.057, p 0.975, `a81_h1_verdict.csv`). A81b frame check 0.481 vs guard 0.8 (`a81b_frame_check.csv`), overridden post hoc; 40 cue + 40 locked genes | locked set only | current | candidate classes, not a validated reader | selection too strict at half size; frame statistic badly chosen |
| A89 (reg.) | hierarchy of cue-like genes over organ groups | GTEx 26 tissues | FAIL: within-group R 0.915 vs random 0.870, p 0.269 (`gtex_cues/a89_verdict.csv`) | random tissue sets | current | no organ grouping of cues | any gene rhythmic in two tissues is in phase |
| A90 (reg., descriptive) | p53 targets locked or cue-like? | GTEx | all 17 undetermined (`gtex_cues/a90_summary.csv`) | clock / cue references | current | gene-level lock is too noisy to call single genes | - |
| A100 (reg.) | are drug targets more cue-like? | GTEx rhythmic pairs | FAIL: 0.285 vs 0.346, p 0.087 (`gtex_cues/a100_summary.csv`); debiased wall-time retention in ventilated donors 0 (A102) | other rhythmic genes | current | - | no difference; retention was noise floor |
| A101 / A101b (reg.) | age dampening of rhythms? | GTEx, host-referenced vs CHIRAL counts | both FAIL; CHIRAL's count loss is mostly death-type composition (`gtex_cues/a101_summary.csv`, `a101b_summary.csv`) | core genes | current | - | composition |
| A85 -> A122 (reg.) | do immune cells follow cues rather than the local clock? | 9 single-cell atlases | A85 PASS on corrected code (0.097, `results/figure_5/output/cancer/a85_verdict.csv`); A122 FAIL: cue minus locked -0.017, p 0.645; matched sets q95 0.119 (`exploratory/A122_cue_controls/output/h122.csv`) | locked set; 200 matched sets | current | none: the excess is low host coupling, not cue-following | not cue-specific |
| A112 (reg.) | do GTEx cue genes lose lock under arrhythmic feeding (mouse liver)? | Greenwell 2019, single tissue | FAIL, wrong sign: cue minus locked +0.230 [-0.096, 0.819], p 0.91 (`exploratory/A111_feeding_controls/output/a112_h112_1.csv`) | locked set | current | the GTEx cue class is not feeding-driven in mouse | single tissue (own-liver host); 24 vs 26 genes |
| A111 H111.3 / A129 (reg.) | does a liver without its clock decouple from the body? | Manella 2021, Bmal1 hepatocyte KO | PASS, specific: liver d* -0.287 [-0.411, -0.118], p 0.001; WAT -0.009, lung +0.019 (RP WP3) | WAT and lung | current | the only causal L handle on disk (mouse) | - |
| A124 (reg.) | does the cue-core offset detect a 12 h behaviour shift (muscle)? | GSE106800, 11 subjects | FAIL: 3.35 h, p 0.211; MDE 7.3 h (`exploratory/A124_night_shift/output/tests.csv`, `mde.csv`) | locked set | current | - | gene-set state too noisy; one-axis readout |
| A257 (reg., not blind) | does muscle follow a 12 h behaviour shift? | GSE106800, 14-gene projection | shift 2.68 h (90 % 1.81 to 3.59), excursion ratio 0.233 (`exploratory/A257_nightshift_muscle/output/score_tests.csv`) | M, R | current | H vs L: the muscle clock lags behaviour and damps | - |
| A136 (reg.) | is the brain-death loss a glucocorticoid effect? | GTEx Hardy 0 | gate H136.1 FAIL: +0.103, p 0.098; H136.3 0.458, p 0.007 (`exploratory/A136_glucocorticoid/output/tests.csv`) | hypoxia score | current | GC exposure is a feature of all ventilated deaths | gate (not robust: without long-term steroid users p 0.021) |
| A248 (reg.) | blood T-cell night peak; GC rhythm by hour | GTEx fast deaths | positive control FAIL (H248.1 p 0.91, fitted peak 15.3 h); H248.2 p 0.10 (`exploratory/A248_glucocorticoid_tcells/output/score_tests.csv`) | 200 matched sets | current | GTEx bulk cannot read T-cell timing; a GC rhythm is more visible on the hour than on the host (descriptive) | RNA shares at death are not a resting baseline |

Archive numbers with no table behind them (not to be cited): A27 D8 part A and D9, A27 D2 / D4 summaries, the A81b
amplitude-weighted shift (-0.50 h) and E1 medians, the A89 post hoc tissue-share table, the A90 PER1 / DBP medians,
the A101b same-29-tissue values. Docs that quote superseded values: the gtex_cues README check row for A46 (+0.0167;
table 0.0199), RP A122 (pre-relearn +0.007, p 0.44; table -0.017, p 0.645), the archive A47 T2D d* (-0.027 is the legacy
module mean; table -0.043 / module minus companion -0.062).

## 2. What the failures have in common, and what a design must do

1. **Gene sets read by an unvalidated state.** A81 / A81b, A85, A112 and A124 build a "cue state" from 40 genes and read
   it like the clock. The clock has a fixed, validated waveform; a cue set has none (A81b's frame check failed), so its
   state is noisy (A124: MDE 7.3 h) and its index rises for any gene set where host loading is low (A122). *Design:*
   no new state readers; per-gene regressions pooled into programs, with a phase reference from timed data (T1) or the
   host itself (T2).
2. **No specificity companion.** A91, A27 D5 / D6 and A92 / A93 have permutation floors but no matched or
   stress-responsive sets; A85 passed until A122 added them. *Design:* every partition count against M and R sets run
   through the same pipeline.
3. **Two times confounded by design.** In GTEx "internal time" includes the person's behaviour (A91), the tissue's
   clock and the body time are measured by the same 14 genes (A91 cannot separate L from B), and ventilation removes
   the day and adds drugs at once (A27, A136, A248). *Design:* use the designs where two times are moved apart by a
   known amount (night shift: H by 12 h, L by 2.7 h; repeat biopsies: X by a known number of hours; mouse KO: L removed).
4. **Small designs read at the gene level.** A90 (every p53 target undetermined), A112 (24 vs 26 genes, single tissue),
   A124 (11 subjects, one readout axis). *Design:* power by simulation at the program level before registering (G1).
5. **Calibration missing.** No analysis of the line was checked on a system where the answer is known; A112 was meant
   to be that check and failed. *Design:* the mouse liver knockout (Manella 2021, timed, three tissues) as calibration
   of the partition before any human design (T0c, G2 (ii)).
6. **Identification of L.** Nothing on disk estimates how large a person-specific tissue deviation from the body time
   is, which bounds what any "tissue clock" reader can detect in humans. *Design:* measure it first (A1 (c), A259).

What phase T reuses: A91 as the B-vs-X backbone (unchanged on current code), extended with companions and the L term
only where A259 shows a detectable deviation; A257's muscle shift and damping as the L reference of T1; the WP3 mouse
data (A111 / A129) as calibration; the GC score (A27 D6b) as the starting reader of T4. Dropped: the A81b cue state as
a reader, organ-group cue hierarchies (A89), cell-type cue indices (A85 / A122).
