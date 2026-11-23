# Prompt for the next agent: mechanism analyses (written 2026-11-23)

You add mechanism to the paper's main claim ("loudness, not time"). Context: `PAPER_PLAN_2026-11-19.md` section 11 and
`STORY_AUDIT_2026-11-23.md` sections 3 and 7.

Bulk data cannot separate dephasing among cells from per-cell damping. A334 / A334b closed that route with public
single-cell data, and A212 could not separate damping from rewiring with two summary statistics. The two tasks below
ask what bulk data can answer:
- is the clock intact but quieter, or rewired?
- what input drives the one local phase shift?

Follow the rules of `prompts/audit_followups_2026-11-23.md` ("Rules for this round"):
- register each task before scoring, power first, one scoring run;
- 14 genes and core loop;
- matched and responsive companions on the same permutations;
- CLAUDE.md statistics and covariate rules;
- machine budget;
- aggregates only for GTEx.

Run after, or alongside, the audit follow-ups. Tasks 2 there and 2 here share the skin series: register them together
or state the order, so that neither scores data the other has opened.

## Task 1: waveform screen, "intact but quieter" against "rewired"

- **Question.** In each condition that loses loudness, is the loss a uniform scaling of the 14-gene waveform (the same
  amplitude factor for every gene, phase relations kept)? Or is it gene-specific?
  - Uniform scaling means dephasing among cells, per-cell damping with the shape kept, or a weaker drive.
  - Gene-specific change means a changed oscillator, or an input acting on some genes.
- **Conditions and data (on disk).**
  - tumours: A304 TCGA pairs, high against low proliferation;
  - chronic lesions: A254;
  - intensive care: A315 groups, fast against ventilated post mortem, age / sex matched as A335;
  - psoriasis: A183 GSE154200, A345 Tsoi spliced and pre-mRNA;
  - as negative references: the A164 / A306 chronic-disease contrasts.
- **Readout.** Per condition and gene: the log amplitude ratio and the phase shift, case against control. Use paired
  estimators where the design is paired, group fits otherwise. Statistic: the across-gene heterogeneity (Q) of each.
- **Nulls.**
  - Uniform-scaling simulation at each condition's design.
  - Genes permuted.
  - Expression-matched and responsive sets: the same heterogeneity computed on their pseudo-waveforms.
- **Calibration first (stop rule).** On the A212 plasmode (real tumours with damping or shared rewiring planted at the
  observed strength), the screen must separate them (registered AUC bar, for example 0.9). If it does not, stop and
  report "not separable at gene level"; that is informative in itself.
- **Readings to fix.**
  - Uniform in tumours, lesions and intensive care: "the clock is intact but quieter".
  - Gene-specific: name the genes (loop arm), and check whether they are PER1, PER2 and NR1D1 (signal layer) or core
    genes.
  - For psoriasis, a delay concentrated on a few genes points to a local input rather than a shifted clock.
- **Caveat to state.** Uniform scaling does not separate dephasing from per-cell damping.

## Task 2: cytokine dose-response of the inflamed-skin delay

- **Question.** Does the lesional delay scale with IL-17 / TNF activity, and less with Th2 activity? That would explain
  why it is largest in psoriasis and smaller in atopic lesions (A345, archive A53).
- **Data.**
  - The psoriasis and atopic lesional / non-lesional series under `$SCCC_DATA/skin_lesional/` (A53 / A54 sets);
  - GSE154200 (A183);
  - Tsoi with pre-mRNA (A345).
- **Readout.**
  - Per lesion, an IL-17 / TNF signature score and a Th2 signature score. The gene lists are fixed in the registration
    from the literature, with no clock gene.
  - Strata are signature tertiles within series. The pooled paired offset per stratum uses `clocksync.paired`, not per
    patient (per-person offsets are too noisy, A259).
  - The trend across strata comes with companions and with a disease-label-permuted null.
  - In Tsoi, per-gene pre-mRNA against mRNA delays, to see which genes move first (descriptive).
- **Blind status.** Not blind for the per-series offsets (A53, A54, A183, A345 seen). The signature strata are new.
  State both.
- **Covariates.** Series as strata. Disease (psoriasis or atopic) is both the contrast and a confounder of signature
  level: report within-disease trends as primary, across diseases as secondary.
- **Readings to fix.** A trend with IL-17 / TNF within psoriasis makes "an inflammatory input delays the clock's
  transcription" the mechanism of the paper's one local phase shift. No trend means the delay is a disease-level
  property, cause open.

## Not in scope
The in-house single-cell time-lapse question (does faster division dephase or damp the clock?) waits for the PI
(PI_QUESTIONS Q44).

## Deliverables
- Amendments committed before scoring; results in `results_preliminary.md`; VERDICT_TABLE rows.
- A short note (proposed, not applied) appended to `PAPER_PLAN_2026-11-19.md` section 11.4.
- Commit with the attribution line, `git push origin main`.
- Reply to the user in at most 15 lines.
