# Story audit (2026-11-23): the results against the story line

Prompt: `prompts/story_audit_2026-11-23.md`. Read, in this order: CLAUDE.md and the memory index; ABSTRACT_DRAFT_2026-11-22;
PAPER_PLAN_2026-11-19; `reviews/plan_review_2026-11-19.md` and REVIEW_2026-11-19; PI_QUESTIONS (open list, Q41-Q43);
VERDICT_TABLE and `triage/registry_2026-11-19_summary.md` (with the registry TSV); the results entries (RP) and
registrations (RE) cited below; METHOD.md (sections 7-9); DATA_ROADMAP_2026-11-20.

No analysis was run, and no plan, abstract, results or PI file was edited. Every number comes from an RP entry or an
output table named with its amendment ID. "Derived" marks arithmetic on those values. Line references are to the files
as of commit 3ea19f7. One Sonnet subagent extracted the statistical-rule checks per contrast (section 3); its numbers
were checked against the RP lines it cited before use.

## Summary

1. **Coverage.** Five newer results (A337 Perrin, A342, A343, A344, A345) are in neither document. Two strong older
   results are missing as well: A251 and A255, the living-person evidence that tissues share one time. Neither
   document cites the positive controls of the time axis (A129, A111).
2. **Claims weaker than written.** Eight abstract phrases need a qualifier or a rewording (section 1.2). The largest
   gaps:
   - the tumour stroma clause rests on a post hoc result with no core-loop run;
   - "cancer cells stay coupled": malignant host coupling failed (A138 H138.1);
   - "organs keep a common body time through intensive care": the registered primary failed (H109.2). The order kept
     is a descriptive per-group readout (A132), re-tested not blind (A293), with offsets compressed toward 0 / 12 h;
   - "psoriatic skin is the exception": atopic lesions run late too;
   - "the only state in which the shared time comes apart": brain-dead donors keep a shared clock state (H317.2).

   In the plan, the CHIRAL age correction has no registered entry.
3. **Coherence.** The adult middle carries the paper. Both ends of the lifespan arc are the least secure:
   - the fetal chapter cannot use clocksync's host-referenced readout;
   - brain death is post hoc and aliased with procurement.

   The abstract's last sentence gives brain death the closing position that Q41 advises against.
4. **Robustness (section 3).** Most headline statistics hold on the core loop; tumour pairing and the psoriatic delay
   sit at the line.
   - Never run on the core loop: tumour stroma pairing, A251 / A255, conserved class gains (A293 H293.4).
   - The intensive-care detection limit was simulated on fast deaths, which have the stronger shared signal, so it is
     optimistic.
   - Brain-death covariates are regressed, not balanced (they cannot be balanced in GTEx).
   - The spread bias of the A164 / A202 / A208 nulls is unmapped.
   - Several companions are not on the same permutations (A183, A254, A332 D1) or not stated (A345, A252, A256,
     A257, A331).
   - RIN is in no chronic-disease balance table.
5. **Additions, by gain over cost.**
   1. Tumour stroma pairing on the core loop.
   2. Inflamed skin against psoriasis, pooled with companions.
   3. CHIRAL age rescore with companions and death-type matching.
   4. The in-house jet-lag data (PI).
   5. The STARNET application (user).

## 1. Coverage

### 1.1 Results with a pass or a powered, informative null that are in neither the plan nor the abstract

**Newer than the plan (2026-11-21 to 2026-11-23).**

| ID | result (registered verdict, key numbers) | blind | where | reason |
|---|---|---|---|---|
| A337 Perrin | H337.1, H337.2 PASS (14 genes): untimed pre-mRNA lead +1.41 h (90 % 1.17-1.66) against timed +1.75 h (1.47-2.05). Core loop H337.2 FAIL: untimed minus timed -1.34 h (-2.19 to -0.70). H337.3 FAIL | blind for unspliced counts | Supplement (lead validation); one line in Figure 1 | the cleanest blind human ground truth of the lead readout; the core FAIL limits any core-loop lead |
| A342 | H342.1-2 PASS: timed skin lead +1.18 / +1.48 h, recovered untimed at both sites. Untimed exposed-minus-protected offset +1.33 h against timed +1.57 h, but the pairing is not beyond its companions (H342.3 FAIL). Untimed log loudness ratio -0.253 against timed log amplitude ratio +0.087 (H342.4-5, descriptive) | rescored on direct introns after the registered subtraction run was seen (erratum, user decision) | Figure 1 (the offset readout against known time in living people); Supplement | first living-human check of the offset itself; H342.4 is the evidence behind METHOD section 9 |
| A343 | H343.1 PASS: acute dexamethasone raises GR-target pre-mRNA far ahead of mRNA (fat D +1.934, companion q95 0.226; liver +2.393, descriptive). PER1 +0.495 log in fat (0.240-0.770) and +0.522 in liver (0.071-0.945), dexamethasone minus none. Core-loop displacement not resolved | registered positive control | Figure 1a, or Supplement | the only human interventional evidence that a systemic signal reaches PER1 within hours: the rationale of the core-loop readout |
| A344 | H344.2 PASS: on 06 / 12 / 18 h samples the untimed lead is within 0.10 h of the full-cycle truth (P +0.10, E -0.09). The primary, H344.1, is undecidable: offset 1.37 h against 1.57 h, difference 90 % -1.68 to 0.88. Per-visit phase error 2.03, 2.59 and 3.27 h (full cycle, three times, two times) | not blind for the full cycle | Figure 1 or the introduction, with H344.1 stated | the only living test of the hospital-hours claim: relative readouts survive, sample times degrade, the offset is not decided |
| A345 | H345.1 PASS: pre-mRNA delay -2.44 h (-4.35 to -1.05), core -1.79 h (-2.88 to -0.71). H345.3 PASS narrowly: -0.345 against matched q05 -0.316 and responsive q05 -0.229 (core -0.299 against -0.293 and -0.092). H345.2 not shown: D -1.05 h (-2.64 to +0.11). AD (descriptive): spliced -1.15 h (p 0.001), pre-mRNA not beyond companions | spliced not blind (A54 discovery series); introns blind | Figure 6d | the first result beyond both companion families for the psoriatic delay (at the pre-mRNA level, narrowly); the AD delay changes the framing (section 2.2) |
| A328 part a | H328.1-2 PASS: del Valle replicates blind (C 0.592, p 0.0021); the core loop covaries on its own, clock-specific (C 0.803; matched q95 0.470, responsive q95 0.605) | blind | Figure 3a-b, wording paused (Q43) | cited in plan 6.2b only; Figure 3(b) ("it rests on the signal-induced genes") contradicts it |
| A340 | H340.1 PASS: the fetal states are not a single slow factor (iota 0.41 / 0.55 / 0.68 against line q05 0.97 / 0.91 / 0.95) | blind | Figure 3 | cited in the abstract's limits table, absent from the plan |
| A328 part b, A338, A339 | no consistent pre-mRNA lead (H328.3 FAIL), uninformative: the running postnatal liver reads the same (0.33 h, -0.30 to 0.78); damping ratio uninformative | blind | Figure 3e "the open test"; Supplement | they close the time-free route on open data |
| A341; scout 2026-11-22 | no GC-responsiveness or maturity ranking of fetal organs (descriptive); no open cohort supports a fetal organ-to-organ offset test | - | Supplement, Discussion | the fetal chapter cannot use the host-referenced readout |

**Older results (A109-A327) not cited.**

| ID | result | blind | where | reason |
|---|---|---|---|---|
| A251 | H251.1 PASS: skin at three sites shares one clock state, r 0.628; H251.2 PASS x2 (layers); H251.4 PASS: untimed dermis minus epidermis -1.58 h against timed -1.98 h; H251.3 FAIL (one site interval beyond ±1.5 h). Calibrated: beyond matched genes (p_cal 0.004) and the broad responsive family (p_cal 0.0075, A320); against the narrow glucocorticoid / immediate-early family not decidable with three to six pairs (p_cal 0.075, A285) | blind (inventory) | main, Figure 1 or 4 | living people, separate biopsies and separate array batches, with a timed check |
| A255 | H255.1 PASS: fat and skin share a clock state, r 0.533 (p_cal matched 0.011, broad responsive 0.0065; narrow family not decidable, 0.057 at a floor of 0.056; A285, A320); H255.2 PASS: the person's own cell line does not share it (0.042); H255.3 PASS: two surgical cohorts, 0.686 (83 persons), beyond both families (A285); H255.4 PASS x2 (co-twins) | blind (inventory) | main, Figure 1 or 4 | the genome control shows that the shared state is body time; no other result has it |
| A129, A111 | hepatocyte Bmal1 knockout: the liver decouples, d* -0.287 (-0.411 to -0.118) at a common spread, companion organs null; Xin daytime feeding: untimed offset changes track timed ones (Fisher-Lee 0.988, exact p 0.042) | not stated | Figure 5 Supplement (positive controls of the time axis) | every "kept" needs a real positive control on the same readout; A336's detection limits are simulated |
| A211 | regenerating mouse liver: synchronous proliferation does not lower clock coherence (T -0.027, MDE 0.064; powered FAIL); amplitude falls moderately (median 0.68 against 0.89 log2, descriptive) | blind (inventory) | Figure 6b Supplement | proliferation alone does not lower coherence; it bears on, but does not settle, whether A304's gain loss is a tumour property |
| A221, A207 | clock genes are lower in tumours than in paired normal tissue; GTEx-normal designs exaggerate the loss and invert it for ARNTL, NPAS2 and NR1D1 (H221.1-3 PASS). Co-expression and order "disruption" scores are present in tumours that follow their normal's clock (dCCD 0.890, 9 of 10 cancers) | blind / partly blind (inventory) | Figure 6 Supplement (corrections) | the literature-cleaning tests the user kept in the paper (A220-A224) are cited only through A220 |
| A165, A180, A167 | single-cell disease cohorts: no phase, amplitude or noise effect (0 of 89 rows; the best rows exclude about 1-1.5 h or 1.3-1.5-fold); TOST: 7 of 8 phase rows equivalent within ±1.5 h; liver fibrosis cell classes 0 of 9 | partly blind / not blind | Figure 6a Supplement | cell-class chronic-disease nulls; bounds, not equivalences |
| A200, A174 | breast surgery: per-patient internal time against recorded resection time, Spearman 0.558 (pseudo q95 0.215); heart AM / PM confirmed on current code | registered re-test / not stated | Supplement (validation) | per-patient times in surgical tissue, with A344's caveat |
| A127, A166, A171, A173 (with A172, A178) | tumour-only clock coherence falls with proliferation or ER-negativity in SCAN-B, METABRIC (ER), lung adenocarcinoma and gastric cancer; not a general tumour property (A172, powered FAIL); the SCAN-B proliferation loss is not robust to gene choice (A178) | blind / partly blind (inventory) | Supplement, as the superseded readout | reviewers know the archive claim; A304 replaces it (Q40), with its mixed record stated |
| A259 | a tissue's departure from its person: upper bound 1 to 3 h (median upper bound 2.20 h), not shown to be clock (H259.1-3 FAIL) | partly blind (inventory) | Supplement | the only per-person bound in hours; it limits what "across a person's organs" can mean (section 5) |
| A230 | inflamed gut epithelium against stroma +0.58 h, MDE 1.71 h | partly blind | drop, or one line in the Figure 6d legend | underpowered for a psoriasis-sized shift |

### 1.2 Claims with no supporting entry, or with a weaker verdict than the text implies

| claim (location) | entry and verdict | gap |
|---|---|---|
| "the surrounding stroma carries it" (abstract l. 15); limits table: "A271 (fibroblasts and endothelium beyond both nulls)" | A271 Part A lists fibroblast 0.686 (responsive q95 0.638) and endothelium 0.677 (0.659) as descriptive "anchors". The clock-specific stromal pairing is post hoc: A149 note 5, fibroblasts 0.699 against pseudo q95 0.542, endothelium 0.701 against 0.527. It could not be confirmed in new cohorts (A226: 11 and 13 new patients). The registered statistic shows a patient-specific pairing only (A151 H151.1 PASS; H151.2 uninterpretable) | post hoc, unconfirmed, and no core-loop run (A336 has no stromal row) |
| "cancer cells stay coupled but grow quieter as they proliferate" (l. 15-16) | A304 is bulk tumour: gain ratio 0.63 (Q_G -0.462, Holm 0.011), departure not raised. Malignant host coupling FAILED (A138 H138.1: 0.076, below the calibrated threshold 0.083). A305 was not powered. A304 itself: "bulk cannot say whether malignant cells carry a weaker clock" | "cancer cells" is not shown; "tumours" is |
| "Tumours follow their own patient's clock" (l. 15) | A138 H138.4 PASS on 14 genes (0.365); at the line on the core loop: 0.252 against matched q95 0.275 and responsive q95 0.313 (A336). Tumour time is a pairing; only a dispersion of 2.4 h or more is detectable (A336 part 3) | "follow" is at the line on the core loop |
| "organs keep a common body time through ... intensive care" (l. 13-14) | The registered primary H109.2 FAILED in ventilated donors (rho 0.261, p 0.165). The order kept in ventilated post-mortem donors comes from a per-group readout registered as descriptive (A132; the between-group test, A132b, is post hoc). A293 H293.2 PASS (rho 0.690) is a registered re-test, not blind. The offsets are compressed toward 0 / 12 h: deviations at 0.58 of fast deaths (A132b); folded distance 0.77 against 1.01-1.02 h (A132). Only a single-organ shift of a median 3.9 h, or a dispersion SD of 1.5 h, is detectable (A336 part 3) | the order is kept, and the compression reads as a shared non-temporal displacement on running clocks (A109 circularity, A132b): "same order, offsets compressed", not "same offsets". "Time" is kept only at the stated limit, and the result is a descriptive readout after a failed primary |
| "What disease changes is how loudly a tissue expresses that time, not the time itself" (l. 14) | Holds for tumours (A304, holds on the core), intensive care (A315 / A335, holds on the core and matched) and lesions (A254 H254.3 at the line; at the line on the core). Chronic systemic disease changes neither (A306; five loudness and departure nulls lose their equivalence on the core, A336). Inflamed skin changes time | a pattern over three conditions stated as a rule |
| "Inflamed psoriatic skin is an exception, running about 1.5 h late" (l. 16-17) | Atopic lesions run late as well: Tsoi AD spliced -1.15 h, p 0.001 (A345, descriptive). The archive A53 per-series table gives AD a median of -0.45 h, 6 of 8 gated series negative, against -1.50 h for psoriasis, 9 of 9 (derived from `figure_5/output/psoriasis_details/a53_skin_per_series.csv`). Psoriasis-specificity was not supported (A183 H183.2, NR1D1-dependent; on the core it falls, not decidable, A336). A254's own reading is "specific to inflamed skin" (hidradenitis -1.07 h, not significant after Holm). Core loop -1.00 h (A335); Tsoi core -1.20 h (A345) | "psoriatic" should be "inflamed skin, largest in psoriasis"; "about 1-1.5 h" on the core loop |
| "brain death ... the only state in which the shared time itself comes apart" (l. 17-18) | A317 H317.2 FAIL: the organs of brain-dead donors keep a shared clock state ("rearranged, not dissolved"). Psoriatic lesions are a local phase offset between a person's compartments, which is internal desynchrony as the project defines it. The order change holds on the core loop (D_rho 0.465, below its MDE 0.622 but significant; A336), aliased with procurement (Q41) | "comes apart" contradicts H317.2; "only" holds for body-wide order only |
| "no method has measured it across a person's organs" (l. 10) | clocksync uses within-person covariance but reports population-level coupling and offsets. A person's own departure is bounded at 1-3 h and not shown to be clock (A259); per-sample times degrade on partial arcs (A344) | reads as a per-person measurement, close to the wording the framing decisions exclude |
| Plan section 1 and claim 5, "the clock layer is built around birth"; Figure 3(b) "it rests on the signal-induced genes" | A328 part a (blind): the core loop covaries on its own, clock-specific (C 0.803). A266 (cortex) and A269 (lung, liver) do not show it | contradicted in one blind cohort; paused under Q43 |
| Plan 2.1: PER1, PER2 and NR1D1 "carry most of the departure's excess (N2, A329) and, in living tissue, most of the person leverage (A335; precision check pending)"; "so it rides mostly on PER1, PER2 and NR1D1" (l. 91-92, 100-101) | The precision check is done (A336). The trio drop is not extreme for the departure (H294b.1: 84th percentile of all drops, 32nd among drops of similar weight) or for the living leverage (A331: 13th / 43rd). A336: no statistic reads "carried by PER1, PER2 and NR1D1" except a mechanical A306 row. Plan line 80-81 already says "not attributable" for the leverage | the plan contradicts itself and A336 |
| Plan 2.3: "lymphocyte content leaves relative gain unchanged (A335: equivalence, post hoc margin)" (l. 131) | A336 tier 1: A295 H295.1 changes verdict on the core loop, PASS: -0.020 per SD, p 0.040 | a small dilution is detected on the core loop (inside the ±0.10 margin) |
| Plan claim 6: "CHIRAL's age dampening, which shrinks with death-type matching (A208, A231 ...)" (l. 245-246) | A208 and A231 test coupling and offsets. A231: "does not contradict reports of rhythm dampening with age". The shrinkage is WP0 item 5: CHIRAL's timed rhythmic-gene count, pooled -0.73 to Hardy-matched -0.13 (CPM); TMM matched minus pooled +0.39, p 0.092. It is not registered and has no companions | no registered entry; the cited IDs do not test it |
| Plan claim 6: "Prognostic clock signatures reduce to proliferation (A191-A194, A220)" (l. 247) | Only A192 is "carried by proliferation". A193 and A194 are "not prognostic" or "not replicated"; A220: "no better than random gene signatures" | narrower than written |
| Plan section 3 and section 10: "amplitude about 0.82", "amplitude 1.02", "amplitude about 0.81", "amplitude 1.28" (derived from variance ratios) | METHOD section 9: a variance ratio measures loudness, not daily amplitude. In timed skin, untimed loudness -0.25 against timed amplitude +0.09 (A342). The Delta values for intensive care and lesions mix loudness with timed amplitude ratios (A257, mouse) | wording, and a scale that is not common |

## 2. Coherence

### 2.1 Does the story line follow from the evidence?

**The middle holds; the two ends are the least secure.**
- **Adult organisation and disease carry the paper.** Living rotation (A252, A256: hold on the core loop); living organ
  sharing (A251, A255, missing from both documents); cell-class participation (A274 and A280 hold against both
  families on the core loop; A283 H283.2 and A288 H288.2 at the line); intensive-care loudness (A315, A335: core and
  age / sex matched); powered chronic-disease nulls; A304; the COFE correction.
- **Before birth** rests on covariation between fetuses within one organ (L1). It is not the host-referenced readout:
  - no open cohort supports an organ-to-organ test (scout 2026-11-22);
  - the time-free route is exhausted (A328 part b, A338, A339 uninformative).

  A328 part a and A340 support "the clock genes move together, in two dimensions", not "driven" or "running".
- **Brain death** is a post hoc contrast (A132b), found after the registered primaries failed (H109.2, H132.1-2). It
  was re-tested on the same donors (A317, A318; not blind) and is aliased with procurement (Q41). It holds on the core loop (A336 part 4).
- **The method's edge is not used at either end.** GTEx deaths cover the day, and there orderers recover the map
  nearly as well (A123b: CHIRAL 0.877 with whole-day sampling). The fetal chapter uses L1. The hospital-hours edge
  (A123b working hours: clocksync 0.945 against orderers at -0.17 or below; A344) belongs to the adult living and
  surgical cohorts: A251-A256, A254, tumours, A345. "Method as enabler" is true of the adult chapter, not of the arc.

**"Loudness, not time" holds as a pattern, not a rule.** Over the human perturbations:
- loudness only: intensive care, tumours, chronic lesions;
- neither: chronic systemic disease (powered);
- time: inflamed skin, and the order between organs in brain death;
- both: the night shift (core loop +4.06 h, excursion about 0.49 kept; A335).

Mouse daytime feeding moves both (A322, A324). "Where human time moves, it moves locally or in one aliased state" is
what the record supports.

**Position of brain death.** The lifespan order puts it last, and "the only state in which the shared time itself
comes apart" makes it the climax. Q41's recommendation is "not the climax without new data".

### 2.2 Contradictions and tensions

| tension | evidence | consequence |
|---|---|---|
| **AD delay against the psoriasis framing** (not against A254) | A345: AD spliced -1.15 h (p 0.001); archive A53: AD median -0.45 h against psoriasis -1.50 h (derived); A183 H183.2: psoriasis minus eczema -1.09 h, not clock-specific, NR1D1-dependent. A254 found no delay in four non-skin lesions (pooled +0.01 h, MDE 0.69 h) and itself reads the shift as "specific to inflamed skin" (hidradenitis -1.07 h). The A345 registration and results text call A254 "no delay in other chronic lesions", which drops that qualifier | the record supports "inflamed skin runs late, psoriasis most", not a psoriasis exception; addition 2 would decide it |
| **Living person loudness falls on the core loop** | A331 0.177 → 0.107, below matched q95 0.130 (A335); the trio drop is at the 13th / 43rd percentile (A336), so the fall may be lost precision. GTEx holds on the core (0.714 against matched q95 0.621), but the matched sets are high there too (median 0.119): a donor-level scale is large in deaths (A336 part 2) | shown as a clock property in deaths only. The living limit is the design (two tissues from one punch biopsy, one array batch), not hospital hours. It does not touch the abstract, but Figure 4f and the model's A_p must say so |
| **Brain-death arrows not stable across gene sets** | D2: lung -1.58 h on 14 genes, -1.47 h (-2.92 to +0.21) on the core; lower-leg skin +1.33 h appears; heart +0.32 → -1.12 h; 3 of 12 intervals exclude 0 on either gene set against 1.2 expected by chance (A332, A335) | "rearranges the order" holds as one statistic; no organ-level arrow can be named |
| **Fetal: driven or running is undecided, and the core loop disagrees between cohorts** | A328 part a: core C 0.803, blind; A266: cortex core C 0.417, p 0.057; A269: lung core p 0.060, liver core not beyond responsive; A328 part b, A338, A339 uninformative; A340 PASS | the abstract's "moves as a whole" stands; the plan's "signal before clock" does not; Q43 open |
| **Intensive care: order kept, offsets compressed** | A109: H109.2 FAIL pooled; the ventilated shared covariance is about half a displacement along one axis (best-matching fraction 0.44), host phases axially concentrated (axial R 0.26 against 0.03 in fast deaths; Rayleigh R 0.014). A132b: deviations at 0.58 of fast deaths | "time kept" means rank order kept with a shared non-temporal displacement; state it, or reviewers will find it in the archive |
| **Tumour pairing at the line, stroma post hoc** | section 1.2 | the tumour sentence of the abstract is the least robust adult sentence |
| **The plan against its own updates** | plan 2.1 l. 91-92 and 100-101 against l. 80-81 and A336; plan 2.3 l. 131 against A336 (A295 H295.1 PASS on the core) | two internal contradictions to fix in the plan |
| **Hospital-hours edge** | real data: A123b PASS (GTEx working hours); A123c H123.5 FAIL (0.679 against the registered 0.7; the timed oracle 0.491); simulation: A145 (no advantage after Holm); living: A344 lead PASS, offset undecidable | the edge rests on one real-data comparison and one partial living test; state it as such in the introduction |
| **Blood as a composition exception** | RP post hoc entry (2026-11-03): blood has no core-clock time signal (cosinor R^2 0.0096, 52nd percentile of genes), but composition removes 6 % of its weak time signal; "the composition explanation is only weakly supported in post-mortem blood" | keep "no readable core clock in blood"; give composition as the likely cause, descriptive |

## 3. Robustness

### 3.1 Headline claims under the two systemic-signal controls

Core loop = the 11 genes without PER1, PER2 and NR1D1, with core-restricted companions (A335, A336). Responsive = the
stress-responsive companion family; "broad" is the calibrated A320 family, "narrow" the glucocorticoid /
immediate-early family.

| claim | core loop | responsive companions | status |
|---|---|---|---|
| untimed readouts recover known time (A123, A270, A273, A279; A342, A344) | A270 holds (human 0.631 → 0.627); A273, A279 offsets hold; A342 core leads within ±0.6 h of timed; A337 Perrin core lead FAIL (-1.34 h) | A270 companions were read with the clock's waveform, which tilts the comparison toward the clock (A335 note) | holds; the core-loop lead is not validated in muscle |
| one time per living person: rotation with elapsed hours (A252, A256) | holds (1.265 → 1.073; 2.970 → 3.902; 2.533 → 3.159; A336) | beyond both (A252 responsive q95 0.738; A256 0.779); adipose narrowly (p 0.045) | holds on both |
| one time per living person: tissues share a state (A251, A255) | **not run** (only without PER1 and NR1D1, PER2 kept: A251 0.637; A255 0.387, no companions) | broad family PASS (A320); narrow family not decidable with few pairs (A285) | gap |
| cell classes carry the clock, lymphocytes do not (A272-A288, A321, A326) | holds against matched genes; lymphocyte bound unchanged | holds in 3 of 6 tests (A283 H283.3, A280 H280.1, A282 post mortem); lost in 3 (A283 H283.2, A288 H288.2 at the line; A282 living, descriptive); cortex not decidable (A336 tier 2) | holds against matched; against responsive in half |
| organ order kept without daily alignment and in intensive care (A293 H293.1-2) | holds (0.784 → 0.831; 0.690 → 0.700) | above all 200 sets on 14 genes (A293) | holds; see 3.3 for the detection limit |
| intensive care lowers loudness (A315, A335) | holds (Q* +0.404 → +0.457, above all matched and responsive sets) | holds | holds; age / sex matched (+0.389 / +0.459) |
| chronic disease leaves time, coupling, gain unchanged (A164, A233, A235, A306) | A164, A233, A235 hold; five A306 loudness and departure nulls lose their equivalence (nothing detected) | none: companions are tissue families (tibial nerve and two skins); offset pseudo-clock q95 1.3-10 h cannot certify a small shift (A208, A231, A233, A235) | holds as nulls on the core, except the A306 equivalences |
| tumours pair with their patient (A138 H138.4) | **at the line** (0.365 → 0.252; matched q95 0.275, responsive q95 0.313) | beyond both on 14 genes, calibrated (A285 p_cal 0.004 / 0.011; A320 broad 0.0055) | at the line |
| tumour stroma carries the patient's time (A149 note 5) | **not run** | descriptive anchors above responsive q95 (A271) | post hoc; gap |
| proliferating tumours lose gain, not departure (A304) | holds (-0.462 → -0.668, below matched q05; departure rise still not detected) | 0.5th percentile on 14 genes | holds; purity-stratified -0.266 (p 0.037) |
| chronic lesions keep time, lose loudness (A254) | time holds (+0.31 h, -0.14 to 0.75); loudness at the line (-0.417 → -0.325) | loudness not beyond responsive on the core | time holds; loudness at the line |
| psoriatic delay (A183, A308, A345) | at the line (-1.00 h, q 0.049, 90 % 0.0003-0.130; A335); pre-mRNA core -1.79 h, H345.3 core passes narrowly (-0.299 against matched q05 -0.293) | core: 5.5th percentile of matched and 1st of responsive sets (A335); pre-mRNA beyond both (A345) | at the line; strengthened at pre-mRNA (discovery series) |
| brain death changes the order (A317, A318) | holds (0.621 → 0.465; adjusted 0.496), below its own MDE (0.622) | above core responsive q95 0.180 | holds; covariates not balanced (3.2) |
| brain death keeps loudness (A315 H315.1, A335 part 3) | equivalent within ±0.357 on the core, pooled and matched | relative to matched genes only; raw variance lower (core Q -0.186) | holds as worded |
| night shift: core loop shifts, loses excursion (A257) | holds (+4.06 h; loss beyond matched q95 0.38 and responsive q95 0.81) | holds | holds (not blind on 14 genes; 11 men) |
| person loudness (A312, A331) | GTEx holds (0.765 → 0.714, above matched q95 0.621, responsive q95 0.692); living at the line (0.107 < 0.130) | living: not beyond responsive on 14 genes (93.5th percentile) | deaths only |
| coupling reads loudness (A291, A294b) | holds (0.717 → 0.490; 0.969 → 0.250); A294b slope not decidable | - | holds |
| composition does not move time (A295) | holds (0.130 → 0.151 h per SD); lymphocyte dilution of gain now passes (-0.020, p 0.040) | - | holds; see the wording fix |
| conserved relative class gains (A293 H293.4) | **not run** (tier 3) | - | gap |
| before birth, clock-shaped covariation (A242-A244, A328 part a, A340) | A328 part a PASS (blind, 0.803); A266 cortex fail (0.417, p 0.057); A269 lung fail (p 0.060), liver not beyond responsive (0.686 against 0.712) | A328 part a beyond both | mixed: one blind cohort positive |

### 3.2 The statistical rules of CLAUDE.md

**A failed gate is uninformative.** No violation was found. Failed gates and imprecise statistics are read as
"uninformative" or "not decidable":
- A328 part b, A338 and A339 ("uninformative, not driven");
- A337 (precision at n = 12);
- A241 (adult pass rate 0.14);
- A109c COVID.

The matched GTEx disease contrasts (A164, A208, A231, A233, A235, A306) state power or MDE instead of a per-arm gate.

**Coupling only at a common spread, or on z2 / psi.** Met in form: A164, A208, A231 and A233 use d*; A306 uses z2 and
psi. In substance, A258 found that d* is biased at unequal spread and that the label permutation is anti-conservative
above |log V ratio| about 0.4.
- A258 flagged rows at anchor noise 0.05 / 0.15: A231 0 / 3 of 54; A233 6 / 33 of 302; A235 5 / 55 of 262.
- A164, A202 and A208 have no spread columns and were not mapped (A258).
- The spread-matched null A258 names as missing was never registered.

The loudness contrasts are variance ratios. They include the person-to-person spread of internal time (METHOD section
9), and the A315 registration records the host spread of ventilated donors at about a tenth of that of fast deaths
(separate fits). Neither the A315 nor the A335 reading reports or uses it.

**A specificity companion on the same permutations.**
- Calibrated placebo forms exist for A251, A255, A138 and A227 (A285, A320), for the pairwise cell-class test (A283)
  and for the psoriatic delay (A308).
- Not tied to the permutations: A254 and A332 D1 ("companions ... not tied to the swaps"); A183 (companions scored
  once per set, a registered deviation).
- Not stated: A345, A252, A256, A257, A331.
- The GTEx disease d* companions are tissue families, not gene sets.
- No companion: A115; A312 cannot have one.

A254 H254.3 ("at the line") and A345 H345.3 ("narrowly") rest on uncalibrated q05 rules. In the one design where
the rule was calibrated (A308, psoriasis), it was conservative.

**Covariates that change rhythm strength are balanced, not regressed.**
- *Brain-death order (A317, A318).* Ischaemic time (258 against 561 min) and RIN are regressed, which removes mean
  effects only (A336), and cannot be balanced across collection modes (A335 part 2). Age (SMD 0.78 before matching)
  and sex were matched for the A315 loudness contrasts only, not for D_rho. The common-range restriction is borderline
  (A318 post hoc: 0.239 on 36 pairs).
- *Intensive-care loudness.* Matched on age and sex, holds. RIN is balanced (SMD 0.14); ischaemic time is not (SMD
  -1.06 after matching) and cannot be.
- *Chronic disease.* Matched within cohort on sex, Hardy, age and ischaemic time; pooled and matched estimates are
  reported where both exist (A164).
  - RIN is in none of the balance tables (A164, A231, A233; A306 uses the A233 matching).
  - Hypertension in ventilated donors keeps an age imbalance (SMD 0.54), with a per-cohort module d* of +0.109
    (p 0.004) against -0.032 in fast deaths (descriptive, A233).
  - T2D and obesity are not separated (BMI SMD 0.68 in the T2D fast design, A235).
  - T2D is unbalanced in the age contrast and not adjusted (SMD 0.46 / 0.79, A231).
- *Before birth.* Age, sex and assay are residualised (mean only). The fetal-against-adult size comparison (A240)
  differs in collection (terminations against deaths). The sampling-hour part was checked (a 4 h adult window does not
  explain the negatives); the termination method is not recorded.
- *A115.* Different persons; site and batch are unknown.
- *A304.* Purity is stratified, and pooled and stratified values are both reported.

### 3.3 Gaps (ranked by what they put at risk)

1. **Tumour stroma** (abstract): post hoc, not confirmed, not on the core loop. Addition 1 closes it.
2. **Living organ sharing (A251, A255)**: the cleanest living evidence for "one time per person" has no core-loop run
   with companions. It is cheap, on the A336 pattern.
3. **The intensive-care "kept".**
   - The detection limit (a single-organ shift of 3.9 h) was simulated on fast-death halvings (122 donors per half;
     A336 registration: "fast deaths only, ... not the ventilated"). Ventilated donors have weaker shared covariance
     (A109 excess 0.310 against 0.469 in fast deaths) and about a tenth of the host spread, so the true limit is
     larger than stated.
   - The relative organ-gain pattern under the overall loss (review question 1) is not tested.
   - The offsets are compressed (A132b).
4. **Brain-death covariates** cannot be balanced in GTEx (R4). This is stated in A336, but the paper must say it
   beside the result, not only the alias.
5. **Spread bias of the chronic-disease nulls.** A164, A202 and A208 are not mapped, and no spread-matched null
   exists. The pooled verdicts are exposed by about 0.03 typically and 0.06-0.15 in flagged rows (A258), in the
   direction "the narrower-spread arm reads lower".
6. **Uncalibrated or unshared companions** in A254, A332 D1, A345, A252, A256, A257 and A331.
7. **Conserved class gains (A293 H293.4)**: not on the core loop.
8. **Provenance.** A336 tier 2 is cited in RP since commit 14d2c79 (2026-11-21), but its scripts are untracked:
   `p1_a227.py`, `p1_a233.py`, `p1_a273.py`, `p1_a279.py`, `p1_a280.py`, `p1_a282.py`, `p1_a283.py`, `p1_a288.py`,
   `p1_a293h3.py`, `p1_a326.py` and `ct_lib.py`. The merged `two_readouts.csv` (129 lines committed, 137 on disk) and
   the VERDICT_TABLE row are uncommitted. This audit leaves them as found.

## 4. Additions (drafts for the user; nothing registered)

Ranked by gain over cost. Data on disk or listed in DATA_ROADMAP_2026-11-20. Each needs its amendment and a power stage
first. All five were chosen after seeing the results (post hoc).

| rank | question | data | statistic | pass / fail means | power (rough) | cost |
|---|---|---|---|---|---|---|
| 1 | Do tumour fibroblasts and endothelium keep the patient's normal stromal clock state on the core loop, beyond both companion families? | A149 / A151 / A227 single-cell tumour-normal cohorts (on disk) | A149 note 5 form (uncorrected complex correlation, tumour against own normal stroma), 11-gene core loop, 200 matched and 200 responsive core-restricted sets on the same permutations; A283-style placebo calibration; Holm over fibroblasts and endothelium | Pass: the abstract keeps "its stroma carries it" with a registered core-loop basis (14-gene values not blind, core values blind). Fail or at the line: drop the stroma clause; the tumour sentence becomes "tumours stay paired with the patient's time" | On 14 genes the margins are 0.157 (fibroblasts) and 0.174 (endothelium) above pseudo q95 (A149 note 5). The trio's share of cross-class sharing in single-cell atlases was +0.04 to +0.14 (A336 tier 2), so decidable is likely but not certain: power stage at the core projector first | hours (the A227 runs took 21-24 min, peak 29-33 GB) |
| 2 | Is the lesional delay specific to psoriasis, or a property of inflamed skin? | 8 gated AD series (298 pairs) and 9 psoriasis series (465 pairs) in `skin_lesional/`; Tsoi AD (A345); GSE154200 eczema (A183) | per-disease pooled offset (`clocksync.paired`), A308-calibrated companions, 14 genes and core; psoriasis minus AD with an equivalence bound | Pass for AD with psoriasis minus AD not different: "inflamed skin runs late" (abstract changes). Psoriasis minus AD beyond its bound: "largest in psoriasis", with AD stated. AD not delayed: the psoriatic exception stands | AD median -0.45 h (derived, archive). A254 had MDE 0.69 h at 184 persons, so about 0.54 h at 298 pairs (derived, square-root scaling): an AD delay of this size is at the edge. The psoriasis-AD difference of about 1 h is decidable | hours; not blind (A53 per-series offsets seen) |
| 3 | Does age lower clock amplitude in GTEx beyond companions once death type is balanced (the CHIRAL correction of claim 6)? | GTEx v10 and the protected phenotypes (on disk) | old against young within cohort, matched as A231 (sex, Hardy, ischaemic time) plus RIN; host-referenced log z2 at a common spread (amplitude-type, METHOD section 9) and the variance readout with matched / responsive Q*; core loop | Pass: CHIRAL's age dampening stands, and the correction is withdrawn. Equivalence: the correction becomes a registered result | A231 matched 1,861 pairs with d* MDE 0.037; amplitude readouts are noisier (power stage) | hours to a day; it corrects the lab's own Science paper, so the PI decides first |
| 4 | Does the untimed readout see organs move apart when they are known to (two jet lags), and do per-organ drive weights predict the order of re-entrainment? | in-house multi-organ jet-lag data (DATA_ROADMAP item 2) | host-referenced offsets per organ against time since the shift; weights from A322 / A324 registered before the data are opened | Pass: a blind positive control of the time axis, and the only route to time constants. Fail: the "time moves rarely" pattern loses its calibration | depends on the design (not known); A322 / A324 read daytime-feeding shifts of 1.6-11.7 h in timed mouse designs | days to weeks; outward-facing (owner's agreement, PI) |
| 5 | Does the order between organs, and the person factor, hold in living people with liver and arteries? | STARNET (DATA_ROADMAP item 1; dbGaP phs001203): about 600 bypass patients, up to seven tissues | the A293 order statistic, the A331 / A336 leverage on the core loop, chronic-disease contrasts with both families | Pass: the organ order and the person factor become living results; this answers the post-mortem objection. Fail: both stay post mortem | about 600 persons with up to seven tissues each, against 680 women with two tissues in A331 (power 0.947) | weeks; outward-facing (the user applies) |

**Not to add.**
- **Lundell 2020 muscle** (delivered): one tissue per person, so no host-referenced readout. A timed single-tissue
  readout would add a second "time moves" point (as A257), but the published report is that time-restricted feeding
  leaves clock-gene rhythms unchanged (to verify; not carded). An expected null at a small n adds little.
- **GSE186805 oral mucosa** (delivered, 32 % median uniquely mapped): one tissue per person. It would only add a third
  timed lead validation after Perrin and skin. The lead matters mainly to the fetal test, which is uninformative.
- **A260's host-loading reader on GTEx.** A260 separates systemic from local genes only where a clock is removed
  (knockout AUC 0.997); a coupled liver has no readable deviation (kappa 0.001). GTEx tissues sit within 1-3 h of
  their host, a departure not shown to be clock (A259), and GTEx has no clockless contrast. The reader would call
  nearly every rhythmic gene host-following, and it would open a new line after the freeze.
- **GSE53702, critical-illness muscle at two times** (DATA_ROADMAP section 4). It is a living counterpart of the
  intensive-care loss, but it has one tissue (L1 loudness only), biopsy hours that are unknown and may differ between
  the two visits (which changes the variance read), arrays and 33 patients. None of A254's single-disease loudness rows
  was clock-specific alone.
- **More fetal or brain-death analyses on public data.** The fetal organ-offset test is not possible on open data
  (scout 2026-11-22), and the lead route is exhausted. B1-B4 are done, and GSE245976 is one organ after procurement.

## 5. Wording

The term rules are kept: internal desynchrony is a phase offset between compartments; blood is a composition
exception (worded as above); brain death stays descriptive; before birth is "a driven state", not "a fetal clock".

| sentence | issue (entry) | proposed |
|---|---|---|
| abstract l. 10: "yet no method has measured it across a person's organs" | reads as a per-person estimate (A259, A344); framing decision l. 28-30 | "yet it has never been measured between the organs of the same people" |
| abstract l. 12-13: "Before birth, the clock network already moves as a whole in selected organs" | "clock network" is close to "a fetal clock"; the core loop is clock-specific in one blind cohort only (A328 part a; A266, A269 not) | "Before birth, the clock genes of selected organs already move together in the clock's pattern, a driven state whose source is open" (subject to Q43) |
| abstract l. 13-14: "organs keep a common body time through chronic disease, intensive care and cancer" | intensive care: order kept as a descriptive readout after a failed primary, offsets compressed (A109, A132, A293); cancer: pairing at the line on the core (A336) | "organs keep their order through chronic disease and intensive care, and tumours stay paired with their patient's time" |
| abstract l. 14: "What disease changes is how loudly a tissue expresses that time, not the time itself" | a rule from three conditions; chronic disease changes neither (A306) | "Tumours, chronic lesions and intensive care change how loudly tissues express that time, not the time itself; chronic systemic disease changes neither" |
| abstract l. 15-16: "the surrounding stroma carries it, and cancer cells stay coupled but grow quieter as they proliferate" | stroma post hoc (A149 note 5, A226); malignant coupling failed (A138 H138.1); A304 is bulk | "as tumours proliferate, they express less of that time without drifting from it" (stroma clause only after addition 1) |
| abstract l. 16-17: "Inflamed psoriatic skin is an exception, running about 1.5 h late, replicated in an independent cohort" | AD also late (A345, A53); core loop -1.00 h (A335) | "Inflamed skin is the one local exception: psoriatic lesions run about 1 to 1.5 h behind the same person's uninvolved skin, replicated in an independent cohort and already in pre-mRNA, and atopic lesions lean the same way" |
| abstract l. 17-18: "brain death rearranges the order between organs, the only state in which the shared time itself comes apart" | H317.2 (shared state kept); psoriasis is a local offset; aliased (Q41) | "in brain-dead organ donors the organs keep a shared clock but in a different order, the only body-wide rearrangement we found, and one this resource cannot separate from organ procurement" |
| limits table, tumour row: "A271 (fibroblasts and endothelium beyond both nulls)" | descriptive anchors; post hoc A149 note 5; A226 | "stromal pairing beyond matched genes: post hoc (A149 note 5), not confirmed (A226); A271 anchors descriptive" |
| limits table, psoriasis row | add A345 (H345.1, H345.3 narrow; H345.2 not shown) and the AD rows | "psoriasis -1.45 h (A183, blind), at the line on the core loop (q 0.049, A335); pre-mRNA delay clock-specific narrowly in the discovery series (A345); AD lesions also late (A345, A53)" |
| plan l. 32-33 and claim 5: "The clock layer is built around birth. Before birth, adrenal cortex and lung carry the signal layer in the clock's pattern without an adult-size clock" | A328 part a (blind core loop) | paused (Q43); at least "an adult-size core loop is excluded in single-cell cortex (A266) but the core loop covaries on its own in whole adrenals (A328 part a)" |
| plan l. 91-92 and 100-101 (departure and person leverage "ride on" PER1, PER2, NR1D1) | A336 gene-drop null | "on the core loop the departure excess falls to 1.28 times matched (A329), as much as dropping three random genes of similar weight (A336): it is not shown to ride on PER1, PER2 and NR1D1" |
| plan l. 131: "lymphocyte content leaves relative gain unchanged" | A336: A295 H295.1 PASS on the core | "lymphocyte content dilutes relative gain slightly: -0.020 per SD on the core loop (p 0.040; A336), within the post hoc ±0.10 margin (A335)" |
| plan section 3 rows and section 10: "amplitude about 0.82 / 1.02 / 0.81 / 1.28" | METHOD section 9 | "loudness about 0.82 ..."; section 10 "variance to loudness"; draw the intensive-care and lesion arrows on a loudness scale, apart from the timed amplitude ratios (A257, mouse) |
| plan section 3, intensive-care row: "order kept (Fisher-Lee rho 0.690)" | A109 H109.2, A132b | "order kept (rho 0.690; a registered re-test, A293, of the descriptive A132 readout; the registered primary H109.2 failed in pooled ventilated donors), offsets compressed toward 0 / 12 h (0.58 of fast deaths)" |
| plan section 3, tumour row: "kept: tumours pair with their patient" | A336 | add "at the line on the core loop (0.252 against matched q95 0.275)" |
| plan section 3, tumour cell classes: "fibroblasts and endothelium carry the patient's time" | A149 note 5 (post hoc), A226, A271 anchors (descriptive) | "fibroblasts and endothelium pair with the patient's normal stroma beyond matched genes (post hoc, not confirmed)" |
| plan l. 77: "every tissue follows the person's amplitude nearly in proportion" | A312 radii carry a donor-level scale that matched genes share in deaths (A336 part 2: matched median 0.119) | "the person's loudness" |
| plan l. 245-246: CHIRAL age correction | WP0 item 5 | "CHIRAL's rhythmic-gene excess in the young shrinks under death-type matching (WP0 item 5, timed count, not registered; TMM p 0.092); coupling and offsets are equal with age (A231)", or drop until addition 3 |
| plan l. 247: "Prognostic clock signatures reduce to proliferation" | A192-A194, A220 | "published clock prognostic signatures add nothing beyond proliferation or random gene sets (A192-A194, A220)" |
| RP A257 header: "falls to about a quarter" | correction note in the same entry (l. 4788) and plan section 8 | "to about a third (14 genes); about half on the core loop (A335)" |

**Understated.**
- **A251 and A255** are absent, though they are the cleanest living evidence for "one time per person".
- **A345** H345.3 is the first pass beyond both families for the psoriatic delay (at pre-mRNA, narrowly). It belongs
  in the text, not yet in the abstract.
- **A343** is the human, interventional demonstration of the "fast signal layer through PER1" in the model (plan
  claim 4).

## 6. For the PI

1. **Inflamed skin.** Should the abstract's exception read "inflamed skin, largest in psoriasis" rather than
   "psoriatic skin", given A345 (AD -1.15 h), the archive AD series and A254's own reading? Register addition 2 first?
2. **Tumour stroma.** Keep the stroma clause only after a registered core-loop re-test (addition 1), or drop it now?
3. **Before birth (Q43).** The plan's term is "a driven state"; the abstract says "the clock network moves as a whole";
   A328 part a supports a covarying core loop in one blind cohort. Which wording?
4. **Brain death (Q41).** The abstract's closing sentence is the climax position. Keep the lifespan order with the
   reworded sentence of section 5, or close the abstract on the main claim?
5. **Intensive care.** Should the abstract keep "through intensive care" once these are stated in the text?
   - its origin: a descriptive readout after H109.2 failed (A132, A293);
   - the compression toward 0 / 12 h;
   - a detection limit simulated on fast deaths.

   Should that limit be recomputed at the ventilated design (power only) before Figure 5a draws it?
6. **Figure 3.** The review's condition for a main figure (S1 scored blind with a timed control) is met formally
   (A328 blind; A337, A342 controls), but the outcome is uninformative. Does that satisfy the condition?
7. **CHIRAL age correction.** Rescore it (addition 3), or remove it from the corrections? It concerns the lab's own
   Science paper.
8. **A251 and A255.** Promote them to a main figure (Figure 1 or 4) as the living-person evidence?
9. **Outward-facing additions.** Ask the owner of the jet-lag data (addition 4)? Apply for STARNET (addition 5)?
10. **A342.** In the registry table, label it "rescored after the registered run was seen (count definition
    corrected)"?

## 7. Addendum after the user's review (2026-11-23)

The user's questions on sections 1-2, and what the record supports. No new analysis; numbers as in sections 1-3.

- **Stroma and cancer cells.**
  - In healthy tissue, stromal-vascular cells carry the clock at high gain: fibroblasts about 1.26 times endothelium
    (A293 H293.4), above lymphocytes (A272-A283).
  - In tumours, the stromal pairing with the patient's normal stroma is post hoc (A149 note 5), unconfirmed (A226),
    and not on the core loop; tumour-stroma gain is not measured.
  - For cancer cells the evidence is indirect: the gain loss is not explained by stromal content. The proliferative
    arm has more stroma, and the purity-stratified loss remains (-0.266, p 0.037; A304).
  - Malignant host coupling was tested and failed (A138 H138.1).
  - Wording: "tumours express less of the patient's time as they proliferate; the loss is not explained by stromal
    content", not "cancer cells stay coupled".
- **Intensive care.** The organs share a clock state per patient (H109.1, excess 0.310) and keep the fast-death order
  (rho 0.690). About half of the shared covariance is a displacement along one axis, mostly glucocorticoid (A109:
  fraction 0.44). This "running clock plus a shared non-temporal factor" reading explains the compression of the
  offsets (0.58 of fast deaths; A132b). The fetal elongation test (A340) is a different analysis. Wording: "same
  order, offsets compressed by a shared non-temporal state". The caveats of section 3.3 item 3 stand.
- **Brain death.** "The healthy order between organs is lost, while the organs still share a clock state" (A317
  H317.1-2), with the procurement alias.
- **Inflamed skin.** The main idea (changes of loudness are more common than shifts of time) is unchanged. What can be
  claimed:
  - psoriatic lesions run late against the same person's skin: replicated blind (A183, -1.45 h), in held-out series
    (A54, -0.93 h) and in pre-mRNA (A345); at the line on the core loop;
  - atopic lesions lean the same way, smaller and less consistent (descriptive);
  - non-skin chronic lesions keep their time (A254, pooled +0.01 h, MDE 0.69 h).

  Psoriasis-specificity cannot be claimed. Wording: "phase shifts are rare and local: among the lesions tested, only
  in inflamed skin, largest in psoriasis; elsewhere lesions lose loudness without a shift". Addition 2 would make the
  atopic part a test.
- **Before birth.** The plan's "signal layer before the clock" is outdated. The core loop covarying on its own is
  shown in one blind cohort (A328 part a). The single-cell cortex (A266), lung and liver (A269) do not show it, and
  driven against running stays open (Q43).
- **PER1, PER2 and NR1D1.**
  - In the bulk results, the core-loop falls are what dropping any three genes of similar weight does (departure:
    32nd percentile; living leverage: 43rd; A336). No statistic is "carried by the trio". The glucocorticoid-bias
    concern is answered there: "not shown to ride on them", not "disproven".
  - In single-cell atlases the trio does carry an extra share of cross-class sharing (1st-7th percentile; 4th-29th
    among similar drops). Those contrasts hold against matched genes, and against responsive genes in half of the
    tests.
- **Person loudness.** Hospital hours are not the limit: MuTHER is untimed clinic sampling, and A331 passes beyond
  matched genes on 14 genes. The limits are these:
  - two tissues from one punch biopsy on one array batch, so shared handling is not separated (not beyond
    responsive genes);
  - a fall to the line on the core loop, likely precision.

  In deaths (27 tissues) it holds on the core loop, beside a large donor-level non-clock scale. More tissues per
  living person (STARNET, addition 5) would decide it.

## 8. Follow-up results (2026-11-23; proposed wording, not applied)

Prompt `prompts/audit_followups_2026-11-23.md`. Each amendment was committed before scoring; numbers from RP.

| gap / addition | result | wording it allows |
|---|---|---|
| Gap 1 / addition 1, tumour stroma (A346) | Fibroblasts: core 0.683 (14 genes 0.686), beyond matched, broad and narrow responsive (narrow descriptive: not decidable at the power stage); endothelium falls on the core (0.510 < matched q95 0.529; trio share 0.17). 14 genes: both pass the calibrated rule (narrow at the line) | "tumour fibroblasts keep the patient's normal fibroblast clock state, clock-specifically, also without PER1, PER2 and NR1D1 (in the discovery cohorts)"; drop endothelium from the clause; "cancer cells" stays out (D346 epithelium not clock-specific) |
| Addition 2, inflamed skin (A347) | Psoriasis -1.51 h / core -0.96 h; AD -0.26 h on 14 genes (p 0.034, not clock-specific), -0.21 h on the core (p 0.097, MDE 0.37 h); psoriasis minus AD -1.25 h / -0.75 h, both intervals below 0; hidradenitis -1.07 h, not clock-specific (descriptive) | the registered reading keeps the psoriatic exception: "psoriatic lesions run about 1 to 1.5 h late (core loop about 1 h), more than atopic lesions, which lean late by a fraction of an hour without a core-loop or clock-specific delay"; replaces section 7's "inflamed skin, largest in psoriasis" |
| Gap 2, living organ sharing (A348) | Core loop PASS in A251 H251.1, H251.2a-b and A255 H255.3 (beyond matched and broad; surgical pairs also beyond narrow); MuTHER fat-skin at the line (matched p_cal 0.057) | "tissues of the same living person share one clock state, also without PER1, PER2 and NR1D1 (skin sites and layers, surgical fat-skin pairs)"; A251 / A255 can enter Figure 1 or 4 |
| Gap 3, intensive care (A349) | Ventilated detection limit: one organ 6.75 h (core 5.6 h), dispersion SD 2.25 h (not 3.9 h / 1.5 h); host spread about half of fast deaths; host-referenced log z2 lower in most tissues, the loss scales with clock strength (rho 0.82); compression 0.58 -> about 0.77 after removing the non-temporal axis | "the order is kept within a detection limit of about 6-7 h for one organ"; "offsets compressed, partly by a shared non-temporal state"; the loudness loss is host-referenced and largest in the strongest clocks (descriptive) |
| Gap 7, class gains (A350) | Core-loop tau within 0.35 in four of five families; mural 0.389 (FAIL; SE x2.4) | "class gains are conserved across organs (on the core loop in four of five classes, not mural cells)" |
| Addition 3, CHIRAL (Task 6) | drafted only (`prompts/chiral_age_rescore_DRAFT_2026-11-23.md`) | none until the user and PI decide |
| Gap 8, provenance | A336 tier-2 scripts and table committed (Task 0) | - |
