# Abstract draft (2026-11-22)

Working abstract agreed with the user on 2026-11-22. It still needs PI review. Journal choice belongs to the user and
the PI; Cell is the current target.
Fill the bracketed counts from output tables (never from memory) before the draft is used anywhere.

## Abstract (about 180 words; revised 2026-11-23 after the story audit)

> Circadian clocks run in nearly every human tissue, and their misalignment between organs, internal desynchrony, is
> thought to contribute to disease, yet it has never been measured between the organs of the same people. We introduce
> an approach that reads each tissue's clock against the same person's other tissues, so it needs no sampling time,
> and apply it across the human lifespan in [N] tissues from [N] donors and [N] disease cohorts. Before birth, the
> clock genes of selected organs already move together in the clock's pattern. In adults, organs keep their order
> through chronic disease and intensive care, and tumours stay paired with their patient's time: the desynchrony that
> mouse models predict for metabolic disease and ageing is absent at our detection limits. What changes is how loudly
> tissues express that time: intensive care, chronic lesions and proliferating tumours quieten tissues without moving
> their time. Time moves locally in inflamed skin, most in psoriasis, whose lesions run 1 to 1.5 h late, replicated in
> an independent cohort. At the end of life, brain-dead organ donors keep a shared clock but lose the healthy order
> between organs.

If it must be shortened, cut "and apply it across the human lifespan in [N] tissues..." first; the counts can go in the
highlights.

**Changes from the 2026-11-22 version** (`STORY_AUDIT_2026-11-23.md` sections 1.2, 5 and 7, discussed with the user).
- **The gap.** Now "never measured between the organs of the same people". The old sentence read as a per-person
  estimate (A259, A344).
- **Before birth.** "Clock genes move together in the clock's pattern" replaces "the clock network moves as a whole".
  The final term waits for Q43; "a driven state whose source is open" is the alternative.
- **Intensive care.** "Keep their order": the order is kept and the offsets are compressed by a shared non-temporal
  state (A109, A132b).
- **Tumours.** "Stay paired" (A138, at the line on the core loop). The stroma clause is dropped until the core-loop
  re-test (`prompts/audit_followups_2026-11-23.md`, task 1).
- **Cancer cells.** "Cancer cells stay coupled" becomes "proliferating tumours quieten": A304 is bulk, and malignant
  coupling failed (A138 H138.1).
- **New sentence: mouse-predicted desynchrony absent** (A235, A164, A306, A231; with detection limits and positive
  controls).
- **Loudness.** "What disease changes" is restricted to the three conditions that move loudness. Chronic systemic
  disease moves neither.
- **Inflamed skin, not psoriasis alone.** Atopic lesions lean the same way (A345, archive A53), and A254 reads the
  shift as specific to inflamed skin. The size is given as 1 to 1.5 h, matching the core-loop values.
- **Brain death.** "Keep a shared clock but lose the healthy order" (A317 H317.1-2). It replaces "the only state in
  which the shared time comes apart" and keeps brain death from reading as the climax (Q41).

## Framing decisions (user, 2026-11-22)

- **Plain language.** A Cell abstract must read for non-specialists. Keep it conceptual and leave the details to the
  introduction.
- **Lifespan order:** method, then before birth, then adult life with disease, then brain death at the end of life.
- **The gap is internal desynchrony measured across one person's organs, not missing timestamps.**
  CYCLOPS, CHIRAL and TEMPO already order untimed samples; BodyTime, TimeSignature and hair-follicle assays give
  single-tissue phases. Do not write "no method measures it at the individual level": it is too broad.
- **The method's practical edge belongs in the introduction, not the abstract.** It works on partial-cycle,
  hospital-hours data (surgical and biopsy cohorts, i.e. most disease and tumour data; GTEx deaths cover the full day).
- **Main claim:** "loudness, not time". State it plainly in the adult section.

## Wording limits (what the record supports)

| sentence | evidence | do not write |
|---|---|---|
| fetal clock genes "move together in the clock's pattern" in selected organs | A328 part a (blind; core loop clock-specific in del Valle), A340 (H340.1), organ selectivity in two species | "running", "24 h", "fetal clock", "maternal" / glucocorticoid entrainment (A338 one-sided, conditional on lifetime; A341 descriptive); the core loop is not shown in the cortex, lung or liver (A266, A269) |
| tumours stay paired with the patient's time | A138 H138.4 (beyond both families on 14 genes; at the line on the core loop, A336) | "follow" without the core-loop caveat; "stroma carries it" (post hoc A149 note 5, not confirmed A226, no core run); "immune cells carry it" (A271) |
| proliferating tumours quieten without drifting | A304 (gain 0.63, departure not raised; holds on the core loop; purity-stratified -0.266) | "cancer cells stay coupled" (A138 H138.1 failed; A304 is bulk) |
| mouse-predicted desynchrony absent | A235 (high-fat-diet liver advance excluded, ±1 h), A164 / A306 (T2D, obesity), A231 (age), A233; positive controls A111, A129, A257 | "human clocks resist desynchrony": the drives stay aligned in these conditions; with conflicting drives time moves (A257) |
| inflamed skin late, most in psoriasis, 1-1.5 h, replicated | A183 (blind, -1.45 h; core -1.00 h, A335), A54 held-out, A345 (pre-mRNA delay; H345.3 narrow); AD A345 and archive A53 (descriptive) | "psoriasis-specific" (A183 H183.2); "clock-specific" in the abstract (at the line on the core loop) |
| brain-dead donors keep a shared clock, lose the healthy order | A317 H317.1-2, A318, A336 part 4 (core loop), descriptive origin A132b | any cause; "time comes apart"; organ-level arrows (not stable across gene sets, A335) |
| intensive care: loudness falls, order kept | A293, A315, A335 (GTEx ventilated) | "amplitude" (variance ratios measure loudness, METHOD section 9); "same offsets" (compressed, A132b); COVID cohorts as evidence (A109c uninformative) |

## Introduction outline (agreed)

1. Why internal synchrony matters: mouse phase relationships between organs; misalignment (shift work, jet lag) and
   disease; only indirect evidence in humans.
2. What existing methods do: population phase ordering of untimed samples (CYCLOPS, CHIRAL, TEMPO); credit the GTEx
   consensus.
3. Why that is not enough:
   - an inferred phase mixes the person's time, the sampling hour and the tissue offset;
   - surgical cohorts cover only part of the cycle;
   - comparing organs within a person needs a reading that does not depend on the hour.
   Evidence: simulations at half a cycle, and agreement with GTEx where both approaches apply.
4. What we do: read each tissue against the same person's other tissues, so shared offsets cancel. A model-free gate
   confirms the shared signal is clock-shaped.
5. What we find: "loudness, not time", leading into the chapters in lifespan order.
6. What it means (added 2026-11-23). The mechanism the mouse work established, drives setting time, holds in human
   tissue. The clock's inertia carries time through acute states, and drive strength, phase coherence and composition
   set loudness.
   - That mechanism predicts robust human time under chronic disease and ageing, where the drives stay aligned. It
     predicts shifts only where drives conflict (night shift, A257) or locally (inflamed skin).
   - The desynchrony extrapolated from mouse models and cultured cells (high-fat feeding, metabolic disease, ageing)
     is therefore absent at powered limits.
   - The field's usual synchrony metrics read mostly loudness (A291, A294b).

## Risks at review (ranked; updated 2026-11-23)

1. Brain death: open confound. State it plainly and show the dissociation between order and loudness.
2. Inflamed skin: the psoriatic delay is at the line on the core loop, and the atopic part is descriptive. Have that
   figure ready (audit follow-up task 2 decides the atopic part).
3. Fetal evidence comes from single-organ measures, not organ-to-organ offsets.
4. **Mechanism.** It is verified at the systems level in human, and the mouse consensus is credited. The cell-level
   cause of lost loudness (dephasing or per-cell damping) is not decidable with public data (A334). Planned: the
   waveform screen and the cytokine dose-response in skin (`prompts/mechanisms_2026-11-23.md`); in-house single-cell
   time-lapse (PI).
5. **Reanalysis only.** The answer is the record:
   - pre-registration, with 70 blind and 42 partly blind of 219 registered tests (registry 2026-11-19);
   - blind replications on unseen data (A256, A274, A280, A183, A328 part a, A321, A337 Perrin);
   - no new human or animal sampling, which is an ethical strength to state, not to apologise for.

   A prospective blind test on the in-house jet-lag data would still strengthen it.

Possible next step, decided by the user and the PI: a Cell presubmission inquiry (this abstract plus a paragraph on the
cancer result); a Cell Press transfer keeps Cell Genomics open.
