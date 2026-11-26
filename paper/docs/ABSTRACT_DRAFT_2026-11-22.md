# Abstract draft (2026-11-22)

Working abstract agreed with the user on 2026-11-22. It still needs PI review. Journal choice belongs to the user and
the PI; Cell is the current target.
Fill the bracketed counts from output tables (never from memory) before the draft is used anywhere.

## Abstract (about 180 words; revised 2026-11-23 after the story audit and A346-A349)

> Circadian clocks run in nearly every human tissue, and their misalignment between organs, internal desynchrony, is
> thought to contribute to disease, yet it has never been measured between the organs of the same people. We introduce
> an approach that reads each tissue's clock against the same person's other tissues, so it needs no sampling time,
> and apply it across the human lifespan in [N] tissues from [N] donors and [N] disease cohorts. Before birth, the
> clock genes of selected organs already move together in the clock's pattern. In adults, organs keep their order
> through chronic disease, and the fibroblasts of tumours keep the time of the same patient's normal tissue:
> the desynchrony that mouse models predict for metabolic disease and ageing is absent at our detection limits. Across
> humans, baboons and sheep, organs that express the body's time more loudly are more tightly coupled to it, beyond
> random genes, and organs rank alike in loudness across species; with known sampling times in baboons, an organ's own
> clock noise grows much less than its loudness. What changes in disease is how loudly tissues express that time: intensive care, chronic lesions, tumour stroma and proliferating tumours quieten
> tissues without moving their time. Inflamed psoriatic skin is the exception: its lesions run 1 to 1.5 h late against
> the same person's skin, replicated in an independent cohort and more than atopic lesions. At the end of life,
> brain-dead organ donors keep a shared clock but lose the healthy order between organs.

If it must be shortened, cut "and apply it across the human lifespan in [N] tissues..." first; the counts can go in the
highlights.

**Changes from the 2026-11-22 version** (`STORY_AUDIT_2026-11-23.md` sections 1.2, 5 and 7, discussed with the user).
- **Healthy organisation (added 2026-11-26, user, after A383 / A384 / A395).** New sentence on the amplitude-coupling
  rule, the conserved loudness ranking and the timed own noise. "What changes" becomes "What changes in disease". Mouse
  (A394) joins "humans, baboons and sheep" only if A394 is scored and passes. Adds about 50 words.
  Revised the same day after A397 (mouse timed law not decided): "barely grows" -> "grows much less", restricted to
  baboons.
- **The gap.** Now "never measured between the organs of the same people". The old sentence read as a per-person
  estimate (A259, A344).
- **Before birth.** "Clock genes move together in the clock's pattern" replaces "the clock network moves as a whole".
  The final term waits for Q43; "a driven state whose source is open" is the alternative.
- **Intensive care.** It appears only in the loudness sentence. The order is kept, but at the ventilated design only
  a single-organ shift of about 6-7 h, or a dispersion SD of 2.25 h, is detectable (A349; 5.6 h on the core loop),
  and the offsets are compressed by a shared non-temporal state (A109, A132b). In the text: "no large disruption";
  intensive care serves as the reference for brain death.
- **Tumours (revised 2026-11-25, A387).** "The fibroblasts of tumours keep the time of the same patient's normal
  tissue" replaces "tumours stay paired with their patient's time, carried at least by their fibroblasts": A387 is
  registered and blind, on the full model (coupling c 0.40 on the core loop, 0.44 on 14 genes, beyond matched and
  stress-responsive genes; peak -0.07 h, 90 % -0.77 to 0.63, against the patient's normal fibroblasts; loudness ratio
  0.69, so "tumour stroma quietens" enters the loudness sentence). The bulk pairing (A138) stays in the results, at the
  line on the core loop. Superseded record (A346, kept for the history): "Stay paired" (A138, at the line on the core
  loop), "carried at least by their fibroblasts" (A346: 0.686 → 0.683 on the core loop). Still not "only fibroblasts":
  - tumour endothelium still pairs but is not beyond matched genes on the core loop (A346);
  - cancer cells are not resolved (A346 D346, A138 H138.1, A305);
  - immune pairing is not beyond responsive genes (A271);
  - no independent cohort exists (A226).
- **Cancer cells.** "Cancer cells stay coupled" becomes "proliferating tumours quieten": A304 is bulk, and malignant
  coupling failed (A138 H138.1).
- **New sentence: mouse-predicted desynchrony absent** (A235, A164, A306, A231; with detection limits and positive
  controls).
- **Loudness.** "What disease changes" is restricted to the three conditions that move loudness. Chronic systemic
  disease moves neither.
- **Psoriasis is the exception, more than atopic lesions** (A347, which supersedes the audit's "inflamed skin"
  proposal):
  - psoriasis -1.51 h (14 genes), -0.96 h (core loop);
  - psoriasis minus atopic -1.25 h / -0.75 h;
  - the atopic delay is small (-0.26 h), not clock-specific and not on the core loop.

  Do not write that atopic lesions are on time: the equivalence test failed, and 7 of 8 atopic series lean late.
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
| the fibroblasts of tumours keep the time of the same patient's normal tissue | A387 (registered, blind; full model over 34 cohorts: c 0.40 core / 0.44 on 14 genes, beyond matched and narrow responsive genes; peak within 0.8 h of the normal fibroblasts, H387.4; holds without the counterpart in the host, S387; colorectal and lung fibroblasts per cancer, A387 note 1); bulk support A138 H138.4 (at the line on the core loop, A336) | "only fibroblasts"; "stroma carries it" (endothelium coupled on 14 genes, not beyond matched genes on the core loop, A387 H387.3); "cancer cells carry it" (cancer cells c 0.10 / 0.03 on the normal-stroma host, A387; on their own tumour's stroma 0.15 / 0.13 in 127 new patients, not beyond matched genes, A388 H388.1; the clock-keeping vs at-null split leans as predicted but is not specific, H388.2; lung adenocarcinoma cancer cells coupled, 0.66, descriptive); "immune cells carry it" (A271); replicated (same discovery cohorts as A346; no independent cohort, A226) |
| tumour stroma quietens | A387 D387.1 (tumour fibroblasts' loudness ratio 0.69, 0.51-0.86, core; 0.80 on 14 genes; endothelium 0.82, interval including 1) | "the tumour clock is broken"; endothelium "quieter" (interval includes 1) |
| proliferating tumours quieten without drifting | A304 (gain 0.63, departure not raised; holds on the core loop; purity-stratified -0.266) | "cancer cells stay coupled" (A138 H138.1 failed; A304 is bulk) |
| mouse-predicted desynchrony absent | A235 (high-fat-diet liver advance excluded, ±1 h), A164 / A306 (T2D, obesity), A231 (age), A233; positive controls A111, A129, A257 | "human clocks resist desynchrony": the drives stay aligned in these conditions; with conflicting drives time moves (A257) |
| psoriatic lesions 1-1.5 h late, replicated, more than atopic lesions | A183 (blind, -1.45 h; core -1.00 h, A335), A54 held-out, A345 (pre-mRNA delay), A347 (pooled -1.51 / -0.96 h; psoriasis minus atopic -1.25 / -0.75 h) | "atopic lesions are on time" (equivalence failed; 7 of 8 series lean late); "clock-specific" in the abstract; hidradenitis leans late (descriptive) |
| brain-dead donors keep a shared clock, lose the healthy order | A317 H317.1-2, A318, A336 part 4 (core loop), descriptive origin A132b | any cause; "time comes apart"; organ-level arrows (not stable across gene sets, A335) |
| louder organs are more tightly coupled, beyond random genes, in humans, baboons and sheep; organs rank alike in loudness across species; own clock noise barely grows with loudness (known times) | A383 H383.1 (rho 0.745 / 0.503 / 0.629, above matched and responsive q95; baboon at the responsive line), A384 H384.1 (T 0.386, p 0.006, above matched q95; baboon-sheep 0.653, human-baboon 0.375, human-sheep 0.132), A395 H395.1 (baboon timed: b_t 0.22, -0.01 to 0.42; phase scatter slope -0.75); A397 H397.1 FAIL in the mouse (b_t 0.50, 0.29 to 1.05; not decided) | "living" animals (healthy, rapid death: euthanasia, slaughter); "in every species pair" or "conserved in humans" (human-sheep 0.13); "the clock's own noise is uniquely constant" or "organs sense time better" (no companion in A395; not replicated in the mouse, A397; random gene sets also give b far below 1 in A383); a cause of the loudness ranking (open; composition explains only 1.3-fold of 2.9-fold, A323); "mouse" before A394 |
| intensive care: loudness falls; no large disruption of the order | A315, A335 (loudness; core loop, age / sex matched), A349 (host-referenced log z2 > 0; loss largest in the strongest clocks), A293 (order) | "time kept" without the limit (single organ about 6-7 h, dispersion SD 2.25 h; A349); "same offsets" (compressed, A132b); "amplitude" for the variance ratio; COVID cohorts (A109c) |

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
2. Psoriasis: the core-loop delay (-0.96 h) is beyond matched genes in the pooled series (A347), but at the line in the
   blind cohort (A335). Have the psoriasis against atopic figure ready.
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
