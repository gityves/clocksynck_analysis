# Abstract draft (2026-11-22)

Working abstract agreed with the user on 2026-11-22. It still needs PI review. Journal choice belongs to the user and
the PI; Cell is the current target.
Fill the bracketed counts from output tables (never from memory) before the draft is used anywhere.

## Abstract (about 170 words)

> Circadian clocks run in nearly every human tissue, and their misalignment between organs, internal desynchrony, is
> thought to contribute to disease, yet no method has measured it across a person's organs. We introduce an approach
> that reads each tissue's clock against the same person's other tissues, so it needs no sampling time, and apply it
> across the human lifespan in [N] tissues from [N] donors and [N] disease cohorts. Before birth, the clock network
> already moves as a whole in selected organs. In adults, organs keep a common body time through chronic disease,
> intensive care and cancer. What disease changes is how loudly a tissue expresses that time, not the time itself.
> Tumours follow their own patient's clock: the surrounding stroma carries it, and cancer cells stay coupled but grow
> quieter as they proliferate. Inflamed psoriatic skin is an exception, running about 1.5 h late, replicated in an
> independent cohort. At the end of life, brain death rearranges the order between organs, the only state in which the
> shared time itself comes apart.

If it must be shortened, cut "and apply it across the human lifespan in [N] tissues..." first; the counts can go in the
highlights.

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
| fetal "moves as a whole in selected organs" | A340 (H340.1 pass), gate, organ selectivity in two species | "running", "24 h", "maternal" / glucocorticoid entrainment (A338 one-sided, conditional on lifetime; A341 descriptive) |
| tumours follow the patient's clock; stroma carries it | A138, A227, A271 (fibroblasts and endothelium beyond both nulls), A304 | "immune cells carry it" (A271: not beyond the stress-responsive null) |
| psoriasis about 1.5 h late, replicated | A183, blind (-1.45 h, 81 persons) | "clock-specific" (A183 matched null at the line; A308 not blind; A335 at the line) |
| brain death rearranges order | A132b, A314-A317 (descriptive) | any cause; "only" is safe only for the conditions tested |
| intensive care: loudness falls, order kept | A293, A335 (GTEx ventilated) | COVID cohorts as evidence (A109c uninformative) |

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

## Risks at review (ranked)

1. Brain death: open confound; state it plainly and show the dissociation between order and amplitude.
2. Psoriasis clock-specificity is at the line; have that figure ready.
3. Fetal evidence comes from single-organ measures, not organ-to-organ offsets.
4. Reanalysis only. A small prospective test (blind run on the in-house jet-lag data) would help most.

Possible next step, decided by the user and the PI: a Cell presubmission inquiry (this abstract plus a paragraph on the
cancer result); a Cell Press transfer keeps Cell Genomics open.
