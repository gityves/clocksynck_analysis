# DRAFT amendment (not registered, not run): CHIRAL age rescore

Written 2026-11-23 for audit addition 3 (`STORY_AUDIT_2026-11-23.md` sections 1.2, 4 and 6 item 7), prompt
`audit_followups_2026-11-23.md` Task 6. It bears on the lab's own published age result (CHIRAL, Science), so the user and
the PI decide whether it is registered and run. Nothing below has been computed. To register: copy into
`registered_experiments.md` with the next free number, date it, commit, then run the power stage.

## A3xx (DRAFT; staged: power, then one scoring run; **NOT BLIND**: WP0 item 5 (CHIRAL timed rhythmic-gene count, old minus young, pooled -0.73 against Hardy-matched -0.13, CPM; TMM matched minus pooled +0.39, p 0.092; not registered, no companions), A208 / A231 (coupling and offsets equal with age; A231 matched 1,861 pairs, d* MDE 0.037), A315 / A335 (loudness readout and its companions); **BLIND**: every host-referenced log z2 and every variance readout by age group, pooled or matched, on either gene set): does age lower clock amplitude in GTEx beyond companions once death type and RNA quality are balanced? **GTEx v10 with protected phenotypes; old against young within cohort**

**Why.** CHIRAL reported dampened rhythms in older donors. Hardy (death type) changes rhythm strength and is unbalanced
with age; WP0 item 5 found the age effect on the timed rhythmic-gene count shrinking under Hardy matching, but that
check is not registered and has no companions. The paper plan cites the correction (claim 6) with IDs that do not test
it (A208, A231).

**Design.** Groups as CHIRAL (old against young; the CHIRAL age cut, to be fixed from the paper before registration),
within cohort (fast deaths, ventilated post-mortem; organ donors separately), matched as A231 (sex, Hardy, ischaemic
time) plus RIN (caliper 0.5); balance table (SMD) before and after matching, pooled and matched estimates both reported.
Tissues: the A231 set.

**Readouts** (14 genes and core loop, `a336_common.core`; companions 200 matched, 200 narrow responsive, 200 broad
responsive, restricted to the core positions for the core readout):
1. Host-referenced log z2 ratio old over young per tissue at a common spread (`Dataset.contrast`, `log_z2_ratio`,
   `d_star`), median over tissues; donor bootstrap 90 %.
2. Loudness: the A315 debiased clock-plane variance ratio, Q = median over tissues of log(W_old / W_young), with Q* =
   Q minus the matched companions' median, and the responsive Q*.

**Hypotheses (draft).**
- **H.1 (primary, matched, core loop):** amplitude lower in old donors: log z2 ratio median < 0 with 90 % interval below
  0, and Q* below the matched q05 and the responsive q05.
- **H.2 (equivalence):** 90 % interval of the matched log z2 median inside +-0.10 (about 10 % in amplitude; margin to
  be set from the power stage and fixed before scoring).
- **D (descriptive):** pooled against matched difference; 14 genes against core; the CHIRAL rhythmic-gene count on
  the matched donors for comparison with WP0 item 5.

**Readings (draft).** H.1 passes on the matched core readout: CHIRAL's age dampening stands, and the paper's correction
is withdrawn. H.2 passes: the age effect is not detected once death type and RIN are balanced, and the correction
becomes a registered result (with its companions). Neither: undecidable at this design; the correction is dropped from
the paper.

**Power first.** Swap null of age labels within matched strata: MDE of the log z2 median and of Q*; H.2 decidable if
the expected 90 % half-width is below the margin.

**Limits.** GTEx sampling hour per donor (time of death) is controlled access and stays aggregated; matched designs
lose the oldest ventilated donors; age and comorbidity are entangled (T2D SMD 0.46-0.79 in A231), stated.

**For the user / PI before registering.** (1) Run it at all, given it corrects the lab's own paper? (2) The age cut and
whether organ donors are included. (3) The equivalence margin. (4) Whether CHIRAL's authors should see the draft first.
