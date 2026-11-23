# Prompt for the next agent: the analyses missing after the story audit (written 2026-11-23)

You close the analysis gaps found by `paper/docs/STORY_AUDIT_2026-11-23.md` (sections 3.3 and 4, and the section 7
addendum). Every task is a new hypothesis or a new readout, so each is an amendment in
`paper/docs/registered_experiments.md`, written and committed **before** its data are scored. The amendment states its
blind status, its power stage and its stop rule (one scoring run). Results go to `paper/docs/results_preliminary.md`,
with a VERDICT_TABLE row.

Start with `CLAUDE.md`, the memory index, the audit (sections 1.2, 3 and 7), METHOD.md sections 2-4 and 9, and the
results entries named under each task.

## Rules for this round
- **Read-only documents.** Do not edit the plan, the abstract, PI_QUESTIONS, or earlier results entries; report new
  results only. Phase 8 waits for the PI.
- **Two readouts.** Every statistic is computed on 14 genes and on the core loop (no PER1, PER2, NR1D1;
  core-restricted companions), as in A336.
- **Companions.**
  - 200 expression-matched and 200 stress-responsive sets, on the same permutations or bootstraps as the clock.
  - Use a calibrated placebo form (A283 / A285 / A320) where the design has few pairs.
  - Report the narrow and the broad responsive families.
- **Statistics.** CLAUDE.md applies:
  - a failed gate is uninformative;
  - coupling only at a common spread or on z2 / psi;
  - covariates that change rhythm strength are balanced (matched or stratified), with pooled and balanced estimates
    both reported;
  - loudness, not amplitude, for variance ratios (METHOD section 9).
- **Compute.**
  - Machine budget (memory: at most about 21 workers and about 83 GB summed caps; at most two Sonnet subagents).
  - Runs over 20 min go through `paper/common/run.sh`, with `watch_jobs.sh`.
  - Verify subagent numbers against the output files before writing them.
- **Data access.** GTEx per-donor data and protected phenotypes: aggregates only, never printed or committed.

## Tasks, in order

**Task 0, housekeeping (ask first).** The A336 tier-2 scripts (`paper/exploratory/A336_two_readouts/scripts/p1_*.py`,
`ct_lib.py`, modified `merge_part1.py`) and the merged `output/two_readouts.csv` are uncommitted, although RP cites
them. Ask the user whether to commit them as they are, after checking that `two_readouts.csv` reproduces the RP
tier-2 numbers.

**Task 1, tumour stroma on the core loop** (audit addition 1; gap 1).
- *Question.* Do tumour fibroblasts and endothelium keep the patient's normal stromal clock state beyond both
  companion families, on 14 genes and on the core loop?
- *Data.* The A149 / A151 / A227 cohorts.
- *Statistic.* The A149 note 5 form: the uncorrected complex correlation between tumour and own-normal stroma.
  - Calibrate it with the A283-style placebo.
  - Holm over the two compartments.
  - Also give the malignant-normal epithelium row as a companion readout.
- *Blind status.* The 14-gene values are known (A149 note 5, A227, A271), so not blind. The core values are blind.
- *Order.* Power stage at the core projector first.
- *Readings to fix.* Pass: "tumour stroma carries the patient's time". At the line or fail: the stroma clause stays
  post hoc.

**Task 2, inflamed skin against psoriasis** (addition 2).
- *Question.* Is the lesional delay specific to psoriasis, or shared by atopic dermatitis?
- *Data.* The gated psoriasis series (A54 / A53 set, 465 pairs) and atopic series (298 pairs) under
  `$SCCC_DATA/skin_lesional/`; the Tsoi atopic pairs (A345 counts); GSE154200 eczema (A183).
- *Statistic.*
  - Pooled paired offset per disease (`clocksync.paired`), with A308-calibrated companions.
  - Psoriasis minus atopic, with an equivalence bound.
  - 14 genes and core loop.
  - Hidradenitis (A254) descriptive.
- *Blind status.* Not blind: the A53 per-series offsets were seen. State it.
- *Power.* Report the MDE for the atopic delay before scoring.
- *Also.* Compute the A345 registered descriptive that was never run: the pre-mRNA lead within lesional and within
  non-lesional samples (A342 lead statistic).

**Task 3, living organ sharing on the core loop** (gap 2).
- *Statistic.* A251 H251.1-2 and A255 H255.1 / H255.3 on the core loop, with the A285 / A320 calibrated companions.
- *Blind status.* The 14-gene values are reproduced first and not blind; the core values are blind.
- *Order.* One run per amendment part, on the A336 pattern.

**Task 4, intensive-care detection limit and host spread** (gap 3; power and descriptive only, no new claim).
- *Detection limit.* Redo the A336 part 3 order detection limit at the ventilated design: 79 ventilated post-mortem
  donors, their own pair coherence, against the fast-death reference.
- *Descriptives.* Report, beside the A315 Q*:
  - the host spread V per group;
  - the host-referenced log z2 (fast against ventilated, at a common spread);
  - the relative organ-gain pattern under the overall loss (review question 1);
  - the A132b compression ratio after removing the A109 non-temporal axis, if registered as a readout.

**Task 5, conserved class gains on the core loop** (gap 7). A293 H293.4 on the core loop: the tier-3 atlas refit.
Scout the cost first; register only if it is under a day of compute.

**Task 6, CHIRAL age rescore** (addition 3). **Draft the amendment only.** It corrects the lab's own Science paper, so
the user and the PI decide before anything is run.
- Old against young within cohort, matched as A231 plus RIN.
- Host-referenced log z2 at a common spread, and the variance readout with matched / responsive Q*.
- Core loop.

**Not in scope.** These are outward-facing and need the user or the PI:
- the jet-lag data;
- the STARNET application;
- Lundell, oral mucosa, and A260 on GTEx (audit section 4: not to add).

## Deliverables
- One amendment per task, each committed before scoring. Each scored result goes in `results_preliminary.md`, with
  a VERDICT_TABLE row.
- A short note appended to `paper/docs/STORY_AUDIT_2026-11-23.md` (section 8): which gap each result closes and which
  wording it allows. Proposed, not applied.
- PI-level points as questions in a reply to the user, not in PI_QUESTIONS.md.
- Commit with the attribution line and `git push origin main`. If the push is refused, say so.
- Reply in at most 15 lines.
