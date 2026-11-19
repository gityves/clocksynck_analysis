# Independent review of the paper plan (2026-11-19)

You are an independent reviewer of the sccc paper plan. You are not the agent that wrote it. Judge it as a critical
senior reviewer for Cell would, and as a co-author deciding whether to invest the next months in it. Your report is
the only output; do not edit any project file.

## Context (read in this order)
1. `CLAUDE.md`: project, method, rules (pre-registration, companions, covariate balancing, controlled-access data).
2. `paper/docs/PAPER_PLAN_2026-11-19.md`: the plan under review: frame "one time, many volumes", claims (section 4),
   figures (section 5), checks and their status (section 6, 6.2b).
3. `paper/docs/MODEL_signal_clock_2026-11-19.md`: the conceptual model (a slow clock behind a fast signal layer).
4. `paper/docs/SYNTHESIS_2026-11-19.md` and `paper/docs/REVIEW_2026-11-19.md`: the evidence synthesis and an earlier
   review.
5. `paper/docs/METHOD.md`: what clocksync estimates (coupling, relative amplitude, offset; companions; gate).
6. `paper/docs/VERDICT_TABLE.md` and `paper/docs/triage/registry_2026-11-19_summary.md`: every registered test with its
   verdict and blind status.
7. `paper/docs/results_preliminary.md`: open the entries the plan cites (search `## A<ID>`); each ends with a `src`
   comment naming its output CSVs under `paper/exploratory/` or `paper/results/`.
8. `paper/docs/literature_positioning.md` and `paper/docs/literature/`: what was searched and found.
9. `paper/docs/PI_QUESTIONS.md`, Q41 and Q42: the open decisions.

## What to judge
For each claim of PAPER_PLAN section 4, and each figure of section 5:
- **Support.** Which registered results carry it, were they blind, powered and clock-specific (beyond matched and
  stress-responsive companions)? Distinguish confirmatory, not-blind, post hoc and descriptive evidence. Is the claim
  worded at the strength the evidence allows?
- **Numbers.** Spot-check at least 10 numbers quoted in the plan against the results entries and, for at least 5, the
  output CSVs themselves. Report every mismatch.
- **Alternative explanations** the plan has not excluded, for example: death-hour distribution and window synchrony,
  composition, ischaemic time and RNA quality, ventilation and Hardy scale, post-mortem versus living tissue, the
  human-versus-mouse waveform, array or batch scales shared by two tissues, regression to the mean in two-step designs.
- **Statistics.** Coupling scales with host spread; failed gates are uninformative; group contrasts need balanced
  covariates (CLAUDE.md). Flag any place where the plan reads a failed gate as a finding, compares coupling across
  unequal spreads, or relies on an unbalanced covariate.
- **Novelty.** Against the literature files and your own knowledge: what is new, what is known, what a reviewer would
  say was shown before. Name papers; mark any you cite from memory as "to verify".
- **Coherence of the frame.** Does "one time, many volumes" (time set per organ and shared by its cell classes;
  loudness factorising into person, organ and class; perturbations moving loudness rather than time) follow from the
  evidence, or is a simpler frame better? Where is the frame weakest?
- **Biological interpretation.** The model has a slow clock with inertia while output genes follow drives quickly.
  Does the plan state clearly what peripheral clocks are for in this picture, and is its reading consistent with
  tissue-specific clock knockouts and liver / skin clock-reconstitution studies?
- **Figures.** Is each figure's message one sentence, carried by its panels, and are the panels the strongest
  available evidence? What would you cut, merge or move to the supplement?
- **Route.** Q41 (Cell Article versus Resource) and the Q42 decisions: give your own recommendation with reasons, without
  deferring to the recommendations written there.

## Rules
- Read-only. No new analyses beyond reading output CSVs; no downloads; no registrations; no commits.
- GTEx per-donor rows and times of death are controlled access: quote aggregates only.
- Be specific: every criticism names the claim, the file and line or amendment, and what would fix it (a rewording, a
  check on data already on disk, or new data).

## Report
Write `paper/docs/reviews/plan_review_2026-11-19.md` (create the folder), at most about 3,000 words:
1. Verdict in five lines: is the plan Cell-ready in substance, and what are the three changes that matter most?
2. Claim-by-claim table: claim, evidence (amendment IDs, blind / powered / specific), strength, problem, fix.
3. Number checks: each checked number, its source, match or mismatch.
4. Alternative explanations and statistical issues, ranked by how much they could change a claim.
5. Novelty and literature.
6. Figures.
7. Route (Q41) and the Q42 decisions.
8. What you would ask the authors if you were the Cell reviewer (at most 10 questions).
