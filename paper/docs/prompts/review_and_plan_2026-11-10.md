# Prompt: fresh review of the project, the method code, the results and the paper plan (2026-11-10)

You join the clocksync project (Yves Paychere, EPFL, Naef lab) with fresh eyes. Previous agents ran about 235
registered analyses and know the material well, but they may also have grown blind to its weak points. Your job is to
form your own view: understand the biology and the literature, review the method code, judge what the results really
support, and propose a paper plan and the questions still worth asking. You decide how deep to go in each part and in
what order; the list below is orientation, not a checklist.

## Context to build

- **Start with `CLAUDE.md`** (repository root): project, rules, environment, data policy. Its rules apply to you.
- **Method:** `paper/docs/METHOD.md`, then the code in `clocksync/src/clocksync/` (entry points `dataset.py`,
  `results.py`, `_core/`), its tests (`clocksync/tests/`), `clocksync/README.md`, `CHANGELOG.md` and `to_add.txt`.
- **Literature:** `paper/manuscript/literature/synthesis.md` (summary), one card per paper in
  `paper/manuscript/literature/notes/` (148, quotes checked against full texts), `paper/docs/literature_positioning.md`.
- **Results:** `paper/docs/VERDICT_TABLE.md` (one row per registered test), `paper/docs/results_preliminary.md` (each
  entry ends with a `<!-- src: -->` comment naming its output tables), `paper/docs/registered_experiments.md` (what was
  registered and when; A109 onwards, earlier in the archive `/home/maxine/Documents/paychere/scCircadianCancer`).
- **Synthesis documents, to read critically rather than adopt:** `CLAIMS_LEDGER_2026-11-03.md` and its audit,
  `FINDINGS_MAP_2026-11-03.md`, `NOVELTY_ASSESSMENT_2026-11-06.md`, `METHOD_EDGE_2026-11-09.md`, `NOTABLE_RESULTS.md`.
- **Plan and decisions:** `paper/docs/provisory_plan_paper.md` (WP0-WP8) and `paper/docs/PI_QUESTIONS.md` (answered
  questions hold the PI's decisions, with dates; open ones are pending).

Some decisions already taken, so you do not reopen them without a new reason:
- **One paper, including the cancer chapter.** Cell first, Nature Methods as fallback.
- **Level of inference:** the method gives population-level estimates of within-person relations, the same level as
  CHIRAL. It does not make individual-level claims (novelty assessment, §2 item 1).
- **Brain death:** the loss of inter-organ order in brain-dead donors is reported descriptively, with no mechanism
  amendments. It is also the candidate first human observation consistent with the mouse SCN-lesion work.
- **Ageing:** A234 shows an advance in GTEx; the living cohorts (A130) point the other way. This is PI question Q26.
- **Phase 8** (naming, scorecard, figures, code hardening, tutorials) is never started without the PI.

## What we would like from you

1. **The method code.** Look for errors (statistical or numerical), fragile assumptions, untested paths, places where
   the public API and its documentation disagree, and improvements a future user would need. Weigh what matters for
   the paper's claims over style.
   - **Known open issues:** gate false positives in no-clock simulations (A145 H145.3, A209); bootstrap intervals that
     under-cover (A154, A162, A163; permutation nulls are used instead).
   - **Decided:** permutations that refit the host per replicate in the analysis scripts are slow but will not be optimised.
   - **Small clear bug with a test:** you may fix it and commit it, then rerun the clocksync tests.
   - **Anything that could change a published number:** report it with the analyses it touches. Do not change it.
2. **The results.** Which claims are solid, which rest on post hoc choices, low power or a single cohort, and where
   documents disagree with their output tables. Spot-check numbers against the CSVs named in the `src` comments. If
   you think a verdict is wrong, say why.
3. **A paper plan.** A storyline (title idea, abstract-level claims, figure sequence with the analyses behind each
   panel) that a Cell reviewer would find important and provocative, and the Nature Methods fallback. Say what is
   missing for each.
4. **New hypotheses or research directions.** Anything the method could answer that has not been asked, including
   leads from the literature, gaps you see in the story, or datasets not yet used. For each: the question, why it
   matters for the paper, data (public, or what access it would need), a rough feasibility and power argument, and
   what each outcome would mean. Be frank when you think nothing more is needed.

## How to work

- **Read-heavy work:** delegate searches to subagents (Sonnet 5.5, at most 2 at a time) when useful. Verify any number
  they report before you use it.
- **New tests:** propose them, do not run them. Analyses are pre-registered (`CLAUDE.md`, rule "Pre-register"). A quick
  exploratory look to judge feasibility is fine if labelled post hoc and kept out of the docs' results.
- **Controlled access:** GTEx phenotypes and time of death are aggregates only, never per-donor rows.
- **Outward-facing steps** (data requests, emails, anything public): suggest them, do not take them.
- **Machine budget** (shared machine): at most 16 workers and about 120 GB in total. Long jobs go through
  `paper/common/run.sh`.

## Deliverable

One report, `paper/docs/REVIEW_2026-11-10.md`, concise and publication-grade, with:
- **Method code:** findings ranked by impact, each with file:line.
- **Results:** assessment, plus the discrepancies found.
- **Paper plan.**
- **Proposed hypotheses and directions,** ranked.
- **Questions for the PI,** if any (also added to `PI_QUESTIONS.md`).

Commit it locally with the attribution lines from your session; there is no remote. Then give the user a short summary.
