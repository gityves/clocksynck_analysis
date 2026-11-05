# Prompt: presentation, continuing with the introduction

You continue building a working presentation of the clocksync project with the user. The user assembles the PowerPoint;
you write the content and make the figures. Work section by section and check each slide with the user before moving
on.

## Project and framing

`clocksync` infers circadian relations between compartments (organs, cell types) of the same donor from untimed
clock-gene pseudobulks: coupling c, relative amplitude |z|, phase offset arg z, read against the clock state the donor's
compartments share (the host). Read `CLAUDE.md` at the repository root first.

**Framing (user, 2026-11-05):**
- **Purpose:** biological. Find internal desynchrony across disease, and describe circadian organisation in health.
- **Why a method is needed:** most human tissue comes from hospital hours (surgery, biopsies). Existing phase-inference
  methods need samples spread over the 24 h cycle (orderers) or timed training data (supervised).
- **The method:** a means, not the purpose. Lead with the biology, then the data constraint, then clocksync as the
  enabling tool.

## State

- **Done:** the model animation, `paper/presentation/scripts/model_animation.py`. Its output is
  `paper/presentation/figures/model_animation/` (`model_01.mp4` … `model_16.png`, one clip per transition and one image
  per pause, plus `model_animation.mp4`).
- **Opening-slide image:** the user is making it with an image-generation prompt already given.
- **Overall plan:** `paper/docs/HANDOFF_2026-11-04.md` gives the structure of the presentation (sections 1-11) and its
  rules. Every number must come from an output table, and GTEx-derived data are shown as aggregates only.
- **Literature:** fully read.
  - `paper/manuscript/literature/notes/<key>.md` has one card per paper (148), with quotes checked against the full text.
  - `paper/manuscript/literature/synthesis.md` §1 summarises the hospital-hours evidence and §2 the methods table.
  - Bibliographic data are in `paper/manuscript/references.bib`.
  - Full texts are in `/mnt/max8tb/Paychere/literature/zotero_export/clocksync/text/` (read only; copyrighted; never
    copy into the repo).
  - The old `refs_*.md` tables were written from abstracts and are superseded for claims: about 45 of their claims were
    wrong.
- **Scripts and figures:** scripts go in `paper/presentation/scripts/`, figures in `paper/presentation/figures/<name>/`.
  Use paperfig (`/home/maxine/Documents/paychere/paperfig`), and load the `dataviz` skill before the first plot.
  paperfig places labels without collisions and reports failures (`pf.place_labels`, `pf.place_all`,
  `pf.check.collisions`). Extend paperfig when something is missing instead of working around it. Example slide style:
  `paper/presentation/example_slides/`.

## Task 1: `paper/presentation/introduction.md`

Write the introduction section of the talk as a slide-by-slide document, reusable for the paper's introduction. The
user wants it bullet-point style, with citation keys from `references.bib` (as `\cite{key}`) so that it can later become
`paper/manuscript/introduction.tex`.

- **Subsections (agreed):**
  1. Circadian cycle
  2. Internal desynchrony
  3. Cancer interaction
  4. Technical challenges, ending on hospital hours and why the model is needed
  - Technical challenges may be split in two.
- **For each slide:** the message (one line); the bullets, each with its citation keys; the figure or visual
  (literature schematic, data plot, or an animation step); and speaker notes if useful.
- **Every literature claim must trace to a card.** Use the card's quote and numbers, and re-check the number in the
  full text when it goes on a slide. Do not use the `refs_*.md` one-liners.
- **Citation rules the user asked for:**
  - **Phase inference:** the literature (orderers, supervised predictors, single-cell methods), including CoPhaser
    (`paychere2025cophaser`, our previous work; clocksync's fixed gene acrophases with fitted amplitude and mean shifts
    come from it) and scRitmo (`salati2026inferring`).
  - **Naef-lab circadian and cell-cycle work:** `bieler2014robust` and `droin2019low` (the cell cycle accelerates the
    clock), plus `talamanca2023sexdimorphic`, `yeung2018transcription`, `atger2015circadian`, `weger2021systematic`
    and `phillips2020circadian`.
  - **Attribution checks:** `feillet2014phase` is not Naef lab, `mure2018diurnal` is Panda/Cooper, and
    `sinturel2021circadian` is Schibler/Dibner. Check each paper's authors on its card before calling it lab work.
- **Key points the reading established** (check the cards; don't overstate):
  - **The full-cycle requirement is stated by the methods themselves:** `anafi2017cyclops` ("It requires data from the
    entire periodic cycle to form an ellipse. Biopsies are almost exclusively obtained during the day."),
    `shilts2018evidence` and `grabe2024quantification`.
  - **Hospital hours in data:**
    - recorded surgery times of 302 glioblastoma patients, analysed 10:00-18:00 (`gonzalezaponte2025circadian`);
    - TCGA inferred phases concentrated in the working day (`li2024tumor`).
  - **GTEx** deaths cover the whole day (`wucher2023day`: 351 / 315 / 266 donors in the day / night / twilight
    windows), so GTEx is the reference where hospital hours are simulated (A123). Call it a non-cancer post-mortem
    reference, not "healthy" (`lonsdale2013genotype`).
  - **Published tumour "uncoupling"** comes from methods that need full coverage applied to working-hours cohorts
    (`ananthasubramaniam2025rhythm`, `vasilyev2026robust`). `li2024tumor` reports agreement in strong-clock tumours.
    Present this as a reason for the method, not as our result.
  - **Closest precedents:** `shilts2018evidence` (CCD: time-free, group-level, blind to phase) and
    `grabe2024quantification` (tumour ordered by its paired normal; assumes full coverage).
  - **Synchrony in health:** `talamanca2023sexdimorphic` (tissues within a few hours) vs `hughey2016differential` (brain
    about 8.5 h behind blood, between people). Desynchrony taxonomy: `galinde2023different`, `mello2026rhythms`.
  - **Chronotherapy is motivation, not proven benefit.** In `levi1997randomised`, response improved but survival did
    not. `levi2010circadian` and `levi2024circadian` report a sex interaction, and only 3 of 18 randomised trials
    improved.
- **Then:** propose the figures for the intro (schematics, or simple plots of published aggregates such as the GTEx
  death-time windows) and build them after the user agrees on the text.

## Task 2: continue the presentation

After the introduction, continue with the method section (projection slide, EM slide, simulations), then the results
sections of `HANDOFF_2026-11-04.md`. Apply the per-slide rules there: claim, provenance and status, plot, source CSV
path, limit; and use the corrections in `AUDIT_CLAIMS_LEDGER_2026-11-03.md`.

## Working rules

- Use the sccc environment's Python. Commit locally with the attribution line; there is no remote. Commit scripts and
  markdown; figures stay git-ignored if their folder is.
- Ask the user before deciding anything about claims, framing or wording that changes what the paper says.
- Keep slides readable: a few bullets each, one message per slide.
