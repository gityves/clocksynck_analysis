# Prompt for the slide agent: finish the results section of the lab presentation (written 2026-11-25)

You continue building the results section of an internal lab presentation, slide by slide, with the user (Yves).
Sections 1-5 and slide 6.1 are built. Your job: 6.2, 6.3 and 6.4, then sections 7 and 8, then a pass over the
backup slides. A second agent (the lead) runs analyses in parallel; see "Coordination" below.

## Read first
- `CLAUDE.md` and the memory index (`MEMORY.md`).
- The `paper-figures` skill (`.claude/skills/paper-figures/SKILL.md`): house style, palette, procedure. Load it before
  writing any plotting code.
- `paper/presentation/results_outline_2026-11-23.md`: the outline you extend. Read its header (phase convention,
  framing, tags) and at least slides 3.5, 4.3, 5.1 and 6.1 to match their format: title as a message sentence,
  **Bullets** (two to four, numbers with amendment IDs), **Figure** (path, script, layout, what each panel shows),
  **Tags**, **Caveats (speaker notes)**, **Source** (output tables).
- `paper/docs/prompts/presentation_results_2026-11-23.md`: the original plan for sections 6-8 and the backup list.
  It predates several results; check every item against `paper/docs/results_preliminary.md` before using it.
- `paper/docs/ABSTRACT_DRAFT_2026-11-22.md` (wording limits) and `paper/docs/STORY_AUDIT_2026-11-23.md` sections 7-8.
- Example scripts: `paper/presentation/scripts/tumours.py` (6.1), `chronic_disease.py` (5.1),
  `holds_organisation.py` (4.1). Split compute from plotting when a figure needs computation (`*_compute.py` writes
  aggregate CSVs to `presentation/output/`, the plot script reads only those).

## What to build
1. **6.2 Corrections of published tumour-clock claims** (the user wants it kept for this presentation): CoFE
   uncoupling as an artefact (A265, A267); prognostic clock signatures add nothing beyond proliferation or random sets
   (A192-A194, A220); clock-gene loss exaggerated by GTEx-normal designs (A221). Check each verdict and number in
   `results_preliminary.md` first (look for errata and later amendments that revise them).
2. **6.3 Chronic lesions keep time and lose a little loudness** (A254; core loop at the line, A336 row 6c).
3. **6.4 Psoriasis, the one local phase shift** (A183, A308, A345, A347, A115; A332 D1; core-loop row 6d of A336:
   the offset is at the line against matched sets, psoriasis minus eczema not decidable).
4. **Section 7 (end of life and the limits).** Brain death is now covered in 4.3 (organ arrangement, A317, A374,
   A377). Ask the user whether section 7 keeps a brain-death slide or becomes a limits slide only, before drafting.
5. **Section 8 (synthesis, 1-2 slides).** Propose the message sentences from the built slides (not from the old plan)
   and let the user choose before drawing anything. Section 5 now includes A379-A382: chronic disease moves neither
   the arrangement nor the absolute loudness at the sizes bounded there.
6. **Small edit to 3.3 speaker notes**, offered to the user and not yet answered: one line on A378 (core vs shell
   gain is not consistent across species; the squamous ratio is partly conserved in the baboon, the fat ratio in
   neither; the organ factor stays open). Ask before adding it.

## How to work (user's preferences, learned over this presentation)
- **One slide at a time, with the user.** For each slide:
  - read the results entries and output tables;
  - propose the title, bullets and panel design in a few lines;
  - build only after the user approves;
  - render the figure and look at the PNG (Read the image), and check that `pf.collisions(fig) == []`;
  - show the user the figure's full absolute path;
  - after approval, commit the script and the outline section together.
- **Numbers come only from output tables or results entries**, never from memory, this prompt or earlier slides.
  Name the table in **Source**.
- **Label results:** blind / not blind / descriptive / post hoc; anything chosen after seeing results is post hoc.
- **No internal identifiers on figures** (no amendment numbers, hypothesis codes or script names). IDs go in bullets,
  tags, notes and docstrings.
- **Phase sign on figures:** peak time, + = later. clocksync offsets are > 0 for earlier: negate them. Offsets
  between two compartments say their direction in words.
- **Wording:**
  - "loudness" for amplitude and variance ratios;
  - internal desynchrony = a phase offset between compartments;
  - use units and words consistent with the neighbouring slides (the user corrected panel labels to match 3.5 and 5.1).
- **Avoid repeats:** do not repeat a panel already shown on another slide. Slides that add little go to backup with a
  note; ask the user.
- **Layout:** keep panel heights consistent across a slide, and keep titles clear of the plot (`pad`).
- **Gene drops:** follow the gene-drop rule in CLAUDE.md. Never write that PER1, PER2 and NR1D1 "carry" an effect
  because it shrinks on the core loop; quote the A336 verdict ("holds", "at the line", "falls") as written there.
- **Honesty:** if a slide's claim is weaker than the old plan says (at the line, not decidable, superseded), say so on
  the slide and tell the user.

## Limits of your role
- **No new tests.** If a slide seems to need an analysis that has not been run:
  - stop and tell the user. New hypotheses must be pre-registered (`docs/registered_experiments.md`), which is the
    lead's job;
  - plots that only redraw existing output CSVs are fine.
- **Do not edit** `paper/docs/registered_experiments.md` or `paper/docs/results_preliminary.md`. If you find a number in
  a results entry that disagrees with its output table, report it to the user.
- **Long jobs** (> 20 min): only through `RUN_CAP_GB=.. paper/common/run.sh <script>` launched from the repo root,
  detached. At most 4 workers: the lead is using about 10 of the 16-worker budget.
- **Controlled-access data:** GTEx time of death and protected phenotypes are aggregates only. Never print, plot or
  commit per-donor values; plotted bins need ≥ 5 donors.

## Coordination with the lead agent
- **Leave the lead's work alone.** The lead is running the power stages of A383/A384
  (`exploratory/A383_species_organisation/`) and A385/A386 (`exploratory/A385_sheep_person/`). Do not touch those
  folders, their logs, or any running job.
- **Commit only your own files**, by path (`git add paper/presentation/scripts/<file>.py
  paper/presentation/results_outline_2026-11-23.md`). Never `git add -A` or `git add .`; the tree has many untracked
  files that are not yours.
- **Do not push.** The user pushes.
- **Room for later slides.** Results from A383-A386 (cross-species organisation, sheep animal loudness and
  glucocorticoid link) may later add a section 3 slide. Leave room for it in the section 3 table as "pending (lead)";
  do not draft it.

## Done when
- 6.2-6.4 and sections 7-8 are in the outline with figures approved by the user and committed.
- The backup list at the end of the outline is up to date.
- You report to the user:
  - the slides built (with their figure paths);
  - what was moved to backup;
  - any number or verdict you found inconsistent.
