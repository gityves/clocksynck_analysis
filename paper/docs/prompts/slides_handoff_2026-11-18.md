# Prompt for a new agent: continue the results slides of the lab presentation (written 2026-11-18)

You take over the slide work of the presentation results section. The audience is the lab (internal talk), not a
journal. The user assembles the PowerPoint; you write the slide text in `paper/presentation/results.md` and make the
figures, one slide at a time. Wait for the user's go before each new slide.

## Read first
1. `CLAUDE.md`, `MEMORY.md`; the project skill `paper-figures` (load it before any plotting; it is a living style
   guide: when the user corrects a figure, add the rule there with the date).
2. `paper/presentation/results.md`: R1 (validation) is finished; R2.1 (organs) and R2.2 (cell types) are drafted with
   speaker notes. Read them to match the format (message, bullets, figure description, speaker notes, sources).
3. `paper/docs/STATUS_2026-11-18.md`; RP entries named below. `git log --oneline -30` (other sessions commit on `main`).
4. Existing scripts in `paper/presentation/scripts/` (`gtex_sweep.py`, `organs_compute.py`, `organs_health.py`,
   `celltypes_health.py`, `_systems.py`): reuse their palette, helpers and layout conventions.

## Chapter order agreed with the user
Method (description + validation, short) -> Organs -> Cells -> What holds the organisation (loss of daily cues, central
drive) -> (Development, only if re-scored) -> Internal desynchrony in disease -> Synthesis. Limits go inside each
chapter plus one final slide.

## Pending changes to the drafted slides
- **R2.1 organs:** reframe the message as gain over a constant departure: each organ expresses the body's time with its
  own gain over a clock-specific departure of nearly constant size; coupling = S / (S + D) (A291, A291 post hoc,
  A294b, A302). Panel D (departure SD, A259) becomes the D of that relation. Bullet 2 and the notes still say "A290
  running": A290 (calibration of the amplitude-coupling relation against stress-responsive genes) runs over the
  weekend of 2026-11-21 by the user, scored with `placebo_a290.py --score` (deviation registered); update when its
  `tests.csv` exists, not before.
- **R2.2 cells:** the user finds the class-gain result the most interesting: "the organ sets the time; each cell class
  sets how strongly it expresses it". Redraw around relative gain per class (A293 Part 3 forest, H293.4), lymphocytes
  as the zero-gain end of the same scale (A281, A287, A294; "no readable clock on the full gene set", core-loop caveat
  post hoc), phase set by the organ (A273, A279, A295; endothelium lead +0.46 h only a lead, H293.3 FAIL). The organ
  factor of gain is inferred, not tested (A298 closed).
- **New, "What holds the organisation":** slow deaths and ventilated living-brain donors keep the arrangement
  clock-specifically, brain-dead donors lose it (A293 Part 1, descriptive for brain death; GTEx is the only data, PI
  Q7); mouse feeding removal keeps liver amplitude (A297).
- **New, disease chapter, four rows:** none (seven chronic conditions: A164 / A202, A233, A235, A306), gain loss
  (proliferating tumours A304, non-skin lesions A254), phase shift (psoriasis: A54, A183, A308), loss of order (brain
  death). Frame it for the lab's prior (mouse jet lag -> obesity, liver disease, diabetes; a reciprocal
  disease-desynchrony loop): cross-sectional data test only disease -> desynchrony; the method is blind to
  misalignment with the external day (constant offsets cancel), so the human null does not contradict the mouse work.
  Psoriasis wording: "about 1.4 h behind the same patient's skin, replicated blind, clock-specific against random genes
  once their noise is removed (A308, not blind); in the local clock (cue genes on time, A78 / A81b archive);
  composition not directly controlled (state mix explains at most a quarter, A69)".
- **Synthesis:** "one body time, broken only locally: tissues lose gain or shift phase; body-wide order fails only
  without central drive". If the modelling agent (`prompts/clock_model_2026-11-18.md`) has reported, check its verdict
  before using any model on a slide; never present the model as evidence.

## Rules carried from the user
- Fair comparisons: score every orderer after its best reflection and rotation (memory `fair-orderer-scoring`).
- Do not suggest a function the data do not show (R2.2 lymphocytes); do not present purely technical explanations as
  findings; companion gene sets drawn as grey ranges or open markers.
- Every number on a slide comes from an output table; GTEx: aggregates only; new statistics are registered amendments
  before they are computed, anything chosen after seeing data is labelled post hoc.
- Look at every rendered PNG; `pf.collisions(fig)` must be `[]`. Commit scripts and `results.md` only (`figures/` and
  `output/` are gitignored).
