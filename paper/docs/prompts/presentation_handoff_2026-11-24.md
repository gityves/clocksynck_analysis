# Prompt: take over the interactive slide making (results section), 2026-11-24

You continue building the results section of the internal lab presentation with the user, one slide at a time. The
task, framing and section list are in `paper/docs/prompts/presentation_results_2026-11-23.md` (read it in full, it
still holds). This note says where the work stands and how the user wants it done.

## Read first

1. `CLAUDE.md` (rules: pre-registration, numbers from output tables, controlled access, commits, no push).
2. The memory index, especially "Presentation step by step": one figure at a time with the user, no batch edits, no
   extraction agents for slides.
3. `.claude/skills/paper-figures/SKILL.md`: load it (Skill tool) before touching any plotting script. It records the
   user's figure rules, including the phase-sign convention and "no internal identifiers on figures" (2026-11-24).
4. `paper/presentation/results_outline_2026-11-23.md`: the slide-by-slide outline (header convention, then per slide:
   title as a message, bullets, figure, tags, speaker notes, source). This is the main deliverable; keep its format.
5. `paper/docs/COUPLING_SOURCES_2026-11-24.md` and `paper/docs/PAPER_PLAN_2026-11-19.md` sections 2.1, 3, 4, 11 (the
   current framing; section 2.1 has the "Where coupling comes from" block).

## State (2026-11-24)

| slide | status | figure (`paper/presentation/figures/...`) | script (`paper/presentation/scripts/...`) |
|---|---|---|---|
| 1.1-1.3 | kept from `results.md` | gtex_organ_map, living_patients, baboon_offsets | same names |
| 1.4, 1.5a, 1.5b, 1.6 | done | living_timed, celltype_controls, positive_controls, living_tissues | same names |
| 2.1-2.3 | done | development, fetal_clock, fetal_signal | same names |
| 3.1 reference organisation | done | reference_map | reference_map.py |
| 3.2 cell classes | done | celltypes_health | celltypes_health.py |
| 3.3 organs differ in loudness | done | organ_gain | organ_gain.py |
| 3.4a louder organs are more coupled (human, baboon; right axis: phase SD [h]) | done | organs_health | organs_health.py |
| 3.4b why: departure of fixed size, not fixed angle | done | organs_departure | organs_departure.py (+ `_compute.py`) |
| 3.5 each person has a clock volume every organ follows | done | person_loudness | person_loudness.py |
| **3.6 blood as the exception** | **next: to build** | - | - |
| 3.last organ arrangement is internal | done | organ_arrangement | organ_arrangement.py |
| 4.1-4.5 drives | skeleton in the outline (section 4 table) | - | - |
| 5-8 | not started; follow `presentation_results_2026-11-23.md` | - | - |

For 3.6: the user's standing reading is that blood has no coherent core clock; its time signal is diurnal cell
composition, an expected exception, not a method failure (memory "Blood has no coherent core clock"). Find the
entries (grep `results_preliminary.md` for whole blood, A262 night-shift blood, the A104 / A294b blood rows) before
proposing the slide.

## How the user works (learned on sections 1-3)

- **One slide at a time.** Re-render the existing figure if there is one, look at the PNG, propose changes in a short
  list, wait for "go ahead". After each change: render, Read the PNG, check `pf.collisions(fig) == []`, fix overlaps
  the check misses (text crossing lines or error bars), then commit the script and the outline text.
- **Observation first, explanation after.** Split a slide in two when the "why" needs its own panels (3.4a / 3.4b).
- **Section 3 stays free of drives.** Drives, cues, glucocorticoid reading and mouse feeding experiments go to section
  4. If a section-3 panel touches them (e.g. 3.5 panel C), flag it.
- **Same visual form across species and panels** (e.g. the coupling curve a^2 / (a^2 + D) on human and baboon).
- **Hours are more intuitive** than unitless coupling: convert where the model allows (phase SD [h] =
  (24 / 2 pi) sqrt((1 - c) / (2 c))), and say it is implied, small-angle, approximate beyond about 3 h.
- **Name estimates by method**, not "measured" vs "implied" when both are estimates on the same donors (3.4b: "model
  estimate (D, S)" vs "gene-split estimate").
- **No internal identifiers on figures** (amendment numbers, hypothesis codes): they go to captions, docstrings and
  notes.
- **Phase convention:** peak time, + = later (CHIRAL Fig. 1B); clocksync offsets are shown negated.
- **Be exact about what a result supports.** The user asks pointed questions ("what was tested exactly?", "how was
  elasticity defined?"). Answer from the registration and the code, and correct earlier statements when they were
  wrong. Examples from this session: A311 does not show tissues are driven below their self-sustained regime; A331's
  MZ / DZ values are not heritability; the A312 elasticity level is about 1.1 for proportional scaling (calibration
  note, post hoc).
- **New questions become registrations.** If the user's question needs a new statistic, propose an amendment; write
  and commit it before scoring (next free number: **A363**; A361 / A362 are reserved by
  `prompts/person_volume_desynchrony_2026-11-24.md`, run by another agent). Calibrations and displays chosen after
  seeing results are labelled post hoc.
- **The user may be on another device**: send rendered PNGs with SendUserFile when they ask to see slides.

## paperfig notes

- `pf.Layout.slide(..., kind="half" | "wide")`; `wgap` is ignored by slide layouts: add space with a blank spacer
  column (`"AxB"`, `blank=("x",)`).
- Two legends on one axis: create the first, `ax.add_artist(leg)`, then the second.
- `ax.secondary_yaxis("right", functions=(fwd, inv))` for derived axes; labels on both panels need the spacer.
- Organ-system colours and legends from `_systems.py` (`COLOR`, `gtex_system`, `baboon_system`, `legend_in`).

## Open items to mention to the user when relevant (do not act unasked)

- About 90 local commits are not pushed; the user runs `! git push origin main`.
- A312 full null calibration (all tissues, about 20 draws) would decide whether H312.1 stands: proposed, not decided.
- A360 on the core loop: not checked.
- The in-house jet-lag data (species, organs, time points, feeding control, location) are needed for the drive-weight
  test of section 4.5.
- The person-volume programme (heritability, exposures) is with another agent; its outcomes feed slide 3.5 and
  section 4.
