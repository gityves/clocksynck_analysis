# Prompt for the next agent: results section of the internal presentation (written 2026-11-23)

Build the outline and slide content of the **results section** of an internal lab presentation. Lifespan order, which
also follows the method:
- before birth, only the model-free gate and the between-individual coherence (L1) are used;
- from adult health on, the full host-referenced method is used.

The framing is that of `ABSTRACT_DRAFT_2026-11-22.md` (revised 2026-11-23) and `PAPER_PLAN_2026-11-19.md` section 11:
"loudness, not time". The mechanism found in mouse (drives set time) holds in human; the desynchrony that mouse models
predict is absent; true internal desynchrony is rare and local.

Read first: CLAUDE.md, the memory index, the existing `paper/presentation/{introduction,method,results}.md` (reuse what
still holds, mark what the new framing replaces), the abstract draft (abstract, change list, wording limits), plan section 11,
`STORY_AUDIT_2026-11-23.md` sections 7-8, and the results entries named below.

## Rules
- **Numbers.** Every number on a slide comes from an output table or a `results_preliminary.md` entry, with its
  amendment ID in the slide notes. Never from this prompt, from memory, or from earlier slides without checking.
- **Wording.** Follow the wording-limits table of the abstract draft. Keep the term rules:
  - internal desynchrony = a phase offset between compartments;
  - "loudness" for variance ratios (METHOD section 9);
  - before birth: "clock genes move together", not "a fetal clock";
  - brain death descriptive.
- **Labels.** Mark each result blind / not blind / descriptive / post hoc on the slide (small tag).
- **Figures.** Use the `paper-figures` skill and the existing notebooks and figures in `paper/presentation/` and
  `paper/results/*/output` where they exist. New plots only from existing output CSVs, as presentation drafts. This is
  not Phase 8: do not regenerate paper figures, rename, or harden code.
- **Data access.** GTEx per-donor data and protected phenotypes: aggregates only.
- **Output.**
  - `paper/presentation/results_outline_2026-11-23.md`: per slide, a title as a message sentence, two to four bullets,
    the panel to show (source file), the numbers with IDs, and the caveat.
  - Draft slides or figure notebooks in `paper/presentation/` if time allows.
  - Commit; push with `git push origin main` (if refused, say so).

## Outline (about 25-30 slides)

**0. The method in one slide** (if not already in the talk's methods part).
- Each tissue read against the same person's other tissues: shared sampling-time offsets cancel.
- Model-free gate, companions (matched and stress-responsive), and the core loop without PER1, PER2 and NR1D1.
- The two layers: a fast signal layer that every cell reads (A275, A279; A343: acute dexamethasone reaches PER1
  pre-mRNA within hours), and the clock layer.

**1. Validation: the readout sees time, and sees it move** (2-3 slides).
- Against known time:
  - GTEx organ map from working-hours deaths (A123; A123b against orderers);
  - dermis against epidermis (A251 H251.4);
  - two skin sites (A342);
  - pre-mRNA lead in timed human muscle and skin (A337 Perrin, A342).
- Hospital hours: A344 (relative readouts survive, sample times degrade; the offset itself undecidable).
- Positive controls of time moving: mouse daytime feeding (A111, A322, A324), hepatocyte and forebrain Bmal1 knockouts
  (A129, A319), cell-class offsets (A273, A279), human night shift (A257).
- The record: pre-registration, blind replications, registry counts (`docs/triage/registry_2026-11-19_summary.md`).
  No new human or animal sampling.
- Caveats: the core-loop lead fails in muscle (A337 Perrin); no advantage in simulation (A145).

**2. Before birth: gates and coherence only** (3 slides).
- Clock-shaped covariation between fetuses in adrenal cortex, lung, liver; not in heart, kidney, adrenal medulla;
  brain: none at 13-22 pcw, covarying before 12 wpc in one cohort (D4, descriptive). Mouse late embryo: hepatocytes
  (A309). Sources: A237-A244, A243, A309, A240 (bounds).
- Blind replication in del Valle; the core loop covaries on its own (A328 part a). Not in single-cell cortex (A266),
  lung or liver (A269).
- Not a single slow factor (A340). The pre-mRNA lead is undecidable (A328 part b, A338, A339).
- Caveats:
  - single organ at a time; no organ-to-organ offsets possible on open data (scout 2026-11-22);
  - driven or running open (Q43);
  - postnatal appearance (A237, A239, D4).

**3. Adult health: one time per person, many volumes** (4-5 slides; full method from here).
- **Living people.** The clock state rotates with elapsed hours between biopsies (A252, A256). Tissues share one time
  (A251, A255; core loop A348; the person's own cell line does not share it, A255 H255.2).
- **Organs.** Order kept without alignment to the day (A293 H293.1); CHIRAL concordance (A123). Composition does not
  move time (A295).
- **Cell classes.** Structural cells carry the clock, lymphocytes not readably: bound about 0.7 of fibroblasts
  (A272-A288, A321, A326). Conserved relative class gains (A293 H293.4; A350 on the core loop, mural fails).
- **Loudness.**
  - Organ gains 9.2-fold, 3.1-fold without blood and testis (A104 table, derived, plan section 10).
  - Coupling mostly reads loudness (A291, A294b).
  - Person loudness shown in deaths (A312, A336 part 2), at the line in living women (A331, A348).
  - Blood is the expected exception: no readable core clock (post hoc entry 2026-11-03).

**4. Drives and inertia** (2-3 slides).
- **Intensive care: drives withdrawn.** Loudness falls (A315, A335: core loop, age / sex matched; A349 host-referenced
  log z2). The loss is largest in the strongest clocks (A349, descriptive). The order shows no large disruption: one
  organ about 6-7 h, or dispersion SD 2.25 h, would be seen (A349); offsets compressed by a shared non-temporal state
  (A109, A132b).
- **Conflicting drives.** Night shift: core loop +4.06 h, excursion about half (A257, A335); the program follows
  behaviour, 11.65 h (A261). Mouse conflict and removal: A297, A322, A324.
- **Inertia.** The clock reads the hour of death (A270, held-out waveform A335).
- **Mechanism table** (plan 11.1): what holds and what failed (A310, A311, A312, A323).

**5. Chronic disease: the desynchrony predicted from mouse is absent** (2 slides).
- T2D and obesity (A164, A235, A306), hypertension, IHD, COPD, depression, smoking (A233, A306), age and sex (A231,
  A208).
- Every null with its detection limit (A336 part 3; A235 MDE). Mouse expectation beside the human result
  (plan 11.2 table).
- Why: drives stay aligned in human chronic disease; where they conflict, time moves (section 4).
- Caveats: spread bias (A258); A306 equivalences lost on the core loop (nothing detected).

**6. Local disease: quieter, rarely shifted** (4-5 slides).
- **Tumours.**
  - They stay paired with their patient's time (A138; at the line on the core loop, A336).
  - The fibroblasts carry it (A346); endothelium and cancer cells are not resolved.
  - As tumours proliferate, they express less of that time without drifting (A304; purity strata).
- **Corrections of published claims.** COFE uncoupling is an artefact (A265, A267). Prognostic clock signatures add
  nothing beyond proliferation or random sets (A192-A194, A220). Clock-gene loss exaggerated by GTEx-normal designs
  (A221).
- **Chronic lesions** keep time (pooled +0.01 h, MDE 0.69 h) and lose loudness (at the line) (A254).
- **Psoriasis, the one local phase shift.**
  - Replicated blind (A183), and present in pre-mRNA (A345).
  - Larger than in atopic lesions (A347); the atopic delay is small, not clock-specific, and not shown to be zero.
  - Uninvolved skin also runs late (A115, descriptive).

**7. End of life, and the limits** (2 slides).
- **Brain death.** Organs keep a shared clock but lose the healthy order (A317 H317.1-2, A318, core loop A336 part 4).
  It cannot be separated from procurement in GTEx (A318 H318.2, A325, A327, B4). Descriptive.
- **Limits.**
  - Cell-level cause of lost loudness undecidable with public data (A334, A351).
  - The skin-delay mechanism is not powered (A352).
  - Coarse intensive-care resolution (A349).
  - Reanalysis only.

**8. Synthesis: loudness, not time** (1-2 slides).
- In human tissue true internal desynchrony is rare and local (psoriasis, brain death); loudness is what moves (intensive
  care, tumours, lesions); chronic disease and ageing move neither.
- The mouse mechanism holds; its desynchrony predictions do not transfer to human disease.
- The field's synchrony metrics mostly read loudness.
- Open questions: dephasing or damping (in-house time-lapse, Q44); what peripheral clocks are for; the fetal source
  (Q43); brain death with procurement samples (Q41).

## Backup slides
- Age advance of internal time (A234).
- Brain-body coupling not clock-specific (A263, A289).
- Adrenal +4.4 h offset read as an effect of death (A296-A303).
- Night-shift blood (A262).
- Tumour coherence history (A127-A173, superseded by A304).
- Chronotherapy window not readable (A249, A250).
- Method limits: spread bias (A258), c intervals approximate.
- CHIRAL age reanalysis (WP0 item 5, unregistered; decision pending).
