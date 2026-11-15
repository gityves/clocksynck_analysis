# Prompt for a lead agent: the four checks that decide the biology paper (written 2026-11-15)

You lead four checks for the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab). The target
line (plan 1, "what internal time controls") ran out of public designs on 2026-11-15 without a replicated headline
(RP A259 to A262; `PI_QUESTIONS.md` Q33, Q34). The user's working decision: write the paper as a **biology paper**
around one thesis, and let these checks decide the venue.

**Thesis: the human clock network is robust, and the places where it breaks are specific.** A person's tissues share
one clock state; it holds across sex, age and chronic disease (group offsets within 0.2 to 0.6 h); it is absent or
loosened in a short list of situations: brain, blood, the fetus except the adrenal cortex, malignant cells against
their microenvironment, brain-dead donors.

Each exception currently rests on one design. Your four checks test whether they survive. If they do, the biology
paper is the stronger one; if two or more fall, the same material is reordered into a Nature Methods paper. Both
outcomes are useful, so a check that removes an exception is a result, not a failure: do not tune a design to keep
one. Work autonomously; report as listed under "Reports and stops".

## Read first
1. `CLAUDE.md` and the memory index (`MEMORY.md`): all entries apply, in particular paper target plans, PI questions
   file, machine budget, blood has no coherent core clock, brain death stays descriptive, estimable bulk designs,
   fetal follow-ups A243/A244, Phase 8 needs PI review.
2. `paper/docs/STATUS_ROADMAP_2026-11-14.md`: section 1 (status, method limits) and the items A2, A3, A4, A7 and T4.
   Its progress note lists what the previous lead agent did.
3. `paper/docs/PI_QUESTIONS.md`: Q32, Q33, Q34. `paper/docs/TIME_CONTROL_INVENTORY_2026-11.md` (for check 4).
4. `paper/docs/prompts/target_line_what_internal_time_controls_2026-11-15.md`: the sections "Rules that bite here"
   and "Pitfalls" apply here unchanged and are not repeated.
5. The entries named under each check below, in `paper/docs/results_preliminary.md` (RP) and
   `paper/docs/registered_experiments.md`.

Before starting, run `git log` and `git status`: another agent may still be working on the same branch (the previous
lead planned T4). If a check below is already registered or in progress, do not duplicate it; take the others and
say so in your first report.

## Check 1. Fetal: is the adrenal signal a running clock, and is the adrenal first? (roadmap A2, A7)
Where it stands: the fetal adrenal cortex carries clock-shaped covariation between fetuses (A244: cortex C 0.608,
n 34, medulla -0.048; A242 replication: 0.552 late, 0.651, 0.732 early;
`exploratory/A244_fetal_adrenal_classes/output/verdicts_a244.csv`,
`exploratory/A242_fetal_adrenal_replication/output/verdicts.csv`). Two weaknesses: a static axis between fetuses
would pass the same test, and the signal is strongest before 12 weeks; and the other organs have not been tested at
the adrenal's size.

- **1a, new amendment (roadmap A2).** Does the cortex statistic survive (i) adjustment for gestational age and
  cohort beyond what A242 / A244 already do, (ii) adjustment for a steroidogenic activity score (candidates STAR,
  CYP17A1, NR4A1; fix the genes and the score from the literature and from adult adrenal, not from the fetal clock
  statistic), (iii) removal of PER1 and NR1D1, the two clock genes most directly induced by glucocorticoid and
  steroidogenic signalling? Not blind: every fetal adrenal statistic of A189, A242 and A244 has been seen; label it
  so. Fix the reading of each outcome in advance, including "the covariation is explained by steroidogenic
  activity", which would move the fetal result to a hypothesis in the supplement. Adult adrenal power for each
  adjusted statistic first (the A242 feasibility scripts do this).
- **1b, score A243 as registered (roadmap A7).** Staged registration: write the design addendum (cohorts, adult
  power), commit it, then score once. Full instructions and pitfalls:
  `paper/docs/prompts/fetal_pooled_organs_and_adrenal_classes_2026-11-11.md` (the A243 part). It may add at most
  60 GB; 434 GB of the 500 GB cap were used on 2026-11-15, so report sizes before downloading.
- Code: `paper/common/development.py`, `exploratory/A242_fetal_adrenal_replication/scripts/`,
  `exploratory/A244_fetal_adrenal_classes/scripts/`.

## Check 2. Tumour: head-to-head with COFE on the same pairs (roadmap A3)
Where it stands: tumour and adjacent normal tissue pair on the clock state (TCGA, 517 pairs: 0.365 against a null
mean of 0.122, pseudo-clock q95 0.260; `exploratory/A138_cancer_estimability/output/score/tests.csv`), and tumour
myeloid and T/NK cells keep the patient's state (A227). COFE (Ananthasubramaniam and colleagues; reference key
`ananthasubramaniam2025rhythm` in `docs/NOVELTY_ASSESSMENT_2026-11-06.md`) orders surgical cohorts without times
and is cited for the opposite reading, tumour clocks uncoupled from the host. A reviewer will ask which is right.

- First read the COFE paper's actual claim and design and state it exactly, with the cohorts used; verify it at the
  source, not from our notes.
- COFE is already wrapped in the project: `exploratory/A123_working_hours/scripts/compute_methods.py`,
  `exploratory/A145_two_step/scripts/lib.py`, `exploratory/A138_cancer_estimability/scripts/q2_gate*.py`. RP A123b,
  A145 and A204 already compare the methods on organs and simulations (in hospital hours every orderer loses the
  organ map; A204: coupled against uncoupled AUC 1.0 against COFE 0.666).
- New amendment: on the TCGA tumour-normal pairs, (i) reproduce COFE's readout as published; (ii) ask whether COFE's
  phases of a patient's tumour and normal sample agree beyond permuted pairs, the question A138 answers with
  `clocksync`; (iii) a simulation on the real design with tumour and normal coupled, and with them uncoupled,
  sampled in surgical hours: which method reports what. The simulation is the decisive part: it says whether the
  disagreement comes from the narrow sampling window. Fix the readings in advance, including the one where COFE is
  right and A138 is carried by something other than the clock (then check the A138 companions again).
- Loader: `common.data.load("tcga_paired")`.

## Check 3. Companions for the organ map and the brain exception (roadmap A4)
Where it stands: A253 showed that pairing between brain regions in controls is matched by stress-responsive gene
sets, and A259 that matched non-clock sets carry 0.30 to 0.74 of a tissue's apparent deviation. The brain exception
(GTEx brain block r\* 0.639 [0.530, 0.728] against a reference median of 0.978, with sheep and baboon;
`exploratory/A142_organ_blocks/output/tests.csv`, `results/figure_2/output/organs/r2_readout.csv`) and several atlas
fits were scored before the stress-responsive companion existed.

- Inventory first, no statistic: for every result planned for a main figure (roadmap section 1.2 and 1.3; the
  atlases of `results/figure_2/scripts/atlas_fits/`), which companions it has (expression-matched pseudo-clock,
  stress-responsive: `clocksync.paired.pseudo_waveforms(..., kind=...)`, `clocksync/_core/controls.py`), and which
  it lacks. Write the table into the amendment.
- New amendment: add the missing companions on the same permutations, with one rule for all: a result is
  clock-specific if it exceeds the q95 of both families. For the brain block the question is directional: is the
  **lower** agreement of brain with the body specific to the clock genes, or do matched and responsive sets show
  the same drop? State in advance what each outcome does to the claim "brain is the exception in three species".
- These are reruns of scored data: not blind, labelled so; original verdicts are not changed, the companion column
  is added, and a result that loses specificity gets an erratum in RP and moves in the roadmap's status tables.

## Check 4. The systemic signal (roadmap T4)
The question: what holds a person's tissues on one clock state? Candidates are glucocorticoids and body temperature.
A named signal is the mechanism the paper lacks; a clean negative says the public data cannot name it.

- Where it stands: the glucocorticoid score of the archive explains little of the host state in fast deaths (R^2
  0.014) and more in ventilated donors (0.234), and nothing of the hour of death (0.015, 0.011;
  `results/figure_3/output/gtex_cues/diag6b_gc.csv`). A136 (donor-management glucocorticoids): gate fails. A248: no
  T-cell night peak in GTEx blood. The inventory keeps the score of A27 D6b as the starting reader.
- Scout first (`paper/exploratory/scouts_2026-11-15/T4_systemic_signal/`), no hypothesis scored: build per-sample
  readers from genes fixed outside the human data. Glucocorticoid targets (FKBP5, TSC22D3, KLF9 as candidates) and
  temperature targets (CIRBP, RBM3, heat-shock genes), selected among genes that stay rhythmic in mouse liver
  without a local clock (`manella2021`, `jouffe2022_bmal1ko`, the A260 classes). Check in the mouse data that each
  reader tracks its signal with the expected phase. Then power: can a fixed phase relation between a reader and the
  host state be told from none, given that agonal stress drives both in GTEx? The 14 clock genes are excluded from
  the readers; PER1 is a glucocorticoid target, so report every result with and without it in the host.
- Register only what the scout supports. The registered question: does the reader keep a fixed phase relation to the
  host state across tissues and people, clock-specifically (against readers built from matched and from
  stress-responsive genes), and is the relation looser in brain? A living confirmation is required for a positive
  to count: MuTHER fat and skin (`$SCCC_DATA/living_multitissue/`) or the repeat-biopsy cohorts
  (`$SCCC_DATA/repeat_biopsy/`), where the reader should advance with the elapsed hours as the clock state does
  (A252 machinery).
- Extension only after the main result, and only after counting the users: GTEx long-term steroid users (protected
  phenotype, aggregates only).

## Order and lanes
Checks are independent. Suggested: start 1b's scout and 3's inventory through at most 2 Sonnet 5.5 subagents (clearly
defined, metadata and tables only; verify what they return), while you write the amendments for 1a and 2. Check 4's
scout follows when a lane is free. Every amendment is committed before its data are scored; next free number at
hand-over is A263 (check the end of `registered_experiments.md`).

## After each scoring
RP entry (headline with verdicts, power first, each readout, reading as registered, limits, code line,
`<!-- src: -->`); a row per readout in `VERDICT_TABLE.md`; `METHOD.md` when a tool or calibration is new; the
roadmap's progress note and, when a standing claim changes, its status tables (numbers from output tables only).
Commit your own files only, with the attribution lines.

## Reports and stops
- Report without waiting: after the first `git` check (which checks you take); after each inventory or scout with
  what you intend to register; after each scoring.
- Stop and wait: (1) when all four checks are scored, or when two exceptions have fallen: write the tally as a new
  question in `PI_QUESTIONS.md` (after Q34): which exceptions stand, which fell, and your recommendation between
  the biology paper and the Nature Methods ordering. The choice of plan and journal is the user's and the PI's.
  (2) Anything outward-facing (requests to authors or data owners, controlled-access data): drafts only.
- A report states verdicts with their numbers, the power behind every negative, what is not shown, and what you do
  next.

## Not in scope
Further plan 1 scouts; brain-death mechanism (replication only, and not part of this prompt); the in-house jet-lag
data (the user's and the PI's decision); Phase 8 (naming, figures, tutorials); rewording roadmap section 2.

## State at hand-over (2026-11-15)
Last commit `e0377aa` (A262 scored, Q34). Registered and unscored: A243. No background job. Tracked files clean
apart from the commit that adds this prompt.
