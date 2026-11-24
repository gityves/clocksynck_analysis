# Prompt: is the person's clock volume a body-wide loss of loudness from conflicting cues? (2026-11-24)

You are working in `/home/maxine/Documents/paychere/scCircadianOverview` (read `CLAUDE.md` first: pre-registration,
companions, covariate balancing, controlled-access rules, `paper/common/run.sh` for long jobs, commit with the
attribution lines, never push). Work with the user step by step: register, run the power stage, report, then score.
Ask before anything outward-facing.

## 1. The idea

Persons differ in how loudly all their organs express the shared clock time: a "volume" shared by every organ, while
the organs stay on one time (each organ's personal departure from the body's time is small and of fixed size). The
paper's rule is that disagreement among the units a measurement averages (cells, cues converging on an organ) appears
as lost loudness. Reading: a quieter person is a candidate for body-wide internal desynchrony below the resolution
(cells or cues), not organs drifting apart. The user's framing: a genetic basis (to test), the mouse mechanism
(conflicting cues cost amplitude), and measurements in human deaths and living cohorts.

Wording ceiling (keep it in every text): even if both tests pass, the claim is "a body-wide loss of loudness, as
conflicting cues produce"; "desynchrony of cells" needs single-cell evidence that public data cannot give (A334).

## 2. What is established (read these entries in `paper/docs/results_preliminary.md`)

| link | entries | status |
|---|---|---|
| persons differ in a volume shared by their organs | A312 (leverage 0.75-0.86; elasticity), A336 part 2 (core loop 0.714 above matched q95 0.621 and responsive q95 0.692), A331 (living women, fat and skin 0.177, above matched, at the line with responsive genes; the size A331's simulation predicted for its design) | have |
| organs stay on one time | A294b, A311, slide 3.4b (`paper/presentation/results_outline_2026-11-23.md`) | have |
| conflicting cues cost amplitude (mechanism) | A322, A324 (mouse feeding inversion; most where organ weights are balanced), A297 (removing a cue keeps amplitude) | have, mouse |
| within a person, in humans | A257 / A335 (night shift: muscle core loop shifts +4.06 h, excursion about 0.49) | have, one tissue |
| exposure and systemic state | A315 (intensive care lowers amplitude; within organ donors lower with ventilation time, descriptive), A360 (volume falls with the glucocorticoid-target state, partial R^2 0.061; RIN and ischaemic time do not explain it) | partial |
| the elasticity level | A312 calibration note (post hoc): the estimator reads about 1.1 for proportional scaling; the H312.1 gradient is not shown beyond the estimator | caveat |
| genetic basis | not tested. NOTE: A331's "monozygotic 0.195, dizygotic 0.152" is the fat-skin leverage computed within each group of women, not a co-twin resemblance; it is not heritability | missing |
| cells desynchronised rather than weaker | A334 (public 10x cannot decide) | missing, out of reach |

Context documents: `paper/docs/COUPLING_SOURCES_2026-11-24.md`, `paper/docs/PAPER_PLAN_2026-11-19.md` section 2.1
(person volume A_p, "Where coupling comes from"), `paper/docs/MODEL_federation_2026-11-19.md` part 5.

## 3. Tasks

### T1. Heritability of the clock volume (MuTHER twins) — register as A361

- Data: MuTHER (E-TABM-1140), fat and skin of the same women, twin pairs with zygosity. Loader:
  `exploratory/A255_living_tissues/scripts/compute_a255.py` `load_muther` (used by
  `exploratory/P1_person_loudness/scripts/score_a331.py`, which already reads `cotwin` and `zygosity`).
- Statistic to design and fix before scoring: the co-twin resemblance of the volume, separately in monozygotic and
  dizygotic pairs, and h^2 = 2 (r_MZ - r_DZ) with a pair bootstrap. Prefer cross-tissue, cross-twin covariances (fat of
  twin 1 with skin of twin 2, and the reverse) via `clocksync.elasticity.radius_covariance`, so that one biopsy and one
  array per woman do not share noise with the co-twin's measurement; radii residualised on the phase harmonics as in
  `radius_covariance`.
- Companions: the same statistic on expression-matched and stress-responsive pseudo-clock sets (A331's construction);
  heritability of generic expression scale would show up there too.
- Check before writing the registration: whether co-twins were sampled on the same visit or array batch (a shared
  batch inflates MZ and DZ resemblance alike; report it as a limit and look for batch fields).
- Power stage first (simulation at A331's design and person spread, pair counts from the data); no clock statistic
  before the registration is committed.
- Readings to state in advance: high h^2 with no exposure effect (T2) = a constitution (trait), not desynchrony;
  exposure effects without h^2 = a state; both = inherited susceptibility with exposure-driven loss of loudness.

### T2. Exposures and the volume (GTEx) — register as A362

- Person volume: host radius per donor (`Dataset.host`, every peripheral tissue an anchor, as A360:
  `exploratory/A358_A360_gain_sources/scripts/a360.py`), residualised on time of death (cos, sin) and the number of
  tissues; companions from matched pseudo-clock hosts (A336 part 2 built them: see
  `exploratory/A336_two_readouts/output/person_gtex*.csv` and its script).
- Exposures (protected phenotypes, `common.gtex.phenotypes`; dbGaP variables): time on the ventilator (`DTHVNTD`,
  `DTHVNTDU`) as a dose-response within ventilated donors (organ donors and Hardy 0; A315 / A349 define the groups);
  long-term steroid use (`MHSTRDLT`; a chronic glucocorticoid signal, predicted lower volume); major depression
  (`MHDPRSSN`; predicted lower volume); season of death (`DTHSEASON`; descriptive). No shift-work or sleep variable
  exists in GTEx.
- Balance (CLAUDE.md covariate rule): match or stratify every contrast on age bracket, sex, Hardy class, cohort and
  number of tissues; report pooled and matched estimates. Every contrast with its companion on the same permutations.
- Blind status to declare: ventilation time is NOT BLIND (A315's descriptive tertile contrast, Q* -0.236; A332 D3);
  steroids, depression and season against the person volume have not been computed (check with grep before claiming
  blind).
- Power stage first (permutation nulls and matched-set spread; detectable effect per exposure given its prevalence).

### T3. Planning note

Write `paper/docs/PERSON_VOLUME_2026-11-24.md`: the chain of section 2 with the outcomes of T1 and T2, what each
outcome means for the reading (trait, state, both), the wording for the paper and for slide 3.5 / section 4, and what
remains out of reach (cells; a within-person repeat design). Propose plan edits to the user; do not apply them
unasked.

## 4. Rules for this work

- GTEx time of death and phenotypes are controlled access: aggregates only (coefficients, correlations, group sizes
  >= 5); never print or commit per-donor rows. MuTHER is open but report aggregates.
- Numbers in text come from output tables. Post hoc choices are labelled. A failed gate is uninformative.
- Long runs through `RUN_CAP_GB=30 paper/common/run.sh ...`; budget at most 16 workers.
- Results go to `results_preliminary.md` (power note, then the scored entry), code under
  `paper/exploratory/A361_volume_heritability/` and `paper/exploratory/A362_volume_exposures/`.
- Stop and ask the user after each power stage, and before changing the paper plan or the slides.
