# Prompt for a lead agent: the two-axis checks of the paper plan (written 2026-11-19)

You continue the paper of the sccc project (method `clocksync`; see `CLAUDE.md`). The user adopted the frame of
`paper/docs/PAPER_PLAN_2026-11-19.md`, "one time, many volumes":
- every compartment is a point in its person's clock plane, time as the angle and loudness as the radius;
- time has one index (the organ), loudness has three factors (person, organ, cell class), and the departure has one
  size;
- the clock layer is built around birth, after the signal layer;
- perturbations move loudness often and time rarely.

Your job: run the assembly tasks and the checks of the plan's section 6, which settle the wording of Figures 1 and
3-7, and report. A check that removes a wording is a result: do not tune a design to keep it. The claim set is
otherwise frozen (REVIEW_2026-11-19 section 2.4): no exploratory amendment beyond this list.

## Read first
1. `CLAUDE.md` and `MEMORY.md`. All entries apply, in particular `check-done-before-proposing`, `one-cell-paper`,
   `brain-death-framing`, `compute-budget` and `subagent-concurrency`.
2. In `paper/docs/`:
   - `PAPER_PLAN_2026-11-19.md` (all);
   - `MODEL_signal_clock_2026-11-19.md` (sections 3, 4 and 6);
   - `MODEL_federation_2026-11-19.md`;
   - `REVIEW_2026-11-19.md` (sections 2.3, 2.4 and 5);
   - `SYNTHESIS_2026-11-19.md`;
   - `STATUS_2026-11-18.md` (sections 2 and 4: development and open checks);
   - PI_QUESTIONS Q32;
   - `prompts/server_mapping_fetal_2026-11-18.md`.
3. RP (`results_preliminary.md`) entries:
   - gain and departure: A291, A293, A294b, A295, A311, A312, A313, A323;
   - ventilation and brain death: A132, A315, A317, A318;
   - lesions and psoriasis: A115, A160, A183, A254, A308;
   - MuTHER loaders: A255, A259, A268;
   - gene classes and the core loop: A275, A321;
   - night shift: A257;
   - baboon: A270, and A36 in the archive;
   - development: A181, A184, A236, A237-A239, A240, A241, A242, A243, A244, A264, A266, A269, A307, A309, and the
     scout of 2026-11-18 on spliced / unspliced counts.
4. `git log --oneline -30`, `git status` and `paper/common/jobs.sh`.
   - Another session may work on the same branch: before each task, check that it is not already done or registered.
   - Commit your own files only. Three uncommitted files are inputs of this work: `PAPER_PLAN_2026-11-19.md`,
     `MODEL_signal_clock_2026-11-19.md` and this prompt. Commit them in one docs commit before the first registration.
   - Leave other untracked files alone.

## Task 1. Assembly (no new statistic; skip what is already done)
- **T3, text fixes** (PAPER_PLAN section 8).
  - Add a dated correction note under the A257 entry of RP (do not rewrite the scored text):
    - R = 0.233 is the misaligned change projected on the aligned chord;
    - the excursion ratio is R / cos s = 0.31, derived from
      `exploratory/A257_nightshift_muscle/output/score_tests.csv`.
  - Replace "amplitude to a quarter" with "excursion to about a third" in SYNTHESIS lines 21, 41 and 136 and in
    MODEL_federation line 63 (line numbers as of 2026-11-19).
- **T1, the two-axis table** (`paper/exploratory/two_axis_map/`, script `assemble_two_axis.py`, output
  `two_axis_table.csv`).
  - One row per condition of PAPER_PLAN section 3.
  - Columns:
    - time (shift in hours, or the order statistic);
    - loudness (amplitude ratio);
    - departure;
    - intervals where the outputs give them;
    - Delta and its loudness share;
    - the entry and its source CSV.
  - Conversions as in PAPER_PLAN section 10. Leave an interval blank where the outputs have none: no new resampling
    here.
- **T2, polar-plot data** (same folder, `polar_points.csv`; aggregates only):
  - organ gains and offsets (A291 and A104 tables);
  - class radii relative to endothelium (A293 Part 3);
  - person spread (A312);
  - departure size (A294b).
- **T4, registry table.** Extend `docs/triage/inventory_2026-11-18.tsv` to every entry from A109 onwards, with:
  - registration and scoring dates (git log of `registered_experiments.md` and RP);
  - blind status, as the entry states it;
  - verdict;
  - the PAPER_PLAN section 5 figure where the entry appears (or "Supplement" / "out").

  Output: `docs/triage/registry_<date>.tsv` and a counts summary. A Sonnet subagent may build it; spot-check 20 rows.
- **T5, developmental loudness table** (same folder, `development_table.csv`; PAPER_PLAN section 2.4 in numbers).
  - Per organ, cell class and stage, give:
    - the clock-shaped covariation as an estimate, or the bound relative to the adult size;
    - its layer (signal or clock);
    - the core-loop result;
    - the source CSV.
  - Sources: the A237-A244, A264-A269, A307 and A309 outputs. No new statistic.
- No final figures; quick-look plots only (the paper-figures skill comes later).

## Task 2. S1, the time-free test of the fetal state (register now, blind; score when the counts arrive)
Why: Figure 3e. Before birth, PER1, PER2 and NR1D1 covary with the core genes in the clock's pattern (A307). This is
either a small entrained loop or an input moving the whole network in that pattern; without sampling times, spliced and
unspliced counts decide (STATUS_2026-11-18 section 2). The server job (`prompts/server_mapping_fetal_2026-11-18.md`)
delivers those counts and the unseen del Valle 2023 adrenal cohort; as of 2026-11-19 they are not on disk.
- **Do not open any server count, or any del Valle 2023 file, before the S1 registration commit.**
- Check that S1 is not already registered (the server prompt calls the analyses pre-registered; nothing was found in
  `registered_experiments.md` on 2026-11-19).
- **Part a, blind replication.** The A242 / A244 statistic (signed coherence C, with expression-matched and responsive
  companions and the adult pass rate) on the 32 fetal adrenals of del Valle 2023 (7, 8, 9 and 11.5 wpc).
  - Add the core-loop variant without PER1, PER2 and NR1D1 (A264, A266).
  - Readings as in A242. A FAIL counts only with power.
- **Part b, the running-loop test** on spliced and unspliced counts of the delivered fetal cohorts: adrenal cortex,
  lung and liver where available, with heart and kidney as the expected negatives.
  - Prediction for a running loop: across individuals, the clock state built from unspliced counts leads the one built
    from spliced counts by a consistent rotation, in the clock's pattern.
  - Prediction for a steady difference in signal level between fetuses: no consistent lead.
  - Fix the statistic, the companions, the expected lead and the readings before opening any count.
  - The expected lead comes from a positive control processed by the same pipeline: the timed skin set the pipeline
    was built on (`config_DUAN_2023_SKIN.yaml`), or another timed set the user approves. If none can be processed under
    the data rules, register with the expected lead from the literature, mark the control as missing, and tell the
    user, who decides whether to add a control to the server job.
  - The lab's spliced / unspliced phase methods (card paychere2025cophaser) where they apply.
- Power first: the design (fetuses per organ, depth) against the adult or timed control.

## Task 3. G1: is the organ factor biological? (register; precondition, power first)
Why: Figure 4g. GTEx organ gains are 9.2-fold apart, and composition with conserved class gains does not explain them
(A323). If organs rank alike in living timed animals, the organ factor is an organ property. If not, a post-mortem
technical scale in GTEx stays a candidate.
- Check first whether A36, A123, A270 or A303 already compared organ amplitudes across species. A123 compared with
  CHIRAL, which is GTEx itself: not a referee.
- **Referees.**
  - Baboon: Mure 2018 (`$SCCC_DATA/baboon`; card mure2018diurnal; A36 scripts in
    `paper/results/figure_1/scripts/species/`).
  - Mouse: Zhang 2014 (GSE54652, 12 organs; card zhang2014circadian). It is not on disk: download the processed series
    matrices under the data rules. `izumo2014_elife` on disk holds figure-summary tables only.
  - Optional, descriptive, human living referee: `muscle_timed`, `adipose_timed`, `blood_timed`,
    `skin_wu2020_gse139301`.
- **Statistic.**
  - Per organ, the timed clock-plane amplitude on the GTEx waveform genes: the A36 / A270 projection, with the same
    gene set as the GTEx gains.
  - One-sided Spearman correlation with the GTEx organ gains (sqrt S, the A294b / A323 values), over organs matched by
    a mapping fixed before scoring.
  - Companion: the same with expression-matched gene sets' timed amplitudes (an organ loud for any gene set).
- **Precondition.** At least 8 matched organs whose timed clock passes (the atlas reads ZT). Power from a simulation at
  the observed organ count.
- **Readings.**
  - rho above 0 and above the companion: organ loudness is an organ property shared across species.
  - rho near 0, with power: GTEx-specific (a technical candidate, or a species difference); "organ" stays "open" in
    Figure 4g.
- **G1b, scout only.** Absolute class gains across the organs of the same donors in the GTEx single-nucleus atlas (16
  donors). Power first; register only if powered.

## Task 4. P1: is person loudness shared across tissues in living people? (register; power first)
Why: Figure 4f. A312 found that persons differ in shared amplitude and that every tissue follows it nearly in
proportion. That was in deaths, where a donor-level technical scale cannot be excluded.
- **Data.** MuTHER (`$SCCC_DATA/living_multitissue/E-TABM-1140`; A255 / A259 / A268 loaders): subcutaneous fat and skin
  of 680 women, from one punch biopsy. The lymphoblastoid lines are not a clock tissue.
- **Statistic.** The leverage of `clocksync.elasticity` with fat and skin as the two halves: corr(r_fat, r_skin) after
  both radii are residualised on the first two harmonics of the shared phase (the module does this).
- **Companions.** 200 expression-matched and 200 responsive sets through the same estimator; report the clock's
  percentile.
  - One biopsy and one array batch for both tissues: the matched sets are the decisive companion.
  - Twins: if twin-pair identifiers are public, bootstrap by pair. Zygosity contrasts only if zygosity is public;
    requesting TwinsUK phenotypes is outward-facing: ask.
- **Power.** Insert the person-amplitude spread estimated in A312 into the MuTHER design with `clocksync.simulate`
  (arrays carry a larger technical floor, A259). Register only if power is at least 0.8 at that size; otherwise record
  "not decidable at this design".
- **Readings.**
  - Clock above the matched q95: persons differ in a loudness shared across tissues in living people, beyond a shared
    technical scale.
  - Within the matched range: not separable from a shared technical scale.
  - Clearly below the A312 size with power: the A312 leverage is specific to deaths.
- **Covariates.** Age, and BMI if public, as descriptive strata, not regressors.

## Task 5. N2: is the departure the signal layer read through PER1, PER2 and NR1D1? (register; one ~25 min run)
Why: model part 6 and Figure 4e (MODEL_signal_clock section 6, N2).
- **Blind check.** The core-loop departure of the A294b tissues must not have been computed before. A296 used the core
  loop for the adrenal offset only, A321 for lymphocytes. If it was computed, the check is not blind: say so in the
  registration.
- **Design.** A294b (`exploratory/A294_floor/scripts/score_a294.py`; GTEx fast deaths, 27 peripheral tissues), refitted
  on the 11-gene core waveform without PER1, PER2 and NR1D1 (A321's set). Each gene set is compared against 100 matched
  sets of its own size.
- **Statistics.**
  1. The H294b.1 excess on the core loop against the full set (0.969), paired over tissues.
  2. Coupling on the core loop against the full set, paired over tissues.
  3. Within tissue, the projection of each sample's departure on its responsive state (A275 glucocorticoid and
     immediate-early classes), against matched sets.
- **Readings (registered).**
  - Leakage: the core-loop excess falls toward 0, coupling rises and the projection is positive. Part 6 then becomes
    "the acute signal state seen through signal-inducible clock genes", and clocksync gets a core-loop readout for
    coupling.
  - Otherwise: the departure is clock-specific, with its cause open.
- Run through `paper/common/run.sh` with 4 workers.

## Task 6. N1: without behaviour, do outputs follow the clock? (archive first; register only what is new)
Why: the one human test of the two-layer claim (Figures 1 and 7).
- **Overlap with the archive.**
  - A91 (`paper/results/figure_3/scripts/gtex_cues/compute_internal_external.py`; archive RP
    `scCircadianCancer/docs/RESULTS_preliminary.md`, A91 near line 2565) found that rhythmic genes following internal
    time in fast deaths keep a much higher ventilated lock than wall-time genes (H91, median difference +1.23,
    p 0.0005).
  - A27 D5 counted internal-only genes in ventilated donors.
  - Read both, and the A81b / A84 / A89 lock tables. If they answer N1, write the reading from them (no new test) and
    close this task.
- **Otherwise, register the per-tissue retention ratio** (MODEL_signal_clock section 6):
  Phi_o = (L_out ventilated / L_out fast) / (L_clock ventilated / L_clock fast).
  - Output program: A91's clock-following genes, excluding clock and responsive genes, selected on a held-out half of
    the fast deaths.
  - Companions: matched non-rhythmic sets and the responsive family.
  - Prediction: Phi near 1, except in meal-driven tissues (liver lowest; archive A33).
- **Power first** at the real design: 79 ventilated post-mortem donors.
- **Covariates.** Balance death type, age, RIN and ischaemic time between fast and ventilated deaths by matching or
  stratifying; report pooled and matched.
- **Registered limits.** The central drive persists in ventilated donors, so Phi counts centrally driven outputs with
  the clock. Continuous feeding is not recorded.

## Task 7. Descriptive companions (register as descriptive, NOT BLIND; one run each)
- **D1, psoriasis and eczema loudness.**
  - The A254 H254.3 statistic (clock-plane variance, lesion over non-lesion, log ratio r;
    `exploratory/A254_paired_lesions/scripts/compute_a254.py`) on GSE154200, with the A183 preparation (psoriasis,
    nummular eczema, atopic dermatitis).
  - Matched and responsive companions on the same per-person swaps.
  - Replication: the four A115 psoriasis series, lesional against uninvolved (`skin_lesional/processed_a115/`).
  - A183 reported r +0.50 for psoriasis without companions: not blind.
- **D2, brain death per organ.** For each of the 12 A317 tissues:
  - the phase change: circular mean over the organ's pairs of the organ-donor minus ventilated post-mortem pair phase,
    with a donor bootstrap on the A317 draws;
  - the loudness ratio: A315 `tissues.csv`, L_t against the matched spread.

  These are the arrows of Figure 5; there is no test.
- **D3, ventilation time.** The order by tertile (A132) beside loudness by tertile (A315 statistic), in organ donors
  and, descriptively, in ventilated post-mortem donors: "loudness spent, order held". GTEx phenotypes: aggregates only.
- **D4, the developmental loudness curve.**
  - Per organ and developmental bin (prenatal windows; infancy and childhood; adolescence and adult): the clock-plane
    variance against matched sets, beside C with its adult pass rate.
  - Data: Cardoso-Moreira 2019 (A237 loaders: seven organs, prenatal and postnatal series) and Kang 2011 (A239: brain).
  - Bins are fixed before the run.
  - Caveat, as in A240: the spread of sampling times differs between bins.
  - Not blind for A237 and A239. Shows when loudness appears (Figure 3c).

## Reports and stops
- **When to report:** after the git check, after Task 1, after the S1 registration, and after each power stage and each
  scoring.
- **Power below 0.8** at the registered size: do not score; record "not decidable at this design" in RP.
- **Data that are not public and on disk, or a download beyond the processed-data rule:** stop and ask. The server job
  of S1 is the one approved exception (user, 2026-11-18).
- **At the end:**
  - update PAPER_PLAN sections 2.4, 3, 5 and 6 (status of each check, wordings that change);
  - write a PI question (next free Q number) with one table: check, verdict, key numbers with their source CSV, and the
    figure and wording each changes;
  - then stop.
- **Not in scope:** new exploratory amendments (claim freeze), final figures, the outward-facing items of PAPER_PLAN
  section 6.4.

## Rules that bite
- **Pre-registration.**
  - Hypothesis, statistic, readings, seeds and what has been seen are committed before scoring.
  - Use staged registrations when the design depends on counts.
  - A test that is not blind is a "registered re-test"; post hoc work is labelled as such.
- **Controlled access.** GTEx time of death and phenotypes: aggregates only, never per-donor rows. Never commit
  controlled-access data.
- **Companions and gates.**
  - Companions run on the same permutations.
  - The matched family is calibrated in bulk (A285). The responsive family is reported, using A320's calibrated broad
    family where living-tissue specificity is claimed.
  - A failed gate is uninformative.
  - Compare coupling across groups only at a common spread (d*) or on z2 / psi.
- **Covariates** that change rhythm strength (death type, ventilation, cohort, RIN, age) are balanced by matching or
  stratifying, not regressed out; report pooled against matched.
- **Compute.**
  - Long runs only through `paper/common/run.sh`, detached, with `watch_jobs.sh` running.
  - At most 16 workers and about 120 GB of summed RUN_CAP_GB; keep 20 GB free.
  - Sonnet subagents for clearly defined tasks, at most two active; verify their numbers against the output files.
- **Downloads are untrusted:** own folder, PROVENANCE.md (URL, date, sha256), parse with `python -I`, never execute.
- **Code.** Use the clocksync public API; add to it only tools useful beyond this project.
- **Numbering and records.**
  - Next free amendment number: A328 (check the end of `registered_experiments.md`).
  - After each scoring: an RP entry (headline with verdicts, design, power, readouts, reading as registered, limits,
    code line, `<!-- src: -->`) and a row in `VERDICT_TABLE.md`.
  - Commit with the attribution lines.
- **Wording.**
  - "Loudness" (or gain) for amplitude relative to the person.
  - "Time" for relative timing between compartments.
  - "Desynchrony" only for phase relations.
  - "Quieter", not "disrupted".
  - Brain death as in the SYNTHESIS wording rules.
  - Before birth: "a driven state carried by signal-induced genes", and "an adult-size core loop excluded down to x of
    the adult size; a smaller entrained loop is not excluded".
