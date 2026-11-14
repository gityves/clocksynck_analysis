# Prompt for a fresh lead agent: the target line, "what internal time controls" (written 2026-11-15)

You lead the next stretch of work on the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab).
The method, `clocksync`, reads the clock state that a person's tissues share (the **host state**) from 14 clock genes
without any time stamp. The user set the paper's target on 2026-11-15: answer **what internal time controls in a
human tissue**: which daily processes follow the body's time, which follow the tissue's own clock, which follow
behaviour or the external day, and which follow none. If that is answered it is a tier-one paper (plan 1). If not,
the paper falls back on the biology as it stands (plan 2: one internal time per person, and where it is absent), and
if that has no headline, on the method (plan 3, Nature Methods).

Your job is the first block of the roadmap: an inventory of what the project already knows on this question, the
scouts that decide which designs are worth registering, the bound in hours that decides whether plan 2 has a
headline, and then the registered analyses that pass. Work autonomously. Report to the user at the points listed
under "Reports and stops"; do not wait for an answer unless that section says so.

A negative outcome is a result here. An underpowered design that cannot separate the times is to be found by the
scout and dropped, not registered and scored. Do not bend a design to make plan 1 succeed: the fallbacks exist so
that you never have to.

## Read first
1. `CLAUDE.md` (repository root): rules, budget, data policy, long-job procedure, git (local commits only).
2. Memory index (`MEMORY.md`), in particular: paper target plans, current handoff, PI questions file, machine budget,
   internal desynchrony term, blood has no coherent core clock, brain death stays descriptive, estimable bulk
   designs, journal framing, Phase 8 needs PI review.
3. `paper/docs/STATUS_ROADMAP_2026-11-14.md` **in full**. It is the plan: status with source tables (section 1), the
   three plans and gates G1 to G3 (section 2), the items A1 to A7, T0 to T6, B1 to B4, C1 to C7 and the order of work
   (section 3). Item names below refer to it. Its literature references are from memory: verify any reference before
   it enters another document.
4. `paper/docs/METHOD.md`, then the public API in `clocksync/src/clocksync/` (`dataset.py`, `paired.py`,
   `rhythms.py`, `decompose.py`, `power.py`, `simulate.py`, `_core/controls.py`).
5. `paper/results/figure_3/scripts/gtex_cues/README.md` **in full**: the existing cue line (scripts, amendments,
   outputs, checks against the legacy).
6. `paper/docs/results_preliminary.md` (RP), these entries: WP3 (A111, A112, A129), A122, A124, A136, A208, A231,
   A232, A233, A234, A235, A248, A249, A250, A252, A255, A256, A257, A258. For archive amendments (A1 to A108c):
   registrations in `/home/maxine/Documents/paychere/scCircadianCancer/docs/PLAN_internal_desync.md`, results in
   `.../docs/RESULTS_preliminary.md`, and `paper/docs/AMENDMENT_INDEX.md` maps each to its migrated script.
7. `paper/docs/registered_experiments.md`: A252, A256, A257 and A258 as models of a registration, a power note and a
   calibration; the header format of an amendment.

## The question, made operational
Four times can be told apart within one person:

| time | read from | what separates it from the others |
|---|---|---|
| the body's internal time | host state shared by the person's compartments | differs from the recorded hour by chronotype and by a shifted day (GTEx: host state against hour of death) |
| the tissue's own clock | the tissue's clock state where it departs from the host state | person-to-person offsets; compartments with a weak or absent clock (tumour, blood, cultured cells) |
| behaviour (sleep, meals) | protocol | night-shift protocol: behaviour moves 12 h, the muscle clock 2.7 h, central time about 0 h |
| the external day | recorded hour | hours elapsed between two biopsies; hour of death |

An answer is a **partition with intervals**: per tissue and per program (not per single gene), the share of daily
variation that follows each time. Plan 1 is taken only if gate G2 holds: the partition is (i) specific: the body and
tissue readers beat expression-matched pseudo-clock and stress-responsive readers; (ii) calibrated: the same
procedure, run untimed, recovers the known classes in mouse liver with and without a local clock; (iii) replicated
in at least two independent designs or tissues; (iv) quantitative. "Mostly imposed by the body" and "mostly local"
are both answers.

## Where it stands (numbers re-read in the tables named; paths relative to `paper/`)
- **Internal against external time, per gene (A91, archive).** `results/figure_3/output/gtex_cues/a91_verdict.csv`:
  H91 passes over 12 tissues (median difference 1.232, p 0.0005). `a91_tissues.csv`, skeletal muscle (220 donors,
  9,313 genes): 1,480 clock-following, 113 time-following, 676 both, against noise floors of 37, 31.5 and 8.5. The
  output is dated 25 October, older than several changes to `clocksync`, and has no companion sets. Script:
  `results/figure_3/scripts/gtex_cues/compute_internal_external.py` (per tissue and gene, log(CPM + 0.5) on the D2
  covariates, the host state of the other tissues, and cos/sin of the hour of death; partial F for each).
- **Behaviour against the muscle clock (A257).** `exploratory/A257_nightshift_muscle/output/score_tests.csv`: after
  a 12 h shift of behaviour the muscle clock moves 2.68 h (90 % 1.81 to 3.59), excursion ratio 0.233; 11 subjects.
  Data: GSE106800, arrays at 07:00 and 19:00 in each condition (44 arrays), in `$SCCC_DATA/geo/GSE106800/`. The
  gene-set attempt on the same data (A124) had a detection limit of 7.3 h.
- **Rotation with elapsed hours (A252, A256).** The clock state advances between two biopsies of one day, clock
  specifically, and replicates in fresh studies (`exploratory/A252_elapsed_time/output/score_tests.csv`,
  `exploratory/A256_elapsed_replication/output/score_tests.csv`). Data: `$SCCC_DATA/repeat_biopsy/` and
  `$SCCC_DATA/living_multitissue/`.
- **A compartment without the donor's state (A255).** Fat and skin of a living woman pair at r 0.533 (680 women);
  her cultured cell line does not (fat 0.042, skin 0.066): `exploratory/A255_living_tissues/output/score_tests.csv`.
- **The earlier cue line is mostly negative or unspecific.** Cue reference A81 / A81b: 40 cue and 40 locked genes,
  frame check 0.481 against a guard of 0.8, used with a post hoc override. A89 (cue hierarchy over organs) fails
  (0.915 against 0.870, p 0.269). A85 (immune cells follow cues) is not cue-specific (A122). A112 (cue against
  locked genes under arrhythmic feeding in mouse liver) fails with the wrong sign. Glucocorticoid score
  (`results/figure_3/output/gtex_cues/diag6b_gc.csv`): R^2 on the host state 0.014 in fast deaths and 0.234 in
  ventilated donors; on the hour of death 0.015 and 0.011. A136 (donor-management glucocorticoids): gate fails.
  A248: GTEx blood does not show the published night peak of T cells.
- **Tumours: not a finding.** Clock-controlled genes of the tumour followed the body time of the adjacent normal
  sample in TCGA (A249, descriptive) but not clock-specifically in CPTAC, where NFIL3 alone replicated (A250).
- **Offsets in hours across groups (already scored).** Old minus young -0.05 h (90 % -0.27 to 0.18, detection limit
  0.20 h; `exploratory/A231_age_coupling/output/score_tests.csv`); women minus men 0.04 h (-0.10 to 0.19; 0.17 h;
  `exploratory/A208_sex_coupling/output/tests.csv`); type 2 diabetes in fast deaths 0.02 h (-0.25 to 0.30; 0.37 h)
  and obesity 0.20 h (-0.23 to 0.64; 0.64 h; `exploratory/A235_metabolic_offset/output/score_tests.csv`); all
  equivalent within 1 h. Five chronic conditions: no reordering of organs (RP A233). The pseudo-clock companions of
  these offsets are wide (|offset| q95 1.3 to 10 h over the four amendments).
- **Method limits you must design around** (roadmap 1.4): absolute coupling is attenuated and its intervals
  undercover; group contrasts of coupling are biased by unequal host spread and label permutation is
  anti-conservative (A258; `clocksync` has `contrast(host_prior="arm")`); the gate passes in a clock-knockout aorta,
  so a claim rests on pairing plus matched and stress-responsive companions, not on the gate alone.

## Order of work

### Step 1. Inventory of the existing cue line (T0a). No new statistic.
Read the registrations, results and migrated scripts of A27 D5 to D8, A46, A47, A81 / A81b, A85, A89, A90, A91,
A92 / A93, A100 to A102b, A112, A122, A124, A136, A248. Write `paper/docs/TIME_CONTROL_INVENTORY_2026-11.md`: one row
per analysis (question, design, verdict with the number and its table, companions used or missing, whether the
output predates current code, what it contributes to the partition, why it failed if it did), then one page on what
the failures have in common and what a design must do to avoid them. This document decides what phase T reuses.
Subagents may read and summarise; you verify every number against its table.

### Step 2. A91 on current code (item A5). A rerun, not a new hypothesis.
Run `compute_internal_external.py` as registered, with `--subsample` first, then in full through `run.sh`. Compare
with the archived output (class counts per tissue, noise floors, H91). The verdict is whatever the registered rule
gives on current code; a change is an erratum in RP, not something to tune away. Keep the archived tables (write
the rerun beside them or keep a dated copy) so the comparison stays reproducible.

### Step 3. Scouts for T1 to T3 and the calibration set (T0b, T0c). Simulation and design only.
Scouts live in `paper/exploratory/scouts_2026-11-15/` with a short `README.md` of findings. A scout may read
designs, sample sheets, clock-gene reliabilities, null distributions and simulated data. It must not compute the
statistic a later amendment will register on the real labels.

- **Calibration (T0c).** Mouse liver with and without a local clock, fed in and out of phase:
  `common.data.load("manella2021")` (GSE159135: hepatocyte Bmal1 knockout against control, ad libitum against
  daytime feeding, liver, fat and lung of the same mice, 12 time points, 2 mice each; read `build.py` for the
  assumption on mouse identity), `greenwell2019`, `jouffe2022_bmal1ko`, `xin2021`, `vipko_gse296764`. The recorded
  times give the truth: genes that stay rhythmic in the knockout liver are body- or feeding-driven; genes that lose
  their rhythm depend on the local clock. Hide the times, take the other tissues as host, and ask whether the
  procedure you intend for humans puts the genes in the right class. This is G2 (ii), and it is the cheapest way to
  learn whether the idea works at all, so do it before any human design. Check first what A111 and A129 already did
  with these data.
- **T2 scout (GTEx, three-way regression).** The difficulty is identification, not power. (1) The tissue's deviation
  is measured in the same RNA sample as the gene, so library-wide technical and agonal variation makes genes
  "follow" it; the floor is the same regression with deviations built from matched pseudo-clock and
  stress-responsive sets. (2) The host state from the other tissues is a noisy measure of body time, and the
  tissue's own 14 clock genes measure body time too, so a gene that purely follows body time will load on the
  deviation: the regression overcalls local control unless the host's reliability is modelled. Quantify both on
  simulated data with the real design (`clocksync.simulate`), and estimate how much of the deviation is real (for
  example from two disjoint halves of the clock genes) before deciding whether T2 is registrable and in which
  tissues. Clock genes and their companions are excluded from the genes scored.
- **T1 scout (night-shift muscle, genome-wide).** With two clock times 12 h apart, a gene following behaviour
  reverses its evening-minus-morning difference and a gene following central time keeps it; a 2.7 h shift lies in
  between and cannot be told from 0 h for a gene of unknown phase. Gene phases and amplitudes from a timed human
  muscle reference (`common.data.load("muscle_timed_gse108539")`) turn the three hypotheses into three predicted
  differences per gene. Simulate on the real arrays (noise from the within-subject residuals, 11 subjects) which
  contrasts are separable when genes are pooled into programs, and how many genes a program needs. A design that
  separates only "behaviour" from "clock or central" is still worth registering if it says so.
- **T3 scout (repeat biopsies).** Which studies have a second biopsy without an intervention in between (exercise,
  clamp, meal), how many persons, and the detection limit for a program rotating with the clock state against one
  that stays flat or follows the intervention.

Gate G1: an item is registered only if its scout shows it can separate the times at the program level. Write the
G1 outcome per item in the scout README and report it.

### Step 4. The bound in hours (item A1), beside steps 1 to 3.
- (a) One table in hours from A208, A231, A233 and A235, from their output tables, with detection limits and the
  companion columns. No new statistic.
- (b) Simulation in the manner of A258 (`exploratory/A258_spread_bias/scripts/{simulate_a258.py,evaluate_a258.py}`):
  is the offset difference between two groups biased when their host spreads differ, and does label permutation
  hold its size? A258 covered d\*, log z2 and log V, not offsets. State which of the four equivalences stand.
- (c) New, registered: the between-tissue offset spread within a person, in hours, with its detection limit, in
  GTEx and in the living cohorts with two tissues. Read first what A133 (`exploratory/A133_axial_decomposition/`),
  A142 and the psi estimates already give; the open question is how much of the residual is true person-specific
  offset and how much is measurement noise. Fix in the amendment the bound that makes it a headline (gate G3).

### Step 5. Register and score what passed G1, then T4.
One amendment per design, written and committed before its data are scored, each with a power note from its scout,
the calibration result it relies on, the two companions, and the reading of every outcome fixed in advance. T2
extends `compute_internal_external.py`; per-gene regressions on the host should use or extend
`clocksync.rhythms.gene_rhythms` (it takes `covariates` and `n_null`), and a tool of general use (a reader of the
tissue deviation, a three-way partition with companions) goes into `clocksync` with a test. T4 (the systemic
signal) is registered after the inventory, starting from the glucocorticoid score of A27 D6 to D8: its main risk is
that agonal stress drives the readers in GTEx, so it needs a living confirmation (MuTHER fat and skin, the
repeat-biopsy cohorts). Items A2 (fetal adrenal covariates) and A3 (COFE head-to-head on TCGA pairs) are
independent of all this and can run beside it when a lane is free. A243 is registered and unscored: score it as
registered when convenient (prompt: `paper/docs/prompts/fetal_pooled_organs_and_adrenal_classes_2026-11-11.md`).

Then gate G2 on everything scored. T5 and T6 follow only if G2 passes; phase B only if it fails.

## After each scoring
- `paper/docs/results_preliminary.md`: one entry per amendment (headline with the verdicts, power first, each
  registered readout, descriptives, reading as registered, limits, code line with runtimes, `<!-- src: -->`).
- `paper/docs/VERDICT_TABLE.md`: one row per registered readout. `paper/docs/METHOD.md` when a tool or a
  calibration is new. `paper/docs/STATUS_ROADMAP_2026-11-14.md`: mark the item done with a pointer to RP; change the
  status tables only with numbers from output tables.
- Commit your own files only, with the attribution lines. Then report.

## Reports and stops
- Report without waiting: after the inventory and the A91 rerun; after the scouts, with the G1 outcome per item and
  the amendments you intend to register; after each scoring; after A1.
- Stop and wait: (1) the outcome of G2 or G3. Which plan the paper takes, and its format and journal, are the
  user's and the PI's decisions: write the evidence and a recommendation as a new question in
  `paper/docs/PI_QUESTIONS.md` (highest number at hand-over: Q32) and pause what depends on it. (2) Anything
  outward-facing: the STARNET application (C2), TwinsUK visit times (C3), requests to authors (C4, C5), the
  in-house jet-lag data (C1, optional, only with the owner's and the PI's agreement). Drafts only. (3) A result
  that changes a standing claim of roadmap section 1.2.
- A report states the verdicts with their numbers, the power behind every negative, what is not shown, and what
  you do next.

## Rules that bite here
- **Pre-registration.** Amendment in `paper/docs/registered_experiments.md`, committed before scoring; next free
  number at hand-over is A259 (check the end of the file). Anything chosen after seeing a result is labelled post
  hoc. A follow-up cut is a new amendment.
- **Numbers** in docs come from output tables, never from memory or from a subagent's message.
- **Statistics.** A failed gate is uninformative, never "decoupled". Compare coupling across groups only at a common
  spread (d\*) or on z2 / psi. Every contrast has a specificity companion on the same permutations. "Internal
  desynchrony" means a phase offset between compartments, nothing else.
- **GTEx hour of death and protected phenotypes** (`data/Gtex_phenotypes/`) are controlled access: aggregates only
  in outputs, logs, docs and messages; never a per-donor row or phase.
- **Blood** has no coherent core clock (diurnal composition): an expected exception. **Brain death** stays
  descriptive: replication only, no mechanism amendment. **Cancer** stays a chapter of this one paper.
- **Bulk designs** need at least two compartments per person; otherwise only signed coherence can be read.
  Single-cell data are read at pseudobulk for the clock.
- **Data.** Processed matrices only, no FASTQ mapping; 434 GB of the 500 GB cap are used, so report before adding
  more than 20 GB. Each download in its own new directory with a `PROVENANCE.md` (URL, date, sha256); downloads are
  untrusted: parse with `python -I`, never execute.
- **Machine.** At most 16 workers and about 100 GB summed over concurrent runs, 20 GB left free; at most 2
  subagents at a time, Sonnet 5.5, for clearly defined tasks. Runs over 20 min only through
  `RUN_CAP_GB=<gb> paper/common/run.sh <script.py> [args]`, detached, with `paper/common/watch_jobs.sh` in the
  background; read `paper/common/jobs.sh` before any output. Scripts log at least every 15 min.
- **Code.** Scripts start with `common.cli.setup`; outputs through `common.paths`; `clocksync` public API first.
  Python: `/home/maxine/miniconda3/envs/sccc/bin/python`.
- **Git.** Local commits only, no remote, no push. The `clocksync/` repository has uncommitted changes that are not
  yours (`.gitignore`, `pyproject.toml`, `tests/test_additions.py`, `src/clocksync/utils.py`, `tests/__init__.py`),
  and the main repository has untracked outputs of other sessions (`paper/exploratory/C4_single_cell/A187/output/`,
  `paper/exploratory/A158_bmal1_bias/sim_2026-10-29/`): leave both alone.
- **Phase 8** (naming, scorecard, final figures, hardening, tutorials) is not started without the PI.

## Pitfalls
- Output table names differ between amendments (`tests.csv`, `score_tests.csv`, `output/score/tests.csv`,
  `output/a123c/tests.csv`): list the folder before citing.
- An RP headline can be stronger than its table: A234 passes but is not clock-specific in its own table; A249 reads
  as a finding and A250 withdrew it. Read the companion columns.
- `pd.read_csv` reads a label "null" as NaN.
- `limit_threads` does not reach loky workers: export `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`
  for joblib runs (the gtex_cues README says the same).
- The shell working directory persists after `cd`: use absolute paths.

## State at hand-over (2026-11-15)
Roadmap written and corrected, nothing of phases A, T, B or C started. Registered and unscored: A243 only. No
background job, no subagent. Tracked files are clean after the commit that adds this prompt.
