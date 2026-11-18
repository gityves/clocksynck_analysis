# Prompt for a new agent: a minimal physical model of the clock's organisation in the human body (written 2026-11-18)

You join the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab; method `clocksync`). The user
wants to know whether one small physical model can unify the organisation results (organs, cell classes, disease) and
propose a new conceptual view of the clock, or whether it would only re-describe them. Your job is to answer that with
a model, a fit and one blind prediction. **"The model is not worth it" is an acceptable outcome**: if it needs about one
parameter per observation, say so and stop. Do not tune the model to rescue it.

## Read first
1. `CLAUDE.md` and the memory index (`MEMORY.md`): every rule applies (pre-registration, numbers from output tables,
   GTEx aggregates only, compute budget, long jobs through `paper/common/run.sh`, explicit `git add`).
2. `paper/docs/STATUS_2026-11-18.md` (claims by tier, what fell).
3. RP (`paper/docs/results_preliminary.md`), entries:
   - organs: A291 and "A291 post hoc", A294b and its correction note, A302, A292, A297;
   - cell classes: A293 (Parts 2 and 3), A273, A279, A295, A281, A284, A287, A294 (lymphocytes);
   - disease: A304, A212, A254, A183, A308, A164, A233, A235, A306;
   - loss of drive: A293 Part 1, A132 (including H132.2, the ventilation-duration tertiles).
4. `paper/docs/PI_QUESTIONS.md` Q7 (no data other than GTEx on brain death; the SCN reading is a hypothesis).
Run `git log` and `git status` first: other sessions work on the same branch. The next free amendment number is the
one after the last `## A` heading in `registered_experiments.md` at the time you register (check again just before
committing; numbers have collided before, see A294 / A294b).

## The candidate model (from the discussion with the user; yours to formalise, change or reject)
Each organ is a population of self-sustained cellular clocks driven by the body's systemic time. The bulk clock state of
cell class c in organ o of person p:

    x_ocp(t) = a_c * R_o * exp(i (theta_p(t) + phi_o + delta_op)) + measurement noise

- a_c: single-cell clock amplitude, a property of the cell class (lymphocytes a ~ 0 in tissue);
- R_o: synchrony of the organ's cells (order parameter), set by drive and coupling against period dispersion;
- phi_o: the organ's locked phase lag;
- delta_op: the person's departure, an additive perturbation of fixed size in the clock plane of the organ's mean
  field, so its angle scales as 1 / (a R) and coupling = S / (S + D) rises with gain.

Observations it must reproduce (numbers: read them from the tables below, not from this list):
1. Shared variance S spans about 85-fold over 27 GTEx peripheral tissues while own variance 2 psi spans about 5-fold and
   barely scales with S (Theil-Sen slopes 0.07 for 2 psi, 0.22 for the clock-specific excess; A291 post hoc, A294b);
   departures are not shared between meal-driven tissues (A302).
2. Measurement does not produce the amplitude-coupling relation (A291).
3. Relative gain of a cell class to endothelium is conserved across organs (A293 Part 3, tau <= 0.31), while the same
   class differs between GTEx organs by a factor of about 1.3 to 2.7.
4. Cell classes share the organ phase (A273, A279 calibrations; no conserved class order, A293 H293.3; pooled
   endothelium lead +0.46 h, a lead only); composition moves gain, not phase (A295).
5. Peripheral clocks are self-sustained: removing the feeding rhythm keeps 0.72 of liver amplitude, against 0.15 for
   additive cues (A297).
6. Disease: proliferating tumours lose gain without a departure rise (A304); psoriatic plaques keep a phase delay of
   about 1.4 h (A183, A308; local clock, cue genes on time); seven chronic conditions change nothing (A164-A306); brain
   death loses the organ arrangement (A293 Part 1, descriptive).
Known tension: in A132 H132.2 (descriptive), organ donors with the longest ventilation keep more order than those with
the shortest, against a free-running drift that grows with time since the drive was lost (ventilation time is not time
since brain death). Report it; do not explain it away.

## Inputs (aggregates; GTEx per-donor rows are controlled access and never leave the compute scripts)
- Organ table: `paper/results/figure_2/output/backbone/a104_gtex_tissues.csv` (S = z2 V, own = 2 psi per tissue);
  `paper/exploratory/A291_attenuation/output/tissues.csv`; `paper/exploratory/A294_floor/output/{tissues,sets}.csv`.
- Cell classes: `paper/exploratory/A293_organisation_levels/output/part3/{families,atlases}.csv`,
  `.../part2/{tests,atlases}.csv`; `paper/exploratory/A295_composition/output/tissues.csv`; A273 / A279 outputs.
- Feeding: `paper/exploratory/A297_feeding_model/output/{arms,tests}.csv`.
- Disease: `paper/exploratory/A304_cancer_gain_departure/output/score/`; `A308_psoriasis_calibration/output/`.
- Brain death: `paper/exploratory/A132_brain_death/output/`, A293 Part 1 `groups.csv`.
- CoPhaser (the lab's single-cell phase inference) is importable in the `sccc` env (`import cophaser`, source
  `/home/maxine/Documents/paychere/CoPhaser`). Read its README and API before using it; prefer its public API.

## Stages (each ends with a short report to the user; stop at every gate)
**Stage 0, literature (half a day; Sonnet subagent allowed, verify every citation).** What is established for (i) bulk
amplitude = single-cell amplitude x synchrony (Welsh 2004, Nagoshi 2004), (ii) amplitude-dependent entrainment and
locking phase (Granada 2009, Abraham 2010, Bordyugov 2015), (iii) Kuramoto / population models of peripheral or tissue
clocks, (iv) renewal / turnover models of epidermal or intestinal clocks. These are from memory: confirm or correct.
Write `paper/docs/literature/clock_population_models_2026-11-18.md` (claim, source, what is new for us). Stop if a
published model already makes our factorisation in human tissue.

**Stage 1, formalisation and a registered adequacy criterion (before any fit).** Write the model (equations, which
observation constrains which parameter), count parameters against independent observations, and state which
alternatives it must beat: (a) linear damped oscillators driven by systemic cues (A297 argues against them in the liver only),
(b) composition alone for organ gain, (c) per-organ free gains without structure. Register as an amendment, before
fitting: the fit statistics, the adequacy criterion (e.g. parameters < observations / 3 and residuals within the
observed bootstrap intervals), and what each outcome means. All data are seen: label the fit "not blind,
descriptive".

**Stage 2, simulation and fit (a few days of work, light compute).** A population simulation (Stuart-Landau or phase
oscillators with amplitude, N cells per organ, period dispersion, coupling, systemic drive, person-level additive
perturbations, Poisson read-out at real depths) that produces the same estimators as the data (`clocksync` on simulated
pseudobulks, not idealised quantities). Fit to observations 1 to 5; check 6 as parameter changes (gain loss = lower a or
R; phase shift = changed phi_o or turnover; loss of drive = free run). Report the adequacy verdict and which
parameters are identifiable. Stop if the criterion fails.

**Stage 3, the blind prediction (register before scoring; power stage first).** The model's sharpest new prediction:
amplitude and synchrony separate in single cells. Within a cell class, single-cell clock amplitude is the same across
organs, and the organ factor of gain is the cell-to-cell phase spread. The same readout decides whether tumour gain loss
(A304) is damping or dephasing, which bulk cannot (A212). First a feasibility and power stage on the body atlases of
A272 / A274 with CoPhaser (per-cell clock amplitude and phase dispersion per class x organ; malignant clocks were
unreadable in pseudobulk, A195): if per-cell clock counts cannot resolve the prediction, write a power note and stop.
Then register the test with companions (expression-matched and stress-responsive gene sets, A272 `pseudo_sets`) and
score it once.

**Optional, only if Stage 2 passes:** a brain-death drift prediction (organ disorder against time since the drive was
lost). Check first, in aggregate counts only, whether GTEx records anything closer to time since brain-death declaration
than ventilation duration (`A132_brain_death/output/scout_variables.csv`). A132 H132.2 has been seen: any test is not
blind.

## Rules specific to this task
- No word "resonance" in the model's name or text unless the fit shows a resonance (gain peaking at a detuning); the
  current data show "each class responds to the organ's drive with its own gain at the organ's phase".
- The model is a synthesis, not evidence: never state a parameter value as a biological finding unless the parameter
  is identifiable and its interval is reported.
- Budget: at most 8 workers and 40 GB for this task (other sessions run; check `paper/common/jobs.sh`).
- Commit scripts and docs only; results to RP with `<!-- src: ... -->` lines; numbers from output files.

## Reports and stops
After Stage 0, after Stage 1 (registration committed), after Stage 2 (adequacy verdict), after the Stage 3 power stage.
Each report: what was done, the verdict in one line, the numbers with their source file, and what you propose next.
