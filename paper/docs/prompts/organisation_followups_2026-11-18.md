# Prompt for a lead agent: follow-ups of the clock-organisation framework (written 2026-11-18)

You continue the paper "The human body keeps internal time" (Yves Paychere, EPFL, Naef lab; method `clocksync`). The
last days produced a framework for how the clock is organised in the body; your job is to apply it where the user
asked: the cancer chapter (now main text), disease, and the fetal section. Work autonomously; report as listed under
"Reports and stops". A negative that removes a claim is a result: do not tune a design to keep one.

## Read first
1. `CLAUDE.md` and the memory index (`MEMORY.md`): all entries apply.
2. `paper/docs/STATUS_2026-11-18.md` (claims by tier, what fell, open decisions). Then `STATUS_2026-11-17.md` section 2.
3. RP (`paper/docs/results_preliminary.md`): A291, A294b (and its correction note), A295, A297, A293, A302 for the
   framework; A138, A212, A227, A271, A127, A166, A171, A173, A215, A219, A224 for cancer; A153, A208, A231, A235, A258
   for disease contrasts; A243, A244, A264, A266, A269 for the fetal section; and `docs/literature/fetal_clocks_animal_2026-11-18.md`.
4. `paper/docs/PI_QUESTIONS.md`: Q24 (cancer chapter), Q32 and Q35 (fetal, tally), Q36 (calibrated pairwise test).
Run `git log` and `git status` first: another session works on the same branch (presentation, A290). Commit your own
files only; never touch its uncommitted files or jobs.

## The framework you apply (in one paragraph)
A tissue's clock state in donor d is: organ term x body's time(d) + departure(d, tissue). The organ term has a **gain**
|z| and a fixed **phase**; the departure has variance 2 psi, random direction between donors. Coupling = S / (S + 2 psi)
with S = z2 V. In GTEx, S spans 85-fold and 2 psi is nearly constant and clock-specific (A291, A294b): differences of
coupling are differences of gain over a constant departure. Group comparisons are allowed on z2 (gain) and psi
(departure), not on coupling (A258: coupling contrasts are biased at unequal host spread). Within organs, structural
cells and macrophages carry the time at the organ's phase; lymphoid cells carry none; composition changes gain, not
phase.

## Task A. Cancer in the gain / departure view (first; the user wants the chapter in the main text)
**A1. Tumour gain and departure by tumour state (TCGA paired tumour-normal, `common.data.load("tcga_paired")`).** Within
each cancer type (never pooled across types without the correction of A3), is the proliferation-associated loss of
tumour clock coherence (A127-A173) a loss of **gain** (z2 of the tumour on the normal-tissue host falls: malignant cells
with weaker or no clock dilute the rhythm, as lymphocytes do in normal tissue) or a rise of **departure** (psi rises:
malignant cells drift from the body's time)? Same for ER status, TP53 and the luminal A contrast where A215-A224 found
effects. Fix the contrasts, the proliferation score (as in A127) and the readings before scoring; companions on the same
permutations (matched; responsive reported, A285: the responsive q95 rule is anti-conservative with few pairs).
**A2. The host problem (the user's caveat).** With one normal sample per patient, the host is that sample's state,
including its own departure: tumour estimates are attenuated by the normal sample's reliability, the same within a
cancer type (contrasts valid) but different between types. Simulate at the real design (`clocksync.simulate.like`)
to show what a pure gain loss and a pure departure rise look like on this host before reading A1.
**A3. Between-cancer comparisons.** Disattenuate tumour estimates with the matching GTEx tissue's reliability
(S / (S + 2 psi) from the A104 table) and say which cancers keep the body's time most and least; descriptive unless
registered with power.
**A4. Single cell (only if powered).** Do malignant cells share the tumour's own stromal-vascular cells' clock state
(pairwise readout, `clocksync.pairwise`, A283 calibrated form)? A149 found non-clock sharing dominates within a tumour
sample; the scout of 2026-11-17 (RP "Scouts of the cell-type program", C2) left power uncomputed. Compute power first;
stop if it fails.

## Task B. Disease in the gain / departure view (GTEx, many tissues per donor: trustworthy host)
Do chronic conditions change the gain or the departure of the organ they affect? Candidates from the protected
phenotypes (aggregates only, never per-donor rows): type 2 diabetes (liver, adipose, muscle, pancreas if present),
hypertension (arteries, heart), smoking (lung), long-term steroids (all tissues; A136 context). Power first (MDE of
z2 and psi differences at the arm sizes); specificity companion on the same permutations; offsets were already shown
equal within about 1 h (A208, A235): do not repeat them. Register only powered contrasts.

## Task C. Fetal section, reframed
The human result (signal-induced covariation in adrenal cortex, lung, liver; heart and kidney negative; no core loop)
matches the animal literature's "maternally driven before autonomous" (literature file above). Write the reframed
fetal paragraph for the status file (no new statistic needed). Optional registered check if cheap: is the fetal
covariation's gene pattern the glucocorticoid-induction pattern (PER1, PER2, NR1D1, TSC22D3, FKBP5 loading together)
rather than the waveform's phase relations? **Do not map FASTQ for spliced / unspliced counts**: CLAUDE.md allows
processed matrices only; this test waits for the user's decision (rule change or an author request, draft only).

## Task D. Housekeeping
- A290 (other session; interrupted): check whether it was restarted; if not and its owner is idle, report to the user
  rather than running it yourself (it is the other session's amendment).
- Write a new PI question after the last one: the before / after tally of 2026-11-15 to 2026-11-18 (STATUS_2026-11-18
  sections 1-3), the claim tiers, and the open decisions (venue, spliced / unspliced, Q37).

## Rules that bite (learned in the last days)
- **Pre-register.** Hypothesis, statistic, readings and seeds committed before data are scored; staged registrations
  (design addendum with arm sizes before scoring) when the design depends on data counts. State what has been seen.
  Smoke tests print design counts only; if a smoke or debug run prints a statistic, say so in the entry.
- **Companions:** the matched family is calibrated in bulk (A285); the responsive q95 rule is anti-conservative with
  few pairs: report it, decide on the matched family or on a calibrated percentile (A283 form).
- **Host leakage:** regressions on a leave-one-out host leak host error into extra regressors in proportion to the
  tissue's loading (A292 post hoc); residual correlations between tissues share host error (correct as A302).
- **Bootstraps must respect the design** (A303: one library per time point per arm makes within-time resampling
  degenerate; check the interval width is not zero).
- **Plumbing pitfalls:** worker processes need file-loaded modules shipped by value (`cloudpickle.register_pickle_by_value`,
  A274); `Dataset(min_donors=...)` drops compartments silently below the default; donor ids must not contain "#"
  (`timed_loadings` splits resections at "#"); sheep tables need `mixed_strata=True`.
- **Long runs** (> 20 min) only through `RUN_CAP_GB=<gb> paper/common/run.sh ...` detached, with
  `paper/common/watch_jobs.sh` running. Machine budget: at most 21 threads, 83 GB summed caps, 2 subagents; the user
  allowed 4 cores for these lines of work.
- GTEx time of death and phenotypes: aggregates only. Numbers in docs from output tables only.

## Registration and writing
Next free amendment number A304 (check the end of `registered_experiments.md`). After each scoring: RP entry (headline
with verdicts, design, power, each readout, reading as registered, limits, code line, `<!-- src: -->`), rows in
`VERDICT_TABLE.md`, and the status file when a claim changes tier. Commit with the attribution lines.

## Reports and stops
- Report after the `git` check, after each power stage, after each scoring.
- Stop and wait: anything outward-facing (author requests, access, journal); a rule change (FASTQ mapping); when tasks A
  and B are scored: write the tally into the PI question and recommend how the cancer chapter reads in the framework.
- Not in scope: Phase 8 (figures, naming, tutorials), brain-death mechanism, new atlases beyond those on disk, pushing.
