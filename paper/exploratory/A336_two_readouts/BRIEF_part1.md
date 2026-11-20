# A336 part 1: brief for wrapper agents

You implement and run, for a list of inventory statistics, the registered A336 part-1 readouts. Repo:
`/home/maxine/Documents/paychere/scCircadianOverview` (read `CLAUDE.md` first). The registration is the last entry of
`paper/docs/registered_experiments.md` (A336, "Part 1"); it is committed and binding. Inventory:
`paper/exploratory/A336_two_readouts/output/inventory.csv` (script, output, criterion, stored value, swap point, notes per
statistic id). Python: `$(conda info --base)/envs/sccc/bin/python`.

## Use the harness
`paper/exploratory/A336_two_readouts/scripts/a336_common.py`: `core()`, `trio()`, `restrict()` (companion rule),
`drops()` (364, or `n=100` when the inventory says gene_drop 100), `drop_percentiles()`, `reading()`, `write_rows()`.
One wrapper per entry: `paper/exploratory/A336_two_readouts/scripts/p1_<entry>.py` (e.g. `p1_a252.py`). Import the
original script as a module (as `A335_review_batch/scripts/a335_1a_nightshift.py` does) and swap the waveform; do not
edit the original scripts. If a swap needs a small change of the original, copy the needed function into the wrapper.
Template for the core/companion pattern: `paper/exploratory/A335_review_batch/scripts/a335_1a_nightshift.py`.

## For each statistic (one row each, via `write_rows(entry, rows)`)
1. **14 genes**: recompute and compare with the stored value (tolerance 1e-6; `reproduced` True/False). If it does not
   reproduce, STOP for that statistic: write the row with `reproduced=False`, no core value, and say why in `notes`.
   Do not compute or print any core value before the 14-gene reproduction check has passed in the same run.
2. **Core**: same data, seeds, statistic; waveform `core()`; companions drawn exactly as originally (same seeds,
   14-gene waveform) then `restrict()`ed to the core positions. Apply the statistic's registered criterion (from the
   inventory / registration of that entry) to the core values: `core_status` = "holds" | "falls" | "at the line" (at the
   line = passes some but not all of its conditions). Write the criterion used in `core_criterion`.
3. **Trio contribution**: `diff_14_minus_core` with a 90 % interval from the statistic's own bootstrap draws applied to
   both waveforms on the same resamples; if the statistic has no bootstrap, 200 unit bootstraps (the statistic's unit:
   persons / donors / pairs) with `np.random.default_rng([336, 1])`. `diff_draws` says which.
4. **Gene-drop null**: the clock statistic only (no companions, no nulls) on every drop from `drops()` (or `drops(n=100)`),
   undefined drops as NaN; `drop_percentiles(values, trio_label, waveforms, direction)`. direction = +1 if the
   original effect is positive in the registered sense (larger = stronger), -1 if negative. The trio drop must equal the
   core value (assert).
5. **Trio alone**: `trio()`; gate: for donor x compartment data `cs.Dataset(...).gate(n_perm=1000, seed=336)` with
   p_pairing and p_acrophase < 0.05; for paired designs (A252, A256, A257) the statistic's own clock-plane test on the
   trio (p < 0.05). Report `trio_alone_value` only if the gate passes (else NaN); `trio_gate_passes` True/False/"not
   computable".
6. `reading = reading(core_status, trio_drop_below_q10)`; `blind` = "core blind" or "core not blind (A335 / A329)" as
   per the registration; `script` = wrapper path.

## Rules
- One scoring run per wrapper. Develop with a `--check` flag that runs on tiny subsets and prints NOTHING about core,
  drop or trio values (assertions only). Then run once for real. No rerun to change a result; a crash before writing
  may be fixed and rerun (say so in `notes`).
- Runs > 20 min: `RUN_CAP_GB=30 paper/common/run.sh <script> [args]`, detached (`setsid nohup ... > /dev/null 2>&1 &`),
  checked with `paper/common/jobs.sh`. At most 4 workers per agent. Call `clocksync.config.limit_threads(...)`.
- GTEx: aggregates only; never write or print per-donor rows, times of death or phenotypes.
- Copy numbers from files, never from memory. Do not edit any doc in `paper/docs/`; do not commit.
- Return: per statistic one line (14 stored / recomputed, core, status, diff and interval, drop percentile all /
  similar (n), trio gate, reading), the wrapper paths, runtimes, and every deviation or problem.
