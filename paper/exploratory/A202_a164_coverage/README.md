# A202: coverage and equivalence size of the A164 interval at its own design

Registration: `paper/docs/registered_experiments.md`, section A202 (commit 408222d). Scripts: `scripts/compute_a202.py`
(simulation, `--arm P` then `--arm B`), `scripts/score_a202.py` (V202.1-2 and descriptive tables). Outputs (git-ignored):
`output/` (`output/subsample/` for the pipeline check). Controlled access: matched sets and labels stay in memory; the
outputs hold per-replicate module estimates per cohort and per-tissue design aggregates only.

## Choices not in the registration
- **A164 code reused** (`compute_a164.py`: matching, `labels`, `dstar`); the cohorts, anchors and matched sets are
  rebuilt exactly as A164 (same function, same minimum of 10 matched cases).
- **One damped copy per tissue.** Some module tissues are anchors of the others (subcutaneous and visceral adipose,
  skeletal muscle). Each tissue is therefore thinned in its own copy of the table and read on its own host, so no host
  contains a damped compartment (the premise stated in the registration). The module statistic is the mean of the
  per-tissue d*, as in A164.
- **Host and c0 fitted once** per cohort x design x tissue, on the donors of that tissue's matched sets: relabelling
  does not change the donor set, so the host (leave-one-out for anchors) and the pooled coupling c0 are label-free.
- **Damping** `clocksync.power.damped` with the host passed in; the arm's loading z is estimated on the arm's donors. If
  the arm's regression is rank-deficient the tissue is left out of that replicate (`n_thin_failed`); the module truth
  is the mean of the realised changes over the tissues kept.
- **Labels.** One uniform key per donor and replicate (shared across tissues); in every set the member with the largest
  key is the case, as the A164 permutations. Seeds `SeedSequence([202, design, cell, replicate, seed, cohort])`.
- **Arm B bootstrap.** As A164 (`boot_chunk`): the resampling unit is the set's original case (joint over tissues), each
  draw brings its set in every tissue with the permuted labels; donors copied, Dataset and host refitted; SE = SD of
  the module mean over 200 replicates. Pooled SE with the fixed A164 weights: sqrt(sum w^2 se^2).
- **Pooling** with the fixed A164 weights (`a164_power.csv`, module, pooled); a replicate enters the pooled tables only
  when both cohorts give a finite module estimate.
- **Coverage** of an interval = |D - truth| <= 1.645 x half-width SD; TOST "equivalent" = D -/+ 1.645 x SE_A164 inside
  +/-0.15 (strict), as A164.
- **Arm B replicate table** (fix 2026-11-04, after the first arm B run stopped in the obesity cells with "need at least 5
  donors"): the copied donors are those of every tissue's sets, as A164 `boot_chunk` (the first version used only the
  target tissue's sets). The contrast still uses only the target's labelled donors, so the estimate is unchanged; the
  chunks already written (t2d cells, 32 obesity delta-0 chunks) were computed with the first version and are kept: for
  them the host is fitted on the same labelled donors, so only the validation of other compartments differed.
