# A164: GTEx T2D and obesity, matched on death circumstances, equivalence and power

Registration: `paper/docs/registered_experiments.md`, section A164. Script: `scripts/compute_a164.py`. Outputs (git-ignored):
`output/` (`output/subsample/` for the check). Controlled access: only aggregates are printed or written; matched sets
stay in memory. The script writes matching, balance and power before it computes any observed estimate.

## Choices not in the registration
- **Matching per tissue.** Sets are formed among donors with data in that tissue (availability differs: liver 74 of 380
  donors in `vent`), so matched counts are per cohort x tissue x design. Cases are processed in rounds (round r gives
  every case its r-th control, cases with fewest candidates first, ties by donor id): no case is favoured by order.
- **Nearest** is lexicographic: absolute age difference (exact age, years), then absolute ischaemic-time difference.
  Ischaemic time is the donor-level value `ischaemic_h` of `common.gtex.load` (first sample of the donor).
  Donors with missing exact age or ischaemic time are not matched.
- Obesity pool: T2D = 1 donors excluded (unknown T2D kept); normal weight 18.5 <= BMI < 25; 1:1, same rules.
- T2D + BMI sensitivity: additional BMI caliper 3 (both BMI known), T2D 1:3. No-steroid sensitivities exclude
  `steroid` = 1 from the pool before matching (both designs).
- **Eligibility** (>= 10 matched cases) is applied per cohort, design and tissue, to companion tissues as well.
- **Balance** SMDs: matched cases vs matched distinct controls (unweighted); before = all pool cases vs controls of the tissue.
- **d\*** = c\*_case - c\*_control (`Dataset.contrast` levels control, case; `d_star`), one `contrast` call per tissue
  (`targets=[t]`, `blocks="none"`, `n_perm=0, n_boot=0`) on the labelled donors only; the host is the pooled anchors of
  those donors. Module = mean of eligible module tissues; module minus companion = module mean minus mean of eligible
  companion tissues.
- **Permutations (2,000).** `contrast` cannot give draws shared across tissues when sets differ per tissue, so the
  script draws them itself: every donor gets one uniform random key per permutation and in each set the member with the
  largest key becomes the case (a uniform label swap within the set, identical donors share their key across tissues).
  One draw feeds the module and the module-minus-companion statistic.
- **Bootstrap (500), matched sets.** `contrast`'s bootstrap is within groups with the host fixed, so it is not used. Per
  replicate the case donors are resampled with replacement (jointly over tissues), each draw brings its set in every
  tissue, donors are copied (`id@draw`), a new `Dataset` is built and the host is refitted in each `contrast` call.
  SE = SD of the replicate statistic.
- **Pooling:** inverse variance with the bootstrap SE of each cohort; the pooled null is the same weights applied to
  the cohort nulls (cohorts have disjoint donors). Weights are fixed at the observed bootstrap SEs. A cohort with no
  eligible module tissue is dropped.
- **90 % interval** of the pooled module d\* = D +/- 1.645 x pooled bootstrap SE (normal). Two-sided add-one
  permutation p on |statistic|. Holm over the two module p values (T2D, obesity); "organ-specific" uses the raw
  module-minus-companion p < 0.05.
- **Power:** MDE = 2.49 x SD of the pooled permutation null; expected half-width = 1.645 x pooled bootstrap SE
  (`a164_power.csv`, per cohort and pooled, for module and module-minus-companion). Powered = half-width <= 0.15.
  Equivalence verdict only when powered, otherwise "not powered (no verdict)".
- **Sensitivities** are point estimates only (no null, no bootstrap); pooled with the weights of the primary design of
  the same family. The A47 side-by-side takes d\* = -loss_star from `figure_3/output/gtex_cues/a47_t2d_*.csv`; the BMI
  entry is the A47 slope per SD (different units, listed for reference).
- Subsample check: `--subsample N` keeps N random donors per cohort, >= 3 matched cases, 20 permutations, 10 replicates.
- Cost: every permutation refits the host through `contrast` (~0.2-0.9 s per tissue), about 3-5 s per permutation.

## Files
`a164_matching.csv`, `a164_balance.csv`, `a164_power.csv` (written first); `a164_tissues.csv` (per cohort x tissue d\*,
bootstrap SE, null SD, p); `a164_tests.csv` (pooled and per cohort module statistics); `a164_tests_hypotheses.csv`
(H164.1-3); `a164_sensitivity.csv`; `a164_a47_side_by_side.csv`; `draws/` (module-free per-tissue d\* of the
permutations and replicates, no donor information).
