# A111 / A112 / A129: feeding and genetic controls (WP3)

Registered in `paper/docs/registered_experiments.md` (section "A111, A112, A129 (WP3 ...)" and the WP3 note, written
before scoring). `scripts/compute_feeding_controls.py` scores:

- **A111** timed truth and informativeness (liver phase shift >= 3 h relative to the other tissues), H111.1 (offsets,
  exact tissue permutation), H111.2 (liver coupling to the non-liver host, DF vs AL) on Manella 2021 (Alb-Cre controls;
  host WAT + lung) and Xin 2021 (NRF vs DRF; host heart + kidney + VAT).
- **H111.3 = A129 primary**: Manella BLKO vs Alb-Cre, ad libitum, liver coupling to the WAT + lung host at the common
  spread d*, one-sided, labels permuted within ZT; companions WAT on lung and lung on WAT on the same permutations.
- **A112** (Greenwell 2019): lock ratio arrhythmic / ad libitum of the mouse orthologues of the A81b cue and locked sets.
- **Descriptive clock-loss readouts**: BLKO liver offset to the host (untimed, timed), WAT-lung coherence, VIP-Cre
  Bmal1 KO (liver-heart), whole-body Bmal1 KO pair counts.

```
run.sh --subsample 1     quick end-to-end check (50 permutations, 20 bootstrap draws), outputs in output/subsample/
run.sh                   full scoring (1,000 permutations within ZT, 200 bootstrap draws, 10,000 gene permutations)
run.sh --only a112       subset of {a111, a112, clockloss}
```
`scripts/_common.py` holds what clocksync does not cover (`Dataset.contrast` needs >= 2 anchors): the contrast engine for
a single tissue as host, with the permutation sequence of the public contrast (asserted identical for the liver:
column `check_null_identical`), the exact tissue-permutation statistic and the expression QC.

Choices beyond the registered note (all flagged in the outputs): see the report of the scoring run. Outputs
(`output/`, git-ignored): `a111_*.csv`, `a112_*.csv`; log `paper/results/logs/compute_feeding_controls.log`.
