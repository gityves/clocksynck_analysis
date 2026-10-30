# A168: aromatase-inhibitor ER blockade and breast-tumour clock coherence

Registration: `paper/docs/registered_experiments.md`, A168 (frozen). Outcome goes to `paper/docs/results_preliminary.md`.

| file | role |
|---|---|
| `scripts/prep_a168.py` | parses the three GEO series matrices and platform annotations into `output/{poetic,letrozole,z1031b}_{meta,expr}.parquet`, `sample_counts.csv`, `prep_log.csv` (no coherence computed) |
| `scripts/a168_lib.py` | loaders, paired designs (`Pairs`: difference-in-differences and time-swap), pseudo-clock sets for a variable gene list, residualised coherence |
| `scripts/run_a168.py` | `--stage power` (null SDs and MDE, no observed statistic) and `--stage tests` (tests, Holm, pseudo-clock companion, sensitivities, descriptive tables) |

Outputs (git-ignored, `output/`): `power.csv`, `null_draws.csv`; after `--stage tests`: `tests.csv`, `pseudo_T_draws.csv`,
`pseudo_sets_{poetic,letrozole}.csv`, `sensitivity.csv`, `coherence_groups.csv`, `treatment_effect_scores.csv`. `output/quick/` is a
smoke test (50 permutations, 10 pseudo sets), not for reading.

Run: `python -I prep_a168.py`, then `RUN_CAP_GB=8 paper/common/run.sh .../run_a168.py --stage power --jobs 2`, then `--stage tests`.

## Choices not in the registration

1. **Probe mapping:** the GEO platform annotation (GPL10558, GPL96, GPL6480), column "Gene symbol". A probe annotated to several symbols
   ("A /// B") counts towards each (needed for the registered 13 genes of GSE20181: NR1D1 sits on a THRA /// NR1D1 probe); probes without
   symbol dropped. Symbols passed through `clocksync.genes.canonical` for the 14 clock genes (aliases such as C1orf51 to CIART).
2. **Log scale:** only GSE20181 (MAS5, maximum 108,791) is linear and is log2-transformed (non-positive values, minimum 0, raised to the
   smallest positive value) before probes are averaged; GSE105777 (maximum 16.9) and GSE87411 (log2 ratios, maximum 12.4) are used as deposited.
3. **Missing values (GSE87411 only):** genes missing in more than 20 % of samples dropped; remaining NaN cells (13 in clock genes) set to the
   gene mean (A166 convention). Counts in `prep_log.csv`.
4. **Scores:** proliferation, purity proxy and oestrogen-response score (Hallmark ESTROGEN_RESPONSE_EARLY and LATE union, from the same
   MSigDB file `h.all.v2024.1.Hs.symbols.gmt` as `sh.hallmark_g2m`, clock genes removed) are the A97 `_zmean` over the whole series, z-scored
   over all samples of the series (including day-90 and month samples).
5. **Pseudo-clock pool:** `a166_lib.pool_genes` minus the oestrogen early and late genes; mean and SD of the pool and of the clock genes taken
   over the samples of the tests (POETIC all 426; GSE20181 pre and day 10-14 only). `a166_lib.pseudo_sets` is hard-wired to the 14 genes, so
   `a168_lib.pseudo_sets` is the same construction for a gene list (reproduces it exactly for the 14, checked). A pseudo set carries the
   amplitudes and phases of the clock genes present (13 for GSE20181).
6. **Pseudo-clock companion statistic:** the pseudo T is the registered statistic of the same design with the observed labels. "Pass" is
   Holm p < 0.05 over H168.1-3; the specificity label follows A166 (clock-specific if T < pseudo 5th percentile).
7. **Seeds and nulls:** SEED = crc32("A168") mod 2^31; permutation nulls in chunks of 250 with `SeedSequence([SEED, test, chunk])`
   (independent of `--jobs`); pseudo sets `[SEED, 10 + series]`; sensitivities `[SEED, 100 + j]`. The tests stage reloads `null_draws.csv`
   from the power stage. Permutation p: `clocksync.stats.perm_p(..., "less", tie=0.0)` as A166.
8. **H168.1 null:** arm labels permuted over the 213 patients (157 / 56), both samples of a patient moving together (as registered).
   H168.2 and H168.3: time labels swapped with probability 1/2 per patient.
9. **Sensitivities** each with its own 2,000-permutation null (same null type as the matching test), p reported descriptively:
   "without PER1" and "without PER1 and NR1D1" (two variants; the registration wording is read this way) for all three tests;
   residualisation on [proliferation, purity proxy] (OLS with intercept, per gene, within each group of samples entering a coherence, so
   within arm x time point under the observed labels and recomputed within the permuted groups in the null); POETIC treated patients in
   thirds of the change in proliferation score (rank thirds, H168.2 statistic within each; no untreated difference-in-differences by third);
   GSE20181 responders (37 patients) and non-responders (15) from the title label ("not assessable" patients out), pre vs day 90 (A vs C
   pairs); GSE87411 BL vs 2-4 weeks.
10. **GSE87411 time labels:** baseline = BL, BL Redo; 2-4 weeks = w2, W2, w2 Redo (a patient has one array per class); "re-analyzed" arrays
    form the second variant as extra baselines (no patient has both). The "M" (month) arrays are not used. Redo arrays are treated as the
    patient's array of that class (no patient has both a Redo and a plain array).
11. **Descriptive tables:** coherence of every arm x time point over all its samples (`coherence_groups.csv`, includes day 90 and months);
    paired change in the three scores per arm with a paired t p-value (`treatment_effect_scores.csv`).
12. **MDE:** 2.49 x SD of the permutation null as registered; METABRIC and SCAN-B effects are the registration's round numbers (-0.165, -0.315),
    stored in `power.csv` for comparison.
