# D1b feasibility: which fetal clock readouts can the Cao 2020 design detect? (2026-11-02, before any registration)

Question (user, 2026-11-02): A181 showed the cross-organ gate cannot pass at the Cao design. Can the same atlas say which
fetal organs or cell types carry a clock, with readouts that do not need fetuses to share organs?

Script `scripts/feasibility.py` (run 2026-11-02, 297 s, 8 workers, `paper/results/logs/feasibility.log`). Outputs in
`output/` (git-ignored): `design.csv`, `passrates.csv`, `subsamples.csv`, `markers.csv`, `markers_expected.csv`.
**No fetal clock statistic (coherence, gate) was computed.** Fetal gene counts were read for option 3 (see below).

Rule: adult material cut to the fetal design (fetal n per organ or class: >= 1,000 cells per organ, >= 200 per class;
binomial thinning to the Cao median depth), 200 subsamples; a readout is feasible if it passes in >= 0.8 of them (A181).

## 1. L1 coherence within one compartment

Coherence R (`decompose.components`, GTEx waveform) against a per-gene donor permutation null (200), pass p < 0.05.

| calibrator | fetal n | pass rate (thinned) |
|---|---|---|
| GTEx Liver (08-18 h deaths) | 7 | 0.805 |
| GTEx Lung | 11 | 0.715 |
| GTEx Adrenal | 7 | 0.705 |
| GTEx Heart LV | 9 | 0.51 |
| GTEx Muscle | 7 | 0.405 |
| GTEx Kidney cortex | 9 | 0.37 |
| GTEx Cerebellum | 5 | 0.25 |
| GTEx Cerebrum (BA9) | 8 | 0.07 |
| heart nuclei: fibroblast (stromal) | 9 | 0.935 |
| heart nuclei: endothelial | 9 | 0.73 |
| heart nuclei: cardiomyocyte | 9 | 0.48 |
| lung HLCA: fibroblast (stromal) | 11 | 0.815 |
| lung HLCA: AT2 (epithelial) | 11 | 0.62 |
| kidney nuclei: proximal tubule | 9 | 0.55 |
| kidney nuclei: endothelial | 6 | 0.245 |
| kidney nuclei: interstitial (stromal) | 9 | 0.085 |

Window check (`--all-day`, `passrates_allday.csv`): the fetal clock phase at collection is unknown. Collection is in
clinic hours, but a fetal clock not entrained to the mother would sit at any phase. GTEx deaths at any hour (wider phase
spread) give nearly the same rates: lung 0.695, liver 0.69, adrenal 0.765, heart 0.525, muscle 0.48, kidney 0.44,
cerebellum 0.38, cerebrum 0.13. The window is not the limit; n is.
Lung endothelial not calibrated (adult pool 15 < 1.5 x 11). Thinning changes little: n is the limit, not depth.
Fetal liver is haematopoietic at this age, so adult bulk liver is a poor calibrator for it.
Feasible: heart stromal (0.935), lung stromal (0.815). The null breaks all gene covariance, so a pass shows clock-shaped
covariance, but not that it is clock-specific: a fetal test would need expression-matched pseudo-clock sets, which need
genome-wide pseudobulks (rebuild from `raw/GSE156793_S3_gene_count.loom.gz`, 13 GB on disk; the 21 GB decompressed loom
would be temporary: 376 + 21 = 397 GB of the 400 GB cap, deleted after the build).

## 2. L2 across cell classes within one organ (`Dataset.gate`, 1,000 permutations)

| calibrator | classes | fetal n (all classes) | pass rate (thinned) |
|---|---|---|---|
| heart nuclei (Reichart 2022 stratum) | cardiomyocyte, fibroblast, endothelial | 9 | 0.94 [0.90, 0.97] |
| kidney nuclei (KPMP stratum) | proximal tubule, interstitial | 9 | 0.76 [0.69, 0.82] |
| lung HLCA | AT2, fibroblast, EC capillary | 11 | not calibrated (adult pool 9) |

Caveat: adult sampling times are unknown and probably spread more widely than hospital hours; GTEx (timed) has no
cell classes. The heart rate is therefore an upper bound for a narrow fetal window.
Feasible: heart (0.94). Kidney is just below the rule.
Adrenal (`--adrenal`, added for A190): no adult adrenal nucleus atlas; heart and kidney classes cut to the adrenal
design (7 fetuses, parenchyma, stromal, endothelial depths) pass 0.90 [0.85, 0.94] and 0.875 [0.82, 0.92].

## 3. Clock maturation markers (descriptive)

`markers.csv`: detection (>= 20 counts) of core and output genes per fetal organ and class. `markers_expected.csv`: counts
expected at the fetal depth if the class expressed the gene at the adult nucleus-atlas level.
- Detection follows depth: every gene is detected in the deep pseudobulks (cerebrum, cerebellum, adrenal, heart, lung,
  kidney parenchyma); output genes are mostly absent in the shallow classes (endothelial, immune, placenta).
- At adult levels most output genes would give < 50 counts in fetal stromal and endothelial classes, so an absence there
  is uninformative. Only the deep parenchymal classes can be read.
- Output-gene levels in a single sample depend on the time of sampling (DBP several-fold over the day), so a level is
  not a maturation readout unless averaged over many fetuses at known or spread times.
- Fetal median counts were seen next to the adult expectation while checking this (`markers_expected.csv`): in heart,
  lung and kidney parenchyma, DBP, NR1D2, PER3 and TEF are at or above the adult expectation. This makes option 3 a seen
  descriptive, not a test.

## Reading

- The atlas can support a small registered test: L2 in fetal heart (cardiomyocyte, fibroblast, endothelial; 9 fetuses,
  89-122 days) and L1 in heart and lung stromal cells. A pass would place a running, coupled clock in the fetal heart by
  13-17 post-conception weeks. A fail stays weak evidence (sensitivity 0.94 at best, adult times unknown).
- It cannot rank organs by onset: most organs and classes are below 0.8, brain cortex lowest (0.07), so "which organ
  first" is not answerable here.
- Age trend within 72-129 days: not estimable (9-11 fetuses per organ).
