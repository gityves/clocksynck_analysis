# nikhat2026transcriptional
**Ref:** Nikhat, Mandal, Veerasubramanian, Chakrabarti, iScience 29:115394 (2026). DOI 10.1016/j.isci.2026.115394. Type: method + data. Lab: Shaon Chakrabarti (NCBS Bangalore).
**Source:** full text

## Question and design
Can the circadian phase of a single cell be read from its clock-gene RNA counts? Multiplexed SABER-FISH of Bmal1, Tef, Nr1d1, Nr1d2 (later also Cry1, Nifl3) in Dex-synchronised NIH3T3 (and lung fibroblasts), about 200-250 cells per time point. A mouse skin-fibroblast scRNA-seq set (12,177 cells after QC, 8,096 genes, ZT2-22 every 4 h) tests non-core genes. **Sampling:** timed, in vitro, a time course after Dex covering about one full cycle; training times 8-32 h post Dex, three centre times as truth.

## Method (if any)
PRECISE: supervised Gaussian-process phase decoder trained on the time course (true phase from wavelets), likelihood over phase per sample; compared with unsupervised SVD/eigengene phase angle and an autoencoder. Needs a timed training series; infers from groups of cells averaged by bootstrap.

## Key results
- Counts are over-dispersed (negative binomial) but time still carries signal: "suggesting that the average levels of gene expression across time are likely to be well separated despite the high intrinsic noise" (Results, Fig. 1E).
- Single cells fail: "of the deviations using single-cell data was close to the diagonal line" of random assignment (Results, Fig. 4B).
- Averaging rescues: "was able to predict the circadian phase with a median" error of about 1 h with 20 cells (0.5 h with 130) from three genes; SVD needs about 70 cells (Fig. 4A, D).
- Gene count needed at single-cell level: 16 clock-like genes still mix; "clusters with accuracy greater than 90% emerged even at single-cell resolution" only at about 50 simulated genes at Nr1d2-level noise (COV 0.81) (Fig. 4G).
- In scRNA-seq, extra bulk-rhythmic genes hurt: "the inclusion of low-noise non-core-clock genes led to a decrease in cluster purity" (Fig. 5D); core genes alone gave about 40% purity; about 240 low-noise genes were added.
- Averaging over about 70 spatially close cells also exposes local phase coherence in desynchronised fibroblasts.

## Not addressed / limits
- In vitro, one synchronised cell line and one mouse skin set; no tissue, no human donors, no disease.
- Supervised decoding needs timed training data; unsupervised SVD only separates sample times. Phase truth for test times limited to three centre points.
- "Like any supervised algorithm, PRECISE makes assumptions regarding the underlying model and noise structure" (Results).
- Noise limit shown for single cells; does not give donor-level pseudobulk size rules in tissue.

## Link to us
- **Intro claim(s):** the pseudobulking-justification claim is supported, with numbers: about 20-70 cells per group for 3 core genes in a culture; single cells are near random.
- **Relation to clocksync:** supports (why donor x compartment pseudobulks); competitor only for the single-cell question.

## Ideas
- Cells per pseudobulk versus precision of coupling/offset: registered as A203 (HLCA, not yet scored); compare its saturation point with their 20-70 cells.
- Cite as the counterpart of phillips2020circadian: same conclusion with 3-6 genes and a supervised decoder; check that tauFisher and Tempo (both cited by the paper as scRNA-seq phase methods that need pseudobulks) are in our bib.
