# phillips2020circadian
**Ref:** Phillips, Hugues, Yeung, Durandau, Nicolas, Naef, Mol Syst Biol 17:e10135 (2021; cited 2020 preprint). DOI 10.15252/msb.202010135. Type: data + method. Lab: Felix Naef (EPFL).
**Source:** full text

## Question and design
How noisy are core-clock mRNA counts in single cells and how much is due to time? smFISH of Bmal1/Cry1 and Nr1d1/Cry1 pairs in confluent, non-dividing NIH3T3 fibroblasts after dexamethasone. **Sampling:** timed in vitro; "measured transcript numbers every 4 h for 24 h (7 time points)" starting 17 h post Dex; about 450 single cells per time point from three replicates; an independent 4-time-point repeat.

## Method (if any)
Bayesian mixture of negative binomials (telegraph/bursting) with a periodic burst frequency (two-harmonic), cell-size scaling of burst size (M2), partial synchronisation by cell-specific phase offsets (M3), extra cell-specific burst variation (M4, preferred by PSIS-LOO). Variance decomposition into clock time, bursting, Poisson, cell area, other extrinsic. No phase inference algorithm.

## Key results
- Means oscillate: "The mean number of Bmal1, Cry1 and Nr1d1 transcripts per cell oscillated along the circadian cycle and ranged from 10 to 35 molecules" (Results, Fig. 1B). Peak/trough fold change 1.5 (Cry1), 1.6 (Bmal1), 2.1 (Nr1d1).
- Large cell-to-cell spread: "The measured coefficients of variation (CV = standard deviation/mean) ranged from 0.5 to 1." (Results) versus 0.17-0.32 expected for Poisson.
- Cell size matters: "The relationship was sublinear, with exponents ranging from 0.42 for Nr1d1 to 0.70 for" Bmal1 (Fig. 2A); cell area explains 4.2% (Nr1d1) to 17.6% (Bmal1) of variance.
- "For all three genes, we found that the variance was dominated by intrinsic noise" (Results, Fig. 4A), bursting above Poisson.
- Time is minor: "Somewhat surprisingly, the contribution of the 24-h cycle to mRNA variance was low across all genes" (Results); with pooling, "the percentage of variance explained by time increases to almost 100% with 1,000 cells" (Fig. 4B).
- Incomplete synchronisation inflates observed correlations; allowing it raises the burst-frequency amplitude by about 5% (M3).

## Not addressed / limits
- One cell line, one synchronisation protocol, two or three genes per cell; no tissue, no cell-cycle resolved analysis beyond size.
- The authors conclude single-cell phase "is encoded post-transcriptionally"; that is a hypothesis, not tested.
- No proof that unsynchronised cells in tissue behave the same; smFISH, not scRNA-seq with dropouts.

## Link to us
- **Intro claim(s):** our entry (10-35 mRNA per cell oscillate; count noise and size heterogeneity make single-cell phase imprecise) is supported; "10-35" is the range of the mean across time. Time share of variance is low per cell and rises with pooling, supporting pseudobulks as the unit.
- **Relation to clocksync:** supports (justifies donor x compartment pseudobulks); comparator: nikhat2026transcriptional gives decoder-level numbers.

## Ideas
- Quote the 1,000-cell figure beside our A203 (cells per pseudobulk, registered, not scored) result.
- Model our within-compartment cell-size/UMI effect as a covariate: not registered; run only if A203 shows cell-level heterogeneity beyond counting noise.
