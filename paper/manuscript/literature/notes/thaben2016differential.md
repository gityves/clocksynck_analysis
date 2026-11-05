# thaben2016differential
**Ref:** Thaben, Westermark, Bioinformatics 32:2800-2808 (2016). DOI 10.1093/bioinformatics/btw309. Type: method. Lab: Westermark (Charite Berlin).
**Source:** full text

## Question and design
Statistics for comparing rhythms between two groups (genotype, light condition). Simulations of 500,000 series pairs with 6, 12, 16, 24 samples over two periods; application to mouse liver under light-dark (LD, Menet/Hughes 2012) versus constant darkness (DD, Hughes 2009), 14,253 non-overlapping transcripts, validated on Zhang 2014 (DD, 24 samples every 2 h) and Menet 2012 (LD, 12 samples every 4 h). **Sampling:** both groups timed on a common clock, evenly spaced over at least one full period, with replicates.

## Method (if any)
DODR (R package). Harmonic regression of both groups, then a joint F-type test (HANOVA) of equal versus different phase/amplitude, a scale test for signal-to-noise differences, and a robust rank-based version (Rfit, Fligner-Killeen). Outputs p-values for change in amplitude, phase and noise. The preceding common approach, a p-value-threshold (PVT) rule, is benchmarked as a baseline. Not limma-based (limma appears in compareRhythms, pelikan2021venn).

## Key results
- PVT fails: "a single P value cutoff generally works very poorly" (Results, Fig. 1); it works "only if interpreted with amplitude cutoffs around" 0.1-0.2, and detection of arrhythmia at small amplitude cutoffs is not feasible (85% of liver amplitudes exceed 0).
- Power depends on n: for a 2-fold amplitude change (0.2 to 0.4) "the recall had sunk to only around 0.5 if the number of samples was 12" (Results, Fig. 3C); about 80% at 16 and 95% at 24 samples (FDR 0.1, prevalence 5%). A 4-fold change is detected with recall near 1.
- Phase shifts above 8 h are found nearly perfectly at any n; amplitude losses are harder than gains (noise 0.1).
- Outliers cut HANOVA power (AUC drops to about 0.95 for phase); robust DODR is marginally worse without outliers.
- Liver LD vs DD: "Using DODR with a FDR of 0.1, we found that 262 of these transcripts exhibit differential rhythmicity" (Results, Fig. 4): 183 amplitude gains in LD, only 10 pure phase shifts; reproducible in the validation sets; IRF7 motifs enriched.
- "Phase differences between LD and DD conditions are only meaningful when there are well-defined rhythms in both conditions" (Results): phase is tested only for rhythmic transcripts.

## Not addressed / limits
- Period changes: "The present methods does not address possible period changes." (Discussion). Not waveform changes either.
- Requires timed samples with the same time grid and at least 12 samples; authors "would not recommend fewer samples" (Discussion). No count-based or cohort data, no donor random effects.
- Absolute amplitude differences across assays are spurious, so validation across technologies is qualitative.

## Link to us
- **Intro claim(s):** our entry says limma-based, which is wrong: it is harmonic ANOVA and rank-based tests (DODR). Otherwise it supports the view that differential rhythmicity is a statistics problem with specific power limits. Our rule that a failed gate is not decoupling matches their warning that low-amplitude series cannot be called arrhythmic.
- **Relation to clocksync:** background; competitor in the timed setting only (needs times, we do not).

## Ideas
- Power note for the paper: their recall curves (n of 12-24, amplitude ratio 2) can calibrate how many donors a coupling contrast needs; our own power tables already exist per contrast (see results_preliminary.md), no new run needed.
- Run DODR on a timed cohort (e.g. a mouse LD/DD tissue pair) and compare its amplitude change with our |z| ratio: related to V52 (paired damping, METHOD.md), not the same design.
