# moller-levet2022diurnal
**Ref:** Möller-Levet, Laing, Archer, Dijk, BMC Biology 2022 (20:63). DOI 10.1186/s12915-022-01258-7. Type: data (re-analysis of public human blood and baboon time series). Lab: Dijk / Archer (Surrey).
**Source:** full text

## Question and design
Overlap and timing of rhythmic transcripts of human whole blood with the 64-tissue baboon atlas (Mure 2018, GSE98965). Human blood (microarray, GSE48113 / GSE39445 / constant-routine set): in-phase diurnal, 19 participants; sleep out of phase, same 19; constant routine, 19. **Sampling:** all timed, full cycle. "Seven 4-hourly blood transcriptome samples (spanning > 24 h) were collected in each condition" (diurnal); constant routine 10 samples 3-hourly, ~40 h. Baboon: 12 samples 2 h apart, no replicates, diurnal 12:12 light-dark; human light period 18 h 40 min.

## Method (if any)
Set overlap of rhythmic genes, circular cross-correlation of acrophase profiles (lag in h) with light onset as reference, circular clustering, GO enrichment. Proof of principle: PLSR from blood probes (excluding clock genes) to the phase of one marker gene in baboon PVN or stomach antrum; trained on 9 participants, validated on 10 (in-phase) and 10 (out-of-phase). Assumes cohort-level acrophases of orthologous genes; needs timed blood samples.

## Key results
- 759 human-blood rhythmic orthologues, 86% rhythmic in at least one baboon tissue. "Of the 759 orthologous genes that were rhythmic in human blood, 652 (86%) were also rhythmic in at least one baboon tissue"
- Only 152/652 (23.3%) stay rhythmic without a sleep-wake cycle. "only 152 genes (23.3%) were also rhythmic in the absence of a sleep-wake cycle"
- Timing correlated significantly in 25/64 tissues, blood on average earlier. "an average earlier timing of 3.21 h (SD 2.47 h) in human blood"
- Core clock genes rhythmic in blood are rhythmic in only 24-46 tissues each; ARNTL is a clear outlier. "Timing of Arntl profiles was on average 10 h phase advanced in human blood compared to baboon tissues/organs" (lag range -10 to +12 h; only 1 of 36 ARNTL probes rhythmic).
- Blood-to-tissue phase predictor: R2 0.63 and 0.67, MAE 2.94 and 2.81 h (ANT, PVN), one blood sample. "predicted the ANT and PVN phase with an R2 of 0.63 and 0.67, and a mean absolute error of 2.94 and 2.81 h"
- Predictor collapses under misalignment. "None of the 45 probes (targeting 41 genes) that were part of the predictor set of the ANT that were rhythmic when sleeping in-phase remained rhythmic when sleeping out of phase."

## Not addressed / limits
- Human and baboon differ in species, photoperiod, lab versus sacrifice; "The baboon data set contained no biological replicates" (12 points per tissue).
- "The human data were obtained in healthy young participants." No shift workers, no disease, no tissue-resolved human data; blood is the only human tissue.
- Overlap measured in rhythmic-gene lists (about 6% per tissue), not in clock-gene coupling; no test of whether blood rhythm tracks other organs within a person.
- Predictor needs timed samples and is not validated outside the training lab protocol.

## Link to us
- **Intro claim(s):** refs_circadian_desynchrony says blood rhythmic genes overlap baboon organ rhythms: supported (86%), but with caveats: only 23% are circadian, and the clock genes themselves are least aligned (ARNTL 10 h).
- **Relation to clocksync:** background / competitor for blood-as-proxy; data (baboon, done: A36/A40; GSE48113, done: A156). Supports that organ rhythms desynchronise when sleep is mistimed (predictor failure).

## Ideas
- Done: baboon organ offsets (A36/A40) and GSE48113 misalignment (A156); not repeated.
- Compare their blood-to-tissue MAE (about 2.9 h, absolute phase, timed samples) with clocksync offset error (0.45 h in baboon, A36) as two different tasks (absolute versus relative to host).
- Their note on organ peak times reconstructed from unlabelled samples (ref 33, bimodality absent) bears on GTEx: check the acrophase distribution of our untimed GTEx organ offsets for bimodality.
