# lellupitiyagedon2019circadian
**Ref:** Lellupitiyage Don et al., Cell Cycle 18:2447-2453 (2019). DOI 10.1080/15384101.2019.1648957. Type: biology (cell lines). Lab: Farkas (UMass Amherst).
**Source:** full text

## Question and design
Do BMAL1 and PER2 oscillate in breast cancer cells when read by continuous luminescence instead of sparse RT-qPCR?
Two cell lines: MCF7 (luminal A, low grade) and MDA-MB-231 (basal, high grade), each with BMAL1:luc and PER2:luc
lentiviral reporters; N = 4 shown per reporter, 12 recordings per line and reporter in total; RT-qPCR in 3 biological
replicates; western blots.
**Sampling:** in vitro, serum-shock synchronised; RT-qPCR every 4 h for 48 h; LumiCycle recordings continuous after the first 24 h transient is removed. No human tissue, no normal cell line in this paper.

## Method (if any)
FFT-based rhythmicity test after quadratic detrending, continuous wavelet period estimate, MetaCycle (Lomb-Scargle,
ARSER, JTK) for qPCR.

## Key results
- RT-qPCR misses it: both transcripts present but "oscillations were not rhythmic in the circadian range (p > 0.05 using Metacycle)" in both lines.
- MCF7: rhythms with typical antiphase; FFT test rhythmic for all BMAL1:luc and 11 of 12 PER2:luc recordings; periods 22.5-23.2 h (BMAL1:luc) and 22.1-23.6 h (PER2:luc) in the main figure; damped amplitude.
- MDA-MB-231: "Very few statistically significant oscillations were observed" (p < 0.05 for 1 of 12 BMAL1:luc and 3 of 12 PER2:luc).
- Reading: clocks are "further disrupted" in higher-grade cells; authors speculate that ER/PR signalling (absent in MDA-MB-231) and EGFR/Ras may explain the difference.

## Not addressed / limits
- Only two cell lines that differ in subtype, ER/PR, EGFR and everything else; grade or "aggressiveness" is not isolated, and no non-malignant control (MCF10A) is run.
- Cell-autonomous, serum-synchronised oscillations in 2D culture say nothing about in vivo tumour clocks or host-tumour synchrony; arrhythmic reporter does not mean an arrhythmic endogenous clock in the tumour.
- Reporter-expressing lines, passage 11-16; no loss-of-function test of the proposed mechanism.

## Link to us
- **Intro claim(s):** tumour clocks are weakened, in a grade-dependent way; also warns that sparse RT-qPCR can wrongly call arrhythmia. Our refs row is SUPPORTED, with two lines only (it should not be cited for a general grade trend).
- **Relation to clocksync:** background.

## Ideas
- Compare with the single-cell/bulk prediction: do malignant epithelium clock states in ER+ (luminal A) vs triple-negative cohorts differ in coupling at a common spread d*? Partly covered by A134/A138.5 (malignant vs adjacent-normal epithelium, 0.583 vs null 0.443, p 0.069, below MDE); an ER-stratified version is not registered.
- Relevant missing references: Gutierrez-Monreal 2016 (MCF10A vs MCF7 series reanalysed in ye2018genomic).
