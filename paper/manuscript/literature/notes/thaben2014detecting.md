# thaben2014detecting
**Ref:** Thaben, Westermark, J Biol Rhythms 29:391-400 (2014). DOI 10.1177/0748730414553029. Type: method. Lab: Westermark (Charite Berlin).
**Source:** full text

## Question and design
JTK_CYCLE only detects symmetric waveforms; RAIN extends rank-based detection to arbitrary waveforms. Synthetic series (noise sd 0.1 relative, from mouse liver microarray residuals, CV constant about 0.10; 100,000 series per setting; sine and sawtooth) and real mouse liver microarray (Hughes 2012: 24 samples, every 2 h, CT0-46) and proteomics (Mauvoisin, Robles 2014: 16 samples every 3 h, pooled mice). **Sampling:** timed, evenly or unevenly spaced, at least one full period of the pre-specified period.

## Method (if any)
RAIN splits each period into a rising and a falling part and tests each with the Mack-Wolfe umbrella test on ranks, so slopes need not mirror each other; scans pre-specified periods, phases and umbrella peak positions, with adaptive Benjamini-Hochberg across them. Exact null distribution. "Missing values and non-even sampling intervals are accommodated for". Run time 1-2 min for 25,000 series with 24 time points.

## Key results
- Sawtooth detection at amplitude-to-noise 2, 16 samples every 3 h, 25% prevalence, FDR 0.1: "RAIN recovered almost 60% of the sawtooth wave forms compared with less than 10% for JTK_CYCLE" (Results, Fig. 2B); equal power for sine waves.
- Mouse liver (BH 0.01): "RAIN detected 2014 rhythmically transcribed genes, compared with 810 genes detected by JTK_CYCLE" (Results, Fig. 3); 1204 RAIN-only.
- Validation: "It turned out that 65% of these 946 transcripts were flagged by the harmonic regression method as being rhythmic in the control data" (Results) versus 40% base rate in the 1 h control series.
- Down-sampled to every 4 h: 1590 (RAIN, cutoff 0.1) vs 701 (JTK, 0.25) genes.
- Proteome: "We applied RAIN to these data and detected 316 and 384 proteins with circadian rhythms in their abundance" (cutoff 0.25), against 195 and 186 originally; JTK found 35 and 25.

## Not addressed / limits
- Tests prespecified periods (24 h used); no phase or amplitude comparison between groups.
- Needs the sampling times; noise model is liver microarray (CV 0.10), not RNA-seq counts or single cells.
- Many RAIN-only hits are asymmetric rhythms in a strongly rhythmic tissue; power for weak rhythms in noisy cell lines is not studied.
- FDR is control-adaptive, high cutoffs (0.25) in the proteome.

## Link to us
- **Intro claim(s):** rank-based nonparametric detection of rhythms of arbitrary waveform in timed data is supported. Also supports the point that a symmetric-cosine assumption is a limit of widely used detectors.
- **Relation to clocksync:** background. Our clock-gene harmonic model is also first-harmonic; waveform asymmetry is untested there.

## Ideas
- On timed GTEx/mouse-liver data, check whether the clock-gene first-harmonic assumption loses asymmetric genes: run RAIN versus cosinor on the clock-gene panel (no registered amendment found).
