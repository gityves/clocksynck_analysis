# hughes2010jtk
**Ref:** Hughes, Hogenesch, Kornacker, J Biol Rhythms 25:372-380 (2010). DOI 10.1177/0748730410379711. Type: method. Lab: Hogenesch (Penn).
**Source:** full text

## Question and design
Nonparametric detection of cycling variables in genome-scale time series. Benchmarks on simulated transcripts (1024; half null) with 48 points at 1 h over 2 days, then 4 h sampling and 24 replicate-pair designs; application to mouse liver, pituitary, NIH3T3 and U2OS microarrays. **Sampling:** fully timed, evenly spaced, at least one (here two) full 24 h cycles; groups defined by sampling time.

## Method (if any)
JTK_CYCLE applies the Jonckheere-Terpstra/Kendall tau test to a family of cosine-shaped orderings over user-set periods and phases; the minimum exact p (Harding look-up table depending only on the replicates per time point) is Bonferroni-adjusted. Outputs: period, phase, amplitude (median sign-adjusted deviation). Rank-based, so outlier-resistant. Assumes known sampling times, balanced groups, and a symmetric (cosine-like) waveform.

## Key results
- Power at low amplitude: "at amplitudes greater than ~1.5, JTK_CYCLE can unambiguously identify all cycling transcripts" (Results, Fig. 1); COSOPT and Fisher G fail below amplitude 2 (noise sd 1, outliers 20 at 1% of points).
- Errors at 48 points (Table 1, p<0.05): JTK 5 false positives and 3 false negatives, COSOPT 4 and 135, Fisher G 10 and 131. "Consequently, the number of false negatives is considerably greater for COSOPT" (Results).
- Accuracy: period R2 0.926 vs 0.732 (COSOPT) and 0.053 (Fisher G); phase R2 0.766; amplitude R2 0.912.
- Real data, BH q<0.05 (Table 2): liver 5425, pituitary 392, NIH3T3 25, U2OS 21 cycling transcripts; the three major periods (about 8, 12, 24 h) are recovered in liver.
- Speed: 48 time points and about 45 k transcripts in 15-20 min versus days for COSOPT.
- Also applied to a human cell-division series (Whitfield 2002), so it is not circadian-specific.

## Not addressed / limits
- Needs the sampling time of each sample and several time points; no untimed or single-snapshot use.
- No missing-data handling or non-even sampling discussion (RAIN adds it); no comparison between conditions: it yields one p per series, so differential rhythmicity by threshold lists is exposed to the problem in pelikan2021venn.
- Waveform symmetry is assumed (see thaben2014detecting for the sawtooth loss).
- Amplitude is a rank statistic, not on a comparable relative scale.

## Link to us
- **Intro claim(s):** nonparametric rhythm detection in timed time series is supported. Marks the boundary of the field: it needs timed samples, which human tissue banks lack.
- **Relation to clocksync:** background (what a timed-data detector does and needs).

## Ideas
- Positive control only: on a timed GTEx or mouse design, JTK/RAIN rhythmicity of clock genes should agree with our compartment amplitude |z| (partly done: V52 same-time variance vs timed amplitude, METHOD.md).
