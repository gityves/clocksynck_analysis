# solanas2017aged
**Ref:** Solanas et al., Cell 2017 (170:678). DOI 10.1016/j.cell.2017.07.035. Type: biology (mouse). Lab: Aznar Benitah / Munoz-Canoves.
**Source:** full text

## Question and design
Sorted epidermal (integrin a6-bright CD34-) and muscle (a7-bright CD34+) stem cells from young adult vs aged (18 months or older) mice, microarray; plus Bmal1 KO, Per1/Per2 double KO, high-fat diet, caloric restriction. **Sampling:** timed, six time points every 4 h (ZT0 lights on, ZT4 ... ZT20) in 12:12 light-dark; JTK p <= 0.05 (epidermal) or <= 0.005 (muscle). Activity monitored in light-dark and constant darkness. Per-time-point animal numbers were not recovered from the text extraction.

## Method (if any)
No new method; JTK rhythm calls, amplitude distributions, GO enrichment, then functional follow-up (S-phase timing, pRPA/gH2Ax, autophagy).

## Key results
- "We detected robust rhythmic expression in 3,144 transcripts from adult epSCs and 1,979 transcripts from adult muSCs, of which 76% and 72%, respectively, ceased to be oscillatory in aged mice". Aged cells gained new rhythmic genes (2,309 aged epSCs; 2,221 aged muSCs); 23-28% remain rhythmic.
- "all core clock genes oscillated with equal amplitude and period in both aged and adult SCs"; average amplitude of rhythmic genes is the same, and aged mice stay behaviourally circadian (same amplitude and phase in constant darkness).
- Aged epSCs lose homeostatic (differentiation) rhythms and gain DNA damage, inflammation (IL-6) rhythms; S-phase peak moves from ZT16 to ZT4 (4-fold fewer cells in S-phase at ZT16); M-phase peak stays at ZT20.
- Bmal1 KO and Per1/2 dKO do not reproduce the age reprogramming; caloric restriction prevents it, high-fat diet does not mimic it.
- Companion paper (sato2017circadian) shows the same output rewiring in liver, with tissue-specific phase profiles.

## Not addressed / limits
- Mouse only; stem cells only (rare sorted cells), no neighbouring niche cells or SCN.
- Whether aged stem cells keep the same phase relative to liver or SCN is not measured: tissues are sorted from separate animals.
- Rhythm counts depend on p cut-offs (0.05 and 0.005), so the 72-76% loss partly reflects power; the amplitude result is the more robust one.

## Link to us
- **Intro claim(s):** refs row age reprograms rhythmic gene expression: supported, but the primary message is that the core clock is intact and only the output programme changes, so clock-gene-based coupling and phase should be age-robust; it does not support age abolishes rhythms.
- **Relation to clocksync:** background; mouse data (no single-cell input).

## Ideas
- Test whether the clocksync clock-gene readout is age-stable in human tissue with the same asymmetry (clock genes stable, output genes rewired): the age contrast on z and psi is A19 (age effects on relearned anchors, MDEs reported in results_preliminary.md) and A130 (internal time vs age); not re-proposed.
- Missing reference to add: Nakamura 2011 (SCN ageing, cited in the paper) for the central-clock side of ageing.
