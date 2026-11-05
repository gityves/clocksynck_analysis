# sancar2021clocks
**Ref:** Sancar A, Van Gelder RN, Science 371:eabb0738 (2021). DOI 10.1126/science.abb0738. Type: review. Lab: Sancar (UNC) / Van Gelder (UW).
**Source:** full text

## Question and design
Review: does clock disruption cause cancer and does time-of-day chemotherapy work? Covers jet lag and metabolic syndrome as disease of desynchrony, shift-work epidemiology, clock gene knockouts, oncogenes (Myc, p53, Ras) acting on the clock, DNA excision-repair rhythms, chronochemotherapy trials.
**Sampling:** not applicable. Mouse repair data (kidneys, cisplatin, LD12:12) and liver transcriptome in constant darkness (44 h) shown in Fig. 1.

## Key results
- "Desynchrony between the internal circadian and external geophysical clock causes at least three pathological conditions: jet lag, sleep phase syndromes, and metabolic syndrome" (Circadian clock disruption and disease); shift work is a significant risk for metabolic syndrome (citations 22-24), with food timing relative to clock phase as proposed mechanism.
- Epidemiology is split: Large studies (155,000 women in Denmark; the Generations study, over 100,000 women in the UK) "revealed no increased cancer risk in the shift-work population".
- Clock genes are not general tumour suppressors: knockouts raise tumorigenesis for some gene-tumour pairs, lower it for others (Cry-less mice are resistant to p53 mutation tumours).
- Excision repair: about 10% of genes in a tissue are circadian; non-transcribed strand repair peaks in early evening, transcribed-strand repair follows each gene's transcription phase, so there is no single best time for DNA damage response.
- Chronochemotherapy: "the empiric application of chronochemotherapy to cancer treatment has not produced sufficiently consistent results to justify routine clinical application"; the small positive ovarian study was not reproduced in a larger one; the EORTC colorectal trial found no overall benefit, a trend toward benefit in men and a 38% mortality increase in women.
- Tumour rhythms in the literature range from arrhythmic, to rhythmic but out of phase with normal tissue, to rhythmic and in phase; "A major challenge toward this goal is the development of noninvasive methods to assess the circadian parameters of cancers."

## Not addressed / limits
- No systematic or quantitative synthesis; no effect sizes with CI for trials; clinical statements rely on the colorectal and ovarian trials only. Does not discuss immunotherapy timing.
- No analysis of patient circadian marker as modifier of timing benefit.

## Link to us
- **Intro claim(s):** clocks, cancer and chronochemotherapy (refs_*.md). Supported, and it is the skeptical counterpoint: our wording should not state that disruption causes cancer or that timing improves outcome in general.
- **Relation to clocksync:** background; names the open question our method can address in humans: relation of tumour and normal tissue clock phase and amplitude without timed samples. Tumour-normal clock match: H138.4 (TCGA), single-cell stromal analogue (results_preliminary.md l.1074).

## Ideas
- The three tumour scenarios (arrhythmic, out of phase, in phase) map to c, |z| and arg z in clocksync; done in part as tumour-vs-own-normal tests (H138.4, A175, A185-A186), no explicit three-class summary.
- Gene-specific repair phase suggests clock readout from clock genes alone cannot predict damage response; cite Yang 2018 PNAS (cisplatin repair map) and Giacchetti 2006.
