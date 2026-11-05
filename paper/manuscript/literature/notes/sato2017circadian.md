# sato2017circadian
**Ref:** Sato et al., Cell 2017 (170:664). DOI 10.1016/j.cell.2017.07.042. Type: biology (mouse). Lab: Sassone-Corsi / Aznar Benitah.
**Source:** full text

## Question and design
Mouse liver, young vs old, normal diet (YND, OND) vs 30% caloric restriction for 6 months (YCR, OCR); stem cells (epidermal, muscle) from the companion paper (solanas2017aged) for tissue comparison. Microarray, JTK_cycle at p < 0.01. **Sampling:** timed, full cycle: sacrifice every 4 h (ZT0-20), 3-6 mice per time point per group (YND 4-5, OND 3-4, YCR 5-6, OCR 3-5), 12:12 light-dark; one mouse per sample (pseudo-series from different animals).

## Method (if any)
No new method; JTK_cycle rhythmic calls, Venn overlap of rhythmic gene sets, amplitude and phase distributions, GO enrichment. Plus acetylation proteomics and NAD+ metabolomics (not needed here).

## Key results
- Young vs old normal diet: 2,626 genes (44.8%) rhythmic only in YND, 1,664 (28.4%) only in OND, 1,575 (26.9%) in both. Ageing therefore rewires rather than abolishes the hepatic oscillating transcriptome: "Aging reprograms the circadian transcriptome in a highly tissue-specific manner" (In Brief).
- "The expression of core clock genes and clock-controlled genes (CCGs) as well as the expression of core clock proteins is virtually unaltered when comparing young and old animals fed a normal diet" (Bmal1, Clock, Rev-erb, Per2, Cry1 by qPCR).
- CR raises amplitude and broadens phase: 1,867 (22.6%) rhythmic genes in YND vs 4,070 (49.2%) in YCR, 2,334 in both; 1,722 vs 3,715 genes in old (OND vs OCR), 1,517 in both; new phase peaks in CR liver are attributed to feeding anticipation.
- Only 74 genes are rhythmic in all three tissues (liver, epidermal and muscle stem cells) in young normal mice: the output programs are tissue-specific. Liver and muscle stem cell phase profiles are alike; the epidermal profile changes with age and reverts under CR.
- Global protein acetylation is rhythmic in young liver and lost in old; CR restores it.

## Not addressed / limits
- Mouse, one sex-and-strain design, inbred animals under controlled light and feeding; no human data. Rhythm calls at p < 0.01 with 3-6 animals per time point, so overlap sizes depend on power.
- No test of between-tissue phase offsets inside the same animals: tissues come from different groups of mice.
- Clock genes are examined only by qPCR in liver (five genes); the other tissues' clock-gene stability is in the companion paper.

## Link to us
- **Intro claim(s):** refs row age effect on tissue rhythms: supported; more precisely it shows rewiring of output genes with an intact core clock, which is the property that makes a clock-gene readout (clocksync) age-robust.
- **Relation to clocksync:** background; mouse data (the 'mouse_ageing' resources, A48/A49 `compute_mouse_ageing.py`) already cover Deota, not this dataset.

## Ideas
- Done: A3 (age effect on GTEx amplitude), A48/A49 (mouse ageing in Deota). Not done: apply the clocksync readout to the Sato and Solanas liver/stem-cell series, to test that a clock-gene panel recovers the same phase in old vs young while the transcriptome-wide rhythmic set changes (needs GEO matrices; check GEO accessions in the STAR Methods before registering).
- Use the 74 tissue-independent genes as a data-driven alternative to the fixed 13-gene panel in a robustness check.
