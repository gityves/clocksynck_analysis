# greco2021integration
**Ref:** Greco, Koronowski, Smith et al., Sci Adv 2021 (7:eabi7828). DOI 10.1126/sciadv.abi7828. Type: biology (hepatocyte-only clock, feeding, muscle clock). Lab: Sassone-Corsi / Aznar Benitah.
**Source:** full text

## Question and design
Isolate the interaction of the liver clock and feeding, and find which rhythmic liver genes need other clocks. Mice, 8-10 weeks, LD 12:12: wild type (WT), Bmal1 KO, Bmal1 Liver-RE (BMAL1 reconstituted only in hepatocytes in an otherwise clock-less mouse, Alfp-cre), each ad libitum (AL) or night feeding (NF, food in the 12 h dark phase, 2 weeks). Readouts: bulk liver RNA-seq at 6 time points (ZT0-ZT20, n = 3), metabolomics and lipidomics (n = 4), ATAC-seq and ChIP-seq (BMAL1, CEBPB); then muscle-specific Bmal1 KO (mKO) liver RNA-seq over the daily cycle; cross-tissue correlation with QENIE (Hybrid Mouse Diversity Panel). **Sampling:** timed, 24 h.

## Method (if any)
JTK_CYCLE, BioCycle, LimoRhyde differential rhythmicity; rhythmic genes classed autonomous (liver clock only), integrated (liver clock plus feeding), feeding-driven, network-dependent (not restored with liver clock plus feeding).

## Key results
- Perturbed: all clocks except hepatocyte (Liver-RE) versus none (KO). NF corrects the phase of core clock genes in Liver-RE: "it corrected the phase of core clock genes and modestly enhanced their amplitude in Liver-RE mice" (Fig. 1F). Hepatocyte clock plus feeding gives WT-like liver transcriptomes in PCA.
- Gene classes (of oscillating genes in WT): "13.7% of oscillating genes were autonomous"; of 1236 nonautonomous genes 666 (53.9%) gained oscillation under NF in Liver-RE, 442 (66.4%) integrated and 224 (33.6%) feeding-driven (also oscillate in KO, lower amplitude); "a large portion of nonautonomous genes (441, 34.4%) failed to oscillate in Liver-RE mice under NF" (network-dependent; oxidation-reduction, fatty acid beta-oxidation).
- Metabolites and temporal coherence: hepatocyte clock plus feeding does not restore much of the lipid and redox metabolome; "temporal coherence between metabolic pathways remained disrupted"; 420 (57.6%) metabolite correlations present only in WT livers.
- Other cells measured: skeletal muscle clock. In muscle-specific Bmal1 KO (bulk liver RNA-seq), liver core clock protein and transcript levels unaltered; 60% of oscillating genes in WT maintained; 371 transcripts (16% of oscillating) show differential rhythmicity: 94 (25%) oscillate in both with different amplitude or phase, 84 (23%) lose, 193 (52%) gain oscillation. "the muscle clock is required for their proper phase alignment"; 30% of network-dependent genes tied to muscle BMAL1. Cross-tissue correlation: muscle has most correlations with network-dependent genes; liver metabolites not restored correlate mostly with muscle, then serum. Muscle-derived serum factors (unknown) regulate many of the same genes in hepatocytes.
- Core claim: "Clock function in peripheral tissues is required, but not sufficient, for this task."

## Not addressed / limits
- Cell-type resolution inside the liver is not addressed (hepatocytes only, bulk liver); the liver non-parenchymal clock is not manipulated. Between-tissue interaction (muscle to liver), not cell-cell within tissue. No sc or sorted data, no phase of non-hepatocyte cells.
- Feeding rhythm imposed by NF is "more robust than the endogenous rhythm" (authors); LD, 2 weeks NF; liver clock rhythms are measured in animals with a feeding cue, so tissue persistence without any rhythmic input is not tested. Muscle effect on 30% of network-dependent genes only; muscle cross-talk factor unidentified.
- Only 3 mice per time point.

## Link to us
- **Intro claim(s):** the tissue clock alone is not sufficient: about 14% of liver rhythmic genes are hepatocyte-clock autonomous; most need feeding and other clocks; so cross-organ organisation matters (tissue clock plus network).
- **Relation to clocksync:** supports / background. Directly relevant to A315 (brain death: tissue amplitude kept, organ arrangement lost): a hepatocyte-only clock plus an imposed feeding rhythm keeps core clock amplitude and phase, while redox and lipid outputs lose rhythm and coherence without the muscle clock (cross-organ phase alignment). Time scale: weeks of chronic knockout, not hours of ventilation. For the within-tissue question: hepatocyte clock is autonomous but needs systemic input to be phased; no pacing by a dominant cell population inside the liver; the muscle effect is an inter-organ pacing effect on phase alignment.

## Ideas
- Test in GTEx whether cross-tissue relations involving skeletal muscle are among those lost in brain-dead donors (A315 follow-up, register first).
- Cited related: Koronowski 2019 (liver-RE first paper), Welz 2019 (epidermal clock), see sinturel2021circadian card.
