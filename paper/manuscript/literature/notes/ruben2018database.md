# ruben2018database
**Ref:** Ruben, Wu, Smith, ..., Anafi, Hogenesch, Sci Transl Med 10:eaat8806 (2018). DOI 10.1126/scitranslmed.aat8806.
Type: data atlas. Lab: Hogenesch (Cincinnati) / Anafi.
**Source:** full text (main; SI not in export)

## Question and design
Human population-level circadian atlas from GTEx without time-of-day: CYCLOPS ordering of 13 of 51 tissues
(≥160 samples, conserved clock network, reproduced clock-gene phase relationships), 632 donors.
**Sampling:** post-mortem (TOD not used).

## Method
CYCLOPS per tissue (independent orderings); cosinor on 160 ordered samples per tissue; drug-target annotation.

## Key results
- "nearly half of protein-coding genes were shown to be cycling in at least 1 of the 13 tissues"; "7486 distinct genes cycle in at least 1 of the 13 tissues".
- 54 ubiquitous cyclers (≥8 tissues); activators "peak, on average, 12 hours out of phase with their E-box targets NR1D1, NR1D2, and PER3".
- Tissue specificity: largest two-tissue overlap of top-500 cyclers "was composed of just 28 genes (5%)"; "266 of the top 500 liver cyclers (53%) were exclusive to liver".
- ~1000 cycling genes are drug targets/transporters/metabolisers; 136 cardiovascular drug targets.
- Persistence despite heterogeneity: "clock gene oscillations persist despite wide interindividual variability".

## Not addressed / limits
- Each tissue ordered separately: no within-donor comparison of tissue phases, so no statement on synchrony.
- Only 13 tissues pass quality: the rest fail ordering (brain weak), so coverage is selective.
- Statistics on inferred order optimistic (CYCLOPS caveat).

## Link to us
- **Intro claims:** human tissues have widespread, tissue-specific rhythmic output; many drug targets cycle (chronotherapy motivation). Our refs entry (CYCLOPS on 13 GTEx tissues, 632 donors, unsupervised, cohort-level) is correct.
- **Relation:** orderer application; shows the cohort-per-tissue design that cannot measure intra-donor synchrony.

## Ideas
- None specific beyond citing it as the CYCLOPS-on-GTEx precedent (and the 13-of-51 tissue pass rate as a limitation of orderers).
