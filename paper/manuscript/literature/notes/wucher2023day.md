# wucher2023day
**Ref:** Wucher, Sodaei, Amador, Irimia, Guigó, PLoS Biol (2023). DOI 10.1371/journal.pbio.3001986. Type: data.
Lab: Roderic Guigó / Manuel Irimia (CRG Barcelona). Not Naef lab.
**Source:** full text

## Question and design
Day-night and seasonal variation of gene expression in 46 GTEx tissues, using the recorded time of death (TOD),
not inferred phase. **Sampling:** post-mortem, TOD over the whole day: "between 08:00 and 17:00 (351 individuals)
and between 21:00 and 05:00 (315 individuals)"; 266 "twilight" donors excluded.

## Method
Differential expression day vs night (covariates: sex, age, BMI, post-mortem interval, season); Lomb-Scargle /
MetaCycle on TOD for rhythmicity; seasons by season of death.

## Key results
- TOD coverage (published aggregate, usable by us): 351 day (9 h window), 315 night (8 h), 266 twilight (7 h): about 37-39 donors per hour in each window, i.e. roughly uniform over 24 h.
- Rhythmicity vs TOD weak: only 187 genes "(1%) that were circadian in at least one tissue", attributed to GTEx limits ("there is no information about the location of death").
- Day-night DE: "12,530 (70%) were differentially expressed between day and night in at least one of the 46 tested tissues"; per tissue "5.5% of expressed protein-coding genes were day-night on average" (lenient: non-adjusted P ≤ 0.05).
- Lung most (2,418 genes), testis, colon transverse, salivary gland fewest; brain small.
- Tissue specificity: day-night genes in "2.6 tissues out of 39" on average.
- Core clock: ARNTL and NPAS2 day-night in most tissues, "NPAS2 shows daynight gene expression differences in a much larger number of tissues than CLOCK (31 versus 6 tissues)".

## Not addressed / limits
- TOD ≠ internal phase (cf. Talamanca 2023); lenient thresholds; no phase or synchrony between tissues.

## Link to us
- **Intro claims:** GTEx deaths cover the whole day (≈ uniform; with the exact counts above), so GTEx is not a hospital-hours cohort: it is the healthy reference where hospital hours can be simulated (A123). Supports correcting `paper/docs/literature_positioning.md` §4. Our refs table: correct that it is Guigó/Irimia, not Naef.
- **Relation:** data/background.

## Ideas
- Cite the 351/315/266 counts directly for the GTEx sampling-time distribution instead of our controlled-access TOD (no per-donor release needed).
- NPAS2 vs CLOCK breadth: check that NPAS2 carries clock signal in our gene set across GTEx compartments.
