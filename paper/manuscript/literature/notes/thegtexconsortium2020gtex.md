# thegtexconsortium2020gtex
**Ref:** The GTEx Consortium, Science 369:1318-1330, 2020. DOI 10.1126/science.aaz1776. Type: data. Lab: GTEx Consortium (Ardlie, Lappalainen, Aguet).
**Source:** full text of the main article (supplement and methods not included: donor-level covariates such as Hardy scale, ischaemic time, ventilator status are not in this text)

## Question and design
GTEx v8 atlas of genetic regulatory effects: eQTL and sQTL mapping in cis and trans, tissue specificity, GWAS integration. Bulk RNA-seq and WGS from post-mortem donors. **Sampling:** post-mortem tissue, collection time relative to clock time or to time of death is not discussed in the article; no circadian analysis anywhere in it.

## Method (if any)
Not a clock paper. QTL mapping per tissue (tissues with at least 70 donors), covariate-adjusted (details in the supplement, "(11)").

## Cohort
- "The GTEx v8 dataset, after quality control (11), consists of 838 donors and 17,382 samples from 52 tissues and two cell lines." (QTL discovery)
- 49 tissues with at least 70 individuals analysed: "a total of 15,201 samples from 838 donors". Figure 1A gives tissue sample numbers (e.g. lung 515, thyroid 574, hypothalamus 170, pituitary 237, 11 brain regions).
- Ancestry and sex: "Of the 838 donors, 715 (85.3%) were European American, 103 (12.3%) African American" and "557 (66.4%) donors were male and 281 (33.5%) female".
- Design aim: "from ~50 tissues from up to ~1000 postmortem donors" and to create "standards and protocols for optimizing postmortem tissue collection and donor recruitment" (Introduction).
- RNA-seq median depth 82.6 million reads; WGS median 32x.

## Key results (relevant to us)
- Cell type composition is a key driver of tissue-level expression: "show that cell type composition is a key factor in understanding gene regulatory mechanisms in human tissues" (abstract). Supports analysing compartments (tissue, cell type) as distinct units.
- Authors warn that sex- and ancestry-associated differences may reflect "cell type composition or environmental exposures" (Sex- and population-biased effects).

## Not addressed / limits
- Time of death, time of collection, Hardy scale, ischaemic time, ventilator status, season, cause of death are not described or used in this article; no circadian, diurnal or time-of-day content. (Searched the full text.) The GTEx circadian use rests on controlled-access phenotypes (dbGaP) and the pilot description (lonsdale2013genotype).
- No statement on whether donors are a healthy reference: donors are deceased, with death circumstances varied; see lonsdale2013genotype for eligibility.

## Link to us
- **Intro claim(s):** our refs row (838 donors, 15,201 samples, 49 tissues; time of death recorded, controlled access) is correct for the numbers; time of death recorded is NOT in this paper (it is in dbGaP phenotypes and the data portal documentation). Cite it only for the cohort sizes.
- **Relation to clocksync:** data we use (healthy reference, multi-tissue per donor). v8 here, we use v10 per A109; cohort numbers will differ and should be taken from our own tables.

## Ideas
- None new. Ventilator/Hardy contrasts and donor-referenced analyses in GTEx: done (A27, A33, A109; METHOD.md validation table).
- Use Fig. 1A tissue sample sizes only to compare with our own per-tissue donor counts; our counts must come from our tables.
