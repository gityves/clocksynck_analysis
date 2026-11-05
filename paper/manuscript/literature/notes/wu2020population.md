# wu2020population
**Ref:** Wu G, Ruben MD, Francey LJ, Smith DF, Sherrill JD, Oblong JE, Mills KJ, Hogenesch JB, Genome Med 12:73 (2020). DOI 10.1186/s13073-020-00768-9. Type: method + human data. Lab: John Hogenesch (Cincinnati).
**Source:** full text

## Question and design
Human skin clock: epidermis vs dermis at three body sites and a 12-gene epidermal phase biomarker (SkinPhaser). Longitudinal group n = 20 healthy Caucasian men (21-49 y), forearm biopsy at 6 am, 12 pm, 6 pm, 12 am (laser capture of epidermis and dermis, microarrays). Population group n = 154 women (20-74 y), one biopsy each from forearm, buttock and cheek, "sample collection times were not recorded". Test: 9 longitudinal participants held out (36 samples) and 18 participants in an independent German suction-blister set (9:30, 14:30, 19:30). **Sampling:** longitudinal = timed full day (4 points); population = untimed.

## Method
- Unit: single epidermal sample. Hybrid design: CYCLOPS orders the untimed population samples together with the timed ones (cohort level; needs a spread over the cycle, seeded with 17 clock genes, 519 epidermis samples used), then ZeitZeiger (supervised) learns 12 genes predicting the CYCLOPS-based phase. Gene acrophases fixed by the trained model; platform-specific (arrays). Noise: sparse PCA, no explicit model.
- Within person: the population design has epidermis and dermis (and three sites) from the same participant, but only clock gene correlation strength (Mantel Z with a mouse reference) is compared between compartments, not offsets.

## Key results
- Accuracy: "the average absolute error across different times of day for the 9 participants was 2.5 h (SD = 1.86 h) for epidermis compared to 3.8 h (SD = 3.92 h) for dermis"; independent set: 3.1 h (SD 2.12 h), order reproduced in 14/18.
- Epidermis clock stronger than dermis: Mantel Z 37.6 vs 32.5 (all sites); PER3 relative amplitude 0.59 vs 0.29 (P = 1.9e-6).
- Layers in phase: "the median phase difference between epidermis and dermis samples is 0.4 h for ARNTL and 0.6 h for NR1D2" (longitudinal, n = 19); ARNTL peaks about 7 h before melatonin.
- CYCLOPS order matches sampling time (circular correlation 0.973 epidermis, 0.925 dermis) on the longitudinal subset.

## Not addressed / limits
- Supervised by CYCLOPS pseudo-times (circular reasoning limits; error 2.5 h is relative to the pipeline phases and women are only in the untimed group).
- Single tissue (skin); does not test skin phase against another organ ("One important question is whether skin clock phase is a good predictor of clock phase in other tissues" is left open). Healthy only. Early-morning samples (6 am) predicted worst.

## Link to us
- **Intro claim(s):** supports that a one-sample biomarker (single biopsy) exists for skin and that CYCLOPS ordering needs a population spread over the day (their population is untimed but presumably spread). Our refs line (supervised, timed training) is partly right, but the training phase itself is CYCLOPS-derived, not timed.
- **Relation to clocksync:** competitor at tissue level (supervised phase from one sample) and a useful dataset: three body sites x two layers per participant in 154 women = natural compartments within person, untimed (a direct clocksync test of coupling across layers/sites).

## Ideas
- Apply clocksync to the 154-participant arm/buttock/cheek x epidermis/dermis data (GEO, Wu 2020, Spörl 2012) to check coupling c and offset between epidermis and dermis (expect offset about 0.4-0.6 h).
- Use SkinPhaser as an external phase for a held-out check of the host state in skin.
