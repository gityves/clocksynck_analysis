# talamanca2023sexdimorphic
**Ref:** Talamanca, Gobet, Naef, Science 379:478 (2023). DOI 10.1126/science.add0846. Type: method (CHIRAL) + data
atlas. Lab: Naef (EPFL), our lab.
**Source:** full text of the main paper; methods and supplementary not in export (TIP/DIP details: bioRxiv
10.1101/2022.06.13.494714, read earlier).

## Question and design
24 h rhythms across the human body, by sex and age. GTEx v8: ~16,000 RNA-seq samples, 914 donors, 46 tissues.
**Sampling:** post-mortem; time of death (TOD) available but unreliable as phase; DIPs "distributed fairly uniformly
along the 24-hour cycle". Seed set benchmarked on time-labelled human muscle (Perrin 2018).

## Method
Covariate removal (linear model) → per tissue, CHIRAL infers tissue internal phases (TIPs) from 12 clock seed genes
→ donor internal phase (DIP) = primary mode of a donor's TIPs, referenced to skeletal muscle. Assumption stated:
"assumes that relative circadian phases of tissues are conserved across donors". Rhythmicity by harmonic regression
on DIPs; sex/age contrasts by model selection (Weger 2021). cSVD summarises the clock across tissues.

## Key results
- TOD is a poor phase: "Time of death (TOD) (available from GTEx) may not reflect circadian phase because of the individuals’ varying chronotypes"; PER2, NR1D1 arrhythmic vs TOD.
- Seed set: 12 genes; "performance decreases with the set size".
- Synchrony: "Clocks across tissues were well synchronized, showing relative offsets of only a few hours; the adrenal gland had the earliest phase".
- Amplitude: "Metabolic tissues (adipose tissue, esophagus, cardiovascular tissues) showed the highest amplitudes, whereas brain tissues and testis had the lowest".
- cSVD clock module: "The first mode, which captured >95% of the variance".
- Morning (7 a.m.) and evening (7 p.m.) waves of rhythmic genes body-wide.
- Sex: "The number of rhythmic genes was higher in females by about twofold at all amplitude thresholds".
- Age (<50 vs >60 y): rhythms "overall strongly damped in the older donors"; coronary artery about half the rhythmic genes in old; clock itself conserved.
- Own limits: "sensitivity to data quality, seed gene sets, complex covariate structures, or sampling bias cannot be fully ruled out".

## Not addressed / limits
- Synchrony is a population-level statement (tissue shifts common to donors); per-donor scatter of tissues around the DIP is not quantified, nor how it varies with disease.
- Requires DIPs spread over 24 h (GTEx has this; surgical cohorts do not).
- Single-tissue or two-tissue donors cannot be phased robustly (mode of TIPs needs many tissues).

## Link to us
- **Intro claims:** human tissue clocks are synchronised with offsets of a few hours (health baseline); amplitude highest in metabolic tissues, lowest in brain/testis; rhythms damp with age and differ by sex; TOD is not phase.
- **Relation:** our lab's precursor. clocksync generalises the fixed tissue shift (arg z_k) and adds coupling c and relative amplitude |z|, without assigning phases (works on hospital hours).

## Ideas
- Reproduce on GTEx the CHIRAL tissue offsets from our arg z_k (sanity figure: our offsets vs their cSVD tissue phases).
- Their brain/testis low amplitude predicts low coupling at fixed spread: check d* comparison is not just amplitude.
- Age: test whether coupling (not only amplitude) declines with age in GTEx at common d* (pre-register).
