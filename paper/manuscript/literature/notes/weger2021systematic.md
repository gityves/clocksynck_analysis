# weger2021systematic
**Ref:** Weger, Gobet, David, Atger, Martin, Phillips, Charpagne, Weger, Naef, Gachon, PNAS 118:e2015803118 (2021). DOI 10.1073/pnas.2015803118. Type: method + biology. Lab: Naef (EPFL) with Gachon (Nestle, then Queensland).
**Source:** full text (main text; SI Appendix and Methods details not in the export)

## Question and design
Separate clock-driven from feeding-driven (system-driven) rhythmic gene expression in mouse liver. Genotypes: Bmal1 KO, Cry1/2 KO, Hlf/Dbp/Tef KO and WT; feeding:
ad libitum (AL) vs night-restricted (NRF) under 12 h light : 12 h dark. **Sampling:** timed, liver RNA-seq "these mice using RNA-Seq (sampled every 4 h, two independent" replicates; PARbZip KO under AL only. Number of days covered not stated in the extracted text. Mouse only.

## Method (if any)
dryR (Differential Rhythmicity Analysis in R, github naef-lab/dryR): harmonic regression plus Bayesian model selection (BIC weights) over a
combinatorial set of models in which amplitude and phase are shared across subsets of conditions, and separate models for mean expression; tailored to RNA-seq count noise.
Key idea: "and phase) for each condition can be shared across subsets" of conditions. Outputs: gene classes clock-driven, food-driven (system-driven), clock-modulated, clock- and food-independent;
differential amplitude (log2 fold change peak-to-trough), phase, mean. Benchmark on simulated data (four conditions, "sampling interval of 4 h and two replicates") against other rhythm and differential-rhythm tools (CircaCompare, DODR, LimoRhyde), comparable sensitivity and specificity.
Requires time-stamped samples covering the cycle in each condition.

## Key results
- Core KOs lose feeding rhythms under AL; NRF restores them: "Bmal1 and Cry1/2 KOs lost their feeding rhythms under AL feeding" (Fig. 1C).
- Food-driven rhythms are 22 % (Bmal1 KO) and 31 % (Cry1/2 KO) of rhythmic genes, with bimodal phases "and on average 36% lower amplitudes" than clock-driven or clock-modulated genes. Clock- and food-independent: 32 % (Bmal1) and 14 % (Cry1/2), e.g. Cirbp, Hsp90ab1.
- With NRF KOs: "Of all rhythmic genes in WT mice, 33% were independent of" a functional clock (system-driven); clock-driven genes: "genes, which lost rhythmicity in the KO mice, represented 20%" of rhythmic genes, plus 7 % clock-modulated.
- Summary of amplitude: "circadian clock drives gene expression rhythms with high amplitudes covering all phases." while system-driven rhythms are low amplitude and bimodal.
- Cry1/2 KO: "subset of rhythmic genes exhibited a phase advance of ∼4 h" (PPARα targets; Bmal1 KO does not).
- PARbZip KO: 18 % fewer rhythmic transcripts; 59 % of WT rhythmic genes unaltered, 28 % lost rhythm, 10 % gained; mean levels of PARbZip targets more affected than rhythms.

## Not addressed / limits
- Liver only, mouse, standard light-dark with defined feeding; no multi-tissue or human data, no disease, no estimation of unknown time.
- Two replicates per time point at 4 h resolution; amplitude estimates at low counts are noisy (simulations show higher noise at low mean counts).
- Clock-independent systemic rhythms (temperature, glucocorticoid) are only partly resolved; they are "clock-independent" but not attributed to a specific signal.
- A rhythm lost in a KO is not proof of direct clock targeting; PARbZip-dependent genes are mostly indirect.

## Link to us
- **Intro claim(s):** refs_cancer_phase_inference.md claim: dryR, Bayesian model selection for differentially rhythmic genes, requires timed series. SUPPORTED. refs_cancer_phase_inference.md also lists this key as the feeding vs clock dissection (Bmal1 KO with timed feeding): SUPPORTED for Bmal1 and Cry1/2 KO; Hlf/Dbp/Tef KO was AL only.
- **Relation to clocksync:** competitor in spirit (compare?): dryR compares rhythm parameters between conditions given time; clocksync targets the case where time is unknown. Also background for what systemic (feeding) input does to non-clock genes.

## Ideas
- Where time is known (GTEx time of death, mouse atlases), compare clocksync's relative amplitude and phase offsets with dryR differential phase and amplitude between tissues/conditions on the same samples.
- Use the clock-driven vs food-driven gene lists (Dataset S2) to restrict the gene panel for a donor's clock state, and test whether clocksync's gate rejects Bmal1 KO liver.
- Related reference we lack: Pelikan 2021 (compareRhythms) already in refs_*; add Thaben and Westermark 2014 (rhythmicity detection, DODR) if comparing method classes.
