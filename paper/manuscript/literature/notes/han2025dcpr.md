# han2025dcpr
**Ref:** Han, Cen, Li, Zhou, Ji, BMC Bioinformatics 27:31 (published 2026; accepted 22 Dec 2025; bib says 2025).
DOI 10.1186/s12859-025-06363-2. Type: method. Lab: Zhiwei Ji (Nanjing Agricultural University).
**Source:** full text

## Question and design
Unsupervised phase reconstruction (DCPR) from untimed bulk transcriptomes; benchmark against CYCLOPS, CHIRAL, Cyclum.
- 12 synthetic sets (1,000 genes, 200 rhythmic): 9 full cycle (24/48 h), 3 incomplete cycle (14 h, 18 h, 18 h),
  uniform or non-uniform intervals, 14-48 samples.
- Real timed sets: rat, human (GSE56931 blood, GSE161566), mouse liver GSE11923, Arabidopsis, wheat; clock KO sets
  (GSE196430, GSE135898).
- Application: untimed AD brain (EC 23, HIP 63, FC 71 samples; 85 control, 72 AD). **Sampling:** unknown (post-mortem).

## Method
PCA pre-sorting on seed genes (homologous rhythmic genes, or top 1,000 by p without prior), cosinor "data
augmentation" (genes replaced by their cosinor fit on the pre-sorted order), deep matrix factorisation + multi-head
attention + autoencoder with a circular 2D bottleneck; one phase per sample. Assumes samples lie on one circle:
cohort-level ordering, no within-person structure, no notion of tissue offsets.

## Key results
- Full-cycle synthetic AUC: "DCPR achieved an average AUC of 0.969 across all datasets, outperforming CHIRAL (0.951), Cyclum (0.932), and CYCLPOS (0.687)" (Results).
- Non-uniform sampling: "When dealing with non-uniform datasets, DCPR and CHIRAL exhibit more robust performance" (Results).
- Incomplete cycles: "Although all methods exhibited performance degradation, DCPR maintained relatively high accuracy. This result confirms the challenge of learning from incomplete cycles" (Results). Numbers only in Suppl. Table 2.
- Real timed data: AUC 0.803-0.897; "CHIRAL (~ 0.741) and Cyclum (~ 0.736) are behind".
- Sample size: "CYCLOPS typically requires more than 250 samples to reliably order circadian phases" (Discussion, citing Anafi 2017).
- Clock disruption: "CHIRAL's design fundamentally relies on the core clock genes of organisms with an intact, healthy circadian system" (Discussion; CHIRAL degrades in Bmal1/Cry KO).
- AD: loss of rhythm of CLOCK, CRY2, DBP in EC etc., from cosinor on inferred order; "validated" against mouse APP/PS1 qPCR read off published figures.

## Not addressed / limits
- "Incomplete cycle" = evenly spaced samples over 14-18 h of the cycle, still several time points spanning most of it; not a narrow window (e.g. 08-17 h ≈ 9 h) with phases inferred against a sampled block.
- AD rhythm calls on a circular order inferred from the same data (circularity), no specificity control, groups ordered separately (phase reference not shared between control and AD).
- Real human benchmarks: blood time courses only; no tissue, no hospital-hours cohort.
- Phase is a cohort-level quantity; nothing on desynchrony between tissues of one person.

## Link to us
- **Intro claim:** technical challenges. Supports the claims that unsupervised orderers degrade when the cycle is not fully covered (only benchmark of incomplete coverage) and that CYCLOPS needs >250 samples. Our refs table quote is correct.
- **Relation:** competitor class (orderer). A145 already compares CHIRAL/CYCLOPS/CoFE; DCPR could be added.

## Ideas
- Add DCPR to the orderer comparison on simulated hospital-hours windows (08-18 h), with the same metric (AUC of error CDF) so the result reads against their 14 h / 18 h numbers.
- Their 14/18 h synthetic sets are a ready template for our simulation figure (sampling window width vs accuracy).
- Missing reference: Cyclum (Liang et al. 2020 Genome Biol, "latent space" key liang2020latent already in bib).
