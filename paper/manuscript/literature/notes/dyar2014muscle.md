# dyar2014muscle
**Ref:** Dyar et al., Molecular Metabolism 2014 (3:29-41; erratum 3:857 corrects an author name only). DOI 10.1016/j.molmet.2013.10.005. Type: biology (muscle-specific knockout, glucose metabolism; transcriptome). Lab: Schiaffino / Dyar (Padova).
**Source:** full text

## Question and design
What does the intrinsic skeletal muscle clock do for glucose metabolism? Two muscle-specific Bmal1 knockouts: mKO (Bmal1 flox x Mlc1f-Cre, active in skeletal muscle fibres) and the inducible imKO (Bmal1 flox x HSA-Cre-ER, human skeletal actin promoter, tamoxifen 1 mg/day for 5 days at 2 months, examined 1 or 5 months later); Cre-negative littermates as controls (tamoxifen also in imKO controls). Tibialis anterior (TA, fast) and soleus (SOL, slow) muscles: Affymetrix microarray, qPCR, insulin-stimulated glucose uptake, PDH activity, metabolomics, telemetric locomotor activity. **Sampling:** timed under 12:12 LD, "Muscles were collected at 4-h intervals" (ZT0-ZT20), three biological replicates per time point and group for the arrays; GEO GSE43071.

## Method (if any)
maSigPro for differential time profiles; JTK_CYCLE (Benjamini-Hochberg q < 0.2) for 24 h rhythmic genes; GSEA on KEGG sets. No inference from untimed data.

## Key results
- **Driver and targeting.** Recombination: "showed a disrupted Bmal1 gene in skeletal muscles but not in heart and liver" (mKO, Southern blot, Supplementary Fig S1B); in imKO "Cre-mediated recombination of the Bmal1 gene was seen in skeletal muscles but not in the heart" (Fig S1F). BMAL1 protein "markedly decreased" in muscle. Recombination efficiency per fibre type or in non-fibre cells of muscle (satellite, endothelial, fibroblast, immune) is not quantified in the text read.
- **Targeted tissue's clock.** Whole TA: Bmal1-dependent genes Rev-erb-alpha, Rev-erb-beta, Dbp, Tef, Hlf, Cry2 "showed markedly reduced oscillation in TA muscles from mKO mice" by microarray and qPCR; similar in SOL. Clock, Npas2 and Nfil3 (E4bp4) are up; Per1 is reduced.
- **Remaining rhythm in whole tissue.** "Per1 expression was also reduced, but Per2 levels were only slightly affected", and the authors say Per2 "was only mildly affected in mKO skeletal muscles" and attribute it to oscillating systemic signals, as in liver knockouts. Of 684 cycling genes in control TA and 1359 in control SOL, "487 genes in TA and 666 genes in SOL lost their 24-h oscillation in the absence of Bmal1": so 197 of 684 (29%) in TA and 693 of 1359 (51%) in SOL keep a 24 h rhythm in the knockout (derived counts; whether they are rhythmic in the knockout or just not lost is not tabulated in the text read). The paper does not attribute the remaining rhythmic genes to non-fibre cells or to systemic cues beyond Per2.
- **Other tissues and behaviour.** Cardiac muscle Bmal1 level and rhythm "were unchanged in cardiac muscle" (Fig S1D); the circadian locomotor rhythm is normal in mKO and imKO, dark-phase activity higher in mKO only. Liver and other tissues' clocks are not measured beyond recombination PCR. Muscle phenotype mild (no sarcopenia as in whole-body knockout); insulin-stimulated glucose uptake reduced in both lines, GLUT4 and TBC1D1 lower, PDH activity reduced (Pdk4 up, Pdp1 down).

## Not addressed / limits
- Whole-muscle arrays: myofibre versus non-fibre cell contributions to the remaining rhythm are not separated; the paper does not discuss cell-type composition.
- 4-h sampling, n = 3, in LD: cannot separate autonomous from feeding- or activity-driven rhythms in non-targeted cells (Per2 attributed to systemic signals without a test).
- Rhythmic-gene lists in the knockout (what remains) are not shown in the text read beyond the lost counts; phase of residual genes not reported.
- Tamoxifen model recombination is shown by PCR; no single-cell resolution.

## Link to us
- **Intro claim(s):** a muscle-specific clock knockout leaves behaviour, heart clock and whole-body glucose tolerance largely unchanged: supports organ-level modularity. For within-tissue organisation: **partly informative, mostly about the targeted fibres**; the retained half of cycling genes in SOL and Per2 persistence would be the non-fibre or systemic share but cannot be assigned. Not evidence for coupling or pacing between cell types.
- **Relation to clocksync:** data we could use. GSE43071: control versus mKO TA and SOL, 6 time points, n = 3, i.e. a direct handle on how much rhythmic signal the whole muscle keeps without the fibre clock; compare with human skeletal muscle (GTEx) where the fibre clock dominates the bulk.

## Ideas
- Download GSE43071 and, per muscle, compute clock-gene amplitude and phase in mKO versus control; check whether the residual 24 h genes (not lost) are enriched for endothelial, fibroblast, immune markers (our stromal-vascular leaders), and whether their phase equals control.
- Muscle-specific rescue paper cited here (BMAL1 constitutively expressed in muscle of Bmal1-null mice) is a related design to compare with Koronowski's liver-RE.
