# wu2019pan
**Ref:** Wu et al., Front Immunol 10:2451 (2019). DOI 10.3389/fimmu.2019.02451. Type: data analysis (pan-cancer, bulk). Lab: Fan / Mao (Nantong University).
**Source:** full text

## Question and design
Are the 14 core clock genes altered across cancers, why (copy number, mutation) and how do they relate to prognosis and to
the immune microenvironment? TCGA RNA-seq, CNV and SNV; 11 cancer types with matched normals (THCA, KIRP, LIHC, STAD, BRCA,
COAD, UCEC, BLCA, KIRC, KICH, PRAD; pairs with n > 10); abstract-level "716 samples across 29 cancer types". Immune
deconvolution (xCell), PARADIGM pathway scores. One mouse liver H3K4me3 24 h ChIP-seq series (re-used, 6 time points).
**Sampling:** human tumours single time point, time of day unknown and not discussed; only the mouse ChIP-seq is timed.

## Method (if any)
DESeq2 tumour vs normal (FC > 2, FDR < 0.05); GISTIC CNV vs mRNA correlation; Spearman rank GSEA of each clock gene against all
genes; correlation of clock genes with xCell immune scores and with checkpoint genes; log-rank survival (5 years, primary tumours).

## Key results
- Expression: "most of the negative regulators, such as PER1, PER2, PER3, CRY1, CRY2, and NR1D2, are down-regulated in a wide range of cancers" while "core positive regulators are upregulated (ARNTL2) or unchanged (CLOCK and ARNTL)".
- Origin: heterozygous deletion of PER1, PER2, PER3, CRY2, NR1D2 (CNV correlates with mRNA); ARNTL2 and NPAS2 heterozygous amplification; homozygous events rare.
- Mutations: "In 446 analyzed patients, 383 (85.87%) show at least one mutation"; PER1/2/3 most mutated (functional consequences unknown).
- Prognosis: high negative-arm expression goes with good outcome, ARNTL2 and NPAS2 high with poor outcome (log-rank only, no covariates).
- Pathways: 45 clock-linked KEGG pathways, Ras and autophagy most affected.
- Immune: clock genes correlate positively with Tregs and negatively with Th1, Th2, NKT and CD8 T-cell scores, positively with CD274, PD-L2, PD-1, CTLA-4 and TGFB1; PD-L1 promoter H3K4me3 rhythmic in mouse liver.

## Not addressed / limits
- "The circadian clock loses its daily rhythmic activity cycle" is inferred from mean mRNA levels; no oscillation, phase, amplitude or co-expression test is done on the tumours. Lower negative-arm mean with unchanged positive arm is not evidence of arrhythmia.
- Correlations with xCell scores are confounded by tumour composition and purity (no correction, no non-clock matched genes); clock-gene expression correlating positively with immune inhibitory genes sits uneasily with the "disrupted clock" reading and indicates composition effects.
- CNV deletions are a general feature of tumour genomes, with no control against random gene sets of similar locus or length.
- No stroma or host reference; "immune evasion" is correlation in bulk.

## Link to us
- **Intro claim(s):** negative-arm clock genes are downregulated by copy-number loss across TCGA tumours. SUPPORTED as stated in our refs row; the link to T-cell exhaustion is bulk correlation, and "disrupted rhythm" is an inference, so do not cite it for loss of oscillation.
- **Relation to clocksync:** background; clock-gene mean shifts in tumour vs normal are exactly what clocksync's constant-offset cancellation is built to ignore.

## Ideas
- Done in part: clock-gene copy number covariate in malignant coherence (A134: copy number accounts for about a tenth of the gap); H138.4 on TCGA pairs.
- Add the same pseudo-gene control (expression- and copy-number-matched gene sets) to the immune-correlation claim (not registered).
