# grabe2024quantification
**Ref:** Grabe, Ananthasubramaniam, Herzel, Sci Rep 14:16238 (2024). DOI 10.1038/s41598-024-66694-7. Type: method + data (TCGA re-analysis). Lab: Hanspeter Herzel / Bharath Ananthasubramaniam (Charite, HU Berlin).
**Source:** full text

## Question and design
Compare clock strength in healthy vs cancerous human lung from untimed snapshots. Time series for gene selection: mouse lung Zhang (GSE54652, 24 samples every 2 h), mouse lung Wolff (GSE201207, 12 per age every 4 h, three ages), baboon lung Mure (GSE98965, 12 samples every 2 h). Target: TCGA-LUAD, "It includes 51 normal samples and 58 tumor samples, each normal sample is matched with at least one tumor sample of the same patient". It is TCGA, not GTEx. **Sampling:** gene-selection data timed, full cycle; TCGA untimed surgical samples, collection times unrecorded.

## Method
Unsupervised on the target, supervised on mouse and baboon for the gene list. (1) Spearman correlations of all genes with PER3 in each timed series; top and bottom 100 genes (phase difference 0 h or 12 h) plus a 6 h group (low r to PER3 and ARNTL, harmonic R2 > 0.8); intersect across datasets: 107 murine genes, a 13-gene mammalian list (HLF, TEF, DBP, KLF9, PER3, ... NPAS2, PDIA6, ARNTL, RORC). (2) In TCGA, compare each gene's Spearman r to PER3 with the r-versus-phase-difference relation from the timed data; circadian correlation matrices with cosine similarity. (3) PCA on the 13 genes, phase = arctan(PC2/PC1), rescaled to 24 h. Unit of inference: one sample, but the shared PCA circle needs the cohort to cover the cycle.

## Key results
- r to phase: "The phase differences between circadian expressed genes are linked tightly to the Spearman rank correlation" (PER3-ARNTL r about -0.97, PER3-TEF about 0.98, PER3-RORC about 0.13).
- Healthy TCGA lung keeps the mouse/baboon phase relations for most of the 13 genes; in tumour, "There is a strong alteration in the circadian behaviour of ARNTL and NPAS2, which are both shifted closer towards rS = 0"; RORC moves from about 0.2 to about -0.4; PER3, TEF, HLF, DBP unchanged.
- Paired design: "Next, we used the reconstructed circadian times of the normal samples from the TCGA data to temporally sort their 58 paired adjacent tumour samples."
- Normal reconstruction R2 ARNTL 0.76, PER3 0.61, RORC 0.19; tumour 0.11, 0.11, 0.02; "is increased by 8 h in the cancerous lung tissue" (PER3 to RORC).
- Caveat on normals: "the healthy human lung tissue was collected in vicinity of the tumor tissue, which could also have increased its similarity to the tumor tissue".

## Not addressed / limits
- No ground truth for human lung. Quality is shown on mouse and baboon only (PCA reconstruction, including mixed-age mouse data).
- Cohort coverage is assumed: "the analysis of circadian rhythms in gene expression data is facilitated, if a dataset covers a full circadian cycle". They admit "clinical samples are collected during surgical working hours that are not evenly spread across the day" and hope chronotypes and shift workers fill the circle; incomplete coverage is listed as an untested option ("An interesting option to further test the robustness of correlation-based time allocation with incomplete temporal coverage of the circadian day in different human tissues").
- One cancer, one gene list, no gate that the shared covariance is clock-shaped, no tests, no confidence intervals; amplitude and R2 depend on the normal sample time estimates.
- Adjacent normal is not a healthy donor.

## Link to us
- **Intro claim(s):** single-compartment snapshot quantification of cancer clocks (technical challenges) and clock weakened in lung cancer. Refs entry says "(human GTEx lung)": **CLAIM WRONG**, the human data are TCGA-LUAD paired normal / tumour; GTEx is only cited through Talamanca. Not single-compartment in one respect: the tumour is ordered with the paired normal's time.
- **Relation to clocksync:** closest precedent of host-referenced inference (normal gives the time, tumour read against it), but with a within-person pair, cohort PCA, and cohort coverage assumed. Also competitor: PER3 correlation vs phase difference.

## Ideas
- Apply their 13-gene PER3-r shift and the paired normal-times-order-tumour procedure to TCGA-LUAD and our other paired cohorts, against clocksync coupling. Related, done: A171 (LUAD clock-gene coherence, TCGA LUAD -0.361, replicated in GSE72094, GSE68465); the Grabe procedure itself is not run (no hit).
- Hospital-hours test of their cohort-coverage assumption using the baboon 8 h windows (A40, perfect ordering) as the positive control.
