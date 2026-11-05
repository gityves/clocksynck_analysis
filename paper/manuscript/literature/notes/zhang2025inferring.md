# zhang2025inferring
**Ref:** Zhang Z, bioRxiv 2025 (10.1101/2025.01.18.633712), single author. Type: method + cancer application. Lab: not stated (single author, no affiliation given).
**Source:** full text (preprint)

## Question and design
Timestamp-free test for clock disruption in tumours from clock-gene phase differences. Calibration on mouse multi-tissue time series (Deota, 22 tissues, every 2 h over 24 h, ad libitum and time-restricted feeding); human reference matrices from GTEx with CHIRAL times (Talamanca); application to TCGA tumour vs normal (8 cancer types with GTEx references, 14 in the TCGA-normal analysis) with overall survival. **Sampling:** reference data are full 24 h; TCGA tumours untimed.

## Method
- Unit: group of samples (phase difference variation, cohort level) and, for ssDistance, single sample by leave-one-out.
- Idea: correlation of two sines is the cosine of their phase difference, so co-expression is converted into phase differences between 13 conserved peripheral clock genes (ADAMTS4, ARNTL, BHLHE41, CIART, DBP, HLF, NFIL3, NPAS2, NR1D1, NR1D2, PER2, PER3, TEF); tumour matrix minus tissue-specific reference matrix = delta phi; ssDistance = Euclidean distance to the reference, normalised.
- Supervised by the reference: gene acrophase relations are fixed from healthy tissue. No explicit noise model. Needs a cohort (co-expression), not one sample, but no sample times. Does not need full-cycle coverage of the test samples, but assumes samples are spread on the cycle so that covariance reflects phase differences.
- Does not compare tissues within an individual.

## Key results
- Co-expression to phase difference works in mouse: "across 22 mouse tissues, the consistency between observed and predicted values" for phase differences and co-expression exceeded 0.9; median phase-difference error 1.58 h (mouse tissues), and for the 13 genes in humans median 2.8 h (2.4-3.5 h) (Fig S2).
- TRF case: predicted 4 h vs observed 4.3 h change (Adamts4-Nfil3); median consistency of delta phi 0.77, MAE 0.91 h (Fig S5A).
- TCGA: "tumor samples had a significantly larger Euclidean distance" (Wilcoxon test, p = 0.028); ARNTL, PER3, NPAS2, NR1D2 altered in all eight cancers.
- ssDistance recovers knockouts: AUC 0.80 (Per1/2 double KO) and 0.86 (Bmal1 KO) in mouse liver.
- Survival: ssDistance stratifies OS in 5 of 8 cancers (COAD p = 0.022, LUSC 0.013, PRAD 0.035, THCA 0.0045, UCEC 0.0074) with opposite directions across types.

## Not addressed / limits
- Own caveat: "our current method struggles to distinguish" between the effects of cycle length and phase differences.
- Assumes sinusoids for reference genes, unverified in tumours; TCGA normals are tumour-adjacent; no tumour purity/composition control (covariance also reflects cell-type mix); human error 2.8 h, large relative to the effects.
- Direction of ssDistance and survival is described inconsistently in the text (high ssDistance is defined as less deviation from the reference, yet high ssDistance is also read as a larger deviation with longer survival).

## Link to us
- **Intro claim(s):** supported as a sample-time-free, within-sample gene-pair method for tumours (cancer interaction / technical challenges). Our refs line is accurate; it is a competitor only in spirit, it does not compare compartments of a donor.
- **Relation to clocksync:** competitor / background. Same input type (clock-gene covariance, no times) but gene-gene phase differences within one tissue versus our donor-by-compartment shared state. Under hospital hours it works (no times), but cannot give coupling.

## Ideas
- Compute ssDistance on GTEx/TCGA tumour-normal pairs and correlate with clocksync coupling c on the same donors (need a specificity companion).
- Check whether the gene-pair phase-difference matrix changes with tumour purity (stromal fraction), the confound our gate must also handle.
