# broadberry2018disrupted
**Ref:** Broadberry et al., Breast Cancer Res 20:125 (2018). DOI 10.1186/s13058-018-1053-4. Type: clinical / biology (human tissue). Lab: Meng and Streuli (Manchester).
**Source:** full text

## Question and design
Is the molecular clock of breast epithelium altered in primary human breast tumours compared with normal breast of the same
patient, and does tumour stroma differ in structure and stiffness? 8 women (46-78 y, no metastasis; Table 1: mostly ER/PR/Her2
negative, mixed histologies) after mastectomy: tumour area and normal tissue at least 4 cm from the margin, same breast.
Per2::Luc reporter in primary epithelial cells (3 patients), tissue qPCR of Bmal1, Clock, Nr1d1 (8 patients), collagen
(n = 6) and atomic force microscopy (5 patients).
**Sampling:** one surgical specimen per region; time of day of surgery or tissue collection not reported; luminescence is ex vivo in culture.

## Method (if any)
Paired tumour vs normal comparison inside each patient (paired ratios); lentiviral Per2::Luc rhythm in cultured cells; ddCT qPCR
normalised to beta-actin.

## Key results
- Culture clock: normal epithelial cells "had a robust ~ 24 h Per2::Luc rhythm"; tumour cells from the same individuals "displayed a weakened rhythm, with much lower amplitude, and the rhythm was not sustained", in all 3 patients.
- Tissue: Bmal1 mRNA lower in tumour than in normal tissue from the same patient (p < 0.0005, n = 8); Clock and Nr1d1 "more variable and showed no consistent changes" (Fig. 5).
- Stroma: more organised fibrillar collagen in tumour regions (n = 6); stiffness "average 30% increase" (1.21 vs 0.90 MPa, p < 0.0001, n = 5) next to tumour ducts.
- Interpretation: stiffer stroma may dampen the epithelial clock (their earlier mouse work).

## Not addressed / limits
- Tiny cohort (3 patients for rhythm, 8 for qPCR); heterogeneous subtypes; no non-clock control genes for the qPCR change.
- Bmal1 mRNA is a single time point at unknown clock time: a lower level can arise from phase, composition (epithelial vs stromal vs fat fraction) or a clock change, and the paper does not separate them. Time of resection is exactly the confounder clocksync removes.
- Cultured epithelial rhythms are cell-autonomous and come after enzymatic digestion; no stromal clock, no host reference, no causal test of stiffness; the authors state that "we still do not know whether the reduction of circadian rhythm is just a consequence of tumorigenesis".
- "Normal" tissue comes from the tumour-bearing breast.

## Link to us
- **Intro claim(s):** clocks are weakened in primary human breast tumour epithelium relative to normal tissue of the same patient. Our refs row is SUPPORTED, but paired human samples means 3 (rhythm) and 8 (qPCR) patients.
- **Relation to clocksync:** background / hypothesis source (stroma stiffness vs clock amplitude); the only human matched tumour-normal functional clock readout.

## Ideas
- Done on larger human pairs: H138.4 (TCGA 517 pairs), H138.5 (52 patients single-cell, malignant vs adjacent-normal epithelium 0.583 vs null 0.443, p 0.069) and A151 (tumour fibroblast/endothelium vs same-patient normal stroma). Resection-hour structure: GSE233242 (42 patients, resection 08-18 h) is scored in results_preliminary (A-series near line 3094).
- Check stromal amplitude: test whether tumour-adjacent stroma shows lower relative amplitude |z| than normal stroma in the same donors (stiffness hypothesis), not yet registered.
