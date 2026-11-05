# relogio2014ras
**Ref:** Relógio et al., PLoS Genet 10(5):e1004338 (2014). DOI 10.1371/journal.pgen.1004338. Type: biology (cell lines) + model. Lab: Kramer (Charité), Sers (Charité), Herzel (Berlin).
**Source:** full text

## Question and design
Why is the clock disturbed in cancer, and can one oncogene do it? Human colon cancer lines (HT29, RKO, SW480, LIM1215, CaCo2,
HCT116, plus Colo205, SW620, SW403, HKe3, HKe3 clone8), U2OS as test line, HaCaT keratinocyte progression series (immortal, I7, II4,
metastatic HRAS-transformed A5RT3), rat 208F fibroblasts with inducible H-Ras (IR2, IR4). Bmal1-promoter luciferase imaging (n = 3), microarray
of 8 lines (leave-one-out), text-mining network, ODE clock model. No patient samples, no tumour-vs-normal tissue, no host.
**Sampling:** continuous in vitro luminescence over 3-6 days after synchronisation (dexamethasone, temperature); qPCR every 3 h from 24 h after sync. No clock-hour sampling of patients.

## Method (if any)
Strong/weak oscillator classification (period + amplitude variation >= 20%), 45-gene discriminator (moderated t-test), ODE model in which RAS lowers BMAL1-mediated transcription.

## Key results
- Heterogeneous clocks inside one tumour type: "we found strong and weak circadian oscillators within the same type of cancer" (Abstract; 2 strong, 4 weak of 6 colon lines, Fig 1).
- 45-gene signature classified "seven out of eight cell lines" correctly, binomial p = 0.035 (Fig 2A).
- HaCaT vs HRAS-transformed A5RT3: "cells showed an advanced phase, of approximately 6 hours" under temperature entrainment (Fig 3D); Bmal1 mRNA "oscillate with a phase about 6 hours shorter in HaCaT cells than in HaCaT A5RT3 cells" (Fig 7A).
- Inducible RAS: "overexpression of RAS disrupts the circadian clock and leads to an increase of the circadian period while RAS inhibition causes a shortening of period length, as predicted by our mathematical simulations" (Abstract); MEK inhibitor UO126 shortens period (HaCaT 22.9 vs A5RT3 24.5 h, treated 22.2 h).
- Kras induction in HKe3 clone8: "HKe3 clone8 cells treated with mifepristone show a large period increase" (25.1 to 37.8 h, Fig 5D).
- Conclusion: "Together, our data demonstrate that perturbations induced by a single oncogene are sufficient to deregulate the mammalian circadian clock."

## Not addressed / limits
- Immortalised lines in culture, one synchronisation; no in vivo tumour, no tissue, no host clock: cannot say whether tumour clocks in patients are damped or shifted.
- Strong/weak classification uses n = 8 lines; the 45-gene list is not clock genes (mostly chromatin, signalling).
- Oncogene effect shown on period/phase of a Bmal1 reporter in cell-autonomous settings; in a body, coupling to host cues would also act.
- Does not compare tumour vs adjacent normal.

## Link to us
- **Intro claim(s):** Oncogenic Ras deregulates the circadian clock in cancer cells (cell lines, model) (refs_cancer_phase_inference.md line 49): supported (cell lines plus model; also HaCaT, rat fibroblasts).
- **Relation to clocksync:** background. Gives a mechanism for tumour-cell clock shifts (period lengthening, ~6 h phase offset) that clocksync would read as phase offset / reduced coupling of the malignant compartment against the host.

## Ideas
- Within-patient malignant-vs-normal phase offset in the paired single-cell cohorts could be compared to RAS/BRAF/MAPK mutation status (mutation calls needed). Not found in A138/A149 (A149 pairs malignant vs normal epithelium, not stratified by RAS status).
- Missing refs: Lévi 2010 (chronotherapy, cited here as [25]) and Fu 2002 (Per2 knockout cancer) are background.
