# kinouchi2020metabolic
**Ref:** Kinouchi and Sassone-Corsi, Nat Rev Cancer 20:645-661 (2020). DOI 10.1038/s41568-020-0291-9. Type: review. Lab: Paolo Sassone-Corsi (UC Irvine).
**Source:** full text

## Question and design
Review of the clock and tumour metabolism: clock gene expression in differentiated vs undifferentiated and tumour cells, non-clock transcription factors (MYC, HIF1a, SREBP, TFEB/TFE3) that link to the clock, and tumour-host communication. No new data, no sampling design; evidence is from mouse models (SCN ablation, jet lag, tumour-bearing mice), cell lines and some human tumour expression data.

## Method (if any)
None.

## Key results
- Clock attenuation in tumours: "expression of core clock genes or proteins is remarkably attenuated particularly in tumours of a higher stage or that are more aggressive". Low BMAL1 is linked to reduced differentiation and poor prognosis in pancreatic adenocarcinoma, melanoma and colorectal cancer.
- Link to differentiation: "undifferentiated, pluripotent stem cells in which the core clock does not exhibit overt oscillatory behaviour"; the hypothesis stated is "We hypothesize that the metabolic switch between differentiation or stemness of cancer may be coupled to the molecular clockwork".
- Host: the SCN connects to peripheral tissues by autonomic and endocrine signals, "enabling peripheral clocks to be temporally synchronized among each other".
- Tumour to host: "circadian gene expression in the liver is distally reprogrammed in mice bearing lung adenocarcinoma, breast cancer" (ref. 239; "distinct subsets of oscillatory genes and metabolites were enriched uniquely in the liver", via lactate and IL-6 proposals).
- Patients: "patients with cancer frequently exhibit circadian disruption of the activity and/or plasma levels of melatonin, cortisol and lymphocytes".
- Chemotherapy context: "the psychological stress linked to chemotherapy within a hospital setting, together with adverse effects such as nausea and vomiting, may impact sleep and host circadian rhythms"; sleep-wake cycles "can be dissociated from circadian rhythms".
- Conclusion on the loop: "the cell-autonomous clock in cancer cells appears to be frequently perturbed, facilitating proliferative metabolism in tumours" (Temporal communication section).

## Not addressed / limits
- No quantitative claim of desynchrony between organs of one patient; inter-organ signals are described mechanistically and from mice.
- Most tumour clock data are expression of single clock genes or lines, not coherence of the whole clock; no measure of coupling, amplitude or phase offsets.
- Human data are mostly bulk tumour expression without time of day.
- Directionality (tumour disrupts host vs disrupted host favours tumour) is described as a "vicious cycle", not resolved.

## Link to us
- **Intro claim(s):** cancer interaction: tumour clock attenuated, tumour reprograms host organ rhythms (mouse liver), host disruption feeds back. Our refs entry (review: circadian homeostasis, metabolism and tumorigenesis; not confirmed): supported for these points; it is a review, so it supports background only and not any human internal desynchrony.
- **Relation to clocksync:** background. The tumour-bearing liver result motivates reading tumour and host compartments against each other; the review gives no human measurement of it.

## Ideas
- Cite Masri 2016 (tumour-bearing mice, liver) as the primary source for host reprogramming and test the same in clocksync terms, host liver vs lung/WAT, if a timed tumour-bearing mouse dataset is available (not checked against our docs; no hit for tumour-bearing or Masri in METHOD.md, registered_experiments.md, results_preliminary.md).
- Check stemness and differentiation as a covariate of clock coupling in single-cell tumour data (cf. proliferation covariates in A97 and A171).
