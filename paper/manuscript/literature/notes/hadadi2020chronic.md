# hadadi2020chronic
**Ref:** Hadadi et al., Nat Commun 11:3193 (2020). DOI 10.1038/s41467-020-16890-6. Type: biology (mouse model). Lab: Acloque (INRAE) / Hadadi (Inserm).
**Source:** full text

## Question and design
Does chronic circadian disruption (CRD) change mammary tumour progression, dissemination and the immune microenvironment?
MMTV:PyMT x MMTV:LUC mice (B6*FVB), 12:12 light-dark (LD, n = 25) vs chronic jet lag (JL, n = 21), 10 weeks from week 6,
analysed at 16 weeks. Bulk mRNA-seq of dissociated primary tumours and bone-marrow mononuclear cells (5 LD vs 5 JL),
flow cytometry of tumour-infiltrating cells, mammosphere and transplant assays, CXCR2 inhibitor (SB265610) in JL mice.
**Sampling:** one terminal time point per mouse. LD mice at ZT3-ZT4; for JL "ZT was not followed" (Methods). No
24 h series, no tumour-vs-host clock comparison; JL contrasted with LD only.

## Method (if any)
None (experimental). Contrast JL vs LD, unpaired t-tests / DESeq2; mammosphere efficiency of MCF12A cells at PER2/BMAL1
peaks after serum shock as an in vitro link.

## Key results
- Metastasis: "The proportion of mice with lung metastasis increased from 28% in LD to 52% in JL" (Fig. 2c); disseminated cells in bone marrow almost 2-fold (p = 0.008 flow, p = 0.035 qPCR).
- Primary tumour: tumour burden higher in JL (p = 0.029) but no difference in onset; transcriptome barely moves: it "did not distinguish between LD and JL conditions" in clustering or PCA (Fig. 3a).
- Clock genes in tumour: "We also observed a decrease of Per2 and Cry2 expression in primary tumours of JL mice" while "other core clock genes remained similar between treatments" (Fig. 4d).
- Largest changes are phototransduction genes in bone marrow: "Rhodopsin by 23-fold, Gnat1 by 11-fold, Rbp3 by 7-fold and Prph2 by 6-fold" down in JL.
- Immune microenvironment: fewer CD45+ cells in JL tumours, "JL tumours had a significantly higher proportion of tumour-supporting MHC IIlow TAMs", more Tregs, higher CD4/CD8 ratio; Cxcl5 up (log2FC > 0.8).
- CXCR2 inhibition lowers lung metastasis in JL, whereas in LD there was "no significant difference in the number of metastases detected in LD mice" (Supp. Fig. 10b): a JL-specific effect (the specificity control).

## Not addressed / limits
- Clock state is never measured in tumour, stroma or host over time: only means of Per2/Cry2 at one undefined circadian time; JL samples have no ZT, so a mean difference conflates clock damage and sampling phase.
- No tumour-vs-adjacent-normal or tumour-vs-host phase, amplitude or synchrony. "Desynchrony" is inferred from the light schedule, not from molecular readouts.
- Bulk tumour mRNA-seq (n = 4-5 per arm) is underpowered; whole-tumour composition changes (CD45+ fraction) can alter clock-gene means.
- Mammosphere oscillation is a cell line (MCF12A), not tumour.

## Link to us
- **Intro claim(s):** chronic circadian disruption in mice promotes metastasis and shifts the immune microenvironment (Cancer interaction); our refs row (stemness, immune changes, lung metastasis) is SUPPORTED.
- **Relation to clocksync:** background; motivates a readout of host-vs-tumour clock state but supplies no molecular synchrony data.

## Ideas
- No public series of tumour plus host organs under JL is used here; ask for the mRNA-seq (GEO, Supplementary Data 1-3) and test whether JL tumours and host marrow differ in clock-gene covariance (coupling) without ZT, the exact use case of clocksync (not done: no mouse JL dataset in registered_experiments).
- Companion mouse data with tumour + liver + clock-timed sampling: Masri 2016 / Hojo 2017 (host liver rhythms reprogrammed by distant tumours), see aiello2020circadian card.
