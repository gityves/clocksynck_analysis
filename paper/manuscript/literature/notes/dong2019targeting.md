# dong2019targeting
**Ref:** Dong et al., Cancer Discov 9:1556-1573 (2019). DOI 10.1158/2159-8290.CD-19-0215. Type: biology (patient-derived cells, xenografts, TCGA). Lab: Rich (UCSD) and Kay (USC).
**Source:** full text

## Question and design
Do glioblastoma stem cells (GSC) keep a clock and need it? Patient-derived GSCs (T387, T3565, others), differentiated glioblastoma cells (DGC), epilepsy-surgery
nonmalignant brain cultures (NM263, NM290), neural stem cells (NSC), astrocytes (NHA). BMAL1-Luc luminescence, shRNA/CRISPR of BMAL1 and CLOCK, BMAL1
ChIP-seq (GSC vs NSC), intracranial xenografts in immunocompromised mice, SR9009/SR9011 (REV-ERB), KL001 and SHP656 (CRY), TCGA glioma expression.
No tumour-vs-adjacent-normal pairing within a patient (NM comes from other patients).
**Sampling:** cell cultures synchronised with dexamethasone or forskolin, luminescence over days; mice treated at unspecified clock time; TCGA bulk, untimed.

## Method (if any)
Not a method paper; BMAL1 ChIP-seq peaks vs H3K27ac; GSEA.

## Key results
- Clock is intact in tumour cells: "patient-derived GSCs and their differentiated progeny displayed circadian rhythms with similar properties to nonmalignant brain cultures derived from surgical resections from patients with epilepsy" (Fig 1A-D); "Canonical rhythms observed in normal brain cells and GSCs suggest that cellular transformation maintains circadian rhythms, despite the activation of oncogenes."
- Selective dependence: "targeting BMAL1 or CLOCK minimally reduced cell proliferation in epilepsy-derived brain cultures or neural stem cells" but potently reduced GSC growth, with modest effect on DGCs (Fig 1E-L, S2).
- Mechanism: knockdown gives G2-M accumulation, apoptosis, lower SOX2, OLIG2, MYC; BMAL1 preferentially binds metabolic genes in GSC vs NSC (Fig 3).
- In vivo: "The life span of tumor-bearing mice nearly doubled upon BMAL1 or CLOCK knockdown relative to controls" (Fig 7A, B).
- Pharmacology: SR9011/SR9009 and KL001 reduce GSC proliferation more than DGC, astrocytes, NM; combination SR9011 + KL001 synergistic in vitro (Fig 6); CRY stabiliser SHP656 "prolonged the survival of mice bearing two different patient-derived GSCs compared with nontreatment" (Fig 7H, I).
- TCGA: "BMAL1 mRNA levels were elevated in glioblastomas relative to lower grades of glioma and other glioma histologies"; high BMAL1 poor prognosis, CRY2, REV-ERBa, PER2, PER3 favourable (Fig 7J-M).

## Not addressed / limits
- Tumour clock measured in culture, not in patients; no host-vs-tumour clock comparison and no time-of-day dosing. Dependence is attributed to non-circadian BMAL1 functions ("beyond canonical circadian circuitry"), so clock amplitude is not the readout.
- Synergy shown in vitro only; in vivo uses single agents (SHP656), not the combination. Number of patient lines modest (2 main GSC lines plus a few).
- Prognostic associations are bulk TCGA mRNA, uncorrected for composition.

## Link to us
- **Intro claim(s):** refs line 56 (GSC depend on BMAL1/CLOCK; CRY and REV-ERB agonists reduce growth, synergistic): supported, with the caveat above (synergy in vitro).
- **Relation to clocksync:** background; shows malignant cells retain rhythmic clocks and that clock-gene dependence is a tumour property, so malignant compartments in glioma single-cell atlases are legitimate clocksync inputs.

## Ideas
- Glioma single-cell cohorts: malignant vs non-malignant brain (neurons/microglia/oligodendrocytes) per patient would give compartments for clocksync; no glioma/GBM cohort appears in METHOD.md, registered_experiments.md or results_preliminary.md (not done; check available datasets).
- Missing ref candidates (verify before citing): Chen et al. 2018 (circadian genes and glioblastoma stem cell biology, check) and Slat et al. 2017 (Per2/BMAL1 in glioma, check).
