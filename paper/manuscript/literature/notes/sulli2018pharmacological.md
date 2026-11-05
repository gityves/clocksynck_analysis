# sulli2018pharmacological
**Ref:** Sulli et al., Nature 553:351-355 (2018). DOI 10.1038/nature25170. Type: biology (cell lines, mouse). Lab: Panda (Salk).
**Source:** full text

## Question and design
Are REV-ERB agonists (SR9009, SR9011) selectively toxic to cancer cells? Cell lines from brain, leukaemia, breast, colon, melanoma (HRAS, KRAS, BRAF V600E, PTEN, b-catenin lesions)
versus normal cells; shRNA of REV-ERBs; p53-null/shRNA; hypoxia; NRAS-driven naevi (Tyr-NrasQ61K mice, SR9009 topical 20 uM, 4 mice per arm); orthotopic glioblastoma (005 cells in C57BL/6J, patient-derived GSC 8.11 in NSG), SR9009 100-200 mg/kg; TCGA/Rembrandt expression.
**Sampling:** not timed; drug dosing twice daily i.p., no clock-hour design; clock never measured (no luminescence, no clock-gene time course).

## Method (if any)
None (pharmacology). Lipidomics, autophagy flux, apoptosis (TUNEL, cleaved caspase 3).

## Key results
- Selectivity: SR9009/SR9011 "are specifically lethal to cancer cells and oncogene-induced senescent cells, including melanocytic naevi" (Abstract).
- Breadth: effective across HRAS, KRAS, BRAF, PTEN-deficient and b-catenin lines; "Agonist-induced apoptosis was largely intact in cells with compromised p53 function" (Fig 1, ED2).
- On-target: "The anticancer activity of SR9009 and SR9011 is abolished following the downregulation of" REV-ERBs by shRNA; "Levels of REV-ERB mRNA are comparable between normal cells and their transformed counterparts" (Fig 1c).
- Mechanism: block of de novo lipogenesis (FAS, SCD1) and autophagy; ROS and hypoxia not involved (ED2-3).
- In vivo: apoptosis in NRAS naevi (TUNEL P = 0.0058; Ulk1 P = 0.025, Atg7 P = 0.007); glioblastoma: "SR9009 effectively and significantly improved survival in two glioblastoma models" (vehicle n = 8 vs SR9009 n = 9, log-rank P = 0.0009, Fig 4j), similar to temozolomide without toxicity.
- NR1D2 expression correlates with survival in brain cancer (Rembrandt, log-rank P < 0.0001, Fig 4d).

## Not addressed / limits
- No readout of the clock itself: BMAL1, PER, or rhythms are not measured in treated or untreated cells; the paper does not test clock dependence or independence. Effect is shown to need REV-ERB, which also has non-clock (metabolic) functions.
- Sensitive cell types versus "normal cells" are a few primary or immortalised lines; no time-of-day dosing; no tumour-vs-host clock comparison.
- Later work questions SR9009 specificity (not in this text; verify before citing).

## Link to us
- **Intro claim(s):** refs line 57 says selectively lethal to cancer cells and oncogene-induced senescent cells, independent of the clock (mouse, cell lines). **Overstated/unsupported: independence from the clock is not shown** (clock not measured; REV-ERB is a core clock component and activity is REV-ERB-dependent). Selective lethality in lines and mouse is supported.
- **Relation to clocksync:** background.

## Ideas
- Wording fix for refs: REV-ERB-dependent, p53-independent selective toxicity; clock involvement not tested.
- Compare with Dong 2019 (GSC need BMAL1/CLOCK; REV-ERB agonists also act) for a clock-targeting-in-glioma sentence; Puram 2016 shows the same with SR9011 in AML.
