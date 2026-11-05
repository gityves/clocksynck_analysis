# fortin2024circadian
**Ref:** Fortin et al., Nat Immunol 25:1257-1269 (2024). DOI 10.1038/s41590-024-01859-0. Type: biology (mouse genetic CRC model). Lab: Masri (UC Irvine), Kessenbrock, Lawson.
**Source:** full text

## Question and design
Does the epithelial clock control tumour immunosuppression and the timing of checkpoint blockade? Villin-Cre mice: WT, Bmal1fl/fl (epithelial Bmal1-/-), Apc+/-, Apc+/-;Bmal1-/-; "the clock is only disrupted in intestinal epithelial cells and remains intact in the periphery". scRNA-seq of CD45+ small-intestine cells (15,234 cells, n = 3 mice per genotype), flow cytometry (n = 4-8 per genotype), intestinal monolayers with cytokine arrays, environmental shift (jet lag/night shift) paradigm, anti-PD-L1 (every 2-3 days) in Apc+/-;Bmal1-/- and subcutaneous MC38 (CRC), CMT167 (lung), D4M-S (melanoma).
**Sampling:** scRNA-seq at a single time (ZT4); rhythm assessed by two times only (ZT4 rest, ZT16 early active) by flow; treatment at ZT4 vs ZT16. No full 24 h series.

## Method (if any)
None; scRNA-seq clustering and differential proportions; no phase estimation.

## Key results
- Epithelial clock loss accelerates tumours: more polyps and spleen weight in Apc+/-;Bmal1-/- than Apc+/- (Ext Data Fig 1c, d).
- Immune landscape: "Apc+/− and Apc+/−;Bmal1−/− mice had increased proportions of neutrophils as well as decreased proportions of CD8+ T cells" (Fig 1); also in Bmal1-/- without tumour and after chronic shifted light (SD).
- MDSC rhythm: "In the spleen, where the circadian clock remains intact, we identified a significant increase in Gr1+ and Gr1+PD-L1+ cells at ZT 16"; the intestinal rhythm is lost in Apc+/-;Bmal1-/- (Fig 5a-d); cumulative cells higher at ZT16 than ZT4.
- Mechanism: "IL-5, CCL5, IL-17 and CXCL9 were found to be exclusively clock-dependent" in Bmal1-/- monolayers, driving neutrophil migration and MDSC-like programmes (Fig 4).
- Therapy: "CR was observed for 1 of 12 (8%) tumors with ZT 4 treatment versus 4 of 14 (29%) tumors with ZT 16 treatment" (MC38, Fig 5j, k); anti-PD-L1 at ZT16 reduces intestinal Gr1+ cells and raises CD8+ T cells in Apc+/-;Bmal1-/- while ZT4 does not (Fig 5f-h); D4M-S also better at ZT16.

## Not addressed / limits
- The tumour clock is deleted, not measured: no tumour-cell clock readout in WT tumours; rhythm of MDSCs rests on two time points (ZT4, ZT16), not a 24 h time course; scRNA-seq at one time.
- Small n (MC38 12 vs 14 tumours; complete responders 1 vs 4); both times reduced MC38 growth. No human data; no adjacent normal vs tumour clock comparison.
- Epithelial clock KO is a strong manipulation; partial amplitude loss as in human cancer is not modelled.

## Link to us
- **Intro claim(s):** refs line 83 (epithelial clock disruption raises inflammatory signalling and PD-L1+ MDSC rhythms, and affects ICB efficacy): partly **overstated**. Clock disruption raises neutrophil/MDSC abundance and abolishes (not raises) the intestinal MDSC rhythm; scRNA-seq is a single-time snapshot; rhythm is two-point flow. Time-of-day efficacy is supported (ZT16 > ZT4) but with small numbers.
- **Relation to clocksync:** background; a tissue-compartment desynchrony example: epithelium clock disrupted while spleen/periphery intact and immune rhythm intact there (internal desynchrony between compartments in one animal).

## Ideas
- Compartment contrast (intestinal vs splenic immune clock in the same animal) is the mouse analogue of clocksync's within-donor design; a clocksync-type test on public mouse tissue series with epithelial Bmal1 deletion would be a controlled positive for "desynchrony across compartments (not found in METHOD.md validation table).
- Missing ref: Aiello 2020 (jet lag B16, already a card) and Masri 2016 (Apc;Bmal1 phenotype, cited as ref 39 here; verify).
