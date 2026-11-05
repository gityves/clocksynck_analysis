# puram2016core
**Ref:** Puram et al., Cell 165:303-316 (2016). DOI 10.1016/j.cell.2016.03.015. Type: biology (mouse + human cell lines). Lab: Ebert / Regev (Broad, Harvard). NB: first author is Rishi Puram (AML), not the HNSCC single-cell atlas of Sidharth Puram (2017).
**Source:** full text

## Question and design
Which transcription factors do leukaemia stem cells (LSC) need? Pooled in vivo shRNA screen (152 TFs, 770 shRNAs, 5 mice per subpool) in murine
MLL-AF9 AML (quaternary transplants), then Clock/Bmal1 validation (shRNA, CRISPR, cDNA rescue, Bmal1f/f conditional knockout, germline Bmal1-/-), human AML
lines (SKM-1, THP-1, NOMO-1, MOLM-13), ChIP-seq of CLOCK and BMAL1 in THP-1, TCGA AML (191 cases). Per2::Luc reporter mice for clock readout.
**Sampling:** untimed for the screens; clock read by continuous bioluminescence of spleen explants (mice in LD12:12). No time-of-day experiments.

## Method (if any)
Not a method paper. shRNA depletion relative to input; Per2::Luc explant rhythms.

## Key results
- Screen: "We identified 35 transcriptional regulators targeted by at least three shRNAs depleted >20-fold" including Clock and Bmal1 (several shRNAs up to 1,000-fold).
- Leukaemia clock is intact: "Continuous bioluminescence recordings of leukemic spleen explants from moribund mice revealed oscillations with a period of 24 hr" (Fig 3A); normal Per2::Luc spleen also 24 h (Fig 4A, B).
- Human AML lines are "highly sensitive to CLOCK and BMAL1 knockdown"; BMAL1 ChIP-seq 9,584 peaks, 140 regions shared with CLOCK, 229 genes, enriched for circadian genes (FDR 2.09E-8).
- Normal tissue spared: "Bmal1 (or Per2) in adult hematopoiesis is well tolerated" (Fig 4J-L, S5); germline Bmal1-/- lowers LSK and LT-HSC fractions, conditional Bmal1 loss does not.
- Pharmacology: SR9011 viability EC50 = 1.8 uM (printed as mM in the extracted text), and "murine leukemia cells are more sensitive than HSPCs to various concentrations of drug" (Fig 5C); human AML lines more sensitive than CD34+ cord blood (Fig 5D).
- Genomics: "In the 191 cases in the TCGA AML dataset, we did not identify mutations or homozygous deletions in CLOCK, BMAL1, or PER2" (Fig S5K).

## Not addressed / limits
- Dependency screen, not rhythm phase/amplitude in patients; no human tumour-vs-normal clock comparison, no time-of-day treatment.
- Rhythms shown only in murine explants and cell lines; human AML clock inferred from TF binding and sensitivity, not oscillation in patient cells.
- Haematopoietic malignancy; relevance to epithelial solid tumours unknown. Specificity control for pathway dependence is the normal HSPC arm, not a random-gene-set control.

## Link to us
- **Intro claim(s):** refs line 55 (Clock and Bmal1 required for AML/leukaemia stem cell growth in vitro and in vivo; clock disruption impairs proliferation and promotes differentiation; human AML lines also dependent): supported (human dependence is in cell lines only, not patients).
- **Relation to clocksync:** background. Says a working clock is a tumour vulnerability, so malignant compartments keep clock genes and clocksync can read them; the intact leukaemia clock is also a counter-example to that tumours lose the clock.

## Ideas
- Not a single-cell atlas: no compartments per patient here.
- Check whether BMAL1/CLOCK targets (140-region set, Table S3) could serve as an independent clock-output gene set to test specificity (pseudo-clock companion sets are done in A149).
