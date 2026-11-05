# ye2018genomic
**Ref:** Ye et al., Cell Syst 6:314-328 (2018). DOI 10.1016/j.cels.2018.01.013. Type: data analysis (pan-cancer, bulk). Lab: Han (UTHealth Houston), Takahashi, Mills, Yoo.
**Source:** full text

## Question and design
Pan-cancer landscape of 14 core clock genes plus 37 clock-associated genes (51 "clock genes"): expression, mutation, methylation,
co-expression, survival and drug-response links. TCGA (32 cancer types; 14 types with paired tumour and normal RNA-seq),
CTRP and GDSC (265 drugs, 1,080 cell lines), plus reanalysis of a published serum-shock time series in MCF10A vs MCF7
(8 time points).
**Sampling:** human tumours single time point, time of biopsy unknown; only the cell-line series is timed.

## Method (if any)
Paired tumour vs normal differential expression (FC > 1.5, corrected p < 0.05); correlation of each clock gene with all
genes in tumour and in normal separately, corrected for tumour purity; JTK_CYCLE / meta2d on cell-line series.
Cohort-level, no sample-level clock state.

## Key results
- Expression: 46 of 51 clock genes (90.2%) changed in at least one cancer; "PER1 was downregulated in 11 cancer types, PER2 was downregulated in 7 cancer types, and PER3 was downregulated in 10 cancer types" (Fig. 1A); ARNTL2 up in 9 types; PER2 up in KIRC (FC = 1.8).
- Loss of clock co-expression: "Compared with normal tissue, significantly fewer genes in tumor tissues correlated with clock genes across all cancer types"; "correlation among core clock genes was reduced significantly" (CRY2 vs eight core genes in normal kidney, reduced or none in KIRP, Fig. S4B).
- Purity: "Pearson correlation coefficient was corrected by tumor purity" (only for the tumour arm).
- Cell lines: "We identified 478 genes with a clear circadian oscillation pattern in MCF10A by MetaCycle::meta2d" that are lost in MCF7; 552 genes oscillate in MCF7 only (reprogramming, not only loss).
- Mutations: breast cancer with core clock-gene mutations, worse survival (log rank p = 0.05); PER1/2/3 mutated: p = 5.9e-5.
- Survival and stage: 31 clock genes prognostic, 23 differ by stage, 43 by subtype; clock genes correlate with 145 of 154 actionable genes and with sensitivity to 265 drugs.

## Not addressed / limits
- "There are several limitations in our study. For example, the time points for obtaining biospecimens from cancer patients are not available" (Discussion): all expression differences are mean levels at unknown circadian time.
- Co-expression loss is correlational and could also arise from tumour heterogeneity, copy number, proliferation or composition; no matched non-clock gene control, and no gate that covariance is clock-shaped.
- No host-vs-tumour comparison within a person other than tumour-normal tissue pairs; no phase or amplitude estimates.

## Link to us
- **Intro claim(s):** pan-cancer clock gene dysregulation and loss of clock co-expression in tumours (Cancer interaction). Our refs row is SUPPORTED; the co-expression result is the closest precedent for a covariance-based readout.
- **Relation to clocksync:** background / precedent; the matched-gene control that this paper lacks is what A138.4 supplies.

## Ideas
- Done: TCGA tumour-normal pairing of clock state with expression-matched pseudo-clock sets (A138, H138.4: 517 pairs, 10 gated cancers, clock-specific, pseudo-clock q95 0.261 vs 0.365).
- Compare with the purity-corrected co-expression drop: compute the same Ye-style count of clock-correlated genes after matching on G2M/proliferation, to test whether the loss is clock-specific (not registered).
