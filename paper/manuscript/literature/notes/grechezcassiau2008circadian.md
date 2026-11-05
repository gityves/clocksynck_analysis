# grechezcassiau2008circadian
**Ref:** Grechez-Cassiau A, Rayet B, Guillaumond F, Teboul M, Delaunay F, J Biol Chem 2008;283(8):4535-4542. DOI 10.1074/jbc.M705576200. Type: biology. Lab: Delaunay.
**Source:** full text

## Question and design
Mouse (C57BL/6j, Bmal1 knock-out and WT littermates, 8-12 weeks), liver, heart, skeletal muscle; C2C12 myotubes (serum shock), NIH3T3 reporter assays, primary hepatocytes. Host (whole-body) clock perturbed by germline Bmal1 loss; no tumour model. **Sampling:** timed, full 24 h: livers at ZT0-20 every 4 h under LD12:12 (n = 3-6 per point) and at CT0-20 in constant darkness; C2C12 sampled every 4 h for 52 h.

## Method (if any)
qPCR (36B4 normalisation), luciferase promoter assays, [3H]thymidine incorporation, p21 siRNA; no inference method.

## Key results
- p21 is a clock output: "confirmed that p21 expression is controlled by an endogenous oscillator, with peak and trough values at constant darkness 0 and constant darkness 8" (Fig 1A, WT liver).
- Not liver-specific: "its rhythmic gene expression pattern was not restricted to the liver as the heart and skeletal muscle displayed a similar gene expression pattern" (Fig S1), unusual since few clock-controlled genes cycle in several tissues.
- In Bmal1-/- liver p21 mRNA is elevated 8-12-fold at ZT8 and rhythm is lost (Fig 1B).
- Control by nuclear receptors: p21 has two conserved ROREs; ROR-gamma activates and REV-ERB-alpha/beta repress; "the lack of cognate E-box element in the p21 promoter and first intron makes unlikely a direct inhibitory mechanism" by BMAL1.
- Phase relation: ROR-gamma peaks at ZT20 in WT liver, "hence 4 h before that of p21, consistent with an activating role" (Fig 4B).
- Function: "mature quiescent hepatocytes retain high proliferative potential"; Bmal1-/- hepatocyte proliferation is delayed by about 24 h and rescued by p21 siRNA (Fig 5, n = 3 experiments, t test).
- Wee1, Cdk6 and Ccne2 are not elevated in Bmal1-/- liver (Wee1 constant); p53 is not rhythmic.

## Not addressed / limits
- No cancer, no tumour clock, no host-tumour comparison; relevance to cancer is indirect (G1 progression as clock target in regenerating liver).
- Whole-body Bmal1-/- is arrhythmic with systemic phenotypes; hepatocyte-autonomous versus systemic effects are not separated in vivo.
- Hepatocytes studied in primary culture; n small, SD error bars, no multiple-testing treatment.

## Link to us
- **Intro claim(s):** the refs_*.md row (BMAL1-p21 via REV-ERB/ROR controls hepatocyte proliferation) is supported. It supports clock gating of proliferation gene programmes, not tumour desynchrony.
- **Relation to clocksync:** background. A caution for gene panels: p21 is high-amplitude and clock-controlled in several organs, so clock-gene coherence in proliferating tumours can mix clock and cell-cycle signal.

## Ideas
- Check proliferation confound in the existing tumour coherence results: done (A97, A127 in paper/docs/results_preliminary.md, SCAN-B); consider reporting CDKN1A as a documented cross-tissue clock-controlled gene.
- Missing reference: Matsuo 2003 (Science, Wee1) and Gery 2005 (Mol Cell, p21 and PER1 in tumours).
