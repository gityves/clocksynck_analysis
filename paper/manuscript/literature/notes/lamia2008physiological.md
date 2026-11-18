# lamia2008physiological
**Ref:** Lamia, Storch, Weitz, PNAS 2008 (105:15172-15177). DOI 10.1073/pnas.0806717105. Type: biology (tissue-specific knockout, physiology). Lab: Weitz (Harvard).
**Source:** full text

## Question and design
Does a peripheral clock matter physiologically? Mice with hepatocyte-specific Bmal1 deletion (Bmal1 conditional allele x albumin-Cre, L-Bmal1-/-; paired male littermates, C57BL/6 x 129, 8-12 weeks old for the global knockout comparison) versus albumin-Cre controls; a global Bmal1-/- and a Per1/Per2 double mutant are used as comparators. Liver and skeletal muscle qPCR of clock genes over one day (Fig 2D, circadian time axis with ticks at 0, 8, 16, 24 h, i.e. about three time points per day; replicate number per time point not stated in the text read; mice entrained to 12:12 LD for 2 weeks, free-running versus LD sampling not stated), liver metabolic genes, blood glucose over 24 h, glucose and insulin tolerance, behaviour. **Sampling:** timed, one 24 h series per genotype, whole liver (no cell fractionation).

## Method (if any)
None (no inference from untimed data). Conventional qPCR (normalised to Hprt), Southern blot of the recombined allele, Western blots.

## Key results
- **Driver and targeting.** Albumin-Cre "acts exclusively in hepatocytes". Southern blot shows the disrupted allele in liver but not in kidney or skeletal muscle (Fig 2B). The residual non-recombined allele in liver "was likely derived from nonhepatocyte cells, such as bile duct cells and vascular epithelial cells" (inference by the authors; recombination efficiency in hepatocytes not quantified).
- **Targeted cells' clock.** Whole-liver Bmal1 expression is severely reduced across the cycle (Fig 2D); Dbp and Rev-erb-alpha show a nearly complete loss of expression; Cry1 is up-regulated, as expected without Bmal1 function. BMAL1 protein lost in liver with a trace residual attributed to non-hepatocyte cells.
- **Whole tissue, non-targeted cells.** Bulk liver (hepatocytes plus non-parenchymal cells) is the only readout, so the contribution of non-hepatocytes is not separated. The one explicit remaining rhythm: "Per2 transcript and protein showed persistent circadian regulation in the livers" (data not shown in the paper); the authors read it as "confirming that Per2 rhythms in the liver can be driven by external signals in the absence of intrinsic clock function". They do not test whether non-hepatocytes carry it.
- Most hepatic metabolic genes lose rhythm: "Most, but not all (e.g., Pepck1), lost circadian regulation" (Fig 3); Glut2 transcript and protein fixed at about the normal daily minimum.
- **Other tissues and behaviour.** Skeletal muscle of the same animals "showed a normal rhythm of expression" of Bmal1 (Fig 2D), Dbp and Rev-erb-alpha also normal in muscle. Locomotor and feeding rhythms "did not differ from controls (Fig. S2), demonstrating intact SCN function"; corticosterone and glucagon at two times of day not different (Table S1). "Histopathological examination" found "no apparent differences in tissue architecture between genotypes" (data not shown).
- Physiology: "exhibited hypoglycemia restricted to the fasting phase of the daily feeding cycle" (abstract), exaggerated glucose clearance.

## Not addressed / limits
- Whole-liver readout only: no non-hepatocyte clock measurement, no phase or amplitude of any non-targeted liver cell, no time series with enough resolution to report phase shifts of the residual Per2 rhythm (data not shown).
- Coarse sampling (about three time points per day from the axis ticks); no free-running design stated for the clock-gene series.
- Cre efficiency in hepatocytes not given as a number; the argument for specificity rests on a Southern blot and Western blot.
- Primary result is physiological (glucose homeostasis), not clock organisation.

## Link to us
- **Intro claim(s):** tissue-specific loss of a clock in one compartment leaves behaviour and other tissues intact (muscle clock normal, SCN output intact): supports. On within-tissue organisation: **not informative** (hepatocyte-targeted, bulk liver readout, non-hepatocytes not measured). The near-complete loss of Dbp and Rev-erb-alpha in bulk liver suggests that bulk clock-gene signal is dominated by hepatocytes or that non-parenchymal rhythms are small (our inference, not tested in the paper); in our atlases, endothelial and stromal cells carry the shared state, so a hepatocyte-only bulk readout is a poor test of their contribution.
- **Relation to clocksync:** background. A hepatocyte Bmal1 knockout is the case we scored as A111/A129 (Manella 2021 line); this paper has no data matrices we could reuse (qPCR panels only).

## Ideas
- If a liver-Bmal1 KO set with single-cell or sorted non-parenchymal cells existed, the shared host state (non-hepatocyte cells) could be read against the KO hepatocyte offset: prediction from coupling, residual non-hepatocyte rhythm keeps phase; from pacing by hepatocytes, it shifts.
- Per2 persistence in whole liver is a candidate for a component driven by systemic or non-targeted cells, to be compared with host-referenced estimates.
