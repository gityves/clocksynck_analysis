# bray2008disruption
**Ref:** Bray et al., Am J Physiol Heart Circ Physiol 2008 (294:H1036-H1047). DOI 10.1152/ajpheart.01291.2007. Type: biology (cell-type-specific dominant-negative clock mutant, cardiac physiology). Lab: Young (Baylor / USDA Children's Nutrition Research Center).
**Source:** full text

## Question and design
Does the cardiomyocyte clock shape cardiac function, metabolism and gene expression? CCM mice (FVB/N) express the dominant-negative CLOCK-delta19 protein under the MHC-alpha promoter in cardiomyocytes; littermate wild types as controls (no Cre-lox: a transgene, not a deletion). In vivo echocardiography (ZT6) and telemetry, ex vivo working-heart perfusion at ZT0, 6, 12, 18, isolated mitochondria, and microarray of atria and ventricles. **Sampling:** timed, 12:12 LD (diurnal, not free-running), 8-week-old males for microarray, 11-14 weeks for physiology; "Atrial and ventricular tissue samples were collected from wild-type and CCM animals every 3 h for a period of 24 h" (n = 4 per time point, 128 samples).

## Method (if any)
Cosine fit (nonlinear least squares, Straume/COSOPT-style permutation test; MMC-beta < 0.05, Fisher G < 0.1, period 20-28 h), two-way ANOVA genotype x time to separate mean shifts from phase or amplitude changes.

## Key results
- **Driver and targeting.** Cardiomyocyte-specific expression of the dominant-negative CLOCK under the MHC-alpha promoter (generated and characterised in an earlier paper, their ref 9). Expression levels, fraction of cardiomyocytes expressing it, and expression in non-myocytes are not given in this paper.
- **Targeted cells' clock.** From the earlier characterisation: the heart clock "is severely attenuated" and clock component and output gene oscillations are attenuated "without affecting circadian clocks within extracardiac tissues" (as cited, not shown here). Here, clock genes (bmal1, per3, dec1, e4bp4, dbp) are confirmed by RT-PCR as altered in whole atria and ventricles (Supplementary Fig 5, not in the text read).
- **Whole tissue, non-targeted cells.** Total RNA is from whole atria and ventricles. The authors state "persistence of attenuated circadian clock gene oscillations in CCM hearts not only is expected, but is essential if one is to state that this model is truly cardiomyocyte specific", reasoning that atria and ventricles "are composed of multiple cell types (cardiomyocytes, vascular smooth muscle cells, endothelial cells, fibroblasts), all of which possess cell autonomous clocks". This is an interpretation: non-myocyte clocks were not measured and the residual rhythm is not attributed to any cell type by data. They also note "oscillation of multiple clock components was only attenuated, as opposed to completely abolished".
- **Transcriptome.** 728 genes oscillate in wild-type atria and 296 in ventricles; "548 and 176 exhibit significantly altered circadian oscillations in CCM mouse atria and ventricles", 45 common. "the majority of which are disrupted at the level of amplitude (as opposed to phase; Supplemental Table 3)". Genes peaking at the dark-to-light transition tend to be constitutively induced, those at light-to-dark constitutively repressed in CCM. Hence about 75% (548/728) of atrial and 59% (176/296) of ventricular rhythmic genes (derived) show significantly altered oscillation; the remainder are not significantly altered.
- **Other tissues and behaviour.** Physical activity diurnal rhythm unchanged ("comparable diurnal variations in physical activity"); blood pressure unchanged; heart rate lowered and its diurnal variation attenuated, persisting ex vivo, so the phenotype is intrinsic to the heart. Other tissues' clocks are not measured here (cited as unaffected from ref 9).
- Physiology: wild-type hearts show diurnal responsiveness to workload (cardiac power 1.7-fold higher at ZT18 vs ZT12) absent in CCM; increased oxygen consumption and fatty acid oxidation.

## Not addressed / limits
- No separation of cardiomyocyte from non-myocyte gene expression (whole atria and ventricles); no sorted-cell or single-cell data.
- Dominant-negative CLOCK may act outside CLOCK/BMAL1 targets and is constitutive (secondary adaptation); the authors say so.
- Diurnal design in LD: cannot separate autonomous rhythms from light, feeding or hormone-driven rhythms in the non-targeted cells.
- GEO deposit is mentioned ("Normalized data have been submitted to GEO archive") without accession in the text read.

## Link to us
- **Intro claim(s):** closest design in this set to our question (a tissue-wide timed series with one cell type's clock disabled), but the cell-type attribution of what remains is the authors' assumption: **suggestive of autonomy of non-targeted cell clocks in the heart (stated, not tested), not informative on coupling or pacing**; amplitude rather than phase loss in the majority of altered genes is compatible with the non-targeted cells keeping the tissue phase, consistent with our finding that organ phase is shared within about 1 h without a fixed leader (our inference, not the authors').
- **Relation to clocksync:** data we could use (if the GEO series can be located): 128 samples, wild type and CCM, atria and ventricles, 3 h resolution, to quantify residual clock-gene amplitude and phase versus wild type, a direct estimate of the non-cardiomyocyte share of the tissue clock signal. Human heart single-nucleus atlases could be asked whether non-myocyte (stromal/endothelial) cells hold more of the shared state.

## Ideas
- Locate the GEO series (search author/title; not given here) and score CCM clock-gene phase and amplitude relative to wild type per tissue (atria, ventricles).
- Compare the fraction of altered rhythmic genes (75% atria, 59% ventricles, derived) with the fraction of Bmal1-dependent rhythmic genes lost in a cell-type-specific knockout in other tissues (Koronowski liver: 89% lost; Dyar muscle: derived 71% TA, 49% SOL).
