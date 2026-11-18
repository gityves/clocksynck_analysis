# marcheva2010disruption
**Ref:** Marcheva et al., Nature 2010 (466:627-631). DOI 10.1038/nature09253. Type: biology (global and pancreas-specific clock mutants, islet physiology). Lab: Bass / Takahashi (Northwestern, UT Southwestern).
**Source:** full text (main text only; supplementary figures and tables not in the text read)

## Question and design
Do pancreatic islets carry a self-sustained clock, and does it control insulin secretion? Islets from Per2Luc reporter mice (live bioluminescence of single islets), global ClockD19 mutants and global Bmal1-/- (islet and whole-body phenotypes), and a pancreas-specific knockout PdxCre Bmal1flx/flx (Cre under the PDX1 promoter; controls Bmal1 flx/flx and PdxCre). Islet and liver qPCR, islet microarray in ClockD19 mice, glucose and insulin tolerance, glucose-stimulated insulin secretion (GSIS) of isolated islets. **Sampling:** timed: 24 h qPCR in WT versus ClockD19 islets (n = 4 mice per genotype per time point); for PdxCre Bmal1flx/flx islet and liver, "three sequential 8-h time points" (n = 4 per genotype per time, Fig 4c); ex vivo islet bioluminescence over 4-5 days; mice 2-4 months (knockout) or 8 months (ClockD19).

## Method (if any)
None. Real-time bioluminescence (Lumicycle period and damping), qPCR, microarray with GO enrichment.

## Key results
- **Islet autonomy (global mutant).** "Here we show that pancreatic islets possess self-sustained circadian gene and protein oscillations of the transcription factors CLOCK and BMAL1." Single-islet PER2-LUC rhythm with period 23.58 +/- 0.3 h, damping after about 3 days, re-initiated by forskolin; ClockD19 islets arrhythmic, with Per2 RNA "reduced and rhythmicity abolished" (qPCR).
- **Islet versus liver differ in clock-controlled genes.** The temporal pattern alterations in mutant islets "were distinct from those in" mutant liver, "reflecting partitioning of metabolic functions within these two tissues at different times of day".
- **Driver and targeting.** PdxCre: BMAL1 immunofluorescence lost "specifically within pancreatic islets (Fig. 4a) and not within brain regions such as the suprachiasmatic nucleus" (arcuate, DMH, PVN in Supplementary Fig 8). Recombination efficiency in islet cells is not given as a number in the text read; recombination in non-endocrine pancreas, in duodenum or in other PDX1-lineage cells is not discussed, although the authors call the line pancreas-specific.
- **Targeted cells' clock.** Islet qPCR of Bmal1, Rev-erb-alpha and S100a6 at three time points shows altered profiles in PdxCre Bmal1flx/flx islets; liver of the same mice shows no such change: "revealed islet-specific alterations in gene expression profiles".
- **Other tissues and behaviour.** Liver clock genes of the same mice are the control (as above). Knockout mice "showed normal circadian activity, feeding rhythms and body weight and composition"; the phenotype "must be due to disruption of the clock network within the islet" in the authors' words.
- **Non-targeted cells in the same tissue.** Not measured: islets are read as a unit (no alpha-cell, delta-cell or exocrine clock readout, no whole-pancreas time series).
- Physiology: knockouts show hyperglycaemia, impaired glucose tolerance, hypoinsulinaemia and reduced insulin secretion, a block at exocytosis (GSIS and KCl, exendin 4, forskolin, 8-bromo-cAMP); more severe than the global mutants, "consistent with compensation occurring in the global mutant mice".

## Not addressed / limits
- Islet-level readouts only (isolated islets, a mixture of endocrine cell types); the beta cell is inferred from the phenotype, not isolated.
- Only three time points in the pancreas-specific knockout: no rhythm, phase or amplitude estimate for any non-targeted cell or for the tissue.
- Ex vivo single-islet bioluminescence shows that islet cells sustain rhythms without the rest of the body but not whether the cells within an islet are coupled.
- No clock-gene readout for tissues other than liver in the knockout.

## Link to us
- **Intro claim(s):** an isolated organ module can oscillate in vitro ("self-sustained"); tissue-specific clock loss gives a metabolic phenotype with intact behaviour and an unaffected liver clock. For within-tissue autonomy, coupling or pacing: **not informative** (no non-targeted-cell or whole-tissue rhythm after the knockout).
- **Relation to clocksync:** background. Islet data are not a part of our human organ set (GTEx pancreas is mostly exocrine).

## Ideas
- Islet microarray of ClockD19 versus WT (Supplementary Tables 1-2 in the paper) is the only large dataset here; no GEO accession is given in the main text read.
- Compare with Perelis 2015 (same lab, whole-islet RNA-seq WT time series) for the islet cycling transcriptome background.
