# veltri2025singlecell
**Ref:** Veltri A, Nukaya M, Korac K, Schwartz PB, Ronnekleiv-Kelly SM, bioRxiv 2025 (10.1101/2025.03.24.644906). Type: data + biology. Lab: Sean Ronnekleiv-Kelly (Wisconsin).
**Source:** full text (preprint)

## Question and design
First scRNA-seq (Parse Evercode fixed cells, about 180k reads/cell) time course of mouse liver: 8-week-old male C57BL/6J, 6 time points every 4 h (ZT3, 7, 11, 15, 19, 23), 34,921 cells in 10 populations (4 hepatocyte subtypes, endothelial, Kupffer, fibroblasts, cholangiocytes, ...). **Sampling:** timed, full 24 h, one mouse per time point (n = 6 animals, no biological replicates), ~6000 cells per time point.

## Method (if any)
- Not a phase-inference method. Per cell population: pseudobulk per ZT; clock health by normalised CV (nCV, a relative amplitude proxy), the clock-correlation matrix (Spearman between clock genes vs a mouse atlas reference, z-statistic) and peak phases; rhythmic genes by 3 random pseudoreplicates x 10 trials followed by meta2d (q < 0.05, called in at least half of the trials).
- Compares compartments within the same animal at each time point (cell types of the same liver), but no coupling, amplitude ratio or offset is estimated, and phase alignment is judged visually, not by a tested offset.

## Key results
- Aligned clock: "aligned expression across cell populations (i.e., narrow window of peak phase of expression)"; anti-phase Bmal1 vs Nr1d2/Per3 in every population.
- nCV varies by cell type, e.g. "For instance, nCV ranged from 0.359 – 1.536 for the core clock genes in" fibroblasts vs 1.126-6.015 in periportal hepatocytes; the authors attribute this to normalisation and cell number.
- Peak phase of clock-driven genes agrees with bulk: Pearson r = 0.8619 (top 100 clock-driven genes) and r = 0.9700 (core clock genes only).
- Plain pseudobulk (n = 1 per ZT) finds almost no rhythmic genes: "evidence of an intact circadian clock in each of the cell clusters (Fig 2-4), none of the core clock genes" were rhythmic (q < 0.05) by pseudobulk even though clocks look intact.
- Cell-type-specific rhythmic genes (e.g. Insig2, Elovl3) found with pseudoreplicates.

## Not addressed / limits
- One mouse per time point, males only, entrained and healthy: alignment is shown at the level of peak phase with no hour-level offsets or uncertainty; fibroblasts (few cells) give no rhythmic genes even per cell.
- Pseudoreplicates are partitions of the same animal's cells, so q-values are anticonservative with respect to biological variance.
- The authors note that perturbations (mistimed feeding, misalignment) may shift cell populations differently: "environmental perturbations (e.g., misalignment or mis-timed feeding) on cell population specific clock" phases are untested.
- Not human, not disease.

## Link to us
- **Intro claim(s):** supports that, in a healthy entrained animal, clocks of cell populations inside one organ are aligned (peak phases within a narrow window), the baseline against which desynchrony in disease is defined. Our refs line is accurate; note that the evidence is qualitative and from n = 6 animals.
- **Relation to clocksync:** data we could use / supports (healthy within-organ baseline; each population compared against a shared clock). It offers the pseudobulk compartment structure we use.

## Ideas
- Re-analyse their per-population pseudobulks with clocksync (6 time points as 6 "donors") to test whether coupling is near 1 and offsets near 0 in the healthy liver (a positive control for the single-cell side; ask for the GEO data).
- Their nCV and clock-correlation z-statistic (Shilts 2018) are a model-free gate analogue to cite.
