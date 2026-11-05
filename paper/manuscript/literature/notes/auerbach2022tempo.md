# auerbach2022tempo
**Ref:** Auerbach, FitzGerald, Li, Nat Commun 13:6580 (2022). DOI 10.1038/s41467-022-34185-w. Type: method (+ new mouse aorta scRNA-seq). Lab: Mingyao Li / Garret FitzGerald (Penn).
**Source:** full text

## Question and design
Unsupervised single-cell circadian phase inference for droplet scRNA-seq. Simulations (500-5000 cells, 3000-20,000 UMI) from unsynchronized cells or a 4-timepoint time course; real data: mouse aorta (10X, ZT0/6/12/18; 18,863 SMCs, 3135 fibroblasts, 288 endothelial, 287 macrophages) and mouse liver hepatocytes (Droin, 18,378 cells). **Sampling:** timed light-dark time courses, 4 points over 24 h (aorta); the timing is used only for evaluation. Simulated unsynchronized cultures draw phases from "23 equally spaced phases over the circadian cycle".

## Method
- Unit: single cell (one phase per cell, with a posterior). Unsupervised w.r.t. time, but needs a user-supplied list of core clock genes and prior on their acrophases (simulated prior error SD 2 h); cell phase prior non-informative (von Mises).
- Gene acrophases: learned (variational posterior, shared by all cells) with prior; de novo cyclers added iteratively. Noise: Negative Binomial counts, 24 h sinusoid mean. Coverage of the cycle: no explicit requirement is stated; all evaluations use cells spanning 24 h.
- Does not compare tissues/cell types within an individual; each cell type is run separately, and the number of animals per time point is not discussed in the main text.

## Key results
- Accuracy vs sample time (median phase error, h, Table 1): Tempo 1.88 (SMC), 2.53 (fibroblast), 3.39 (endothelial), 2.99 (macrophage), 3.78 (hepatocyte); Cyclops 4.6-5.7, Cyclum 4.2-5.7, PCA 3.4-6.0.
- Uncertainty calibrated: "Tempo's uncertainty quantifications are well-calibrated for both unsynchronized and time course data".
- Stability: "Tempo's median circular standard deviation was less than 1 h for all simulation settings for which stability was evaluated".
- De novo cyclers: 189/109/87/28/117 called (SMC/fibroblast/endothelial/macrophage/hepatocyte) but "Tempo favored cell phase estimates based on the core clock genes alone rather than those including de novo cyclers for all real light-dark cycle cell type datasets".
- Cyclops/Cyclum unstable (circular SD overlaps uniform), a point we can cite against two-step orderers.

## Not addressed / limits
- Human data not analysed; mouse, healthy, entrained, with timed validation. Endothelial cells (weak clock) unstable: "27% of cells had circular standard deviations <1 h".
- Gene sinusoid, shared by all cells: "Tempo assumes the sinusoidal parameters are shared by all input cells" (limits cell-type or donor-specific amplitude/offset).
- No cross-compartment or cross-donor inference; no hospital-hours scenario.

## Link to us
- **Intro claim(s):** our refs entry (unsupervised, von Mises prior, needs clock gene prior) is supported; the earlier note cells span the cycle: not verified stays correct, it is implicit in the simulations, not stated.
- **Relation to clocksync:** background for single-cell phase inference; competitor only at cell level. Could not estimate coupling between compartments of one donor under hospital hours: it gives a phase per cell, but the donor-level offset between compartments needs a common reference that Tempo does not provide.

## Ideas
- Run Tempo on single-cell compartments of one donor and compare per-compartment mean phase offsets to clocksync arg z (A187 atlases).
- Cite Tempo's per-cell-type instability (weak clock) as a reason for compartment-level pseudobulk.
