# riba2022cell
**Ref:** Riba, Oravecz, Durik et al., Nat Commun 13:2865 (2022). DOI 10.1038/s41467-022-30545-8. Type: method + data. Lab: Nacho Molina (IGBMC).
**Source:** full text

## Question and design
Cell-cycle (not circadian) phase of single cells from scRNA-seq with RNA velocity. Three asynchronously cycling systems: mouse ES cells (>5000 cells, 2i+LIF), mouse pancreatic ductal progenitors (public), human fibroblasts (5367 cells, 3086 in the proliferating subpopulation). Chemistry: 10x. **Sampling:** no time axis at all; one snapshot of an unsynchronised population in which cells are spread over the whole cycle; orthogonal check on FACS-sorted G1/S/G2M bulk.

## Method (if any)
DeepCycle: autoencoder with one circular latent variable (theta, via cos/sin in the first decoder layer) fitted to the spliced and unspliced counts of GOterm:cell_cycle genes (n = 532). Unsupervised, no gene acrophases supplied; assumes cells cover the whole cycle and cycling genes draw closed loops in the unspliced-spliced plane. Output: theta per cell (monotone but arbitrary-speed function of time), later aligned to phases with marker genes.

## Key results
- Depth: median UMI 31977 (mESC) and 27319 (fibroblasts), 8043 in ductal cells. "while the ductal cells are as low as 8 thousand UMIs per cell, the median value of 8043" (Results, Fig. 1D).
- Theta is not time: "so we cannot directly infer temporal dynamics with it" (Results, Fig. 2).
- Internal validity: RNA counts per cell double across theta (fold change 2.1, 2.0, 2.1). "show a positive fold change of 2.1, 2.0, and 2.1 for mESCs, ductal cells" (Fig. 3C). mESC phase fractions 19% G1, 40% S, 41% G2/M vs 22-26/42-51/27-32% by propidium iodide.
- Bulk RNA-seq of sorted cells agrees with a shift: "The cells from the scRNA-seq sorted by transcriptional phase are close to the correct phase with a slight shift." (Results).
- Beats Cyclum, which gave inconsistent circles over repeated runs; agrees with Revelio, whereas Cyclone failed in mESC (Supp. Figs S7-S9).
- Finds a G1 branch point where a subpopulation leaves the cycle (quiescence-like).

## Not addressed / limits
- Needs deep libraries and genes known to cycle: "To train DeepCycle, we used the expression of unspliced and spliced RNAs of the genes in the GOterm:cell_cycle" (Results). Not tested at 10x depths of tissue atlases (only one low-depth set, 8 k UMI).
- scVelo failed on these genes, so the method rests on its own flexible fit; no ground-truth timed data, only FACS and global features.
- Phase boundaries still need markers: the authors say "In the future, the usage of gene markers needs to be replaced by the adoption of methods relying on dynamical features of gene expression" (Discussion).
- Cell cycle only; no circadian clock, no donors, no disease contrasts.

## Link to us
- **Intro claim(s):** the unsupervised RNA-velocity claim is supported. Corrected detail: it is unsupervised for the phase but supervised by a gene list (GO cell_cycle). Supports the technical-challenges contrast: cell-cycle phase needs no timing because cells are asynchronous, whereas circadian phase is shared by all cells of a donor.
- **Relation to clocksync:** background (cell-cycle confounder; comparator for CoPhaser cell-cycle phase, A187).

## Ideas
- Compare DeepCycle with the CoPhaser cell-cycle phase on the A187 studies (done in part: A187 fitted CoPhaser in 18 studies; the DeepCycle comparison itself is not registered).
- Use the RNA-count doubling across theta as a model-free sanity check of any cell-cycle phase on our scRNA-seq (fold change near 2).
