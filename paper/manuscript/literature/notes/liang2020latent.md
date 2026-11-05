# liang2020latent
**Ref:** Liang S, Wang F, Han J, Chen K, Nat Commun 11:1441 (2020). DOI 10.1038/s41467-020-15295-9. Type: method (cell cycle). Lab: Ken Chen (MD Anderson).
**Source:** full text

## Question and design
Cyclum: sinusoidal autoencoder that finds a circular trajectory in scRNA-seq expression space to infer cell-cycle pseudotime and to remove the cell-cycle effect. Validation: 4 labelled cell-cycle datasets (mESC 288 cells, PC3 361, MDA-MB-231 342, H9 227; flow-sorted G0/G1, S, G2/M) and a virtual tumour; melanoma 4645 cells, 19 patients. **Sampling:** asynchronously cycling cells, no circadian data. Not a circadian paper.

## Method
- Unit: cell (one angle per cell). Unsupervised ab initio; no gene list needed; genes are standardised, no gene filtering.
- Model: encoder to a circular latent with sinusoidal activation, linear decoder; least-square error with L2 regularisation ("We use the least square error as the optimization target"). Squared-error (Gaussian) loss, no count model. Gene phases/loadings learned.
- Coverage of cycle: assumes cells "sampled from various stages of a periodic process"; no explicit check. Does not compare tissues or donors.

## Key results
- Cell-cycle classification accuracy beats Slalom, Cyclone, reCAT, Oscope, Cyclops and PC1/2 on all four datasets (about 0.7-0.75 median for the mESC subsampling; figure only, exact values not in text).
- Robust to number of cells and genes (no loss above 10,000 genes), while reCAT fails below 5000 genes.
- GSEA NES for cell-cycle genes 1.57 (Cyclum) vs 1.07/1.06 (PC1/PC2).
- Removing the cell-cycle effect improves separation of two virtual subclones.
- Circadian use is only future work: "We also plan to investigate the potential of applying Cyclum to characterize other periodic processes, such as circadian rhythms".

## Not addressed / limits
- No circadian validation; Tempo (Auerbach 2022) shows Cyclum is unstable on circadian scRNA-seq (median phase error 4.2-5.7 h).
- Limit of the circular assumption: "when the latent process does not fit the circular manifold well, the method may bring limited benefit."
- Only sinusoidal, shared across cells; no uncertainty.

## Link to us
- **Intro claim(s):** our refs line (circular latent, unsupervised, cells cover cycle) is right for the cell cycle; it should not be cited as a circadian phase tool, only as the origin of the circular-autoencoder idea used by CYCLOPS-style methods and Cyclum benchmarks in Tempo.
- **Relation to clocksync:** background / low relevance. It cannot estimate donor-level coupling or offsets; a cell-cycle confound in tumour scRNA-seq could, however, affect clock-gene pseudobulks.

## Ideas
- Cell-cycle effects as a confound for tumour compartments: regress a Cyclum or Seurat cell-cycle score out of pseudobulks before clocksync (A187 cohorts) and check that coupling does not change.
