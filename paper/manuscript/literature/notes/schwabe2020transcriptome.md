# schwabe2020transcriptome
**Ref:** Schwabe, Formichetti, Junker, Falcke, Rajewsky, Mol Syst Biol 16:e9946 (2020). DOI 10.15252/msb.20209946. Type: method + biology. Lab: Falcke (MDC) / Rajewsky / Junker.
**Source:** full text

## Question and design
Shape of the cell cycle in transcriptome space. Deep Drop-seq of asynchronous HeLaS3 (1477 cells, about 11,000 UMI per cell, 1031 highly variable genes); replicated in a second HeLaS3, HEK293 and 3T3. **Sampling:** none; an unsynchronised snapshot, assumed to cover the whole cycle. Cell-cycle phase is only a marker-score label (Whitfield genes).

## Method (if any)
Revelio: PCA, then a linear rotation to "dynamical components" DC1-DC2 in which cycling cells form an annulus (DC3 adds a helix from non-cycling expression). Order by the angle around the origin = cell-cycle pseudotime; removing the plane removes the cycle. Linear and gene-list free apart from the marker labels used to choose the rotation. RNA velocity is used only afterwards to confirm the direction of motion, not to infer the phase.

## Key results
- The cycle lives in a plane: "cells move along a planar circular trajectory in transcriptome space during the cycle" (Abstract); non-cycling expression adds a third axis.
- Reproduced in other lines: "We were able to reproduce the results with other cell lines, including an additional HeLaS3, a HEK293, and a 3T3 data set" (Results; annulus in HeLa, disc-like in 3T3).
- Depth rule: "at least 600 cells with an average depth of 4,500 UMIs or at least 1,000 cells" (Results, Appendix Fig S6); the limiting factor is depth/dropout.
- Pseudotime for free: "ordering the cells by their angle in a clockwise motion around the origin of the DC1-DC2 plane" gives temporal order.
- Each gene is up-regulated once per cycle; two gene groups drive it.

## Not addressed / limits
- The authors state the scope: "The analysis is limited to information on transcriptional regulation along the cell cycle."
- Immortalised proliferating lines only; tissue cells with a few cycling cells are not tested. Angle is a rank order, not time or duration.
- Statistical independence of components is not established (Methods, ICA remark). No circadian signal, no donors.

## Link to us
- **Intro claim(s):** our refs entry says cell-cycle ordering from unspliced/spliced counts. That is wrong for this paper: Revelio is a linear PCA rotation on spliced expression; velocity is a check only (the spliced/unspliced approach is riba2022cell).
- **Relation to clocksync:** background (cell-cycle covariate; why proliferation contaminates clock-gene pseudobulks, A179/A187).

## Ideas
- Regress the Revelio cell-cycle plane out of clock-gene counts before pseudobulking and refit coupling: related to A179 (bulk cell-cycle composition, resolved negative) and A187 (single-cell locking, not estimable); the cell-level removal is not registered.
