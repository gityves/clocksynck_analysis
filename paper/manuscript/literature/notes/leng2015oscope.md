# leng2015oscope
**Ref:** Leng, Chu, Barry, Li, Choi, Li, Jiang, Stewart, Thomson, Kendziorski, Nat Methods 12:947-950 (2015). DOI 10.1038/nmeth.3549. Type: method. Lab: Christina Kendziorski (UW Madison).
**Source:** full text (Brief Communication; online Methods and supplement partly in the export)

## Question and design
Find oscillating genes and the cyclic order of cells in unsynchronised scRNA-seq, without time labels. Tests: a published bulk cell-cycle time series (Whitfield et al. 2002) with shuffled sample order; H1 hESC scRNA-seq (three replicates, n = 213; Fluidigm C1) and H1-FUCCI hESC (n = 247). Cell cycle, not circadian. **Sampling:** snapshot; "unsynchronized cells in a variety of different states are collected at the same calendar time T". No timed series for the scRNA-seq.

## Method
Unsupervised. (1) Fit a paired-sine model to all gene pairs and keep pairs with an elliptical scatter; (2) K-medoids into groups of similar frequency but possibly different phases; (3) per group, extended nearest-insertion (travelling-salesman style) recovers the cyclic order of cells within one base cycle, direction not fixed; phases then from standard time-course tools. Unit of inference: the cell (a group of genes defines its own order). No time scale, no external anchor. Authors state that "unlike previous approaches, ours does not attempt to construct a linear order on the basis of minimizing change among adjacent samples".

## Key results
- Shuffled bulk cell-cycle data: "The top group identified by Oscope had 151 genes, 116 of which were validated as oscillating in earlier work".
- hESC: "One of the top groups identified by the K-medoids algorithm in Oscope contained 29 genes", 21 annotated as cell cycle (GO:0007049); the reconstructed order recovers G1, S, G2/M; "classified 72% of H1-FUCCI hESCs into the correct phase".
- Technical artefact: a second group followed capture position on the C1 chip: "We found that 403 genes showed strong artificial trends" across the three runs, also in public data sets.

## Not addressed / limits
- Coverage of the cycle is implicit, not stated as a requirement: the method "capitalizes on the fact that cells from an unsynchronized population represent distinct states in a system". No explicit statement that cells must cover all phases, no test with partial coverage.
- Static data only give one cycle; no period, no clock time, no direction.
- No circadian data, no tissues, no individuals, no comparison of compartments within a donor. Not tested on a gene set as small as the clock genes.
- The C1 chip-position artefact shows that non-biological cyclic structure can be picked up as an oscillator.

## Link to us
- **Intro claim(s):** single-cell, unsupervised oscillation detection and ordering, assumption that the sample covers the cycle (technical challenges). Our refs entry (unsupervised; ordering assumes cells cover the cycle; no explicit statement verified) is correct: no explicit statement exists, the assumption follows from the design (cells in an unsynchronised population).
- **Relation to clocksync:** background. Different object (cell cycle, cells, within-sample phase) and a different target; relevant as the origin of cyclic-ordering ideas that CYCLOPS and CHIRAL inherit.

## Ideas
- Use the C1 capture-site artefact as a warning for the technical-covariate control of shared-covariance gates: a cyclic or planar structure in clock genes may be a batch trend. Our gate checks clock shape; whether it flags a plate-position trend has not been tested (no hit for Oscope in docs).
- Missing reference: Trapnell 2014 (Monocle), cited here as ref 4, is in the bib.
