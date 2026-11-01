# A187 registration draft (superseded: registered as primary in registered_experiments.md, commit fab1005)

## A187 (draft 2026-11-01, written after the C4 scout (`scouts_2026-11-01/c4_single_cell/SCOUT_C4.md`: cell counts, depth, clock-gene detection) and the CoPhaser cell-cycle refits (`C4_single_cell/pilot`, `C4_single_cell/fits`: phases and QC only), before any clock-gene value is read across cell-cycle phase; the archive C2 result on overlapping cells is known: CoPhaser G2/M − non-cycling offset −0.02 h [−0.55, 0.51], 9 studies, `figure_5/output/cancer/decision.csv`): is the clock of malignant cells locked to their cell cycle? Readout level **L3** (within-patient phase offset between cell-cycle states, `clocksync.paired`), **secondary** (re-test of archive C2 with current CoPhaser fits)

**Why:** coupled-oscillator work predicts that the mammalian clock phase-locks to the cell cycle (Bieler 2014, Feillet
2014, Droin 2019). Locking would make the proliferation coherence loss (A97, A175, A186) a mechanism. The bulk test
(A179) found no linear cell-cycle composition effect; the archive single-cell test (C2) used older CoPhaser fits, two of
them degenerate (GSE183904, SRP327758: 100 % cycling, S→G2M order failed). Both are valid in the current fits.
**Data:** CNV-malignant cells of every study whose current CoPhaser fit is valid (validation plot and S→G2M order,
`fits/output/qc.csv`; fit chosen without clock genes, which are excluded from rhythmic, context and score sets), plus Wu
and Pal breast epithelial tumour cells. Studies whose kept seed is unstable (circular correlation with another usable
seed < 0.8 over cycling cells) are excluded.
**Bins (per patient):** NC = cycling probability < 0.5; among cycling cells, each cell goes to the nearer (circularly) of
the study's S-score peak and G2M-score peak along theta: S bin, G2/M bin. Clock-gene counts are summed per patient × bin;
depth = total UMI of the bin. A patient enters a contrast with ≥ 50 cells in each of its two bins.
**Estimator:** `clocksync.paired` (counts, GTEx waveform, unit = patient, cohort = study, reference = NC), offset =
arg(z_bin / z_NC) in hours; reference gate on the NC pseudobulks per study (`gated_only=True`).
**Power and stop rule (before the test):** MDE = 2.49 × SD of the per-patient label-swap null of the pooled offset. If MDE >
1.0 h the test is not run ("not estimable"). (Locking 1:1 to a cycle of about a day separates G2/M and G1 clocks by several
hours.)
**Test H187.1 (two-sided):** pooled G2/M − NC offset ≠ 0; per-patient label-swap null, 10,000.
**Specificity (same swaps):** 200 expression-matched pseudo-clock sets (`pseudo_means` from the malignant pseudobulks) and
the responsive companion. "Clock-specific" only if |offset| exceeds the pseudo 95th percentile of |offset|.
**Artefact checks (pre-specified, descriptive):** (a) without the short-lived clock mRNAs NR1D1 and PER1 (`drop_genes`,
as in A183; chosen a priori, not from the archive per-gene table), i.e. mitotic transcript loss; (b) the same contrast in cycling T/NK and myeloid cells where
fits allow (no locking expected to show there if it is a malignant-cell property; an equal offset means artefact);
(c) library-size-matched subsampling of the G2/M bin to the NC depth.
**Secondary H187.2 (runs only if MDE ≤ 1.0 h):** |pooled offset| in high minus low within-study tertile of patient cycling
fraction > 0 (tertile labels permuted within study, 10,000).
**Descriptive:** S − NC offset; r (log variance ratio) per contrast; per-study offsets; leave-one-study-out; breast by ER.
**Reading:**
- H187.1 passes clock-specifically and survives (a)-(c): the malignant clock is phase-locked to the cell cycle. That
  makes the proliferation coherence loss a candidate mechanism, and H187.2 tests the dose.
- A fail at MDE ≤ 1.0 h bounds any G2/M–G1 clock offset below about 1 h, so phase locking does not explain the
  coherence loss. This replicates archive C2 on corrected fits.
**Code:** `paper/exploratory/C4_single_cell/` (A187).
