# A183 registration DRAFT (not committed to registered_experiments.md; for the PI, Q22)

## A183 (draft 2026-11-01, written after the GSE154200 prep, gate and label-permutation nulls (`prep_gse154200/output/{gates.csv,power.csv}`) and before any lesional-minus-non-lesional or between-diagnosis value is computed; the archive psoriasis result is seen (A54-A59); HANDOFF_2026-10-31 H1): does the psoriatic lesion delay replicate blind in a new single-centre cohort, and is it specific to psoriasis among inflammatory skin diseases?

**Why:** A54 replicated a lesional clock delay in psoriasis on held-out GEO series (0.93 h, p 0.004; discovery 1.50 h);
A59 found no psoriasis-sized delay in other lesional skin diseases (191 pairs, 7 conditions, pooled −0.39 h, p 0.064).
Both used multi-series pools. GSE154200 is one centre and one pipeline, with psoriasis, nummular eczema and atopic eczema
sampled the same way (lesional and non-lesional punch biopsy at one sitting), so diagnosis is not confounded with series.
Readout level **L3** (within-person phase offset, `clocksync.paired`); d* is not estimable with two compartments.
**Data:** `$SCCC_DATA/skin_gse154200/GSE154200_countMat.txt.gz`, persons with both biopsies in one batch: psoriasis 87,
nummular 52, atopic 46 (one psoriasis person with the two biopsies in different batches dropped). Gate passes in each
diagnosis and pooled (p ≤ 0.002, `gates.csv`). All 14 clock genes present.
**Estimator:** `clocksync.paired` (counts, GTEx waveform, one projector over all persons, cohort = batch as in A54), offset
= arg(z_lesional / z_nonlesional) in hours, negative = lesion delayed.
**Tests (one-sided, Holm over the two):**
- **H183.1:** psoriasis lesional offset < 0 (per-person label-swap null, 10,000). MDE 0.81 h (`power.csv`); the A54
  held-out value is −0.93 h.
- **H183.2:** psoriasis offset − mean offset of the two eczemas < 0 (diagnosis labels permuted within batch, 10,000).
  MDE 0.69 h. Prediction from A59: eczema offsets near 0.
**Specificity companion (both tests, same permutations):** 200 expression-matched pseudo-clock sets (`pseudo_means` from
the genome-wide count matrix) and the responsive-gene companion (`pseudo_responsive=True`). "Clock-specific" only if the
observed offset is beyond the pseudo q05 (H183.1) or the pseudo-difference q05 (H183.2).
**Descriptive (no Holm):** offsets of nummular and atopic eczema with bootstrap intervals; coherence and r (log variance
ratio) per diagnosis; leave-one-gene-out offsets; each test without PER1 and NR1D1; offset by batch.
**Reading:**
- H183.1 passes clock-specifically → the psoriatic delay replicates in a third, blind, single-centre setting.
- H183.2 passes → the delay is a psoriasis property, not a property of eczematous inflammation in the same pipeline.
- A fail at MDE 0.81 h bounds the psoriatic delay below about 0.8 h in this cohort; it does not retract A54 (held-out).
**Code:** `paper/exploratory/A183_skin_gse154200/` (from `scouts_2026-11-01/prep_gse154200/scripts/prep_gse154200.py`).
