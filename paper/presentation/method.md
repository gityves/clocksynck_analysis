# Method section: slide by slide

Order (agreed 2026-11-09): M1 projection, M2 model (the video), M3 gate, M4 fit (host by EM on anchors, then the
readout of every compartment), M5 simulation. The pseudo-clock companion (random expression-matched gene sets) is a
check on results, shown after the readout, not part of the gate. Code references are to `clocksync/src/clocksync/`.

---

## M1. Projection: each pseudobulk becomes one point in the clock plane

**Message:** the 14 clock genes of one donor × compartment are summarised by one complex number β, read on a fixed
reference waveform.

**Bullets**
- Pseudobulk: one per donor × compartment; summed counts of 14 clock genes and the genome-wide depth.
- Log fractions y_g = log(count_g / depth), centred over genes: a common shift (depth, capture, RNA content) cancels.
- Weighted least squares on the reference plane X_g = amp_g (cos ψ_g, sin ψ_g): one β = λ e^{iφ} per pseudobulk, with
  its Poisson covariance V.
- Linear, and the same for every pseudobulk of a compartment: a constant offset (platform, cell-type baseline) moves
  the whole compartment by a constant, which the model absorbs (μ_k, next slide).
- Concept: gene acrophases fixed and shared, amplitude and means free per context, as in CoPhaser
  \cite{paychere2025cophaser}.

**Visual:** `figures/projection/projection.png` (script `scripts/projection.py`, simulated).
- A: the GTEx reference waveform, 14 arrows (angle acrophase, length relative amplitude).
- B: one fibroblast pseudobulk, gene log fractions minus the compartment mean (bars) and the plane fit (dots).
- C: the clock plane, three compartments × 30 donors sampled 08-18 h, 2 SD Poisson ellipses; each cloud is an arc of its
  cycle; a constant per-gene offset moves the hepatocyte cloud rigidly.

**Speaker notes**
- The waveform was fitted once on GTEx with time of death: shared acrophases, a mesor and an amplitude per tissue
  (`waveform.fit_waveform`, `data/human_clock_gtex.npz`). Core genes: NR1D1 3.6 h, DBP 8.5 h, PER2 11.7 h, CRY1 17.9 h,
  ARNTL 21.0 h; largest relative amplitudes ARNTL 0.50, NPAS2 0.33, NR1D1 0.32 (`cs.gtex().table()`). Only the
  relative acrophases matter: the zero cancels with the offsets.
- Why fixed: gene loadings are not identified from between-donor covariance when sampling does not cover the day
  (METHOD.md §1). Refitting gene timing failed in simulation (A139: V139.1 and V139.2 FAIL), so the conserved waveform
  is stated as an assumption.
- The choice of waveform does not drive the GTEx organ map: A123 uses a cross-fitted waveform; with the all-donor
  waveform, offsets 0.945 and amplitudes 0.970 (unchanged).
- The angle of one β is not a time: the compartment's centre is unknown and, in hospital hours, is not the centre of
  the cycle (panel C: each cloud is an 8-10 h arc). This is why a model is needed (M2).

**Limits**
- 14 genes, one human waveform (GTEx global; per-tissue waveforms exist, `cs.gtex("Lung")`).
- The gate tests for a clock, not for the waveform: loops shifted unequally pass it and return an offset that is not
  the average gene shift, with lower |z| and c (A205, V205.1 PASS, V205.2 FAIL).

<!-- src: waveform: clocksync.gtex().table() (data/human_clock_gtex.npz); A123 sensitivity: results_preliminary.md
A123; A139: paper/exploratory/A139_waveform_refit/output/{v139,d139_3,gene_level,failures}.csv; A205:
paper/exploratory/A205_partial_shift_gate/output/{tests,cells,per_rep}.csv; projection: _core/plane.py -->

---

## M2. The model (the video)

`figures/model_animation/` (`scripts/model_animation.py`): β_ik = μ_k + z_k f_i + e_ik, built term by term. To write.
Add on the last still or on M4: c_k = |z_k|² V / (|z_k|² V + 2σ_k²), with V the host spread.

## M3. Gate: is a shared, clock-shaped state observable?

Two nulls on the cross-covariances of compartment pairs projected on the clock plane: shuffled donor pairing and
permuted acrophases (`Dataset.gate`). A failed gate is uninformative. To write.

## M4. Fit: host by EM on anchors, then every compartment read against it

EM (SQUAREM) on the anchor compartments gives each donor's host state; each compartment is then regressed on it (z,
σ², c, |z|, arg z; anchors on their leave-one-out host; donor bootstrap with the host refitted). To write.

## M5. Simulation: known parameters recovered

To write. Then the pseudo-clock companion (A150 as the example) before the results.
