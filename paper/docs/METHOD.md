# Method (as implemented, 2026-10-13)

The code as it stands.

## 1. Pseudobulks and the clock plane (`pseudobulk.py`, `clockplane.py`, `reference.py`)
- **Pseudobulks.** One per donor x compartment (cell class), with these columns:
  - the summed counts of the 14 clock genes of the GTEx waveform;
  - the summed genome-wide depth;
  - n_cells.
  - A stratum (study x assay, region) carries one intercept per compartment.
- **Reference waveform.** log E[y_g] = m_g + lam amp_g cos(phi - psi_g), fitted on GTEx with time of death
  (`data/human_clock_gtex.npz`). Its relative amplitudes and acrophases are taken as conserved.
  - It is used only to project and to gate.
  - It is not re-estimated from the data, because gene loadings are not identified from between-donor covariance
    when sampling does not cover the day.
- **Projection.** y_g = log(max(count_g, 0.5) / depth), centred over genes, projected by weighted least squares onto
  the plane X_g = amp_g (cos psi_g, sin psi_g):
  - the projection gives one complex coefficient beta per pseudobulk;
  - Poisson measurement covariance V = L diag(1/count) L^T;
  - a constant offset of a compartment shifts beta by a constant.

## 2. Gate (`Study.gate`, model-free)
Cross-covariances C_kl of the centred log fractions of compartment pairs across donors carry only shared variation.
- `clock_cov`: the projection of C_kl onto the clock plane, against shuffled donor pairing. Is there shared signal?
- `clock_share`: the same projection as a share, against permuted acrophases. Is the shared signal clock-shaped?
- A fit passes when both p < 0.05. A failing fit is uninformative, never "decoupled".
- Gate power for failures: A17 step 5 (`power.gate_power`).

## 3. Host state and coupling (`factor.py`, `anchored.py`)
- **Factor model.** beta_i = mu_g + M(z_k) f_s + e_i, with e_i ~ N(0, psi_k I + V_i) and f_s ~ N(0, Sigma_f).
  - Exact EM (SQUAREM) gives the donor host state f_s, with its posterior mean and covariance S_u.
- **Anchored host** (A17).
  - The host is fitted on anchor compartments only.
  - Anchor families are learned leave one study out from the model-free pair coherence:
    the |complex correlation| over shared donors, minus a shuffled-pairing null,
    pooled by DerSimonian-Laird.
  - Sensitivity: the fixed stromal-vascular rule, with neurons in CNS tissues.
  - An anchor's own host leaves it out.
- **Coupling** of compartment k: regression of its beta on the host posterior means, with intercepts per stratum.
  - Small-sample corrections: |z|^2 minus tr cov(z); psi raised by its degrees of freedom.
  - c_k = |z|^2 V / (|z|^2 V + 2 psi), where V = tr(cov f + mean S_u) is the host spread.
  - c is a variance share, so it depends on the host spread (sampling times). It is compared within a fit; across
    groups it is compared at a common spread (below).
- **Relative amplitude and offset:** |z_k| and arg z_k, reported reference-free (each minus the mean over
  compartments).
- **Intervals:** donor bootstrap with the host refitted.
- **Boundary** (ML psi at its floor): a lower bound by test inversion at the real design (`anchored.lower_bound`).
- **Identity diagnostics:**
  - free-vs-anchored host agreement and the bootstrap stability of the host, from the complex correlation of host
    states;
  - a fit with stability < 0.8 is flagged.

## 4. Group contrasts (`anchored.contrast`, `desync/common/common_v.py`)
- **Pooled host.** One anchored host over all labelled donors; per group, z, psi and V of k on that host.
  - d = c_case - c_control.
  - p: labels permuted within blocks, host fixed.
  - Interval: bootstrap within groups. Both groups at the boundary: non-estimable.
- **Common host spread (d\*).** c*_g = z2_g V0 / (z2_g V0 + 2 psi_g), with V0 the host variance over all labelled
  donors. It removes most of a sampling-window difference between groups (A18 addendum 3).
  - Also reported: the log case/control ratios of V, psi and z2.
- **Module tests:** the mean d over a family, against shared permutations.
- **Power:** MDE by binomial thinning of k's host-locked part in the case group (`power.mde`).

## 5. Covariates of coupling (`covariate.py`, A19)
- **Model.** beta = a_g + b x + M(z) f + M(w) x f + e, with e ~ N(0, (psi e^{gamma x} + |z + w x|^2 s_u) I + V),
  where s_u = tr S_u / 2.
  - Fitted by ML: GLS profiled over (log psi, gamma).
- **Estimands:**
  - locked amplitude Re(w/z);
  - phase drift Im(w/z) x 24 / 2pi, in h per unit x;
  - own noise gamma;
  - coupling slope.
- **Inference.**
  - Wald statistics against a null simulated at each compartment's design (no effect).
  - Bias correction: minus the median under that null.
  - Per-compartment design check: recovery of planted effects within 25 %.
- **Validation.** GTEx (timed): the untimed slopes rank tissues like the timed ones and are compressed about x0.4.

## 6. Between cells (`cells/vonmises.py`, `desync/common/cell_split.py`, A18 part 2)
**Scope in the paper (user, 2026-11-04, Q23):** only D, only in malignant cells (WP0 item 2; about 4.3 clock UMI per
cell, where D separates coherent from scrambled gene phases: 4.70 vs 0.42 in the calibrated simulation). The healthy-atlas
lambda and R (A18 part 2) are dropped: their per-cell files had ARNTL = 0 (integrity note 2026-11-04) and the readout was
inconclusive (A144 part 2). Resolution depends on clock UMI per cell; eligibility below stays as the minimum.
- **D.** Within-sample in-plane dispersion from fraction-scale cross-gene moments, per donor x compartment (>= 100
  cells), minus a permuted-acrophase null. It needs no time and no host.
- **Coherent amplitude.** rho^2 = c S / (1 - R_W^2), from the anchored coupling, the between-donor clock-plane
  variance S and the sampling-window synchrony R_W.
- **lambda and R.** The per-cell amplitude lambda and the synchrony R solve lambda R = rho with the von Mises
  dispersion model. They rest on Assumption W (V11).
- **Eligibility:** >= 20 groups, median clock UMI per cell >= 2, D above its null.

## 7. Validations (where to find them, all in the archived /home/maxine/Documents/paychere/scCircadianCancer)
| check | result | location |
|---|---|---|
| GTEx hidden times (V2) and benchmark against COFE, CYCLOPS 2.0, CHIRAL, TimeTeller | only sccc recovers between-tissue offsets untimed | `paper/code/fig1_method/` |
| Mouse colon (regenerating, day 17 after DSS), known ZT (V1); heart AM/PM (V3); spike-ins in biopsies (V5) | pass | `docs/RESULTS_preliminary.md` |
| A17 acceptance (identity, coverage, boundary, power) | pass, except coverage 0.88 at c = 0.9 | `outputs/v12_anchored/` |
| Cell-type coupling against timed truth (V13: regenerating mouse colon, skin; A21: undamaged vs regenerating colon) | anchored and corrected estimators track truth; the old free estimator overestimates | `outputs/v13_celltype_truth/` |
| ORPHEUS estimator under subject-phase error (V14) | bias only at low dispersion | `outputs/v14_orpheus_bias/` |
| Age model (A19 validations 1-3) | see the A19 note | `desync/20_age/output/` |
| Internal collection window (A25: V25.1-V25.4) and per-patient internal time against recorded resection time (V25.5, breast GSE233242, post hoc) | V25.2 passes, V25.3 fails; V25.5: patient order Spearman 0.56, median error 1.1 h with the hospital prior (curvature flip wrong in this 7 h, 2-compartment design) | `desync/common_output/a25_*`, `desync/validation/` |
| Untimed per-gene rhythms against GTEx (A27 part 1, post hoc): bulk, degraded towards atlas conditions, single-cell atlases | bulk: untimed >= timed oracle (atlas-like R 0.57 vs 0.40); single-cell atlases fail for non-clock genes (~12-15 h offset, cause open); clock-plane genes align | `desync/validation/gonogo/` |
| Ventilator cases and glucocorticoids (A27 part 2, post hoc); host without PER1 | internal clock decoupled from time of death (host and CHIRAL DIP agree); clock outputs keep their rhythm; sustained GC moves the core clock; surgical GC only induces PER1; A25 and V25.5 unchanged without PER1 | `desync/validation/glucocorticoid/` |
| Same-time variance ratio against timed amplitude ratio (V52, 572 GTEx tissue pairs) | Spearman 0.94 (0.91 at 08-18 h, 92 % sign agreement), attenuated (slope 0.64-0.68) | `desync/validation/paired_damping/` |
| Extended clock waveform (A28): 14 clock genes + tissue-stable clock targets, selection/evaluation split of GTEx | fails: single-sample error -2.5 % (CI includes 0), no harm; error is donor-level, not measurement noise; 14-gene waveform stands. D10 (exploratory): atlases predict the phase of internal-time follower genes (R 0.28-0.52) but not of other rhythmic genes | `desync/validation/extended_clock/` |
| Bulk clock from cell types (A29) | fails (Spearman 0.17); misalignment loss <= 0.13 | `desync/validation/bulk_clock/` |
| Atlas antiphase by gene class (A30) | fasting/feeding hypothesis rejected twice; antiphase is a concentrated axis | `desync/validation/antiphase/` |
| Per-cell-type gene phases vs true ZT, mouse colon (A31) | pass: R 0.98-0.99, all classes | `desync/validation/celltype_truth/` |
| New Census atlases and donor subsampling (A31b, A39) | followers do not replicate; lung HLCA (200 donors) predicts GTEx gene rhythms (R 0.56); needs ~200 donors | `desync/validation/gonogo_new/` |
| Blood composition control and blood's own clock (A32, A35) | GTEx blood uninformative; living blood clock weak (rho 0.31), ARNTL/RORA compositional | `desync/validation/blood_composition/`, `desync/validation/blood_clock/` |
| Ventilated donors by tissue (A33) | liver, pancreas, muscle decouple most; not a food-entrained class effect | `desync/validation/ventilator_tissues/` |
| Baboon multi-tissue (A36) | untimed ordering 0.4 h, tissue offsets 0.5 h; cross-species untestable | `desync/validation/baboon/` |
| Baboon in working-hours windows (A40) | 8 h windows: ordering perfect, direction 11/12; tissue offsets r 0.52 (1.2 h); centre needs species frame | `desync/validation/baboon/` |
| Species-specific vs human reference (A41) | baboon reference no better out of sample (offset r change -0.01; direction 12/12 vs 11/12); sheep results unchanged | `desync/validation/baboon/` |
| Cross-species messages 2 and 3 (A42) | brain couples less to the periphery in sheep and baboon; sheep brain host partly separate (r* 0.63); tissue ranks not conserved | `desync/validation/cross_species/` |
| Living human skin, epidermis vs dermis (A37) | untimed offset -1.6 vs timed -2.2 h; per-sample rho 0.505 | `desync/validation/skin_paired/` |
| MASLD hepatocyte coupling (A38) | 0.68 -> 0.54, n.s. (MDE 0.45) | `desync/validation/masld/` |
| Farm-animal GTEx (A34) | sheep gate passes (343 animals, 45 tissues); tissue phase order not conserved with human (r -0.19) | `desync/validation/farm_gtex/` |
| Gene-specific partial shifts in one compartment (A205; Saini 2013 pattern, loop genes +6 h, others +12 h) | gate passes (100 %); offset not bounded by the gene shifts (12.9-14.3 h); incoherence lowers c and abs(z): a lower amplitude does not separate a damped from an internally incoherent clock | `paper/exploratory/A205_partial_shift_gate/output/` |
| Shift common to the clock genes of a sample (2026-11-10; simulation, 200 data sets per scenario) | without a clock, a common shift of sd 0.2 passes the permutation test of the unsigned coherence R in 0.625 (n 10) and 1.000 (n 40) of the data sets; the signed coherence with a per-sample offset fitted out (`decompose.signed_coherence`) passes in 0.050 and 0.035 and detects a clock in 1.000 at every shift (R: 0.690 at n 10). Earlier L1 results use R (PI question Q31) | `paper/exploratory/scouts_2026-11-10/feasibility_fetal/` (`shift_artefact.csv`) |

## 8. Not identified, not implemented
- **Not identified without times:** absolute amplitude and phase, the host's own ageing, and per-gene waveforms.
- **Not implemented:** composition standardisation within a class.

## 9. Implementation details stated by the code (added 2026-10-17, from the STAR Methods check)
- Poisson covariance uses counts floored at 0.5 (as the projection), not raw 1/count.
- Host spread V is computed from host means centred within strata.
- Relative amplitude and offset come from the full anchored host; an anchor's coupling uses its leave-one-out host.
- Boundary: psi <= 1e-7; the test-inversion grid includes c = 1.0.
- Coherent amplitude uses the host-locked variance |z|^2 V from the anchored regression (equal to c S / (1 - R_W^2)
  under the model, not numerically).
- Per-cell dispersion null: 10 permuted-acrophase planes; "D above its null" is a point comparison, not a test.
- Covariate model: the package accepts >= 10 rows; the A19 rule of >= 30 donors is enforced by the drivers.
- Developmental series (A237-A239, `paper/common/development.py`): log fractions are residualised gene by gene within
  a compartment on an age polynomial, sex and assay before L1 or the gate; pseudo-clock sets receive the same
  adjustment; L1 uses the signed coherence; the gate is compared with pseudo sets on `Gate.pairing_z` (shared
  clock-plane covariance standardised by its pairing null).
- `sccc.api` defaults (200 gate permutations, 200 bootstraps, 40 MDE replicates) are lower than the paper settings
  (1000, 500, 100).
