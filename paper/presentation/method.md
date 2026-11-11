# Method section: slide by slide

Order (agreed 2026-11-09): M1 projection, M2 model (the video), M3 gate, M4 fit (host by EM on anchors, then the
readout of every compartment), M5 simulation. The pseudo-clock companion (random expression-matched gene sets) is a
check on results, shown after the readout, not part of the gate. Code references are to `clocksync/src/clocksync/`.

---

## M1. Projection: each pseudobulk becomes one point in the clock plane

**Message:** the 14 clock genes of one donor × compartment are summarised by one point β, read on a fixed reference
waveform.

**Model of the counts** (donor i, compartment k, gene g; n counts, D genome-wide depth):

$$\log \mathbb{E}[X_{ikg}] = \log N_{ik} + s_{ik} + m_{kg} + \mathrm{amp}_g\,\lambda_{ik}\cos(\varphi_{ik} - \psi_g)$$

| term | what it is | status |
|---|---|---|
| $$N_{ik}$$ | Sequencing depth | **free**, one per pseudobulk |
| $$s_{ik}$$ | common shift of all genes (depth definition, RNA content, capture) | **free**, cancels (centring over genes) |
| $$m_{kg}$$ | gene baseline of the compartment (cell type, platform, any per-gene constant) | **free**, per compartment × stratum; becomes a constant shift of β, absorbed by μ_k |
| $$λ_{ik}, φ_{ik}$$ | clock amplitude and phase of this pseudobulk: $$β_{ik} = λ_{ik} (sin φ_{ik}, cos φ_{ik})$$ | **free**, one per pseudobulk |
| $$ψ_g, amp_g$$ | gene acrophase and relative amplitude | **fixed** (GTEx), assumed the same in every donor, compartment and condition |

Not allowed by the model: genes of one compartment or condition changing timing relative to each other (ψ_g differs)
or their amplitude ratios. Such a change still passes the gate and is read as an offset plus lower |z| and c (A205).

**Bullets** (in the order of the figures)
1. *Figure A, reference waveform.* 14 clock genes, each with an acrophase ψ_g and a relative amplitude amp_g, fitted
   once on GTEx and kept fixed: gene loadings are not identifiable from untimed, hospital-hours data. Fixed shared
   acrophases with free amplitude and means per context, as in CoPhaser \cite{paychere2025cophaser}.
2. *Figure B, one pseudobulk → one β.* Log fractions y_ikg = log(max(n_ikg, 0.5) / D_ik), centred over genes, fitted
   by weighted least squares on the plane X_g = amp_g (sin ψ_g, cos ψ_g):

   $$\hat\beta_{ik} = \arg\min_\beta \sum_g w_{kg}\,\big(\tilde y_{ikg} - \tilde X_g\beta\big)^2 = (\tilde X^\top W_k \tilde X)^{-1}\tilde X^\top W_k\,P\,y_{ik} \equiv L_k\,y_{ik}$$

   with P = I − 11ᵀ/G (centring over genes), X̃ = PX, W_k = diag(w_kg), w_kg = 1 / (1/μ_kg + α)
   (inverse variance of a log count under NB variance μ + αμ²; μ_kg the median count of gene g in compartment k,
   α = 0.05, `tau2` in the code), and Poisson covariance V_ik = L_k diag(1/max(n_ikg, 0.5)) L_kᵀ.
   In the figure, each gene's log fraction is taken relative to the compartment's mean over donors and divided by its
   fixed reference amplitude amp_g: the model then predicts one cosine λ cos(φ − ψ_g) for all genes. Its height λ and
   peak φ are β relative to the compartment's mean (by linearity, L_k (y_ik − ȳ_k) = β̂_ik − mean β̂_k).
3. *Figure C, the clock plane.* L_k 1 = 0 and L_k X = I, so β̂_ik ≈ λ_ik (sin φ_ik, cos φ_ik) + L_k m_k + noise: the common
   shift s_ik is gone and the gene baselines add one constant per compartment. A compartment's donors form an arc
   (hospital hours cover part of the cycle); a platform or cell-type offset moves the whole arc, not its shape.

**Visual:** three one-panel figures, shown in order (`scripts/projection.py`, simulated, 150 × 120 mm each).
- `figures/projection/projection_A.png`: the GTEx waveform, one point per gene (angle ψ_g, radius amp_g).
- `figures/projection/projection_B.png`: one fibroblast pseudobulk: (log fraction − compartment mean) / reference
  amplitude against acrophase, which falls on one cosine; point area = weight in the fit; λ and φ marked.
- `figures/projection/projection_C.png`: β₁, β₂ of 30 donors × 3 compartments sampled 08-18 h, 2 SD Poisson ellipses;
  the hepatocyte arc moved by a constant offset; panel B's donor ringed, with the arrow from the fibroblast mean (+)
  that B's cosine encodes (length λ, direction φ).

**Video (alternative to figure B, or before it):** `figures/projection_animation/` (`scripts/projection_animation.py`,
50 s; `proj_01.mp4` … `proj_16.png`, one clip and one still per step). Text on screen only.
1-2. The 14 genes over the day; a sample at internal phase φ; each gene is high when φ is near its peak ψ_g.
3-4. Remove the baseline m_g, divide by A_g: every gene reads cos(φ − ψ_g). Display only: the package projects y
     itself, so β keeps the constant L_k m_k (nothing is removed and added back; μ_k absorbs it in M2).
5-6. Place each gene at ψ_g: one cosine, peaking at φ; it follows φ around the day.
7-8. A noisy sample: the weighted least-squares cosine (the package's projector); weights w_g A_g² against the
     unweighted fit (simulated sample: weighted 15.2 h, unweighted 16.2 h, true 15.0 h).

**Speaker notes**
- Drawing convention as in the video: phase 0 on top, turning clockwise; β₁ = λ sin φ is horizontal. The package
  stores (cos, sin); the order of the two components is a convention.
- The waveform was fitted on GTEx fast deaths (Hardy 1-2; 244 donors, 5,611 samples, cultured cells excluded; v10 counts,
  v8 time of death) against clock time of death (0 = midnight): shared acrophases, a mesor and an amplitude per tissue
  (`waveform.fit_waveform`; `clocksync/src/clocksync/data/gtex_waveform.npz`, provenance in its `meta` field). Core genes: NR1D1 3.6 h, DBP 8.5 h, PER2 11.7 h, CRY1 17.9 h,
  ARNTL 21.0 h; largest relative amplitudes ARNTL 0.50, NPAS2 0.33, NR1D1 0.32 (`cs.gtex().table()`). Only relative
  acrophases matter: the zero cancels with the offsets. Time of death is a noisy phase (`talamanca2023sexdimorphic`);
  that noise enters the waveform, not the per-donor readout.
- Why fixed: refitting gene timing failed in simulation (A139: V139.1 and V139.2 FAIL), so the conserved waveform is
  stated as an assumption (METHOD.md §1).
- The waveform choice does not drive the GTEx organ map: A123 uses a cross-fitted waveform; with the all-donor
  waveform, offsets 0.945 and amplitudes 0.970 (unchanged).
- Figure B, how to read it: x = the gene's reference acrophase; y = how far this donor's gene sits from the
  compartment's average, in units of the gene's reference amplitude. If the donor's clock is ahead of or behind the
  compartment average, genes peaking at that phase are up and the opposite ones down, all by the same scaled amount:
  one cosine. Low-amplitude genes (PER1, CRY2, NR1D2) scatter most after the division; they also weigh least.
- Why the compartment mean: gene baselines are unknown, but constant within a compartment. Subtracting the
  compartment mean removes them; by linearity this is the same as taking β relative to the cloud's mean in C. That
  mean is not the centre of the cycle in hospital hours (M2).
- The angle of one β is not a time: the compartment's centre is unknown and, in hospital hours, is not the centre of
  the cycle (figure C). This is why the model of M2 is needed.

**Limits**
- 14 genes, one human waveform (GTEx global; per-tissue waveforms exist, `cs.gtex("Lung")`).
- The conserved waveform is an assumption; violations are not flagged by the gate (A205: V205.1 PASS, V205.2 FAIL).

<!-- src: waveform: clocksync.gtex().table() (data/human_clock_gtex.npz); projection: _core/plane.py (projector,
gene_weights, project), _core/study.py (tau2 = 0.05); A123 sensitivity: results_preliminary.md A123; A139:
paper/exploratory/A139_waveform_refit/output/{v139,d139_3,gene_level,failures}.csv; A205:
paper/exploratory/A205_partial_shift_gate/output/{tests,cells,per_rep}.csv -->

---

## M2. The model (the video)

`figures/model_animation/` (`scripts/model_animation.py`): β_ik = μ_k + z_k f_i + e_ik, built term by term. To write.
Add on the last still or on M4: c_k = |z_k|² V / (|z_k|² V + 2σ_k²), with V the host spread.

## M3. Gate: is a shared, clock-shaped state observable?

**Message:** before any model is fitted, a model-free test asks whether the compartments of a donor share variation,
and whether that shared variation lies in the clock plane. Both are needed; a failed gate is uninformative, never
"decoupled".

**Bullets** (in the order of the figures)
1. *Figure A, what the gate looks at.* For each pair of compartments k, l, the gene × gene cross-covariance across the
   donors they share, of the gene-centred log fractions:

   $$C_{kl} = \frac{1}{n_{kl}-1}\sum_i (y_{ik} - \bar y_k)(y_{il} - \bar y_l)^\top$$

   Measurement noise is independent between compartments, so C_kl keeps only what they share. A shared clock gives
   C_kl = X z_k Σ_f z_l* Xᵀ: rank 2, inside the clock plane; genes peaking together co-vary, opposite genes
   anti-co-vary.
2. *Figure B, is there shared signal in the clock plane?* With P the orthogonal projector on the clock plane:

   $$\text{clock\_cov} = \sum_{k<l} \lVert P\,C_{kl}\,P \rVert_F^2$$

   Null 1: donors shuffled within each compartment (each keeps its data, the pairing is lost).
3. *Figure C, is the shared signal clock-shaped?*

   $$\text{clock\_share} = \frac{\text{clock\_cov}}{\sum_{k<l} \lVert C_{kl} \rVert_F^2}$$

   Null 2: planes built from acrophases shuffled across genes (any 2-D gene pattern, not the clock's).
4. Pass = both p < 0.05. Each null alone can be fooled: clocks that are not shared have a clock-shaped (high-share)
   but unpaired covariance; a shared non-clock factor is paired but not clock-shaped.

**Visual:** three one-panel figures (`scripts/gate.py`, simulated, 150 × 120 mm each; 30 donors × 3 compartments,
deaths 08-18 h, 1,000 permutations).
- `figures/gate/gate_A.png`: fibroblast × hepatocyte cross-correlation, genes ordered by acrophase (correlation for
  display; the gate uses the covariance). Blocks, not bands: with 08-18 h deaths the host spread is elongated.
- `figures/gate/gate_B.png`: clock_cov relative to the pairing null (grey), three cohorts; filled = p < 0.05.
  Shared clock p < 0.001; clocks not shared p = 0.469; shared non-clock factor p < 0.001.
- `figures/gate/gate_C.png`: clock_share against the acrophase null. Shared clock 0.88 (null mean 0.10), p < 0.001;
  clocks not shared 0.53, p < 0.001; shared non-clock factor 0.01 (null mean 0.03), p = 0.582.

**Speaker notes**
- Model-free: no host, no EM, no coupling; only the 14-gene covariances and the waveform's acrophases. It runs on the
  gene-centred log fractions, not on β, and P is unweighted.
- Why cross-covariances and not each compartment's own covariance: within a compartment, Poisson noise and the
  compartment's own clock add variance in the plane; between compartments only shared variation survives in
  expectation.
- Specificity (A154, N0, no clock): pass rate 0.040-0.085 over 8 designs (mildly anti-conservative at 100 donors,
  K 8, single-cell depth). With a shared non-clock factor (N1): at most about 0.10-0.14.
- Sensitivity is limited: under strong shared non-clock structure 0.26-0.415 of true clocks pass (A154, P arm); adult
  GTEx at the fetal organ-presence design passes 0.20-0.21 (A181). Hence: a failed gate is never evidence of absent or
  decoupled clocks.
- In A145, the gate passed 4 of 30 no-clock replicates (0.133; H145.3 fail): the gate alone is not a specificity
  guarantee; the pseudo-clock companion (after the readout) is the second check.
- The pseudo-clock companion is not inside `Dataset.gate`: it is run by the paper's drivers (`paired.pseudo_waveforms`:
  expression-matched non-clock genes given the clock genes' amplitudes and acrophases), mostly on the readouts (A150
  coupling; pairing and offsets in A227-A229, A187), and on the gate itself in some amendments (A237-A239: gate
  `pairing_z` against pseudo sets).

**Limits**
- Needs several donors sharing compartments (A181: the design, not depth, limited the fetal data).
- A conserved-waveform violation (genes shifting relative to each other) still passes (A205); the gate tests clock
  shape, not the exact waveform.

<!-- src: gate: _core/study.py (Study.gate, cross_covariances, plane), dataset.py (Dataset.gate), _core/controls.py
(shuffle_pairing), waveform.py (permuted); figure numbers: printed by scripts/gate.py; A154, A181, gate sensitivity
statement, A145 H145.3: results_preliminary.md (exploratory/A154_gate_calibration/output/summary.csv,
v154_1_n0.csv, n1_map.csv; A181_fetal_depth/output/v181_1_passrates.csv) -->

---

## M4. Fit: host by EM on anchors, then every compartment read against it

**Message:** two steps. The anchors define each donor's host clock state; every compartment is then read against
it. The host is never fitted on the compartments being tested, so a disrupted target cannot pull it.

**Model** (M2): β_ik = μ_k + z_k f_i + e_ik, with f_i ~ N(0, Σ_f), e_ik ~ N(0, σ_k² I + V_ik), z_ref = 1.

Split in two slides: M4a, the EM (video below, bullet 1); M4b, the readout (figures B and C, bullets 2-3).

**Video, slide M4a** (`figures/em_animation/`, `scripts/em_animation.py`, 48 s; `em_01.mp4` … `em_16.png`, one
clip and one still per step; text on screen only). Left the data (3 anchors × 6 patients), right the host plane; the
model is written at the bottom with the term of each step highlighted, and the step's formula under it.
1. The data: known β_ik, V_ik; unknown μ_k, z_k, σ_k², f_i, Σ_f.
2. Start: μ_k the compartment mean (cross), z_k = 1 (arrow pointing to 0 h), σ_k² and Σ_f wide.
3. E-step: each compartment's guess of the host, ĝ_ik = z_k⁻¹(β_ik − μ_k) (small markers in the host plane;
   endothelial guesses are off because z_E = 1 is wrong).
4. E-step: the guesses combined, f̂_i = S_i Σ_k w_ik ĝ_ik, S_i = (Σ_f⁻¹ + Σ_k w_ik)⁻¹, w_ik = |z_k|²/(σ_k² + v_ik)
   (exact for isotropic V; the package uses the 2 × 2 form), with S_i and Σ_f as ellipses.
5. M-step: μ_k and z_k by least squares of β on the host (the host pattern, hollow markers, shifted, turned and
   scaled onto each compartment; residual lines).
6. M-step: σ_k² (dotted circles: what is left beyond V_ik) and Σ_f.
7. Twenty iterations: z_E turns to its value, the guesses of each patient agree.
8. Result: f̂_i and S_i, the regressor of M4b.
The E and M steps are the package's (`factor._EM.estep / mstep`), from the guess of step 2; SQUAREM is not shown.
Fitted: |z_E| 0.80 at 3.1 h, |z_S| 0.87 at 1.2 h (true 0.75 at 3 h, 0.90 at 1 h; 6 patients).

**Bullets** (in the order of the figures)
1. *Figure A, step 1: the host, by EM on the anchors only.* f_i is latent, so the fit alternates:
   - E-step: given the parameters, each donor's host posterior is Gaussian, mean f̂_i and covariance S_i (closed
     form, pooling all of that donor's anchors, each weighted by (σ_k² I + V_ik)⁻¹);
   - M-step: given those moments, μ_k and z_k by weighted least squares, then σ_k², and Σ_f = mean(f̂ f̂ᵀ + S);
   - exact likelihood, accelerated by SQUAREM. Anchors are learned (A17: stromal-vascular; neurons in CNS).
2. *Figure B, step 2: the readout, a regression on the host.* Each compartment k, anchor or not:

   $$\beta_{ik} = \mu_{k,q} + z_k\,\hat f_i + e_{ik}$$

   weighted least squares with one intercept per stratum q (σ_k² iterated). An anchor is read on the host fitted
   without it (leave-one-out, needs ≥ 3 anchors), so its own noise does not enter its regressor.
3. *Figure C, the outputs*, with 95 % donor-bootstrap intervals (the host refitted in every replicate):

   $$c_k = \frac{|z_k|^2 V}{|z_k|^2 V + 2\sigma_k^2},\qquad V = \mathrm{tr}\big(\mathrm{cov}(\hat f) + \overline{S}\big)$$

   coupling c (share of k's clock variation carried by host time), relative amplitude |z_k| and offset arg z_k
   (reported relative to the mean over compartments, or to a reference). Small-sample corrections: |z|² minus
   tr cov(ẑ), σ² raised by the regression's degrees of freedom.

**Visual:** three one-panel figures (`scripts/fit.py`, simulated, 150 × 120 mm each; 40 donors × 6 compartments,
deaths 08-18 h; anchors fibroblast, endothelial, smooth muscle; 500 bootstraps).
- `figures/fit/fit_A.png`: the host f̂_i with 2 SD posterior ellipses, coloured by the true sampling hour (not used
  by the method): the host orders the donors along the arc of their hours.
- `figures/fit/fit_B.png`: observed β (centred) against the fitted z_k f̂_i, both coordinates; hepatocyte c 0.93 on
  the diagonal, macrophage 0.64 around it, lymphocyte (own clock, unshared) 0.00: ẑ ≈ 0, all its variation off the
  line.
- `figures/fit/fit_C.png`: c with intervals against the truth. Targets: hepatocyte 0.931 [0.876, 0.967], true 0.938;
  macrophage 0.641 [0.524, 0.759], true 0.634; lymphocyte 0.000 [0.000, 0.076], true 0. Anchors read low: 0.86-0.89
  against 0.92-0.96. Offsets (h, vs fibroblast): endothelial 3.01 (true 3.0), smooth muscle 1.15 (1.0), hepatocyte
  −3.85 (−4.0), macrophage 2.72 (2.0); amplitudes 0.69, 0.90, 1.14, 0.62 (true 0.75, 0.90, 1.15, 0.60).

**Speaker notes**
- Why two steps rather than one joint EM over all compartments: a target that is disrupted (disease, tumour) would
  enter the host it is compared with. Anchors are chosen before the target is scored.
- Why |z|² and σ² together: c is the signal share; |z| alone mixes coupling with amplitude, σ² alone with depth.
- c depends on the host spread V, i.e. on the sampling window: compare c within a fit, or across groups only at a
  common spread (d*) or on z2 / psi (rules in CLAUDE.md).
- c is attenuated: the host is estimated from a few noisy anchors. A162: estimates are about 0.87-0.92 of the truth
  with anchors at c 0.8 (0.80-0.88 at 0.6), and intervals undercover at high c (median coverage 0.18 at c_target 0.9,
  c_anchor 0.8). Anchors on the leave-one-out host read 0.68-0.73 for a true 0.8. Hence absolute c is reported
  without intervals, as ranks, Δc and d*; contrasts are unaffected.
- Boundary: if the ML σ² is at its floor, c is a lower bound (test inversion at the real design); if |z|² minus its
  noise is ≤ 0, c reads 0 because the loading is noise, not because the compartment is uncoupled.
- A155: reading an anchor on a host that contains it inflates its |z| by 13-20 % with 3 anchors; hence the
  leave-one-out readout (the default).

<!-- src: _core/factor.py (fit, _EM), _core/host.py (host, _hosts, on_host, readout_z, estimates), dataset.py
(Dataset.fit); figure numbers: printed by scripts/fit.py; A162, A155: results_preliminary.md
(exploratory/A154_gate_calibration/output/a162/summary.csv; A155_anchor_self_inclusion/output/summary_bias.csv) -->

---

## M4b. Outputs: what each donor and each compartment gets, and what the analyses use

**Message:** one fit gives a host state per donor and four numbers per compartment. The analyses use the
compartments' numbers, and compare groups only on quantities that do not depend on the sampling window.

**Bullets**
1. *Per donor: the host state* f̂_i (with S_i), a point in the clock plane (units of the reference anchor, centred over
   donors). Its position orders the donors along their internal time; its angle (`phase_h`) is measured around the
   cohort's mean, so it reads as a clock time, up to one constant, only when the cohort covers the day (under hospital
   hours the centre lies inside the arc). Its radius is not an amplitude (shrunk by the prior). Used as the regressor
   of every compartment, and in aggregate only: the host spread V, the coverage of the day (Rayleigh R of host
   phases), and validation against a recorded time where one exists.
2. *Per compartment k*, from β_ik = μ_k + z_k f̂_i + e_ik:

   | output | symbol | read as | used for |
   |---|---|---|---|
   | coupling | c_k = \|z_k\|²V / (\|z_k\|²V + 2σ_k²) | share of k's clock variation that follows the host | ranking compartments within a fit; differences between groups |
   | relative amplitude | \|z_k\| | k's clock swing relative to the reference (or the mean) | organ map; loss of amplitude in disease |
   | phase offset | arg z_k (h) | how far ahead of the reference k's clock runs (> 0: advanced, it reads later) | organ map; phase shifts in disease |
   | own noise | σ_k² (`psi`) | k's clock variation not shared with the host | why c changes: amplitude lost or noise gained |
   | host spread | V | how much of the day the donors span | not a property of k: the denominator to control |

3. *Comparing groups* (`Dataset.contrast`): one host on the pooled anchors, then per compartment and group:
   - Δc = c_case − c_ref, and Δc* at a common spread V0 (each group's |z|² and σ² with the same V0), which removes a
     difference in sampling windows between groups;
   - log ratios of |z|² and σ² (does the clock weaken, or does the noise grow?) and the offset difference (h);
   - donor labels permuted within blocks, the same permutations for every compartment and for the specificity
     companion (M6).

**How c and V are obtained** (per compartment k, `host.on_host`; speaker notes or a backup slide)
1. *Regression.* β_ik = μ_{k,q} + z_k f̂_i + e_ik over k's donors, weighted least squares with weights
   (σ_k² I + V_ik)⁻¹, alternating with the ML σ_k² (Fisher scoring) until both are stable. Gives ẑ_k, its 2 × 2
   covariance cov ẑ, and σ²_ML. The regressor is the posterior mean f̂_i, not a noisy point estimate: regressing on
   E[f | anchors] leaves the slope unattenuated (regression calibration), as long as k is not in the host (hence the
   leave-one-out host for anchors).
2. *Signal.* |z_k|² = |ẑ_k|² − tr cov ẑ, truncated at 0. |ẑ|² overestimates |z|² by its sampling variance
   (E|ẑ|² = |z|² + tr cov ẑ); after truncation a 0 means "loading at noise level" (`boundary_type` z2_truncated).
3. *Own noise.* σ_k² = σ²_ML + p/(m − p)·(σ²_ML + v̄), with p the number of regression parameters (2 per stratum
   intercept + 2 for z), m = 2 × donors the number of observations, v̄ the mean measurement variance. The ML residual
   variance is low by p/m of the total residual, and all of it falls on σ² because V_ik is known. σ²_ML at its floor:
   c is a lower bound (`boundary_type` psi_floor).
4. *Host spread over k's donors.* V = tr( cov(f̂_i) + mean_i S_i ), f̂ centred within strata, over the donors that
   have compartment k. Law of total variance: Var f = Var E[f | data] + E Var[f | data]; cov f̂ alone is too small
   because the prior shrinks f̂, and mean S puts back what the shrinkage removed. V is per compartment because each
   compartment has its own donors.
5. *Coupling.* The variance of k's clock state in the plane splits into tr(z_k Σ z_kᵀ) = |z_k|² V (z is a rotation
   times a scale) and tr(σ_k² I) = 2σ_k²; measurement noise V_ik is excluded:
   c_k = |z_k|² V / (|z_k|² V + 2σ_k²).
6. *Gauge.* The host's scale and rotation are fixed only by the reference (z_ref = 1). Replacing f by a·f changes z to
   z/a and V to |a|²V: |z|² V, and so c, do not change. |z| and arg z are reported relative to the reference (or the
   mean over compartments). In `Fit.table`, z2 and V of an anchor are on its leave-one-out host (another gauge):
   only their product is comparable across rows (in the simulation fibroblast z2 1.90, V 0.264; others V ≈ 0.51).
7. *Intervals.* Donor bootstrap, the host refitted in every replicate, all five steps repeated.

**Visual:** `figures/outputs/outputs.png` (`scripts/outputs.py`, wide; the simulated fit of M4, point estimates):
- left, *per donor*: the host states f̂_i with S_i (1 SD), one highlighted; dashed: cov f̂ + mean S at 1 SD; its trace
  is V (0.51);
- middle, *per compartment* (macrophage): the data, the host pattern (grey, z = 1) and the fitted z_k f̂ (blue), the
  host pattern turned by arg z_k (+2.7 h) and scaled by |z_k| (0.62); arrows z = 1 and z_k from μ_k; dotted σ_k;
- right, *coupling*: per compartment, |z_k|²V (colour) and 2σ_k² (grey) stacked; c is the coloured share
  (fibroblast 0.89, endothelial 0.86, smooth muscle 0.86, hepatocyte 0.93, macrophage 0.64, lymphocyte 0.00).
The table of bullet 2 can sit under it, or replace bullet 2 on the slide; bullet 3 as one line.

**Speaker notes**
- Constant offsets are not identified: the absolute phase of the host is set by the reference anchor, so a single
  donor's host phase is relative to the cohort, not a clock time. Phase offsets between compartments need no
  absolute phase.
- c alone is not comparable between cohorts: it grows with V. The same compartment reads lower in a cohort sampled
  over a narrower window. Hence Δc* (common V0), or |z|² and σ² separately.
- Δc says whether something changed; the |z|² and σ² ratios say what changed (amplitude loss vs added noise), which
  is what the biology needs.
- Per-donor host states are never released for GTEx (controlled time of death): only aggregates (V, R, correlations).
- Package columns (`Fit.table`): coupling (+ lo, hi), rel_amplitude, log_amp_rel, offset_h (+ lo, hi), z2, psi, V,
  own_amplitude = √(|z|²V + 2σ²), boundary; `Fit.host`: x, y, phase, radius; `Fit.host_cov`, `Fit.sigma_f`.
  `Contrast.table`: c per group, d, d_star, p, p_star, log_z2_ratio, log_psi_ratio, offset_diff_h, p_offset.

<!-- src: results.py (Fit, Contrast docstrings), dataset.py (fit, contrast: V0 = V of the pooled host, d_star from
each group's z2 and psi at V0) -->

---

## M5. Simulation: known parameters recovered

To write.

## M6. Pseudo-clock companion: is the coupling specific to clock genes?

**Message:** coupling measures what a compartment shares with its donor's host, and a donor's compartments share more
than a clock. The same pipeline is run on expression-matched non-clock genes; a coupling is called clock-specific
only above theirs (Δc = c_clock − c_pseudo).

**Bullets** (in the order of the figures)
1. *Figure A, why it is needed.* Genotype (eQTL), sex, age, agonal state, ambient RNA: anything shared by a donor's
   compartments moves the 14 clock genes too, and its share in the clock plane is read as coupling. In simulation, a
   compartment with its own, unshared clock reads c ≈ 0 without such donor effects, and c 0.32 (median of 5) with a
   donor effect of sd 0.3 on log expression, while the gate still passes (4 of 5 runs at sd 0.3; the anchors share a
   clock). The gate is pooled
   over pairs; the readout is per compartment and is not protected by it.
2. *Figure B, the companion.* Pseudo-clock sets: 14 non-clock genes matched to the clock genes on mean expression
   (quantile bins), given the clock genes' amplitudes and acrophases (`paired.pseudo_waveforms`), each run through the
   same projection, host fit and readout. Registered test (A150, 5 atlases, 53 compartments): mean coupling 0.393
   with clock genes against 200 pseudo sets, mean 0.242, 95th percentile 0.313; p 0.010.
3. *Figure C, per family.* Above the pseudo 95th percentile: fibroblast 0.670 (pseudo median 0.302), endothelial
   0.520 (0.187), epithelial 0.517 (0.265), myeloid 0.328 (0.148). Not distinguishable from pseudo-clocks: lymphocyte
   0.069 (0.071), mural 0.613 (0.401, q95 0.705), parenchyma and neuron (one compartment each, wide q95).

**Visual:** three one-panel figures (`scripts/companion.py`, 150 × 120 mm each).
- `figures/companion/companion_A.png`: simulated (30 donors × 4 compartments, deaths 08-18 h, anchors fibroblast and
  endothelial, 5 seeds per value): coupling against the sd of donor effects shared by all compartments; a shared clock
  (hepatocyte, c ≈ 0.5 throughout) and an unshared clock (lymphocyte, 0.01 at sd 0 → 0.32 at sd 0.3).
- `figures/companion/companion_B.png`: A150 histogram of 200 pseudo-set means, the clock-gene mean, the 95th
  percentile.
- `figures/companion/companion_C.png`: A150 per family: clock-gene median (diamond, filled above the pseudo 95th
  percentile) and the pseudo median-to-95th-percentile range (grey).

**Speaker notes**
- Where it is used: on the readouts mostly (A150 coupling; pairing and offsets in A227-A229, A187); on the gate's
  `pairing_z` in the developmental series (A237-A239). It is run by the paper's drivers, not inside `Dataset.gate`.
- It is conservative: pseudo genes carry no rhythm, so their coupling is all donor sharing; clock genes carry both.
  A compartment below the pseudo q95 is "not clock-specific at this design", not "uncoupled".
- Example in real data: tumour B / plasma cells pair with the patient's normal tissue (pairing 0.420 vs 0.338,
  p 0.010) but pseudo sets pair as much (q95 0.465): pairing without clock specificity (A227 H227.3).
- Two constructions exist: A150 and A159 differ in pseudo-set construction (myeloid Δc 0.180 on A150, about 0.07 on
  A159); quote each with its method.
- A responsive companion (glucocorticoid and immediate-early genes, `kind="responsive"`) asks whether a shared pulse
  would move as much as the clock genes.

<!-- src: A150: paper/exploratory/A150_atlas_specificity/output/{h150_1.csv,pseudo_pooled.csv,per_family.csv}
(class-fixed rerun, results_preliminary.md "Reruns of possibly stale results"); A227: results_preliminary.md H227.3;
A159 vs A150: same rerun entry; simulation: printed by scripts/companion.py; pseudo sets: clocksync paired.py
(pseudo_waveforms), _core/controls.py (pseudo_clock_genes) -->

Then the results.
