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

EM (SQUAREM) on the anchor compartments gives each donor's host state; each compartment is then regressed on it (z,
σ², c, |z|, arg z; anchors on their leave-one-out host; donor bootstrap with the host refitted). To write.

## M5. Simulation: known parameters recovered

To write. Then the pseudo-clock companion (A150 as the example) before the results.
