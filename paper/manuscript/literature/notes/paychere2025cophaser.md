# paychere2025cophaser
**Ref:** Paychère, Salati, Gobet, Naef, bioRxiv (v posted April 9, 2026). DOI 10.64898/2025.12.24.696376.
Type: method. Lab: Naef (EPFL), our lab (our previous work).
**Source:** full text

## Question and design
Learn context-dependent periodic manifolds in scRNA-seq: separate cycle phase from cell context (cell type, state,
disease) while keeping a phase definition comparable across contexts. Applications: cell cycle (RPE1 FUCCI/10x,
mouse embryo VASA-seq, breast cancer 69 specimens / 55 patients, paediatric AML, spatial ovarian cancer), circadian
(mouse aorta, ZT 0/6/12/18, Auerbach 2022), menstrual cycle (endometrium, endometriosis), somite clock.

## Method
VAE with a phase (f-space, Power Spherical posterior) and context (z-space, Gaussian) latent; MINE penalty on their
mutual information. Expected counts: "a periodic component Fg (φc ), a baseline expression µcg , and a context-dependent mean shift ∆µg (zc )", plus "a context-dependent cycling amplitude λ(zc ) that scales the periodic component". Negative binomial likelihood. Fourier profiles (gene acrophases) are shared across contexts; contexts shift means and scale amplitude only. Optional Bernoulli cycling indicator.

## Key results
- Simulation: phases and context-dependent mean/amplitude recovered; without z-space, "severe artifacts in the inferred phases".
- Mouse aorta circadian: "the inferred phases closely tracked the ZT, with a MAD of 1.67 h across all cell types"; lower than Tempo, gains for fibroblasts and macrophages.
- Aorta amplitude: "global circadian rhythmicity has the highest amplitude in smooth muscle cells and fibroblasts, while appearing slightly dampened in endothelial cells and macrophages".
- Cell-type-specific clock-gene means: Npas2, Hlf nearly silent in endothelial cells.
- Breast cancer cell cycle: TNBC highest proliferation; CDC20 phase-specific vs KRT14 constitutive overexpression.
- Somite clock coupled to the cell cycle around G1/S.
- Stated assumption: "CoPhaser assumes the existence of a shared cycling gene program across contexts".

## Not addressed / limits
- Within-sample phases: needs cells covering the cycle (fine for cell cycle; for circadian in vivo, one sample = one time, so phase spread across cells is small and the aorta set needed 4 ZTs).
- No donor-level structure: does not ask whether compartments of one donor share a clock.
- MAD vs ZT is aligned per cell type independently: inter-cell-type phase offsets are not evaluated.

## Link to us
- **Intro claims:** phase-inference toolbox (single cell); concept source for clocksync: fixed gene acrophases shared across samples, with fitted amplitude scaling and mean shifts per context (clocksync: per-compartment μ_k and z_k on a shared clock state f_i). Cell-type differences in clock amplitude (SMC, fibroblast high; endothelial, macrophage lower).
- **Relation:** our previous work; precursor of the model, not a competitor.

## Ideas
- Cite the aorta amplitude ranking as independent support if clocksync |z| ranks stromal/fibroblast > endothelial/immune in human data.
- Combine: CoPhaser cell-cycle phase to control proliferation in tumour compartments before clocksync (cell-cycle gating, cf. Nagoshi 2004).
