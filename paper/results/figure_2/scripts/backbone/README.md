# backbone: the stromal-vascular backbone, per-cell synchrony, strong clocks are coupled

Where the host-locked clock sits across cell families in healthy tissue (stromal and vascular cells and neurons carry it;
immune cells hardly), how much of a compartment's clock is per-cell amplitude and how much within-sample synchrony, and
whether tissues with strong own clocks are the more tightly coupled. The per-dataset fits are in
`atlas_fits`, the anchors and host identity in `host_identity`. Runs on `clocksync.Dataset` plus
`clocksync.cells.vonmises`. No plotting code. Outputs: `paper/results/figure_2/output/backbone/` (`.../subsample/` for checks).

```
compute_oscillator_gtex.py   A104 (kept as written: reproduces legacy, rho 0.734)
compute_family_summary.py    healthy couplings by cell family under three estimators
compute_per_cell.py          A18 part 2: per-cell lambda and R (healthy and disease contrasts)
compute_cancer_split.py      A18 part 2, contrast 3: per-cell lambda and R of the cancer cohorts (tumour, adjacent normal, copy-number malignant)
compute_positive_control.py  A18 part 2 positive control (V10): mouse dermis with known ZT, mixtures of ZT groups
compute_followups_a18.py     A18 part 3 items 1, 3, 4: depth thinning, PD DMV neurons by Braak stage, HBCC against control halves
run.sh [--subsample N]
```
Imports `Target` and helpers from `atlas_fits/_shared.py` and `family` from `host_identity/anchor_rule.py`.

## Analyses
| script | amendment, hypothesis | estimand | legacy script replaced |
|---|---|---|---|
| `compute_oscillator_gtex.py` | A104: H104.1, S1 | Spearman rho(own amplitude, coupling) over GTEx tissues, periphery and brain | `a104_gtex.py` (S2 is in organs) |
| `compute_family_summary.py` | A17 step 8 and the estimator note of 2026-10-13 (descriptive) | per gate-passing healthy fit (34 datasets): point coupling of every compartment under learned anchors, fixed anchors and the free host with the A17 corrections (`Dataset(anchors=None)`); per family the median over compartments and a percentile interval from a 2,000-draw bootstrap over datasets; variants with and without the Tabula Sapiens confirmation fit | `unified.py` (family medians of the PLAN note) |
| `compute_per_cell.py` | A18 part 2 (exploratory, as labelled there) | per target (hlca_fibrotic, lung_hlca, ild_natri, islets_t2d, pd, age_mssm, AD_mssm, AD_radc, AD_seaad_dlpfc), compartment and group: D (within-sample in-plane dispersion over (donor, compartment) cell groups of >= 100 cells, minus its acrophase-permuted null, 10 planes), rho^2 = s / (1 - R_W^2) with s = z2 V, lambda and R from lambda R = rho and the von Mises model; case-minus-control contrasts of log lambda, R, D with donor-bootstrap intervals (100); eligibility: >= 20 cell groups, median clock UMI per cell >= 2, D above null | `cell_split.py`, `x5_vm_human.py` (`split_lr`) |
| `compute_cancer_split.py` | A18 part 2, contrast 3 (exploratory) | per compartment pooled over the 13 cohorts with copy-number calls and barcoded cell tables (21 gate-passing cohorts of the frozen A17 table): D of (patient, compartment) groups of >= 100 cells minus its acrophase-permuted null (10 planes), cell-weighted; rho^2 = s / (1 - R_W^2), s = 2 psi c / (1 - c) of the cohort's learned-anchor fit (copy-number table for `tumour:malignant` / `tumour:epithelial_nonmalignant`, main table for the rest), R_W = sinc(5 omega) (08-18 h), rho = median over cohorts; lambda, R as above; paired contrasts of R and log lambda over the same patients (>= 5), delete-one-cohort jackknife SE; eligibility as above | `cancer_split.py` (+ the groups of `x5_vm_human.py --cnv`) |
| `compute_positive_control.py` | A18 part 2 positive control (V10) | mouse dermis (GSE223109, 18 samples, ZT 2-22): D per (sample, cell type) group (>= 100 cells), rho from the timed design (between-sample variance of the pseudobulk clock coefficients minus the measurement variance, over 1 - R_W^2), lambda and R per cell type; two-sample mixtures Delta h apart (equal cell numbers): R against the prediction R abs(cos(omega Delta / 2)); every row also after binomial thinning to 2.4 clock UMI per cell | `positive_control.py`, `v10_cell_dispersion.py` (data) |
| `compute_followups_a18.py` | A18 part 3, items 1, 3, 4 (exploratory) | (1) case pseudobulks thinned binomially to the control median depth of the compartment, anchored contrast rerun (200 permutations, 200 bootstraps), depth and cells per group; (3) PD DMV compartments: donors, cells, depth per group, the contrast restricted to pseudobulks of >= 50 cells, coupling by Braak Lewy-body stage (0, 1-4, 5-6; each against stage 0 on its own two-group host, bootstrap SE); (4) BD and SCZ (HBCC) against each of 10 random halves of the shared controls (200 permutations) | `followups.py` (`thin_one`, `pd_dmv`, `pd_braak`, `hbcc_one`) |

Not in this folder: A18 part 1 (windows: `windows`), A18 part 3 items 2 (copy-number coupling: `cancer`) and 5
(cortex modules: `lung_breaks`); the validation of the dispersion estimator on mouse skin (V10, the timed waveform fit and the
mixture design of the original) is in `calibration`; here only the A18 positive-control readout.

## Outputs
`a104_gtex.csv`, `a104_gtex_tissues.csv`; `family_couplings.csv` (target, compartment, family, coupling under the three
estimators, anchor flags), `family_medians.csv` (variant, estimator, family, median, interval, counts);
`a18_cell_split.csv`, `a18_cell_split_contrasts.csv` (legacy columns; `boot` replicates are not saved);
`a18_cancer_groups.csv` (D per group), `a18_cancer_split.csv`, `a18_cancer_split_contrasts.csv`;
`a18_positive_control.csv`, `a18_positive_control_groups.csv`; `a18_thinning.csv`, `a18_pd_dmv.csv`, `a18_pd_braak.csv`, `a18_hbcc.csv`
(legacy columns; the Braak table adds `c_stage0`, `d_vs_0` and its interval). Cells: the builders `data/cells_census`
(`$SCCC_DATA/cells/<dataset>/cells.parquet`, per-cell clock-gene counts of 7 Census datasets), `data/cells_cancer` (cancer cohorts) and
`data/cells_mouse_skin` (mouse dermis).

## Checks
`run.sh --subsample 30` runs the six scripts (38 min for the first three on a busy machine, the per-cell step dominating; 30 donors, 100 gate permutations, 5 bootstraps). Full data against
`legacy/desync/common_output/`:

| readout | here | legacy | note |
|---|---|---|---|
| A104 | see the A104 commit (rho 0.734 identical) | | unchanged |
| family medians, learned anchors (all 34 fits with Tabula Sapiens, 267 compartments) | lymphocyte 0.09, myeloid 0.27, endothelial 0.55, epithelial 0.40, fibroblast 0.65, mural 0.67, neuron 0.47, glia 0.21 | identical (PLAN note) | |
| couplings, 212 compartment rows of the 22 non-Tabula fits | learned max difference 0.0, fixed 0.0, free-corrected 3e-6 | `a17_unified_healthy.csv` | |
| free-corrected family medians | e.g. fibroblast 0.71, mural 0.68, neuron 0.78, glia 0.35 | identical | lymphocyte 0.30, epithelial 0.60 here vs 0.35, 0.46: the legacy medians also include the gate-failed fits (274 compartments); only gate-passing fits here |
| per-cell D, 5 targets (hlca_fibrotic 33, pd 63, age_mssm 18, ild_natri 29, islets_t2d 16 rows) | identical (max 4e-15) | `a18_cell_split.csv` | seed 0 |
| per-cell rho, R | median relative difference of rho 1.7-4.6 %, median |R| difference 0.003-0.009 | | host refitted per group (deviation 1) |
| eligible flags, same 159 rows | identical (13/18 age_mssm, 7/33 hlca_fibrotic, 29/63 pd, 2/29 ild_natri, 0/16 islets) | | |
| cancer split, 13 cohorts (2,249 groups): D per group, pooled D | identical (max 2e-13) | `vm_cnv_groups.csv`, `a18_cancer_split.csv` | the groups are computed here from the cell tables, one random stream in the legacy order |
| cancer split, rho, lambda, R, jackknife SEs of R and log lambda (26 compartments), eligible flags (6) | identical (max 4e-15) | | e.g. `tumour:malignant` R 0.333, `adjacent_normal:epithelial` R 0.413 |
| cancer split contrasts (6) | identical (0) | `a18_cancer_split_contrasts.csv` | malignant vs adjacent-normal epithelium d_R -0.249 (SE 0.160), d log lambda -0.161 (0.172); tumour vs normal fibroblasts d_R -0.300 (0.208) |
| positive control, 16 rows (4 real, 12 mixed; unthinned and thinned) | identical (max difference 0) | `a18_positive_control.csv` | fibroblasts R 0.884 (thinned 0.891), immune 0.609 (0.671); R_W 0.059 / 0.100 |
| thinning, 13 contrasts (5 targets) | depth, cells, d and d thinned identical (6e-14); p within Monte Carlo (max 0.11 observed, 0.20 thinned) | `a18_thinning.csv` | DCM mural +0.132 / +0.132, ILD EC_capillary +0.723 / +0.615, HLCA T -0.140 / -0.120 |
| PD DMV contrast (cells >= 50), 8 compartments; describe table | c_control, c_case identical; p within Monte Carlo (max 0.08); donors, cells, depth identical | `a18_pd_dmv.csv` | neuron_other -0.446, p 0.046 (legacy 0.063) |
| PD DMV by Braak stage, neuron_other | 5-6: 0.000 vs stage 0 0.446 (host of stages 0 and 5-6; 0.540 on the host of 0 and 1-4); 1-4: 0.233; 5-6 minus 0: -0.446 (bootstrap interval -0.887, 0.069) | 0.524 / 0.231 / 0.000; -0.524 (jackknife SE 0.285) | deviation 5 (host and SE differ); 5-6 minus 0 negative in 6 of 7 estimable compartments, as legacy (ependymal +0.171) |
| HBCC halves, 18 compartment x disease means over 10 splits | mean d identical (max 0.005); share of halves with p < 0.05 within 0.15 | `a18_hbcc.csv` | SCZ inhibitory_CGE -0.36 (90 % of halves), excitatory -0.28 |

Subsample check of `compute_cancer_split.py` (WP0 item 2, 2026-10-23): the archive's stored subsample `a18_cancer_groups.csv` came from an earlier script version that drew each permuted plane and its splits interleaved; only `D_null` (and the columns derived from it) differ, D is deterministic. The draw order here is the one of the archive's full run: a copy of the current archive script and this one give byte-identical outputs at `--subsample 12 --seed 0`, and the archive full run in this order (snapshot `a18_cancer_groups.csv`) matches the legacy `vm_cnv_groups.csv` (2,249 groups, D_null to 7e-13); the migrated full run is left to M3. `D_null` is a mean over 10 permuted planes and is noisy per group (the two orders differ by up to 144 per group, same mean); read it pooled.

## Deviations
1. **Host in the per-cell rho.** The legacy held the pooled anchored host fixed and took s = 2 psi c / (1 - c) of the group's
   donors. `Dataset.contrast` reports only c and ratios per group, not z2, psi, V, so here the host is refitted on each
   group's donors (and on each bootstrap draw) and s = z2 V of `Dataset.fit`. rho differs by 2-5 % at the median; R by a few
   thousandths. Bootstrap draws also refit the host, as the A17 bootstrap does (legacy: host fixed).
2. **Seeds.** One seed (0) for all targets (this reproduces the legacy D exactly); bootstrap draws differ by Monte Carlo.
3. **Family medians** include the Tabula Sapiens fit in variant `all` (as the legacy PLAN note: 267 compartments) and
   exclude it in `no_tabula`.
4. **Not rerun:** bootstraps of the per-cell contrasts (full 100 draws take hours per target), the full `run.sh`.
5. **Braak groups (A18 part 3 item 3).** `Dataset.contrast` fits one host on the donors of the two groups it compares, so the stage-0
   coupling and the other stage's coupling are read on the host of stages 0 and that stage (legacy: the host of all donors, fixed;
   jackknife SE, a deviation of the original already). The SE here is the donor-bootstrap SD of `Dataset.contrast` (host fixed within
   the bootstrap). The pre-registered Spearman permutation test was replaced in the original by stage groups; the same here.
6. **Cancer split cohorts.** Only the cohorts whose legacy dispersion table carries barcodes can be split by copy-number call (the
   first extraction has none): 13 of the 21, as legacy. `data/cells_cancer/build.py` rebuilds the tables with barcodes for the others.
7. **Positive control.** `Dataset` needs two compartments: the pseudobulks of one cell type are duplicated under another name to read the
   clock-plane coefficients and their variances from `Dataset.study.beta` / `.V` (the projection is per compartment: identical to the
   legacy single-compartment `Study`; rho identical to the digit above).

## Local helpers (not in clocksync)
- `compute_cancer_split.py`: `split_lr` (also used by `compute_positive_control.py`), `study_s` (s = 2 psi c / (1 - c) from `Dataset.fit`),
  `cell_groups` (copy-number split of the tumour epithelium); `compute_positive_control.py`: `coherent_variance` (reads `Dataset.study.beta`
  and `.V`); `compute_followups_a18.py`: `braak_stage`.
- `split_lr` (lambda, R from rho and D): `clocksync.cells.vonmises.solve` takes D and |Gamma|^2 but no given rho.
- `VM._split_projectors`: used from `clocksync.cells.vonmises` (a private name; no public constructor of the projectors).
- `coherent_s`: s = z2 V of one compartment on a donor set (host refitted).
- `map_cells`: cell-to-compartment rules of the builders (`paper/common/datasets/_atlas/{lung,cortex,census}.py`,
  `data/islets`, `data/parkinson`) applied to the per-cell tables.

## API gaps (for clocksync)
- Per-group `z2`, `psi`, `V` in the `contrast` table (and on the fixed pooled host) to replace deviation 1.
- Public `cells.vonmises.split(D, rho, c, LA)` and projector constructor; a `cells` entry point from per-cell counts and a
  compartment map.
