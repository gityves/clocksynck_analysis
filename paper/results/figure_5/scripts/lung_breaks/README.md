# lung_breaks: fibrotic-lung T cells and the other candidate breaks

Is a coupling loss in a cell type specific, and what is lost (amplitude, noise, coherence, alignment)? Everything runs on `clocksync`
(`Dataset.gate / fit / contrast`, `Fit.compare`, `decompose.components`, `stats`) plus `common/` and the builders
(`lung_natri`, `hlca_fibrotic`, `hlca_tstates`, `natri_tstates`, `islets`, `skin_ad`, `parkinson`, `psychad_*`, `blood_onek1k`). No plotting code.
Outputs: `paper/results/figure_5/output/lung_breaks/` (`.../subsample/` for checks). `run.sh [--subsample N]`; `_shared.py` (loaders, `module_test`).

| script | amendment | content | legacy replaced |
|---|---|---|---|
| `compute_modules.py` | A18 addenda 1-3 (module test, common spread) | module mean, module minus stroma, named cell minus rest of d (and d*), one shared label permutation; T2D, ILD, HLCA fibrotic, AD skin lesional / non-lesional | `common/module_breaks.py` |
| `compute_decomposition.py` | A98 | per class and group: amp, psi, coupling, log plane / off-plane, coherence, share, alignment; case minus control, donor bootstrap, reading rule | `validation/breaks_decomposition/scripts/a98_breaks.py` |
| `compute_a5_ild.py` | A5 | gates per group, per-class contrast, lymphoid and fibroblast primaries at the >= 20 and >= 15 donor rules (Holm), IPF only, >= 50 cells | `07_lung_ild/scripts/ild.py` |
| `compute_a6_tstates.py` | A6 | resident-like vs circulating-like T (healthy HLCA; fibrotic HLCA; Natri cases), core score, CD8, >= 50 cells, depth check, CD4 / CD8, non-cycling, composition, Natri contrast per state, OneK1K blood | `08_lung_t_residency/scripts/residency.py` |
| `compute_a6_stress.py` (+ `_stress.py`) | A6 secondary (a2), stress check | R^2 of stress gene groups (immediate-early, heat shock, hypoxia) and clock outputs on the host state per T-state compartment; R^2 of the host state on the stress axis | `08_lung_t_residency/scripts/residency.py` (`stress=True`) |
| `compute_a7_activation.py` | A7 part 2 (tests 1-3) | T activation score (CD38, HLA-DRA, GZMB, IFNG, PDCD1, CXCL13); Natri case - control score; coupling and shared amplitude of high vs low activation (HLCA fibrotic, Natri cases, healthy HLCA) | `08_lung_t_residency/scripts/activation.py` |
| `compute_a18_followups.py` | A18 part 3, items 1, 3, 4 (exploratory) | depth thinning of the fibrotic stromal increases; PD DMV neurons and Braak dose-response; HBCC (BD, SCZ) vs halves of the controls | `common/followups.py` |
| `compute_a13.py` | A13 | cortex DLBD and secondaries (neuron / glia mean d, neuron minus glia companion); islets: healthy fit, beta contrast with companion alpha | `03_cortex_conditions`, `15_islets` |

## Checks (full donors; permutations / bootstraps reduced as stated)
| readout | here | legacy | note |
|---|---|---|---|
| module statistics d and d* (15 + 15 values, all 5 contrasts) | identical (max diff 0) | `a18_module_breaks{,_common_v}.csv` | |
| module p, ILD / HLCA, 200 permutations | ILD 0.020 / 0.015 / 0.021; HLCA 0.050 / 0.49 / 0.21 | 0.011 / 0.012 / 0.010; 0.042 / 0.50 / 0.24 | Monte Carlo (2,000 in legacy, default here) |
| A98, 5 contrasts x 10 quantities, 200 bootstraps | 46 of 50 readings equal; components (log plane / off, plane, off, coherence, share) identical | `a98_differences.csv` | differing: ILD amp (no change vs contributes), HLCA align (contributes vs no change), T2D share (contributes vs no change), T2D amp (no change vs contributes); amp / coupling / psi differ by <= 0.05 (group-specific projector, boundary amplitude is sqrt(z2) = 0) and the host is refitted in each draw |
| A13 T2D beta c (control, case), d, d* | 0.707, 0.310, -0.397, -0.190 | 0.707, 0.310, -0.397 | A17 value |
| A6 gates (healthy HLCA state fits) | clock share 0.22, pass | 0.22 | |
| A6 composition (Natri, 10,000 perms) | resident-like 0.302 vs 0.296, p 0.91; cycling 0.000 vs 0.003, p 0.009 | 0.302 / 0.296, 0.91; 0.011 | |
| A6 (a) resident - circulating coupling, healthy HLCA (100 bootstraps) | +0.01 [-0.06, 0.18]; core -0.03; >= 50 cells +0.01 | -0.39 [-0.63, 0.00] (free fit, circulating at the boundary) | not comparable: anchored host, T states are non-anchors, couplings ~0.1 for both |
| A6 Natri contrast per T state | infeasible (< 15 donors) | infeasible | |
| A6 stress: R^2 of the host state on the stress axis (HLCA healthy, fibrotic, Natri cases) | 0.000, 0.004, 0.003 | 0.000, 0.005, 0.004 | host state is not a stress axis; donors 175 / 68 / 60 vs 198 / 71 / 61 (A17 donor rule); major-axis share 0.73 / 0.54 / 0.54 vs 0.76 / 0.53 / 0.60 |
| A6 stress: clock outputs R^2 on the host state, lung:AT1 (healthy HLCA), p | 0.212, p 0.005 (null 0.030); IEG, heat shock, hypoxia 0.017, 0.012, 0.011 (p 0.67-0.91) | same pattern (fibrotic AT2: clock outputs 0.349, p 0.005; stress groups p 0.73-0.99) | different compartments shown; same reading: only clock outputs follow the host |
| A7 test 1, Natri case - control activation score | -0.300 vs 0.136, diff 0.437, 20 / 44 donors | identical (p 0.001 legacy; here the add-one p is 0.000 by rounding at 10,000 permutations) | |
| A7 test 2, T coupling high - low activation: HLCA fibrotic, Natri cases | -0.020 [-0.21, 0.33] p 0.32; 0.000 p 1.00 (both halves at the boundary) | +0.178 [-0.73, 0.94] p 0.63; -0.024 p 0.94 | not comparable in value (A17 anchored host, boundary couplings; legacy a separate reference fit); none significant in either, H7 not supported |
| A7 test 3, healthy HLCA | -0.127 [-0.34, 0.04] p 0.13 | -0.119 [-0.63, 0.48] p 0.83 | |
| A18-3.1 thinning: DCM mural / endothelial / fibroblast, diff after thinning | 0.128, 0.153, -0.043 | 0.132, 0.154, -0.029 | ILD (Natri) EC capillary, ciliated: 0.674, 0.424 vs 0.615, 0.423; p_thinned 0.010 / 0.035 / 0.657 vs 0.010 / 0.025 / 0.761 (200 permutations); the fibrotic increases survive thinning |
| A18-3.3 / 3.4 PD Braak and HBCC | run (`a18_pd_*`, `a18_hbcc.csv`); HBCC BD inhibitory CGE half A: d -0.258, p 0.025 | d -0.263, p 0.015 | PD: DMV contrasts 1-4 vs 0 all non-estimable for glutamatergic / GABAergic (1-3 donors); OPC d -0.376 (p 0.16) |

## M3 full run (2026-10-25)
`run.sh --jobs 4`, 18 min, peak 0.7 GB. Every table equals the archive paper outputs (`outputs/paper/figure_6_lung_breaks`)
byte for byte or to 1e-9 (A6, A7, A13, A98, A18 follow-ups), except `modules.csv`: the d and d* values are identical, p now
uses 2,000 label permutations (archive 200; ILD module mean / minus stroma / named minus rest, d: p 0.011, 0.012, 0.015; legacy
0.011, 0.012, 0.010; HLCA fibrotic 0.041, 0.51, 0.23; legacy 0.042, 0.50, 0.24), and the islets T2D and AD skin (lesional,
non-lesional) targets are new rows. `a7_scores.csv` (per-donor activation scores) is not written here.
<!-- src: figure_5/output/lung_breaks/modules.csv; archive outputs/paper/figure_6_lung_breaks/; legacy p from the Checks table above -->

## Deviations
- **Estimator.** A5 and A6 (pre-A17) used a separate reference fit per group / donor set; here the A17 learned anchors and `Dataset.contrast`
  (pooled host), or the anchored fit of each donor set. Results are therefore read as A17-style, not as the original numbers.
- A98: group-specific estimates come from the two groups' pseudobulks of a class as two non-anchor compartments of one `Dataset` (same pooled
  host, but one projector per group); the bootstrap refits the host (legacy: fixed); skin:fibroblast (an anchor, 2 other anchors) is not computed
  (`Dataset.fit` needs >= 3 anchors).
- `module_test` refits the host in every permutation (cost: ~1.5 s each); default 2,000 permutations as registered.
- A5: gates per group on the group's donors; primaries are infeasible at the 20-donor rule (as legacy), reported at the 15-donor rule (exploratory).
- Builders `hlca_tstates` / `natri_tstates` have no legacy table: stored as parquet by `load_built` on first use (`$SCCC_DATA/pseudobulk`).
- Braak stages (A98 PD) are read from the raw public donor metadata (`parkinson_regions/donors.csv`).

## Not ported / not rerun
- A12 common-axis check is in `host_identity` (`compute_axis_reading.py`).
- A18 part 3 items 1, 3, 4 are also run by `backbone/compute_followups_a18.py` (same estimands; the two scripts agree, the outputs are duplicates).
- A7 test 2: the A7 panel builders `hlca_tstates_a7`, `natri_tstates_a7` are used as they are (whole T class); no separate pre-A17 reference fit.
- Not rerun in full: module permutations (2,000), A98 with 1,000 bootstraps, A5 variants and A6 fibrotic / mixing parts (subsample only), A13 cortex.

## Deviations of the ported items
- **A6 stress, A7, A18-3:** A17 learned anchors and `Dataset.contrast` (pooled host) replace the pre-A17 per-set reference fit, so couplings (A7, PD) differ in value from the legacy; couplings at the boundary (psi = 0) are upper bounds. The stress R^2 of the host state and the clock-output R^2 use `Dataset.data` states; permutations of the host states over donors (200).
- **A18-3.4:** the MSSM controls comparison and the MSSM SCZ MDE are not run (as in the legacy run note).
- A18-3.3 bootstrap intervals: 500 replicates, host fixed.

## API gaps
- `Contrast` returns a p per compartment only: no module / function-of-d statistic with its shared permutation (`_shared.module_test`, local).
- `Contrast.table` has no per-group z2, psi, V (needed by A98).
- `Dataset.fit` needs >= 3 anchors; leave-one-out with 2 remaining is not possible.
- `Dataset.data` (per-sample states) is used for alignment.
