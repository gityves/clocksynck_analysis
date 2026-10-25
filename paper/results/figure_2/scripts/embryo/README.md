# embryo: when does a shared clock state appear?

The mouse embryo (Qiu 2024, E8-E18.75, clock time known), the human fetal atlases (Suo 2022, Cao 2020) and the
first-trimester placenta (Arutyunyan 2023). Everything runs on `clocksync` (`Dataset.gate` with and without `targets` /
`hosts`, `Dataset.timed`) plus `common/` and the builders `paper/common/datasets/{mouse_qiu,fetal_suo,fetal_cao,placenta_arutyunyan}`.
No plotting code. Outputs: `paper/results/figure_2/output/embryo/` (`.../subsample/` for checks).

```
run.sh [--subsample N]      all scripts in order (full run by default)
_shared.py                  eligibility rules, gate rows and Holm over gates, the mouse-embryo windows, and the local helpers below
                            (also imported by natural_experiments)
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_embryo_timed.py` | A108: H108.1 (+ sensitivity, pseudo-clock controls); A108b: S1-S3, D1 | does the host phase of the late embryo (E15.0-E18.75, 16 embryos, 100 cell types) follow clock time 24 x frac(E-day)? R = \|mean exp(i (phi - w t))\|, clock times permuted (10,000); residualised on E-day; 3 pseudo-clock gene sets; leave one batch out; without blood and liver; CNS only; periphery only; glucocorticoid score (descriptive) | `natural_experiments/scripts/a108_embryo.py` (H108.1), `a108b_embryo.py` |
| `compute_embryo_emergence.py` | A108 (emergence, descriptive); A108c V2 (descriptive) | model-free gate over the cell types of an embryo in four windows; per lineage (>= 3 cell types) and window; first appearance = earliest passing window | `a108_embryo.py` (gates), `a108c_lineages.py` (`lineage_gate`) |
| `compute_embryo_lineages.py` | A108 secondary; A108c V1 (descriptive) | per trajectory, complex cosinor of the lineage clock-plane state on clock time (E-day term), BH over trajectories; the same in the early and late half of the window (8 embryos each, batch as group, t permuted within batch) | `a108_embryo.py` (`trajectory_rhythms`), `a108c_lineages.py` (`timed_half`) |
| `compute_fetal_atlases.py` | A107: H107.1, H107.2; A107b: H107b.1-4 | gate over fetal organs (Suo, >= 1000 cells) and organ macrophages (>= 50 cells); Cao whole organs, stromal cells per organ, lung and heart compartments; Holm over the four A107b gates | `a107_fetal.py` (H107.1-2), `a107b_cao.py`, `a105b_organs.py` (gate wrapper, Holm) |
| `compute_placenta.py` | A106: H106.1, H106.2 and its ambient check | gate of the maternal compartments; cross gate of fetal VCT and SCT against the maternal compartments (2000 permutations); same cross gate on synthetic fetal pseudobulks without a clock of their own (a = 10 x XIST soup share); pass stands if the synthetic data fail and the observed `clock_cov` is > 5x the synthetic | `a106_placenta.py`, `natexp.py` (`cross_gate`) |

Secondaries of A106 / A107 (S1, S2, early / late) are in `natural_experiments`.

## Outputs
`a108_timed.csv`, `a108_gates.csv`, `a108c_lineage_gates.csv`, `a108_trajectories.csv`, `a108c_timed_halves.csv`,
`a107_gates.csv`, `a107b_gates.csv`, `a106_gates.csv`, `a106_soup.csv`. Column names follow the legacy files, except
`phase_h` for `peak_phase_h`, `p_pairing` / `p_acrophase` rows carry `passed`, and the trajectory names have no `emb:` / `fetal:`
prefix. The Cao `fetal_cao` table was not rebuilt from the loom (see `data/INDEX_natexp_cancer.md`); `load` uses `from_legacy`.

## Checks
Every script runs with `--subsample` (donors; minutes). Full point estimates against `legacy/desync/validation/natural_experiments/output/`,
full permutation numbers except the last row:

| readout | result here | legacy | note |
|---|---|---|---|
| H108.1 R, p, Fisher-Lee r, offset, error (16 embryos, 100 cell types) | 0.516, 0.016, 0.312, 15.77 h, 2.02 h | identical | 10,000 permutations, same seeds |
| H108.1 sensitivity and the 3 pseudo-clock controls (matched sets) | identical (R 0.533; 0.191, 0.334, 0.194) | | |
| A108b S1 (4 leave-one-batch-out), S2, S3 x2, D1 | identical (R 0.431-0.677; 0.487; 0.421; 0.472; D1 p 0.260) | | |
| A108 trajectories (25: amp, phase, p, q) | identical | | 10,000 permutations, same seed |
| A108c V1 timed halves (25 x 2: amp, phase, p) | identical (max difference 0) | | |
| A108 window gates, four windows: embryos, cell types, clock share | 32 / 26 / 0.0480; 11 / 75 / 0.0439; 13 / 63 / 0.0505; 16 / 100 / 0.0516 | identical (0.0480, 0.0439, 0.0505, 0.0516) | 100 permutations (the pre-registered 1000 take ~50 min for the 100-cell-type window; the full run of the script hit a 25 min limit); p_pairing at its floor 0.0099 in the first three as in legacy (0.001); p_acrophase 0.545, 0.535, 0.267, 0.099 vs 0.634, 0.552, 0.332, 0.083 (Monte Carlo at 100 permutations) |
| A108c lineage gates (4 of 100 checked: Epithelial E8-E10 and E15-E18.75, Endothelium E15-E18.75, CNS neurons E10.25-E12.75) | clock share, null mean, embryos, cell types identical; p within Monte Carlo (e.g. 0.026 vs 0.033) | | 1000 permutations, other seeds |
| H107.1 whole organs (20 fetuses, 4 organs) | p 0.240 / 0.042, share 0.1001 | identical | |
| H107.2 macrophages | 17 fetuses, 2 organs: p 0.053 / 0.007 | 18 fetuses, 3 organs: 0.063 / 0.020 | builder, see Deviations |
| H107b.1-4 (23, 23, 11, 9 fetuses) | p_pairing 0.746, 0.269, 0.028, 0.093; Holm 0.746, 0.538, 0.112, 0.279 | identical | |
| H106.1 maternal gate (16 pregnancies) | 0.012 / 0.030 | identical | |
| H106.2 cross gate VCT+SCT (13) | p 0.507 / 0.131, clock_cov 8.404 | identical | |
| H106.2 ambient check | p_pairing 0.280, clock_cov 1.824, ratio 4.61 (< 5) | 0.357, 1.992, ratio 4.22 (< 5) | soup share differs, see Deviations; same verdict |

## M3 full run (2026-10-25)
`run.sh --jobs 4`, 1.9 h, peak 0.7 GB. A106, A107, A107b, A108 trajectories and A108c timed halves equal the archive outputs
(`outputs/paper/figure_3_embryo`) byte for byte; `a108_timed.csv` to 5e-6 h. First full runs here: the A108 window gates at
the pre-registered 1,000 permutations (clock share 0.0480, 0.0439, 0.0505, 0.0516; p_acrophase 0.634, 0.552, 0.332, 0.083, as
legacy; no window passes) and the A108c lineage gates (104 rows, 33 testable; identical to the legacy
`a108c_lineage_gates.csv` in every column including p, 0 difference; one passes: Mesoderm E15.0-E18.75).
<!-- src: figure_2/output/embryo/{a108_gates,a108c_lineage_gates}.csv; legacy desync/validation/natural_experiments/output/a108c_lineage_gates.csv -->

## Deviations
- **Majority assay (data).** The standard tables apply the legacy majority-assay rule over all compartments of a donor; the
  legacy `fetal` scripts applied it over the compartments of each call (macrophages alone, fibroblast/stromal alone).
  H107.2 therefore loses spleen macrophages (7 of 8 rows >= 50 cells) and one fetus; H107.1 and A107b are unaffected.
  Fix in the builder (keep the assay column and apply the rule in the analysis) if the legacy numbers are needed.
- **Soup share (A106).** Computed from the standard table (majority-assay rows, `XIST` / `depth`) instead of the legacy
  `ambient_genes.csv` (all assay rows): VCT 0.024 vs 0.047, EVT 0.009 vs 0.026, Hofbauer identical. The ambient check is
  stricter for the legacy numbers; the verdict is the same.
- **A108 window gates:** the script uses the pre-registered 1000 permutations (about 1 h for the four windows, then the lineage gates on `--jobs`); the check above used 100.
- **Trajectory names** carry no prefix and the lineage table is read from the standard table with
  `S.trajectory_of()` (reads the legacy pseudobulk `obs`).

## Local helpers (not in clocksync), all in `_shared.py`
- `complex_cosinor`, `timed_wave`: complex rotating-wave regression of a lineage state with a permutation p
  (`clocksync.timed.cosinor` is real-valued per coordinate); the states are read from `Dataset.data`.
- `residualised`: states residualised on a donor covariate per compartment before the host fit (replaces `Dataset.data`).
- `matched_sets`: pseudo-clock sets matched on log mean expression from the 67-gene panel (the panel is too small for
  `pseudo_clock_genes`); `qiu_panel`: the 67-gene embryo table (the standard table keeps the 14 clock genes).
- `trajectory_of`: cell type to lineage map. `soup_share`, `synthetic_ambient`: A106 ambient check.
