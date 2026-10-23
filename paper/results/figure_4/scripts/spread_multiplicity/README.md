# spread_multiplicity: coupling spread, hierarchical pooling, relative phase, multiplicity

Coupling is a variance share that grows with the host spread, and the contrast table has hundreds of tests. This folder
holds the spread-free readouts of the 22 contrasts of the A20 table (d*, log ratios of slope and own noise, relative phase),
their hierarchical pooling by disease class x cell family, and the multiplicity of the per-contrast tests. Runs on
`clocksync` (`Dataset.contrast`, `stats.bh`) plus `common/` and the builders; the Bayesian model is a local copy
(`_hierarchy.py`, "local: not in clocksync", REFACTOR_PLAN §1.3). No plotting code. Outputs:
`paper/results/figure_4/output/spread_multiplicity/` (`.../subsample/` for checks).

```
run.sh [--subsample N]   spread inputs -> hierarchy -> multiplicity (the last reads resilience/a17_contrast.csv:
                         run that folder's compute_contrasts_a17.py first, same --subsample)
_shared.py               A20 class table, estimands, donor bootstrap of the contrast estimands (local)
_hierarchy.py            verbatim copy of legacy/src/sccc/hierarchy.py (NumPyro NUTS)
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_spread_inputs.py` | A20 inputs; A22 note (own-d* detector); A23 estimand and inference | per contrast and compartment y1 = d*, y2 = log z2 ratio, y3 = log psi ratio, y4 = relative phase (h), bootstrap SEs; anchors confirmed by the detector; R23.3 IBDverse Crohn's on a two-compartment host | `a20_inputs.py`, `common_v.py`, `a22_group_anchors.py` (detector), `a23_phase.py` (inputs, ibd) |
| `compute_hierarchy.py` | A20: R20.1-R20.4, sensitivities; A23: R23.1 | pooled class x family effects, lfsr calls beside BH-adjusted permutation p, resilience interval (+-0.10); y4 pooled with prior scale 2 h | `a20_pool.py`, `a23_phase.py` (pool) |
| `compute_multiplicity.py` | descriptive (A17 reading rules, A20 motivation) | p distribution, BH and Holm over every contrast x compartment test, spread decomposition (V, z2, psi ratios) | `unified.py` (BH / Holm columns) |
| `compute_phase_cancer.py` | A23: R23.4 (cancer arm) | y4 = phase of tumour epithelium minus adjacent-normal epithelium on the host of the learned A17 anchors (h; > 0: tumour later), patient bootstrap (500, host refitted), DerSimonian-Laird over studies | `desync/common/a23_phase.py` (`cancer`) |

## Outputs
`a20_inputs{,_fixed}.csv`, `a22_flags.csv`, `a23_inputs.csv`, `a23_ibd.csv`; `a20_{cells,rows,diag,readouts,sensitivity}.csv`,
`a23_{cells,rows,diag,cancer}.csv`, `multiplicity.csv`, `multiplicity_{summary,p_histogram,per_target,spread}.csv`. Inputs keep the legacy
column names (`y1`-`y4`, `se1`-`se4`, `mde4_h`, `p`, `p_star`, `p_offset`).

## Checks
Every script runs with `--subsample` (donors per contrast, 20 bootstraps, short chains). Full point estimates against `legacy/desync/common_output/`
(3 contrasts, dcm, myocarditis and skin_ad, 22 compartment rows, 200 bootstrap replicates; the full 22-contrast run needs hours and was not rerun):

| readout | result here | legacy | note |
|---|---|---|---|
| A20 inputs y1, y2, y3 (22 rows) | identical (max difference 0) | | |
| A20 SE ratio se1 / se2 / se3 (median) | 1.05 / 1.00 / 1.03 | | host refitted in each replicate here (legacy: host and V0 fixed), other seed |
| A23 y4 (22 rows, confirmed anchors), SE ratio | identical (1e-15); 1.03; median SE 0.76 h | | |
| `_hierarchy.fit` on the legacy `a20_inputs.csv`, y1, 4 chains x 2,000 draws: all class x family medians | max difference 3.5e-9 | | the local copy reproduces the legacy pooled estimates exactly |
| A23 R23.4 cancer: studies estimable here / legacy | 6 estimable, 4 enter the pool (SE > 4 h excluded) / 8 estimable, 8 in the pool | | pooled -1.04 h [-5.11, 3.03] (4 studies) vs +0.54 h [-2.06, 3.14] (8); null stated with its interval either way |
| A23 R23.4 per study (ERP149257, GSE132465, GSE181294, GSE270680, GSE281736, GSE326201) | -0.75, 11.17, 0.00, 3.60, -2.23, -6.12 h | -1.47, 6.67, 0.21, 2.86, -2.23, -6.75 h | GSE281736 identical; the others differ by the bootstrap draws and the host (see deviations) |
| multiplicity (subsample only) | runs; counts depend on the sample | | the full counts need the full figure_5 table (legacy: 307 tests) |

## Deviations
- **Bootstrap SE.** `Dataset.contrast` bootstraps d only. `_shared.boot_se` draws donors within groups, renames repeats `<id>.<j>`
  and rebuilds the Dataset (projection, host) in every replicate; the legacy bootstrap kept host and V0 fixed. SEs agree within 5 % (check table).
- **A22 detector.** Applied on the learned anchors of each contrast with the own-d* rule of the A22 note; the V17.3 positive control and
  the A22 part 2 refit of A20 on confirmed anchors are not run (the A23 inputs use the confirmed anchors; the A20 fit uses the learned anchors, as A20).
- **A23 family "other"** excluded as in A20; rows with SE > 4 h do not enter the pool. IBDverse (R23.3) uses `contrast` with two anchors
  (non-anchor targets only); the cancer arm R23.4 is `compute_phase_cancer.py`.
- **Multiplicity** counts every learned-rule row with a p in figure_5's table (all 36 contrasts); the A20 BH column in `a20_rows.csv`
  is over the 22 pooled contrasts only.
- **Seeds:** one seed per contrast (`--seed` + index).

- **R23.4 cancer arm.** `Dataset.fit` needs 3 anchors after removing the two epithelial compartments (the registration: 2): GSE131907,
  ERP156718 and GSE181919 (2 anchors) are not estimable here (legacy: GSE131907 +0.18 h, SE 1.7). The phase is read only where
  the bias-corrected amplitudes z2 of both epithelial compartments are positive (`Dataset.fit` table; legacy: raw phase), so
  GSE178341 (31 patients, 4 anchors) is not estimable here (legacy +5.10 h, SE 3.2). Both rules reduce the estimable studies from 8 to 6.

## Not in this folder (not ported)
- V16 / V18 simulation validations of the hierarchical model (null calibration, planted recovery) and the V15 null-simulation SEs (sensitivity `se_v15`).
- A22 V17.1-V17.3 and the part 2 re-estimation.
- Full runs: the 22-contrast spread inputs (hours), all sensitivity fits, the hierarchy on the full inputs.

## API gaps
- Bootstrap intervals / SEs of d*, the log ratios and the relative phase from `Dataset.contrast` (with the host fixed).
- Hierarchical pooling across contrasts (paper-only by plan; local copy here).
- A module statistic with its own shared-permutation p.
