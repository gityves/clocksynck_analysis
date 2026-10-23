# natural_experiments: fatal COVID-19, fetal secondaries, placenta sensitivity, GTEx reference

Supplement to figure 2 (development) and figure 5 (fatal COVID-19): the gates and couplings of the COVID-19 lungs and
organs (Melms 2021, Delorey 2021), the length-of-ventilation contrast, the secondary readouts of the fetal and placenta
amendments, and the GTEx power reference of the 14-donor gate. Everything runs on `clocksync` (`Dataset.gate` with and
without `targets` / `hosts`, `fit`, `contrast`, `covariate`) plus `common/` and the builders
`paper/common/datasets/{covid_melms,covid_delorey,fetal_suo,placenta_arutyunyan}`. Shared helpers (eligibility, gate rows, Holm,
soup share, synthetic ambient) are imported from `../figure_3_embryo/_shared.py`. No plotting code. Outputs:
`paper/results/figure_4/output/natural_experiments/` (`.../subsample/` for checks). GTEx: aggregates only (pass rates, per-subset
p values without donor identity).

```
run.sh [--subsample N]      all scripts in order (full run by default; for compute_gtex_reference N = number of subsets)
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_covid_gates.py` | A105: H105.1, H105.2, D1; A105b: H105b.1, S1 | model-free gate of the Melms COVID arm (classes >= 20 cells and >= 10 donors; >= 50 cells sensitivity), the 7 surgical controls (D1), the four Delorey organs (H105.2), the classes within each organ with Holm over the four (H105b.1), lung without the Melms-overlap donors, endothelium across organs (S1) | `a105_covid.py` (gates), `a105b_organs.py` |
| `compute_covid_couplings.py` | A105: H105.2 coupling, D1 (descriptive) | coupling of each Delorey organ on the host of the other organs (leave-one-out, bootstrap, permutation p); Melms lung classes, COVID vs controls on one pooled stromal host (`contrast`) | `a105_covid.py` (`target_coupling`, `per_arm`), `rerun_loo.py` |
| `compute_covid_ventilation.py` | A105: H105.3, D2, D3; A105c (post hoc, no verdict) | epithelial coupling, IMV >= 9 days vs <= 1 day: d, d*, myeloid companion on the same permutations, MDE; Delorey lungs (D3); both studies pooled (A105c); PMI, symptom-to-death and IMV covariates (D2) | `a105_covid.py` (contrast, D2, D3), `a105c_pooled.py` |
| `compute_fetal_secondary.py` | A107: S1, S2, early / late (descriptive) | stromal gate per organ; per-organ coupling on the host of the other organs; the whole-organ gate in early and late fetuses | `a107_fetal.py`, `rerun_loo.py` |
| `compute_placenta_sensitivity.py` | A106: S1, S2 (S3 not run: 4 pregnancies) | per fetal compartment cross gate with its ambient check; fetal-only gate; coupling of the fetal compartments on the maternal host | `a106_placenta.py`, `natexp.py` (`target_coupling`) |
| `compute_gtex_reference.py` | A105: GTEx reference (descriptive) | pass rate of the gate for 50 random 14-donor subsets of GTEx fast deaths and of ventilator cases (4 tissues of the Delorey organs, >= 3 per donor, 200 permutations) | `a105_covid.py` (`gtex_reference`) |

## Outputs
`a105_gates.csv`, `a105b_gates.csv`, `a105_delorey_organ_coupling.csv`, `a105_arms.csv`, `a105_imv_contrast{,_companion}.csv`,
`a105_mde.csv`, `a105_delorey_lung_imv{,_companion}.csv`, `a105c_pooled{,_companion}.csv`, `a105_covariates.csv`,
`a107_secondary_gates.csv`, `a107_organ_coupling.csv`, `a106_sensitivity_gates.csv`, `a106_fetal_coupling.csv`,
`a105_gtex_reference{,_summary}.csv`. Column names are those of `clocksync` (`coupling_lo` / `coupling_hi`, `p_coupling`, `p_less`,
`donors_short` / `donors_long`, `V0`, `log_z2_ratio`), not the legacy ones (`c_lo`, `p_perm`, `p_one_sided`, `n_short`); `p` and `p_star` of
`contrast` are two-sided, `p_less` is the one-sided p of d < 0 (the legacy `p_one_sided`). Compartment names carry no
`covid:` / `lung:` / `placenta:` / `fetal:` prefix in the single-class tables.

## Checks
Every script runs with `--subsample` (donors; minutes; `compute_gtex_reference`: number of subsets). Full point estimates
against `legacy/desync/validation/natural_experiments/output/` with the full permutation numbers (bootstrap 200 as in legacy):

| readout | result here | legacy | note |
|---|---|---|---|
| A105 gates: H105.1 (19), D1 (7), sensitivity, H105.2 (23) | p_pairing / p_acrophase 0.001 / 0.003, 0.131 / 0.069, 0.001 / 0.005, 0.284 / 0.032 | identical | |
| A105b: heart, liver, kidney, lung, lung without overlap, S1 endothelium | all p, Holm (0.054, 0.015, 0.012, 0.054) | identical | |
| H105.2 organ coupling (4 organs: c, psi, n) | 0, 0, 0, 0.601 (kidney); p_coupling 0.013 | 0, 0, 0, 0.601; p 0.011 | intervals differ (e.g. heart upper 0.481 vs 0.435): bootstrap draws; legacy `amplitude` column (\|z\|) not reproduced, `z2` instead |
| H105.3 contrast (8 short, 9 long) | epithelium c 0.713 / 0.639, d -0.075, p_less 0.399, d* 0.147; myeloid c 0.398 / 0.000, d -0.398, p_less 0.133 | identical | the interval is from the host-fixed bootstrap (legacy had none) |
| H105.3 MDE of the epithelium | inf (curve: power 0.06 at s = 1.0, c_pooled 0.600) | inf, 0.600 | |
| D3 Delorey lung (5 short, 7 long) | epithelium c 0.923 / 0.375, d -0.549, p_less 0.077, d* -0.586; myeloid 0.623 / 0.666 | identical | |
| A105c pooled (13 short, 16 long) | epithelium c 0.781 / 0.455, d -0.327, p_less 0.082, d* -0.239; myeloid d -0.025 | identical (p 0.085) | interval [-0.659, 0.011] vs [-0.689, 0.033]: bootstrap draws (500) |
| D2 covariates (PMI, symptom-to-death, IMV; epithelium and myeloid) | coupling slope, locked amplitude, drift identical (e.g. PMI drift -2.038 h) | identical | t-statistics replaced by p and by the null-centred value `_c` |
| D1 arms (Melms lung, per class c and d) | pooled-host contrast (epithelium: control 0.515, COVID 0.691) | each arm on its own host (control 0.547, COVID 0.571) | different estimand, see Deviations |
| A107 S1, early, late gates | p_pairing 0.493, 0.131, 0.921 | 0.500, 0.131, 0.921 | S1: 20 vs 19 fetuses (majority assay, see figure_2 README); early and late identical |
| A107 S2 organ coupling | c 0.049, 0.191, 0.484, 0.391; p 0.306, 0.137, 0.045, 0.104 | c identical; p 0.306, 0.119, 0.043, 0.130 | Monte Carlo (1000 permutations; legacy 1000) |
| A106 S1 cross gates (EVT, Hofbauer, fibroblast), S2 | p identical (0.991, 0.231, 0.254; 0.077), `clock_cov` identical | | |
| A106 S1 ambient checks | `clock_cov` 0.354, 0.405, 0.760 | 0.350, 0.405, 1.346 | soup share differs, see Deviations |
| A106 fetal coupling (5 compartments) | c 0, 0.091, 0.191, 0.342, 0.424; p 1, 0.269, 0.316, 0.021, 0.041 | c identical; p 1, 0.248, 0.307, 0.010, 0.043 | Monte Carlo |
| GTEx reference (50 subsets x 2 cohorts) | pass 1.00 (fast), 0.84 (ventilator); every subset p identical | identical | |

## Deviations
- **D1 per-arm coupling.** The legacy `per_arm` fitted each arm on its own two-anchor host (endothelial, fibroblast).
  `Dataset.fit` needs >= 3 anchors, so `a105_arms.csv` reads both arms on the pooled host with `Dataset.contrast` (c per arm,
  d, d*, spread-free log ratios). The estimand is the one of the contrast amendments (common host), but the numbers
  differ from `a105_arms.csv` of the legacy run (descriptive only).
- **Dataset of the contrasts holds the labelled donors only** (H105.3, A105c): states are centred per compartment over the
  donors of the dataset, as in the legacy `Study` of the labelled donors. With all COVID donors in the dataset d changes
  (-0.120 vs -0.075).
- **Majority assay and soup share** as in `../figure_3_embryo/README.md` (A107 S1: 20 vs 19 fetuses; A106 ambient checks).
- **Leave-one-out host (A105 H105.2, A107 S2).** As in the legacy `rerun_loo.py`; `Dataset.fit` reports `z2`, not the legacy
  `amplitude` (|z| on the leave-one-out host).
- **Seeds** of the bootstraps differ from the legacy streams; intervals agree to about 0.05, p within Monte Carlo.

## Local helpers (not in clocksync)
All in `../figure_3_embryo/_shared.py` (`soup_share`, `synthetic_ambient`, `holm_max`, `eligible`); none is a clock
quantity except the ambient mixture. API gap: `Dataset.fit` with 2 anchors and non-anchor targets (legacy `per_arm`).
