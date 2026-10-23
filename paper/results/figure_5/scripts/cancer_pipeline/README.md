# cancer_pipeline: the single-cell cancer cohort pipeline

The pre-registered per-cohort pipeline of `docs/cancer/ANALYSIS_PLAN.md` (the "21-study pipeline"): does tumour epithelium
carry its own clock, or is it carried by the host time of the patient's other compartments? One script, no plotting code.
Outputs: `paper/results/figure_5/output/cancer_pipeline/studies/` (`.../subsample/studies/` for checks).

```
run.sh [--subsample N]     compute_cohort_pipeline.py (full run: hours, 72 cohort fits with 500 resamples each)
_shared.py                 loads cancer/_shared.py (loaders and the local joint-model helpers)
```

Split with `cancer`: figure 6 holds the anchored estimates and the decomposition of the tumour clock (A52, A83,
A86, A87, A97, A17 cancer, A18.3); this folder holds what the 2026-10-03 plan specifies per cohort with the joint factor
model: gate, E1-E3, H1-H3 contrasts, E5 own rhythm, E4 gene loadings and the microenvironment contrasts.

## Analysis (`compute_cohort_pipeline.py`)
| step | hypothesis / estimand | public call |
|---|---|---|
| eligibility, reference compartment | >= 6 patients, >= 3 compartments; reference = first of tumour fibroblast, endothelial, myeloid, T_NK, mural, B_plasma | `Dataset` |
| gate | shared, clock-shaped state observable (1,000 permutations); a cohort enters the pooled tests only if both p < 0.05 | `Dataset.gate` |
| E1-E3 | coupling, amplitude, offset of every compartment, joint model with the reference compartment, patient bootstrap | local (`_shared.joint_estimates`) |
| H2, H3 (`_contrasts`) | log amplitude ratio and offset difference of the malignant compartment against every other compartment (tumour epithelium in `main`, copy-number-called malignant cells in `cnv`) | local (`joint_contrasts`) |
| microenvironment (`_micro`) | each tumour compartment against the same compartment of adjacent normal tissue; non-malignant tumour epithelium against normal epithelium; BH over comparisons | local (`joint_contrasts`) |
| E5 own rhythm | clock-plane variance of each compartment's residual after its host-coupled part, against permuted-acrophase planes (500) | local (`own_rhythm`) |
| E4 gene loadings (descriptive) | per-gene timing deviation `dev_h` of the malignant compartment from the conserved waveform, pooled per gene over gated cohorts | local (`gene_loadings`) |
| pooling | DerSimonian-Laird over gated cohorts, SE = bootstrap interval width / 3.92, pan-cancer and per cancer type | `stats.meta`, `stats.bh` |

Legacy replaced: `legacy/scripts/run_study.py` (per cohort) and `legacy/scripts/m4_summary.py` (pooling). Data:
`cancer_sc` builder (`from_legacy(study, cnv)` = the legacy pseudobulk tables; `build()` from the counts equals them on
slices). Cohorts: the 43 of `STUDIES` (variant `main`) and the 29 with copy-number calls (`cnv`).

Outputs (per cohort `<study>[_cnv]_*.csv`): `gate`, `estimates`, `e5`, `e4`, `contrasts`, `micro`; pooled:
`summary_{main,cnv}_gate.csv` (with the pass flag), `summary_{main,cnv}_e5.csv` (cohorts with an own rhythm, p < 0.05, and
median `coupling_clock` per tissue x compartment type), `summary_{main,cnv}_e4_genes.csv`,
`summary_{main,cnv,main_micro,cnv_micro}_{contrasts,meta,meta_by_cancer}.csv`. Column names follow the legacy files
(`n_donors` for `n_resections` in the gate table).

## Checks
`run.sh --subsample 12` runs both variants (three small cohorts, 12 patients, 20 resamples, 50 permutations; a cohort that
fails the gate or is singular at that size is reported and does not stop the run). Full point estimates, with 60 resamples instead of 500 (point
estimates do not depend on it), against `outputs/studies/` of the legacy run:

| cohort | readout | result | legacy |
|---|---|---|---|
| SRP273482 (9 patients, 8 compartments) | gate p (pairing, acrophase), clock share | 0.001, 0.006, 0.304 | 0.001, 0.006, 0.304 |
| SRP273482, ERP149257 | E1-E3 (amplitude, offset, coupling, psi), 22 compartments | identical (max difference 3e-11) | |
| both | E5 (rho, p, coupling_clock) | identical (rho: 2e-6, from the permutation mean) | |
| both | H2/H3 contrasts and microenvironment contrasts, point estimates (30 rows) | identical (max difference 4e-11) | |
| both | contrast intervals | differ by up to 4 h (60 vs 500 resamples) | |
| both | E4 gene loadings | computed; no legacy output at the pipeline level (E4 was used by the cycling follow-ups only) | |

The pooled meta-analysis (`summary_*`) was run on these two cohorts only; the 72 full fits (500 resamples) were not rerun.

## Deviations
- **Gate.** `Dataset.gate` (same estimator, checked identical above); the gate table stores `n_donors` instead of
  `n_resections` and the passing flag is added in the summary.
- **Subsample.** With `--subsample N`, at least 12 patients per cohort are kept (singular fits below).
- **Estimator.** E1-E5 use the joint factor model of `Dataset.study` through four thin local wrappers
  (`_shared.joint_estimates`, `joint_contrasts`, `own_rhythm`, `gene_loadings`); `Dataset.fit` (anchored host) is a
  different estimator and is used in `cancer/compute_anchored_cancer.py`, not here.
- **Ported elsewhere (work package of 2026-10-22).** The analyses that this README used to list as not ported are in other folders:
  A24 cancer arm, the cycling follow-ups T1-T3, C1-C4 and C2b, A85, A84 TCGA arm and the A17 anchors of the 21 cancer cohorts in
  `cancer`; A18 part 2 for cancer (per-cell dispersion, `cancer_split.py`) in `backbone`; A23 cancer arm (R23.4) in
  `spread_multiplicity`. Per-cell inputs: `data/cells_cancer`, `data/cells_cancer_cycle`.
- **Full run.** The 72 fits with 500 resamples (hours) were not rerun; the check above is the point estimates of two cohorts. The
  pipeline was run end to end on three cohorts (`run.sh --subsample 12`, 2026-10-22): per-cohort tables for the main and the
  copy-number variants (`gate`, `estimates`, `e5`, `e4`, `contrasts`, `micro`) and every `summary_*` table are written; two of
  the three cohorts (GSE132465, GSE181919) fail the gate at 12 patients and are reported as such, SRP273482 passes (pairing p 0.020,
  clock-shape p 0.039, 9 patients, 7 compartments). With `--subsample` the pooled tables hold one gated cohort only.
