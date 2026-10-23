# resilience: ageing, the disease contrasts with their power, fatty liver and gut

Does the internal clock coupling of a cell type change with age or disease, and how large a change could the data have
seen? Everything runs on `clocksync` (`Dataset.gate`, `contrast` with d, d*, companion and MDE, `covariate`, `fit`,
`timed.agreement`, `stats`) plus `common/` and the builders of `paper/common/datasets/`. No plotting code. Outputs:
`paper/results/figure_4/output/resilience/` (`.../subsample/` for checks). GTEx scripts write aggregates only.

```
run.sh [--subsample N]   all scripts in order (full run by default); the --mde steps (hours) are listed in run.sh, not run by default
_shared.py               A17 registry, contrast loader (labels, blocks), anchors, Holm helper
```

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_contrasts_a17.py` | A17 steps 3, 5, 8 (reading rules; A18 add. 3 for d*) | the 36 contrasts of A1-A16 (321 estimable contrast x compartment rows, 307 with a p) under the learned and the fixed anchor rule: gate, c per group, d (bootstrap interval, label-permutation p), d*, p*, log ratios of V, z2, psi, phase offset, Holm over the original primaries, specificity companion, optional MDE (`--mde`) | `common/run_anchored.py` (contrast, mde), `unified.py`, `common_v.py` |
| `compute_procurement.py` | A19 part 1 (P1a, P1b) | is the fibrotic-lung T-cell loss explained by procurement (donor lungs vs living tissue vs explants)? d, d*, companion lung:fibroblast, lymphocyte module, MDE of lung:T | `common/a19_procurement.py` |
| `compute_age.py` | A19 part 2: H19.1, H19.2 and secondaries | locked amplitude, phase drift, own noise per decade for each compartment of the eligible healthy fits; DerSimonian-Laird pooling over fits (own noise rises with age; immune minus stromal-vascular locked amplitude) | `20_age/scripts/{age_fit,ages,pool}.py` |
| `compute_mouse_ageing.py` | A48: H48.1, H48.2, H48.3 | mouse (Wolff 2023, known CT): untimed tissue offsets against timed phases at each age; peripheral synchrony 27 vs 6 months, timed and untimed (d*); hypothalamus against the peripheral median | `mouse_ageing/scripts/a48_wolff.py` |
| `compute_gtex_age_tiebreak.py` | A94 H94.1 tie-break sensitivity (post hoc) | spread of the observed H94.1 over 100 age-tertile tie-breaks | none |
| `compute_gtex_age.py` | A94: H94.1-H94.3, timed slopes | GTEx: individual amplitude (radius), alignment between people, brain-body coherence by age tertile; timed amplitude and phase slopes per decade | `gtex_covariates/scripts/{gtex_states,a94_age_sex_bmi}.py` (A94 part) |
| `compute_gtex_sex_bmi.py` | A99 (H99, sex); A96 part 1 (H96.1, BMI, T2D, liver disease); A94 secondary cue-vs-core | GTEx: radius (individual amplitude) female minus male at equal n, fast and ventilated (labels permuted within Hardy x age tertile); liver timed log amplitude per SD of BMI (permuted within Hardy class), radius on BMI, T2D, liver disease; timed amplitude of the A81b cue state relative to the core state by age tertile (needs figure_S07 `a81b_cue_reference.csv`) | `a94_age_sex_bmi.py` (A99, A96 part 1, cue-core) |
| `compute_a45_learned.py` | A45 sensitivity: the learned A17 anchor rule | TAURUS inflamed vs non-inflamed gut on anchors learned leave one study out (the rule of figure_S04 `anchor_rule.select`, tables of figure_S04), same estimand, companion and permutations as the A45 primary | `a45_taurus.py` (`a45_learned`) |
| `compute_gtex_chiral_age.py` | A101 (core genes only) | CHIRAL-like age dampening of the 14 core genes on the internal phase: pooled, composition-matched, fast deaths | `paired_damping/scripts/a101_chiral_age.py` |
| `compute_masld_gut.py` | A43-A44 (fatty liver, pooled open cohorts), A45 (TAURUS inflamed vs non-inflamed gut) | hepatocyte / T-cell d with specificity companion (hepatocyte minus cholangiocyte, T minus fibroblast), leave-one-study-out, without living donors, per study + random-effects pool, sensitivities of A45 | `masld/scripts/{masld,a43_gse289173,a44_pooled}.py`, `ibd_taurus/scripts/a45_taurus.py` |

A105 main (fatal COVID-19: the length-of-ventilation contrast H105.3, D2, D3, A105c, MDE) is in
`natural_experiments/compute_covid_ventilation.py` (gates and couplings in the same folder); it is not repeated here.

## Outputs
- `a17_{gates,contrast,contrast_companion,mde}.csv`: gate per target; contrast table (`rule` = learned / fixed; `primary`,
  `p_holm_primary`; columns of `Contrast.table`); companion table (learned rule; companion = the target's fixed-rule
  stromal-vascular anchors, neurons in brain); MDE with `--mde`.
- `a19_procurement{,_companion,_module,_mde}.csv`.
- `a19_age_{registry,compartments,pooled}.csv`.
- `a48_{tissue_qc,gates,timed_age,untimed_age,h481,h482,h482_tissues,h483}.csv`.
- `a94_{summary,groups,timed_slopes,age_by_hardy,cue_core}.csv`; `a99_a96_summary.csv`, `a99_{groups,timed_sex}.csv`, `a96_gtex_metabolic.csv`; `a45_learned_{analyses,gates,contrast,companion,anchors,families}.csv`; `a101_core_{tissues,summary}.csv`.
- `a44_{analyses,gates,contrast,companion,meta,mde}.csv`, `a45_{analyses,gates,contrast,companion,mde}.csv`.
Column names are `clocksync`'s (`d`, `lo`, `hi`, `p`, `p_less`, `d_star`, `p_star`, `c_control`, `c_case`, `donors_*`), not the legacy ones
(`diff`, `p_perm`, `c_low` / `c_high`, `p_one_sided`); `p` is two-sided, `p_less` the one-sided p of d < 0.

## Checks
Every script runs with `--subsample` (donors; minutes; `compute_mouse_ageing`: number of permutations). Full point estimates against the legacy
output files (full permutation numbers, bootstrap 200 where a bootstrap is needed; the legacy used 500):

| readout | result here | legacy | note |
|---|---|---|---|
| A17 contrasts, 5 targets (skin_ad, dcm, islets_t2d, myocarditis, ckd_sc; 41 estimable rows, learned rule): c_control, c_case, d, d*, log V ratio | identical (max difference 0) | | gates p_pairing 0.001, p_acrophase 0.001-0.011 |
| same, permutation p of d / p* | median difference 0.018; p* max 0.047 | | Monte Carlo (other seed); ckd_sc differs more (kidney:podocyte 0.24 vs 0.54, proximal tubule 0.44 vs 0.71, d identical): not explained (permutation blocks of the two-stratum kidney table suspected) |
| same, bootstrap interval width | 0.97 x legacy (median) | | 200 vs 500 replicates |
| A19 P1a / P1b, lung:T d (-0.051 / -0.107), d* (-0.067 / -0.128), p (0.61 / 0.081), lung:fibroblast d (-0.166 / 0.178) | | -0.028 / -0.106; -0.043 / -0.127; 0.79 / 0.085; -0.176 / 0.176 | all 36 compartment x contrast rows: c and d differ by <= 0.06 (median 0.01) because `Dataset.contrast` fits the host on the contrast's donors only (legacy: one host on all 269 donors) |
| A19 part 2, heart fit (7 compartments): locked amplitude, drift, gamma | identical raw estimates; `_c` within 0.03 (drift up to 0.04 h); p within 0.07; bootstrap SE within 10-30 % | | simulated-null draws and bootstrap differ (the host is refitted in each replicate here) |
| A48 gates (clock share 0.729, 0.704, 0.700 at 6, 18, 27 months) | identical | identical | |
| A48 H48.2 timed R^2 27-6 months -0.0259 (p 0.190, MDE 0.074); untimed d* +0.0067 (p 0.590, MDE 0.055); d -0.0135 (p 0.296) | | -0.0259 (0.207, 0.078); +0.0067 (0.602, 0.055); -0.0135 (0.311) | observed identical, p within Monte Carlo (1,000 vs 5,000 permutations) |
| A48 H48.3 (c hypothalamus 0.858 / 0.838 / 0.404; peripheral median 0.897 / 0.934 / 0.925) | identical | identical | |
| A48 H48.1 r (all six tissues: 0.902, 0.975, 0.932; p 0.011, 0.012, 0.003) | | 0.933, 0.989, 0.933 (p 0.010, 0.011, 0.004) | PASS 3/3 in both; offsets here are read on one full host and centred on their mean (`Fit.table.offset_h`), the legacy z argument of each tissue on its own leave-one-out host |
| A44 pooled D1 hepatocyte: c 0.759 -> 0.470, d -0.290, p_less 0.001, d* -0.175 (p* 0.005) | | 0.761 -> 0.462, -0.299, 0.001, -0.174 (0.012 two-sided) | D2 -0.383 and leave-one-study-out within 0.02 |
| A44 specificity (hepatocyte minus cholangiocyte), D1 | -0.135, p 0.215 (two-sided) | -0.128, one-sided 0.135 (two-sided 0.248) | companion fails here as in legacy; D2 -0.290 (p 0.028) vs -0.290 (0.026) |
| A44 no living donors / leave out GSE289173 | d -0.123 (p_less 0.062) / -0.164 (0.062) | -0.135 (0.057) / -0.174 (0.054) | |
| A44 per-study random-effects pool (hepatocyte) | -0.241 [-0.489, 0.007] | -0.203 [-0.386, -0.020] | per-study d differ for GSE244832 (-0.610 vs -0.503; boundary couplings) and because the detector also acts per study (macrophage dropped in GSE289173) |
| A45 T cell, five analyses (d, d*, c_control, c_case; companion d_T - d_fib and p) | identical (e.g. primary 0.112 -> 0.050, d -0.062, d* -0.092; companion -0.139, p 0.267) | identical | |
| A94 (full, 1,000 permutations): tertile cuts fast 22/55/62, ventilated 21/44/55; H94.1 fast 0.011 (p 0.76), H94.2 0.0024 (p 0.62), H94.3 0.036 (p_lower 0.84); timed phase slope 0.248 h/decade (fast), 0.467 (Hardy 2) | | cuts identical; 0.008 (0.74), 0.0062 (0.56), 0.033 (0.80); 0.231, 0.452 | ventilated H94.1 here 0.234 (p 0.003), legacy 0.171 (p 0.030): the pooled H94.1 is 0.122 (p 0.020, PASS) here, 0.090 (p 0.081, FAIL) there. Explained (post hoc, `compute_gtex_age_tiebreak.py`, 100 tie-breaks): the per-sample states are identical to legacy (checked for 4 tissues), but ages are integers, so tertile membership at the cut points depends on the random tie-break, and the stream differs. Over tie-breaks the observed H94.1 is fast -0.024 (sd 0.021), ventilated 0.169 (sd 0.031, legacy 0.171, here 0.234), pooled 0.073 (sd 0.018; 5-95 % 0.046-0.099; legacy 0.090, here 0.122). Neither run is a bug; H94.1 as pre-registered is a random variable and its verdict (p 0.081 vs 0.020) is not stable across tie-breaks |
| A99 sex (1,000 permutations): H99 fast, ventilated, mean (p) | -0.0127 (0.74), -0.0946 (0.24), -0.0537 (0.30); FAIL | -0.0127 (0.78), -0.0946 (0.22), -0.0537 (0.36); FAIL | point estimates identical; p differ by up to 0.06 (random age-tertile tie-breaks in the strata, other streams); radius, C, R by sex and timed sex slopes (0.132 log amp, +0.35 h) identical |
| A96 part 1 liver, log amplitude per SD of BMI (n 124), p_lower | +0.0610, 0.789; FAIL | +0.0610, 0.780; FAIL | per-tissue BMI slopes, radius slopes and T2D radius ratios identical for 7 tissues |
| A94 cue vs core, 13 tissues: median change in log(cue / core) amplitude and cue-core offset (h), oldest minus youngest | +0.048, -0.05 | +0.029, -0.00 | same tissues; per-tissue differences up to 0.2 and 1.3 h (cue state built from the reference file `a81b_cue_reference.csv`, not the legacy `A81.waveforms`) |
| A45 learned anchors (full, 1,000 permutations): anchors; T-cell c, d, d*; companion | fibroblast, endothelial, mural (as fixed rule); 0.112 -> 0.050, d -0.062, d* -0.092; -0.139, p 0.27 | identical anchors; 0.112 -> 0.050, -0.062, -0.092; -0.139 (p 0.25) | identical |
| A101 core genes (500 bootstraps): D_core pooled / matched / fast | 0.027 [-0.081, 0.131] / 0.051 [-0.069, 0.129] / -0.005 [-0.208, 0.083]; pooled - matched -0.024 [-0.108, 0.099] | 0.080 / 0.087 / 0.031 | all intervals include 0 in both; tissue set (18 here, A91 set there) and covariate coding differ |

## Deviations
- **Specificity companion (A17 contrasts).** A17 and the 36 original contrasts had no common companion; the companion
  here is the mean d of the fixed-rule anchors of the table (stromal-vascular, neurons in brain), chosen after the
  legacy results and therefore post hoc and descriptive. A44 (cholangiocyte) and A45 (fibroblast) use their
  pre-registered companions.
- **Host donors.** `Dataset.contrast` fits the pooled host on the contrast's labelled donors; legacy A19 part 1 used
  one host on all donors of the table (see check table). A17 contrasts are unaffected (same donors).
- **Primary p of A44.** Hepatocyte d < 0 is read from `p_less`; `p` is two-sided.
- **A17 pooled coherence confirmation (A44, A45)** is not recomputed: the legacy run kept all anchors in every analysis;
  the A22 own-d* detector is applied (it also drops macrophage in the single-study GSE289173 contrast, where legacy
  used the pooled anchors).
- **A19 part 2.** (i) No design check: the legacy rule that a compartment enters the pool only if it recovers planted
  effects (`_ok`) needs a simulator with planted covariate effects, which the package does not expose, so every
  compartment with >= 30 aged donors enters (the pooled estimates are therefore not comparable with legacy `pooled.csv`).
  (ii) The bootstrap SE refits the projection and the host in each replicate (legacy: host fixed). (iii) The
  pre-specified permutation p is not computed; `_p` is the p against the simulated null (the primary of the A19 note).
  (iv) Fits are the gate-passing healthy tables of the builders; ages are the builders' `age` (exact or bin midpoints).
  In `--subsample` runs the eligibility minimum is min(30, N / 2).
- **A48.** H48.1 uses `Fit.table.offset_h` (one full host, centred), see check table; p from 10,000 random permutations
  (legacy: exact enumeration of the 6! orderings); the hypothalamus rows are the non-anchor rows of the same fit.
  H48.2 recomputes the pooled contrast in each permutation (no shared-permutation API for a module mean).
- **A101.** Core genes only (see below); tissues are the 15 v2_gtex host tissues plus liver, adrenal gland and pancreas
  with >= 60 donors (legacy: tissues with >= 50 A91 output genes); covariates: sex, age bracket midpoint, Hardy dummies,
  ischaemic time and its square.
- **A94.** Equal-n subsampling and the permutation scheme as legacy; per-sample states use `clocksync.paired.projector` /
  `states` (library-size depth); donors with duplicate libraries are reduced before the projector.
- **A99 / A96 part 1 / cue-core.** Same states and estimators as A94 (imported from `compute_gtex_age.py`); the covariate ischaemic time is the tissue's
  own sample (`sample_ischaemic_h`, local), as the legacy. The age-tertile strata of the A99 null use a random tie-break (stream differs from the legacy).
- **A45 learned rule.** TAURUS is not in the learning set; the pair coherences used for the confirmation step come from `Dataset.coherence` on all
  compartments of the eligible samples (the legacy: on the fixed anchors only). Same anchors result. `anchor_rule.py` and the healthy-fit tables are
  those of host_identity (imported / read; run its `compute_coherence.py` first).
- **Seeds.** One fixed seed per script; sequences differ from legacy wherever a legacy script shared one generator.

## Not in this folder (not ported)
- **A19 part 1 secondary and part 3:** the ORPHEUS discrepancy (v14 bias simulation, A18 per-cell dispersion on AD_mssm,
  AD_radc, AD_seaad_dlpfc): needs cell-level tables and the per-cell estimators (`cs.cells`) on cell extracts that no builder produces.
- **A19 part 2 validations** (null and planted-effect acceptance simulations, GTEx timed check `gtex_check.py`) and the permutation p.
- **A48:** none of the H48 hypotheses is dropped; A49 (Deota) belongs to figure_S12.
- **A94:** none (A99, A96 part 1 and the cue-versus-core secondary are `compute_gtex_sex_bmi.py`). A96 part 2 (single-nucleus MASLD and islets, cue-referenced amplitude) is figure_S09 / cancer.
- **A101 / A101b:** H101.1, H101.2 on the output genes (A91 classes) and A101b (CHIRAL's own measure) are in `gtex_cues/compute_chiral_age.py` (the A91 classes are its input).
- **A43 alone** (`a43_gse289173.py`: GSE289173 as the single replication cohort) is the per-study row `a44_study_GSE289173`; the A23 phase columns of A44 / A45 (use `figure_S06`'s estimator).
- **A17 healthy fits, anchors, per-cell, gate power, MDE of every row** (`--mde` is implemented and tested on a subsample; the
  full 307-row MDE run is hours and was not run). The healthy fits belong to figure_4 / figure_S05.
- **Not rerun in full:** A19 part 2 (all fits, 200 x host-refitted bootstraps), A17 contrasts for the 31 other targets,
  every `--mde`, A44 / A45 bootstraps beyond 200.

## API gaps (local helpers marked "local: not in clocksync")
- Local helpers of `compute_gtex_sex_bmi.py`: `sample_ischaemic_h` (per-sample ischaemic time; `gtex.load` keeps one value per donor), `cue_waveform` (the A81b reference file -> `Waveform`), `timed_slope` with an age covariate.
- `Dataset.covariate` has no bootstrap and does not run on `Dataset.resample` replicates (the covariate filter and the
  `<id>#<j>` relabelling disagree): `compute_age.boot_se` rebuilds a Dataset from resampled pseudobulk rows.
- No planted-effect simulation at a real design (A19 design check).
- No per-sample clock-plane state API: `compute_mouse_ageing` reads `Dataset.data`; `compute_gtex_age` uses `paired.projector` / `states`
  and local radius, alignment, coherence and timed-slope estimators; A101 uses local amplitude regressions on `Fit.host`.
- `Dataset.contrast` has no module (mean over compartments) statistic with its own permutation p: A48 H48.2 loops.
- `Dataset.gene_rhythms` has no covariates and no fixed host from all donors with regression on a subset (blocks A101 on output genes).
- `Contrast.mde` is returned per compartment only through `contrast(mde=True)` (one call per compartment here).
