# psoriasis: the psoriatic delay and its chain of checks

The psoriatic lesion keeps a clock-shaped covariance but its clock plane is turned relative to the uninvolved skin of the
same patient and visit (delay, offset < 0) and carries more variance (r > 0). The chain: replication on held-out series,
dependence on IL-17 activity and proliferation, resolution with treatment and response, and where the delay sits (core
clock or timing cues). Everything runs on `clocksync.paired` (and its public building blocks) plus `common/` and the
builders `paper/common/datasets/{skin_paired_discovery,skin_paired_heldout,skin_paired_treatment,skin_timed_gse205155,psoriasis_sc}`.
No plotting code. Outputs: `paper/results/figure_5/output/psoriasis/` (`.../subsample/` for checks).

```
run.sh [--subsample N]      all scripts in order (full run by default); compute_replication.py first
_shared.py                  tables (discovery, held-out, treatment visits), the A54 pipeline, score and visit helpers,
                            timed-skin cue sets; also imported by psoriasis_details
_singlecell.py              single-cell skin helpers (A58; also used by psoriasis_details)
_cues.py                    GTEx cue reference, disease items, A84 matched sets, A81b timed-skin offsets (A81b; also figure_S09)
```
Every script also accepts `--n-perm` / `--n-boot` (overriding the pre-registered resample sizes, used for the checks).

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_replication.py` | A54: H54.1, H54.2, H54.3 + composition check | lesion vs uninvolved skin of the same visit, nine discovery and seven held-out psoriasis series: r = log tr Cov(beta_lesion) - log tr Cov(beta_uninvolved) (Poisson-corrected, one projector), offset = arg ccorr (h, < 0: delayed); gate on the uninvolved samples; pooled over gated series (r: inverse-variance with DerSimonian-Laird reported, offset: weighted resultant, weights 1/se^2), per-pair label swaps (10,000), pseudo-clock q95, leave-one-gene-out and per-gene shifts, offset in pairs below / above the median keratinocyte-minus-dermal score | `paired_damping/scripts/a54_skin.py` (+ `paired.py`) |
| `compute_dose.py` | A56: H56.1; A73: H73.1 | is the delay larger in pairs with higher lesional IL-17 (IFN-gamma secondary) or proliferation (differentiation loss secondary)? pooled offset of the high half minus the low half (split at the within-series median of the lesional-minus-uninvolved score), labels permuted within series (2,000) | `a56_dose.py`, `a73_proliferation.py` |
| `compute_treatment.py` | A56: H56.2 (+ post hoc interaction) | does the delay resolve after treatment? pooled baseline offset minus pooled post offset (treated arms of GSE117468 and GSE201827 week 12), weights 1/se^2 of the baseline r, baseline/post swapped per patient (2,000); secondaries: placebo, week 52, by drug. Post hoc (labelled): treated minus placebo change, arm labels permuted within series | `a56_treatment.py`, `a56_interaction_posthoc.py` |
| `compute_response.py` | A57: H57.1; A74: H74.1; A73 (treatment, descriptive) | does the change of the offset follow PASI75 response (A57), proliferation reduction within responders (A74), IL-17, IFN-gamma, differentiation restoration? change in one group minus the other, labels permuted across patients (2,000) | `a57_response.py`, `a74_response_prolif.py`, `a73_proliferation.py` (treatment part) |
| `compute_cues.py` | A78: H78.1 | is the delay in the local clock or in the systemic timing cues? reference planes from the timed healthy skin (systemic set: rhythmic in both layers, in phase within 0.5 h; epidermal-output set), pooled offsets per set on the gated psoriasis series with identical per-pair swaps (10,000); core minus systemic < 0 | `a78_local_vs_cues.py` (+ `a67_timed_layers.py` for the sets) |
| `compute_sc_status.py` | A58: H58.1, H58.2 | which cells carry the delay? psoriasis single cells / nuclei (GSE278330, GSE228421, GSE173706; builder `psoriasis_sc`), unit = study \| patient \| visit, compartments `L:<class>` / `N:<class>` on one free host: circular mean over classes of offset(L:k) - offset(N:k) < 0, statuses swapped per unit (2,000), gate; keratinocyte minus stromal-vascular and immune means, per-class offsets and log amplitude ratios (unit bootstrap 500) | `a58_sc.py` (+ `build_psoriasis_sc.py`) |
| `compute_cues_gtex.py` | A81b: H81.2, H81.3 | the GTEx cue reference (S07 `a81b_cue_reference.csv`) in healthy timed skin and psoriasis. H81.2 (GSE205155, 11 subjects): per-subject epidermis-minus-dermis offset of core and cue state (A76 regression), core minus cue > 0, subject bootstrap (10,000). H81.3: the A78 estimator with the cue plane as reference, core minus cue < 0, per-pair swaps (10,000), sets core / cue / locked / systemic | `a82_cues_diseases.py` (`timed_skin`, psoriasis arm), `a81b_cues.py` |

Not in this folder: A62 and A77 (layers), A53 and A55 and A59 (other diseases, damping), A67 and A76 (timed skin), A79
(companion of A78) and the single-cell and cue follow-ups of A53 (TAURUS), A60-A61, A63-A66, A68-A69, A82, A84, A88 are in
`psoriasis_details`. Not rewritten, see Open items below: A80 (`cancer`), A72 (`species`).

## Outputs
`a54_{discovery,heldout}_{per_series,pooled,genes}.csv` (per-series table: gate p, r, se, offset, coherence, amplitude-ratio
bounds; `pooled`: pooled r, offset, p values, H54 verdicts, composition split; `genes`: leave-one-gene-out offset and per-gene
shift), `a56_dose{,_series}.csv`, `a73_{proliferation,series,treatment}.csv`, `a56_treatment{,_series}.csv`,
`a56_interaction_posthoc.csv`, `a57_response.csv`, `a74_response_prolif.csv`, `a78_gene_sets.csv`,
`a78_psoriasis_{sets,verdict}.csv`, `a58_{counts,gate,classes,verdict}.csv`, `a81b_h2_verdict.csv`,
`a81b_timed_skin_subjects_summary.csv` (public GEO data, per subject), `a81b_h3_{sets,verdict}.csv`. Column names follow the legacy files, except: the per-series table has `cohort`
(`series|disease`) besides `series` and `ccorr_abs` is the legacy `ccorr_abs` (modulus of the complex correlation, `x_corr`
repeats it); `offset_low/high_composition_diff_h` is the legacy `offset_low/high_keratin_diff_h`; the pooled table has
`p_offset_two_sided`; the A78 verdict has `core_minus_other_h` / `p_two_other` (legacy `core_minus_output_h`) and `p_two`; the A81b H81.3 verdict names the
cue contrast `core_minus_cue_h` (legacy `core_minus_systemic_h`, the A78 key) and `core_minus_locked_h` / `p_two_locked` (legacy `core_minus_other_h` /
`p_two_other`).
All data are public GEO series.

## Checks
Every script runs with `--subsample 12` (12 random patients per series, 50 resamples; 23 min for the folder, the A54 step dominating because of the pseudo-clock control on all genes). `compute_sc_status.py` ran with `--subsample 20 --n-perm 8 --n-boot 8` (18 units; 3.5 min) and `compute_cues_gtex.py` with
`--subsample 12` (13 min, mostly reading the all-gene series). Full point
estimates against `legacy/desync/validation/paired_damping/output/`, with the full resample sizes except where listed.

**M3 full rerun (2026-10-25, `run.sh --jobs 4`, 1.9 h, peak 34.1 GB in `compute_replication.py`):** every point estimate
equals the archive copy; 17 tables sha256-identical, the A54 pooled, A78 and A81b H81.3 tables differ only in swap p
(registered 10,000 swaps instead of 1,000 / 2,000; no verdict changes); A58 now at its registered 2,000 swaps and 500
bootstraps (rows marked M3):

| readout | result here | legacy | note |
|---|---|---|---|
| A54 held-out: n pairs, r, offset, LOGO max, genes negative, composition split | 129, 0.1556, -0.9324 h, -0.6225 h, 9, -0.944 / -0.624 h | identical | H54.1 FAIL, H54.2 PASS, H54.3 FAIL as legacy |
| A54 held-out p (swap r, offset), M3 | 0.107, 0.0042 (two-sided 0.0080) | 0.107, 0.0042 | 10,000 swaps (the subsample check used 1,000: 0.107, 0.0060) |
| A54 discovery: n pairs, r, offset, LOGO max, genes negative | 465, 0.3585, -1.4994 h, -0.9272 h, 12 | identical | H54.1-3 PASS; M3: swap p 1e-4 at the 10,000-swap floor, as legacy |
| A54 per series (16 series: gate p, r, se, offset) | identical | | same random stream (gate, bootstrap and pseudo draws in the legacy order) |
| A56 H56.1 IL-17 high - low (n 594, 15 series), p | -0.0278 h, 0.472 | identical | FAIL; IFN-gamma -0.2444 h, 0.191 also identical |
| A73 H73.1 proliferation, differentiation loss | -0.1839 h, 0.263; -0.0014 h, 0.507 | identical | FAIL |
| A56 H56.2 treated (116 patients): baseline, post, difference, p | -1.8239, -0.4750, -1.3489 h, 0.0030 | identical | PASS; placebo, week 52, brodalumab, ustekinumab rows identical |
| A56 post hoc treated minus placebo | -0.5028 h, 0.283 | identical | |
| A57 H57.1 responders - non-responders (74 / 31), p | -2.2773 h, 0.0425 | identical | PASS; median split, placebo PASI50, GSE201827 rows identical |
| A74 H74.1 proliferation reduction within responders, p | -3.7501 h, 0.0155 | identical | PASS; IL-17, IFN-gamma, differentiation and all-patient rows identical |
| A73 treatment (descriptive) GSE117468 all / treated, GSE201827 | -4.0512 h (p 0.0005), -3.5244 h (0.0080), -0.4431 h (0.335) | -4.0512, -3.5244, -0.4431 h (0.0005, 0.0085, 0.324) | estimates identical; p differ by Monte Carlo (the legacy stream continues from the dose permutations) |
| A78 H78.1 core minus systemic (15 series, 594 pairs), p | -1.9712 h, 0.0002 (M3, 10,000 swaps) | -1.9712 h, 0.0002 | PASS; core minus epidermal output -0.8645 h identical; the 80 reference genes and amplitudes identical |
| A58 gate (31 units, 14 compartments): clock share, p pairing, p acrophase | 0.0957, 0.0050, 0.0199 | identical | `Dataset.gate`, same seed |
| A58 H58.1 statistic (circular mean of L - N offsets over 7 classes), keratinocyte minus stromal-vascular, minus immune | -0.616 h, p 0.313 (M3, 2,000 swaps); -1.273 h [-7.23, 4.22]; +1.520 h [-9.69, 13.02] (500 unit bootstraps) | -0.612 h, p 0.317; -1.357 h; +1.354 h | H58.1 FAIL as legacy; estimator differs in the second digit (Deviations, free-host estimator) |
| A58 per class L - N offset (T, endothelial, fibroblast, keratinocyte, melanocyte, mural, myeloid), h | -3.04, -0.82, +0.80, -0.75, -0.62, +1.60, -1.49 | -2.93, -0.83, +0.78, -0.84, -0.60, +1.59, -1.46 | log amplitude ratios 1.19, 0.17, 0.25, 0.01, -0.26, 0.24, 1.53 vs 1.17, 0.19, 0.29, -0.02, -0.28, 0.23, 1.50 |
| A81b H81.2 (11 subjects, 14 core / 30 cue genes): core, cue, core minus cue, p upper, interval, subjects core > cue | 1.935 h, 0.434 h, 1.5007 h, 0.0838, [-1.169, 3.200], 8 | identical | FAIL; 10,000 subject bootstraps, same stream |
| A81b H81.3 (15 series, 594 pairs): core minus cue, p lower; core minus locked | -1.6891 h, 0.0003 (M3, 10,000 swaps); -1.1669 h, p two-sided 0.093 | -1.6891 h, 0.0003 (10,000); -1.1669 h | PASS; the legacy swap stream is reproduced |

## Deviations
- **Resample sizes in the checks.** A54 used 1,000 (not 10,000) swaps and A78 2,000 swaps; point estimates do not depend
  on them. The code defaults are the pre-registered sizes.
- **Hand-over files.** `compute_dose.py`, `compute_cues.py` read the A54 per-series tables (gated series, weights
  1/se^2) written by `compute_replication.py` in the same output folder (also with `--subsample`); the legacy scripts read
  the same files from `output/`.
- **Subsample mode.** `thin` keeps N random patients per series (treatment: per file); series with fewer than 5 pairs are
  dropped (10 otherwise).
- **Random streams.** A73 treatment descriptive runs use their own stream (seed 73) in `compute_response.py`; the legacy
  script continued the stream of the dose permutations, so only these p values differ (Monte Carlo).
- **Union tables.** The seven held-out and nine discovery series form one table per script; genes absent from a series
  are all-NaN columns of the union and are dropped per series (`split_pairs`), as the legacy did by reading one file at a time.
- **Free-host estimator of the single-cell scripts (A58).** The legacy `Study.fit` took the loadings of the joint factor model on
  all compartments; `Dataset(anchors=None).fit` reads each compartment on the posterior host state (offsets and amplitudes on the
  full host, reference-free; the coupling is the leave-one-out one). All quantities used here are differences of offsets or log
  amplitudes within the same fit and do not depend on the centring, but they differ in the second digit (table above). The
  pre-registered verdicts do not change. Full resample sizes (2,000 swaps, 500 bootstraps, about 9 s per fit and 4 jobs: > 1 h) were
  not run; `--subsample` runs and the full point estimates are the checks.
- **Single-cell builders** read the processed pseudobulk tables of the legacy raw builders (annotation by Leiden clustering and
  markers is not re-implemented; see `data/INDEX_atlases.md`).
- **A81b arms.** The cue reference is read from `paper/results/figure_3/output/gtex_cues/` (documented file), and the
  contrasts are named `core_minus_cue_h` (see Outputs). The frame check of the reference fails its guard (0.48, post hoc override
  of A81b) as in the legacy.

## Local helpers (not in clocksync)
- `_shared.lesional_all_genes`, `treatment_samples` (with all gene columns and the `pasi` column; the builders keep the 14
  clock genes only, the IL-17, proliferation, composition scores and the pseudo-clock control need other genes).
- `_shared.compare_cue_sets`: `paired.compare_sets` on a union table (see API gaps: genes absent in a cohort are all-NaN
  columns and `compare_sets` does not drop them). Same estimator and swap stream (outputs identical to the legacy).
- `_shared.cosinor_f`, `per_layer`, `reference_sets`: per-gene cosinor with an F test on the timed healthy skin (the
  A78 reference sets). `clocksync.timed.cosinor` has no test.
- `_shared.timed_skin` reads `skin_timed_gse205155.build._read()` (the builder's `build()` keeps the 14 clock genes only).
- `_shared.analyse` is the A54 driver (gate, pooling, LOGO, gene shifts, composition) on `cs.paired` and the public blocks
  (`cohort_states`, `ccorr`, `pooled_offset`); it could move into the package as `paired(..., gene_checks=...)`.
- `_singlecell.py` (A58, shared with figure_S09): status-split tables, `fit_states` (`Dataset(anchors=None).fit`), unit bootstrap and status swaps, wrapped intervals.
- `_cues.py` (A81b, shared with figure_S09): cue-set waveforms from the S07 reference, `timed_core_cue` (per-subject epidermis-dermis
  offsets, the A76 regression), TCGA all-gene paired table (`tcga_paired` builder `_pairs`), the A84 per-gene aggregation and
  tissue-matched sets (copy of the S07 rules), GTEx row reader of A88 / A68 (`gtex_rows`).
- `_shared.compare_cue_sets` also accepts a callable cohort -> set dict (tissue-matched sets, A84) and TCGA pairs (levels normal / tumour).

## API gaps (for clocksync)
- `paired()` pools only upper-tail r and reports `pseudo_q95`/`pseudo_median`: a lower-tail test and `pseudo_q05` (A53, A52)
  and the per-cohort pseudo-clock vectors are missing.
- `paired.compare_sets` should drop genes that are all-NaN in a cohort (union tables).
- A baseline-vs-post change (A56, A57, A73, A74: offset change between two visits of the same patients, label swaps within
  patient or labels permuted across patients) has no public call: `compute_treatment.py` and `compute_response.py` build it
  from `cohort_states`, `ccorr`, `pooled_offset` and `boot_se`.
- A free-host fit that returns the joint-model loadings (legacy `Study.fit`: offsets and amplitudes relative to a reference compartment, no
  posterior-mean regression) and a `fit(..., n_boot=0, gate=False)` that is fast; the free-host fit takes about 9 s for 14 compartments x 31 units,
  which makes 2,000-swap nulls hours long (the legacy was as slow).
- Status-split designs (case / uninvolved compartments of the same unit, one host): swap of the status labels per unit and unit bootstrap
  of `Dataset.fit` quantities (here `_singlecell.swapped`, `boot`), and a per-class contrast `L:k` minus `N:k`.
- `compare_sets` with a per-cohort set dict (tissue-matched sets).

## Open items
A80 (treatment, TCGA) belongs to `cancer`; A72 (mouse timed skin) to `species`. The A82, A84 and A88 follow-ups of
the cue reference, the single-cell follow-ups (A60, A61, A63-A66, A69), A53 H53.1 (TAURUS) and A68 (GTEx skin) are in
`psoriasis_details`.
