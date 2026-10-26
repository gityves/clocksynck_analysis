# psoriasis_details: layers, other diseases, timed skin, specificity of the cue comparison

The secondary and detail analyses of the psoriatic delay (`psoriasis`): the same estimator in other diseases and
in the gut, resolved by tissue layer, and in healthy timed skin where the truth is known. Runs on `clocksync.paired`
(and its public building blocks) plus `common/` and the builders `paper/common/datasets/{skin_paired_discovery,
skin_paired_other,gut_bulk_ibd,skin_layers,skin_timed_gse205155,psoriasis_sc,psoriasis_sc_bins,psoriasis_sc_kc_states,healthy_skin_bins,
healthy_skin_sites,taurus,tcga_paired}` and GTEx through `psoriasis/_cues.py`; the A54 pipeline and the table helpers are
`psoriasis/_shared.py` (imported as `_shared`). No plotting code. Outputs: `paper/results/figure_5/output/psoriasis_details/`
(`.../subsample/` for checks).

```
run.sh [--subsample N]      all scripts in order (full run by default); compute_other_diseases.py before compute_cues_companion.py
_details.py                 table loaders of the other series, the A53 per-cohort damping pipeline
_attribution.py            mixing model of A60 / A61 / A69 (NNLS fractions, gains, predicted delay)
(imports `_singlecell.py` and `_cues.py` of psoriasis: single-cell helpers, GTEx cue reference and disease items)
```
Every script accepts `--n-perm` / `--n-boot` (overriding the pre-registered resample sizes, used for the checks).

## Analyses

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_other_diseases.py` | A53: H53.2; A55: H55.1-3; A59: H59.1 (+ conditions) | A54 estimators (gate, r, offset, pseudo-clock, LOGO, composition) on other series. A53: the damping hypothesis r < 0 in psoriasis (primary) and atopic dermatitis (secondary), 17 discovery series, lower-tail swap p and pseudo q05 (data gave r > 0, which led to A54). A55: inflamed vs non-inflamed gut, six series, CD, UC, same bowel segment; composition genes epithelial vs immune. A59: other skin diseases pooled (H59.1: offset < 0), without mycosis fungoides, each condition, linear-scale series descriptive | `a53_skin.py`, `a55_ibd.py`, `a59_other.py` (+ `a54_skin.analyse`) |
| `compute_cues_companion.py` | A79: H79.1 | specificity of A78: the core-minus-systemic comparison (timed-skin cue plane, identical per-pair swaps) in the gated other-skin series (pooled) and in atopic dermatitis; PASS = no core-minus-systemic < 0 at p < 0.05 | `a79_companion.py` (+ `a78_local_vs_cues.compare`) |
| `compute_layers.py` | A62: H62.1; A77: H77.1 | A62: lesion vs uninvolved offset in the epidermis and dermis (Visium GSE206391, GSE202011; GeoMx GSE310145); one projector per series and layer, pooled complex correlation from cross-moments summed over series, gate on the series-centred uninvolved epidermis, swap null (10,000); psoriasis vulgaris only, GeoMx epidermis segment only, without GeoMx, leave one series out. A77: basal layer vs suprabasal epidermis vs dermis of the same GeoMx sections (6 patients): offset, patient bootstrap (1,000), compression (lesional minus uninvolved) against all 64 patient-level status swaps | `a62_layers.py`, `a77_geomx_layers.py` |
| `compute_timed_layers.py` | A67: H67.1; A76: H76.1 | healthy living skin, epidermis and dermis every 4 h (GSE205155, 11 subjects): A67 per-gene cosinor with subject intercepts, epidermis-dermis acrophase shift of the 14 clock genes (circular mean, subject bootstrap 1,000, circular SD) and genome-wide (descriptive); A76 per-subject epidermis-dermis offset of the clock-plane state against the epidermal differentiation index (Spearman, one-sided permutation 10,000) | `a67_timed_layers.py`, `a76_subject_offsets.py` |
| `compute_taurus_amplitude.py` | A53: H53.1 | tissue-wide clock amplitude of inflamed vs non-inflamed gut biopsies of the same endoscopy (TAURUS, builder `taurus`): unit = patient x timepoint, `I:<class>` / `N:<class>` compartments on one free host, mean over classes (>= 10 units) of log amplitude ratio < 0, statuses swapped per unit (2,000) and below the 5th percentile of 200 expression-matched pseudo-clock sets (all-gene sample x class pseudobulks); gate | `a53_taurus.py` (+ `ibd_taurus/scripts/a45_taurus.py` loaders) |
| `compute_sc_attribution.py` | A60, A61 (descriptive; scaled NNLS post hoc) | which cells could carry the bulk lesion delay: single-cell coupling map (builder `psoriasis_sc`), NNLS deconvolution of four bulk series with whole-skin class signatures, mixing model (predicted bulk delay: full, composition only, within-class only, per-class contributions), gene fingerprint against the A54 per-gene shifts; A61 the same over keratinocyte states (`psoriasis_sc_kc_states`) with the reading rule (state mix / dysregulation); unit bootstrap 500 / 200 | `a60_attribution.py`, `a61_kc_states.py` |
| `compute_differentiation.py` | A63: H63.1; A64, A66 (descriptive); A65: H65.1, H65.2 | the keratinocyte clock along differentiation: healthy donors, offset(kc_d4) - offset(kc_d1) > 0 against bin labels permuted within donor (1,000) (`healthy_skin_bins`); offsets relative to the dermal classes and the layer-balanced fit (A64); compression of the gradient G and of the epidermis-dermis offset E in lesions, statuses swapped per unit (1,000) (`psoriasis_sc_bins`); body sites (`healthy_skin_sites`) | `a63_differentiation.py`, `a64_epi_dermis.py`, `a65_compression.py`, `a66_sites.py` |
| `compute_bin_mix.py` | A69 (descriptive) | composition-only predicted delay from the keratinocyte bin shares of lesional and uninvolved skin with healthy-skin bin states (A63 fit; donor bootstrap 500) | `a69_bin_mix.py` |
| `compute_gtex_skin.py` | A68: H68.1, H68.2 | are differentiated-keratinocyte-rich GTEx skin samples (fast deaths, time of death) phase-advanced? pooled complex slope of the clock-plane state on the differentiation / epithelial score, scores permuted within tissue (2,000); aggregates only | `a68_gtex_skin.py` |
| `compute_cues_diseases.py` | A82 | core clock vs the GTEx timing-cue plane (S07 reference) in atopic dermatitis, other skin diseases, inflamed gut and TCGA: the A78 estimator with the cue plane as reference, core minus cue, two-sided, Holm over four units; secondaries by condition, CD / UC, cancer | `a82_cues_diseases.py` |
| `compute_cues_matched.py` | A84 | the same with tissue-matched cue sets (colon for IBD and COAD, per-cancer sets for TCGA, skin for psoriasis and atopic dermatitis) built from the 26-tissue selection table; H84 Holm over IBD and TCGA | `a84_matched_cues.py` |
| `compute_cue_context.py` | A88: H88 | core-cue coupling relative to core-locked (rho) per unit, minus the matched GTEx fast-death tissue: surgical (TCGA normals) and endoscopy (IBD non-inflamed) minus outpatient skin < 0, unit bootstrap (10,000) | `a88_cue_context.py` |

## Outputs
`a53_skin_{per_series,pooled}.csv`, `a55_ibd{,_CD,_UC,_samesegment}_{per_series,pooled,genes}.csv`,
`a59_{all,no_mf,linear_descriptive,<condition>}_{per_series,pooled,genes}.csv`, `a79_{other_diseases,atopic_dermatitis}_{sets,verdict}.csv`,
`a79_verdict.csv`, `a62_{layers,series,verdict}.csv`, `a77_{offsets,compression,verdict}.csv`, `a67_{clock_genes,genomewide_rhythmic,verdict}.csv`,
`a76_{verdict,summary}.csv`, `a53_taurus_{counts,gate,classes,nulls,verdict}.csv`,
`a60_{weights,gene_shares,fingerprint}.csv`, `a60_{all,no_GSE278330}_{coupling,prediction}.csv`,
`a61_{all,no_GSE278330}_{coupling,state_offsets,prediction,state_shares,verdict}.csv`, `a63_{healthy,nonlesional,all}_compartments.csv`,
`a63_verdict.csv`, `a64_{offsets,balanced_fit}.csv`, `a65_{all,no_GSE278330}_offsets.csv`, `a65_verdict.csv`, `a66_sites.csv`,
`a68_{gtex_skin,verdict}.csv` (aggregates), `a69_{bin_mix,bin_shares}.csv`, `a82_<unit>_{sets,verdict}.csv`, `a82_summary.csv`,
`a84_set_<tissue>_{cue,locked}.csv`, `a84_sets_info.csv`, `a84_<unit>_{sets,verdict}.csv`, `a84_summary.csv`, `a88_{gtex,units,verdict}.csv`.
Column names follow the legacy files, except: per-series tables carry `cohort` (`series|disease`
for skin) besides `series`; pooled tables list the series accessions (without the disease suffix of the legacy names);
`a53_skin_pooled.csv` has `p_swap_lower` (legacy `p_swap`), `H53_2` (legacy `verdict`) and is sorted by disease;
`a79_*_verdict.csv` has `core_minus_other_h` and `p_two`, `p_two_other` (legacy `core_minus_output_h`); in the `a59_*`
and `a55_*` pooled tables the hypothesis keys are `H59_1/2/3`, `H55_1/2/3` as written by the shared A54 driver (variance
ratio, offset, genes), so the legacy "H59_2" (offset test) is `H59_2` here too. All data are public GEO series (the legacy
`a76` per-subject table is not written) except the GTEx aggregates of A68 and A88 (controlled access: no per-donor rows). The cue
contrasts are named `core_minus_cue_h`, `core_minus_locked_h`, `p_two_locked` (legacy `core_minus_systemic_h`, `core_minus_other_h`,
`p_two_other`, the A78 keys); the A82 / A84 summaries carry Holm p values (`p_holm`, `H82` / `H84`).

## Checks
Every script runs with `--subsample 12` (26 min for the folder; 12 random patients per series; the timed skin and the GeoMx A77 keep all subjects
and only reduce resamples). Full point estimates against `legacy/desync/validation/paired_damping/output/`:

| readout | result here | legacy | note |
|---|---|---|---|
| A53 pooled psoriasis (9 series, 465 pairs) r, se, p_swap lower, pseudo q05 | 0.3611, 0.0494, 1.0, -0.245 | identical | H53.2 FAIL; AD (8 series, 298 pairs) 0.2882, 0.0659, 1.0, -0.124 identical |
| A53 per series (18 cohorts: gate, r, se, pseudo p) | identical | | same random stream (file x disease order, gate, bootstrap, pseudo) |
| A55 all gut (4 gated series, 563 pairs) r, offset, LOGO max, genes negative | -0.0908, -0.3510 h, -0.0300 h, 9 | identical | H55.1 FAIL, H55.2 PASS, H55.3 FAIL |
| A55 CD / UC / same segment: r, offset | -0.0530, -0.5032 h / -0.1543, -0.4938 h / -0.1564, -0.4106 h | identical | |
| A55 offset one-sided p (all / CD / UC / same segment) | 0.0297 / 0.0391 / 0.0451 / 0.0691 | 0.0297 / 0.0391 / 0.0481 / 0.0691 | 10,000 swaps (same pre-registered size); UC differs by Monte Carlo (the swap stream is the same for the other three) |
| A59 all other skin (13 series, 191 pairs) r, offset, LOGO max | 0.2642, -0.3873 h, -0.2708 h | identical | swap p 0.0010 (1,000 swaps; legacy 0.0009), offset p 0.055 (legacy 0.064) |
| A59 without mycosis fungoides, acne, vitiligo, linear (descriptive): r, offset | identical | | p within Monte Carlo (1,000 vs 10,000) |
| A62 epidermis / dermis (18 / 16 pairs): offset, p; gate | -2.514 h, 0.159 / -1.971 h, 0.093; share 0.3838, p 0.005 | -2.514, 0.162 / -1.971, 0.093; identical | H62.1 FAIL; all 10 sensitivities identical (2,000 swaps here) |
| A77 offsets (3 pairs x 2 statuses), compression, H77.1 | identical (bootstrap intervals to 3 digits) | | FAIL; compression p 0.406, 0.281, 0.141 |
| A67 clock genes: 10 rhythmic in both layers, mean shift [95 % interval], circular SD, genome-wide | 1.147 h [0.605, 1.644], 0.819 h, 1,533 genes | identical | H67.1 PASS |
| A76 Spearman(differentiation, offset), p; mean offset | -0.0909, 0.608; 1.935 h | identical | FAIL |
| A79 other diseases: core minus systemic, p, core minus epidermal output | -0.7523 h, 0.221, -1.9396 h | identical estimates; p 0.221 (legacy 0.221) | PASS (no contrast) |
| A79 atopic dermatitis: core minus systemic, p | -0.3388 h, 0.196 | -0.3388 h, 0.181 | 2,000 vs 10,000 swaps |
| A53 H53.1 TAURUS (29 resections, 91 pairs): gate clock share, p; statistic | 0.5724, p 0.005; 0.2194 | 0.5724, 0.005; 0.2196 | statistic: free-host vs joint-model estimator; swaps and pseudo sets run only on the subsample (full: 40 + 40 was started and stopped, 9 s per fit); legacy H53.1 FAIL (p_swap 0.83) |
| A60 deconvolution (units 31): full, composition, within, h | 2.28, 2.84, 0.54 | 2.16, 2.82, 0.51 | bootstrap 10 here vs 500, so intervals not compared; fingerprint Spearman per class identical; coupling levels (e.g. T L 0.05 vs 0.23) differ: leave-one-out coupling vs joint-model coupling |
| A61 (31 units): L - N offset per class (endothelial, fibroblast, kc_basal, kc_cycling, kc_differentiated, mural), h | -0.475, 0.950, 1.434, 1.856, -0.676, 1.530 | -0.584, 0.897, 1.438, 1.848, -0.777, 1.522 | kc states negative 1, reading undetermined: identical (10 vs 500 bootstraps) |
| A63 d4 - d1 (healthy 22 donors / non-lesional 32 / all 54), h; cycling - d1 | 0.806 / 0.740 / 1.018; -0.459 / -0.969 / -0.504 | 0.809 / 0.759 / 1.024; -0.469 / -0.966 / -0.506 | gate share and p identical; H63.1 FAIL as legacy; 10 permutations / bootstraps here |
| A64 offsets of kc bins versus dermal (healthy) | 0.56, 1.03, 1.62, 1.81, 1.83 h | 0.56, 1.03, 1.62, 1.83, 1.84 h | intervals differ (10 bootstraps) |
| A65 (31 units): G_L, G_N, dG; E_L, E_N | 1.416, 3.444, -2.027; 3.738, 3.533 | 1.398, 3.426, -2.029; 3.724, 3.513 | `no_GSE278330` (16 units) G_L 1.196 vs 1.415, E_N 5.73 vs 3.26: small-sample, 10 bootstraps |
| A66 sites (non-volar, 18 donors) kc_d1 - fibroblast | 2.397 h | 2.364 h | volar (4 donors): all 10 bootstraps failed here (300 in legacy), no estimate |
| A69 composition-only delay, all units, healthy phases / all donors | -0.193 / -0.259 h | -0.196 / -0.266 h | 30 vs 500 bootstraps; intervals not compared |
| A68 GTEx skin: H68.1, H68.2; timing consistency (lower leg, suprapubic) | FAIL, FAIL; 0.1551, 0.2120 | identical | aggregates only; full |
| A82 (core minus cue, h): AD, other skin, IBD, TCGA | -0.267, -0.519, +0.225, -0.237 | -0.267, -0.519, +0.225, -0.183 | TCGA pooled differs (10 cancers, 517 pairs; per-cancer causes not investigated); p within Monte Carlo for the others (1,000 here vs 10,000 swaps); all "no offset" |
| A84 matched sets (core minus cue, h): IBD, TCGA | -0.363, -0.826 | -0.363, -0.759 | TCGA as above |
| A88 surgical minus outpatient (12 / 23 units), diff, p lower | 0.6435, 0.997 (0.190, 1.112) | 0.6435, 0.998 (0.189, 1.123) | FAIL as legacy; aggregates only |

## M3 full run (2026-10-26)
`run.sh --jobs 4` (2026-10-24/25: A53 to A84 up to `compute_cues_diseases`, then killed for memory at 30.6 GB after 3 h 24 min); `compute_cues_matched.py` and
`compute_cue_context.py` rerun with `--jobs 1` (17 and 9 min, peak 6.9 / 6.1 GB; `logs/jobs.tsv`). Against the archive copies
(`outputs/paper/figure_S09_psoriasis_details`, reduced resampling) and the legacy full runs:
- **Point estimates identical** to the archive copies in every table (offsets, correlations, gate shares, couplings; `a60_*_coupling` within 2e-7). Resamples
  now at their registered sizes: A59 and A62 / A79 / A82 / A84 swaps 10,000, A60 / A61 / A65 / A69 bootstraps 500 (A61 without GSE278330 200), A63 / A65
  permutations 1,000, A66 300, A53 H53.1 2,000 swaps and 200 pseudo-clock sets.
- **No registered verdict changes.** H53.2 FAIL (psoriasis r 0.3611, pseudo q05 -0.2450; `a53_skin_pooled.csv`), H53.1 FAIL (stat 0.2194, p_swap 0.8286;
  `a53_taurus_verdict.csv`, first full run, legacy 0.2196 / 0.8276), H55.1 FAIL / H55.2 PASS (offset -0.3510 h, p 0.0297; `a55_ibd_pooled.csv`), H59 FAIL
  (offset -0.3873 h, p 0.0638; `a59_all_pooled.csv`), H62.1 FAIL (p 0.1617; `a62_layers.csv`), H63.1 FAIL (p 0.3097; `a63_verdict.csv`), H65.1 / H65.2 FAIL
  (p_dG 0.3027, p_dE 0.5534; `a65_verdict.csv`), H67.1 PASS, H68 FAIL, H76.1 FAIL, H77.1 FAIL, H79.1 PASS (p_lower 0.2214; `a79_verdict.csv`), H82 and H84 no
  offset (Holm p 1.0; `a82_summary.csv`, `a84_summary.csv`), H88 FAIL (diff 0.6435 h, p 0.9976; `a88_verdict.csv`), A61 reading undetermined (all units,
  composition -0.4866 h) / state mix (without GSE278330, -0.9735 h; `a61_*_verdict.csv`, `a61_*_prediction.csv`).
- **p values** equal the legacy full runs where the swap size is the same (A59, A62, A82, A84, A88: identical; A55 UC 0.0451 vs legacy 0.0481, Monte Carlo)
  and the text of the archive RESULTS (A59 p 0.064, A63 p 0.31, A65 p 0.30 / 0.55).
- **One sensitivity row changes side of 0.05:** A59 without mycosis fungoides, offset p 0.0480 (archive copy, PASS) -> 0.0539 (FAIL; `a59_no_mf_pooled.csv`;
  legacy 0.0539). Descriptive discoid lupus (`a59_cutaneous_lupus_discoid_pooled.csv`) variance p 0.0559 -> 0.0360, H59.1 still FAIL.
- **Bootstrap intervals widen** with 500 resamples (the archive copies used 10 to 30): A65 dG [-3.7936, 1.1781] -> [-10.7782, 6.4399]; A63 d4 - d1 over all
  donors [0.4606, 2.2343] -> [-0.2842, 2.1585]; A60 full predicted delay [0.0340, 5.2128] -> [-2.1043, 5.4183] (composition stays above 0: lower bound
  0.6692); A64 healthy kc_d4 vs dermal lower bound 0.9936 -> -0.9689. The registered readouts of A63 and A65 are the swap p, not these intervals; the legacy
  full-run intervals agree with the new ones.
- **A66 volar** (4 donors) remains not estimable (`a66_sites.csv`, n_boot 0); the legacy run had estimates.
- **TCGA pooled core minus cue** differs from the legacy full run (A82 -0.2369 vs -0.1831 h, A84 -0.8256 vs -0.7589 h; `a82_tcga_verdict.csv`,
  `a84_tcga_verdict.csv`); per-cancer rows are identical to legacy, no decision depends on it (p_two 0.9189 / 0.9925).
- New tables without archive copy: `a53_taurus_{classes,counts,nulls,verdict}.csv`.
<!-- src: figure_5/output/psoriasis_details/*.csv (M3, 2026-10-24/26) vs m3_snapshot_backup/2026-10-24 and legacy/desync/validation/paired_damping/output -->
The Checks table above (subsample figures) and the Open item "not rerun in full" are superseded by this run.

## Deviations
- **Resample sizes in the checks.** A53 and A59 used 1,000 swaps, A55 10,000, A62 and A79 2,000 (code defaults: pre-registered); the A55 offset p at 1,000 swaps (0.045 for all gut) was 2.8 Monte Carlo SE from the 10,000-swap value (0.030), so it was rerun at 10,000.
- **A53 pipeline.** `clocksync.paired` pools all cohorts together and tests r > 0; A53 pools per disease and tests r < 0
  with a pseudo-clock q05, so the per-cohort loop is written here from the public blocks (`reference_gate`, `cohort_states`,
  `log_var_ratio`, `boot_se`, `pseudo_waveforms`, `swap_null`) with the legacy random stream (`_details.damping`).
- **Subsample mode.** Units are thinned per series before pairing; series with fewer than 5 pairs (10 otherwise) drop.
- **Series with different gene sets.** Union tables of series have all-NaN columns for genes a series lacks; they are dropped
  per series (`_shared.split_pairs`, `compare_cue_sets`), as the legacy did by reading one file at a time.
- **Free-host estimator of the single-cell and TAURUS scripts (A53 H53.1, A58-A66, A69).** `Dataset(anchors=None).fit` reads each
  compartment on the posterior host state instead of the joint-model loadings of the legacy `Study.fit`; differences of offsets and log
  amplitudes agree to the second digit, coupling levels do not (leave-one-out coupling). Verdicts are unchanged.
- **Resample sizes of the single-cell checks** (10 to 40 resamples, see the table); full sizes (500 bootstraps, 1,000 permutations, 2,000 swaps)
  are hours at 9 s per free-host fit and were not run.
- **Builders.** `psoriasis_sc`, `psoriasis_sc_bins`, `psoriasis_sc_kc_states`, `healthy_skin_bins`, `healthy_skin_sites`, `skin_dermis_zt`, `taurus`
  (see `data/INDEX_atlases.md`) read the processed pseudobulk tables of the legacy builders.

## Local helpers (not in clocksync)
- `_details.damping` (A53, see above); `compute_layers.pooled_corr` (complex correlation from cross-moments summed over
  series: `paired` pools correlations by weights, A62 sums the moments) and its pair-swap; the offsets of two compartments of
  the same units with patient bootstrap and exact status swaps (A77).
- `_shared.cosinor_f`, `per_layer`, `timed_skin` (per-gene cosinor with F test, the A67 shift; `clocksync.timed.cosinor` has
  no test and is per series), the A76 complex regression per subject and layer (`compute_timed_layers.a76`).
- `_shared.compare_cue_sets` (`paired.compare_sets` with absent genes dropped per cohort).
- Builder access: `skin_timed_gse205155.build._read()` (all-gene array matrix), `lesional_all_genes` (all gene columns).

## API gaps (for clocksync)
- `paired`: lower-tail pooled test (`alternative`) and `pseudo_q05` / per-cohort pseudo vectors; pooling per group (disease)
  in one call; a pooled-correlation option from summed cross-moments (small series, A62).
- `compare_sets`: drop genes absent from a cohort (all-NaN columns).
- `timed`: a cosinor with an F test per gene (vectorised over genes, subject intercepts) and a complex regression of a
  clock-plane state on clock time (A76).

## Open items
Not rerun in full (resamples; point estimates are checked above): the A53 H53.1 swaps (2,000) and pseudo-clock sets (200), the A60 / A61 / A65 bootstraps (500),
A63 / A64 permutations (1,000) and the A66 bootstraps (300), the A69 bootstraps (500), the A82 / A84 TCGA arms at 10,000 swaps. A72 (mouse timed skin) belongs to
`species`. Nothing in this package's list needs `legacy/` any more.
