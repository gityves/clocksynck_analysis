# Single-cell atlas and gene-wide builders

One folder per dataset, `build.py` with `build(donors=None) -> (pb, donors)` (donor x compartment pseudobulks: `donor`,
`compartment`, `stratum`, `n_cells`, `depth`, genes) and, where a legacy table exists, `from_legacy()`. The module docstring of
each builder gives source, amendment, A17 target names and the legacy builder replaced. Shared code is in `_atlas/`
(census, organ, cortex, lung, heart, skin, blood, ibd, kong, immune_aging, colon_chrono, genemap, genome, scored, tstates,
heldout, `legacy.py`; `a17_targets.csv` maps every A17 target of `legacy/desync/common_output/a17_registry.csv` (70 rows) to its
dataset and donor columns).

## Checks
`python _atlas/check.py <dataset> [--n 8] [--twins]` compares `from_legacy()` (full legacy table) and `build(donors=slice)`
(random donors plus one case and one control of every contrast column) with the legacy rows of the same donors.
Run 2026-10-21 (`--n 4`): **all 71 datasets that have a legacy table agree exactly** (same rows, `max_abs_diff` 0 over
n_cells, depth and counts, strata equal, shared donor-table columns equal). Rows only in the slice belong to compartments
dropped by the legacy `min_donors` rule, which needs the whole dataset.

Not an exact check (no stored legacy table):
- `*_genome` (all genes) and `natri_tstates` / `hlca_tstates` (T-state compartments) have no `from_legacy`; their clock-gene
  columns are compared with a twin (`--twins`). `hlca_tstates` against `hlca_fibrotic`: exact. `natri_tstates` (non-cycling)
  against `lung_natri_control`: exact on common rows (4 extra rows are compartments the twin drops).
- **`lung_hlca_genome` vs `lung_hlca` differ**: the genome table follows the A4 rules of the gene-rhythm tests (stratum
  `assay|study`, >= 20 cells, healthy plus tumour-adjacent donors), `lung_hlca` the A17 lung rules (stratum `study|assay`,
  other cell filters). Strata strings differ in order and a few rare compartments have 1-20 more cells. The legacy A4
  genome-wide table was never stored, so this is not an equivalence test; the other `*_genome` builders have no twin check run.
- Cancer and species datasets are listed in `INDEX_natexp_cancer.md`, bulk and paired ones in `INDEX_bulk.md`.

## Datasets
| group | folders | amendments | A17 targets |
|---|---|---|---|
| brain (cortex) | `psychad_hbcc`, `psychad_mssm`, `psychad_radc`, `seaad_dlpfc`, `seaad_mtg`, `parkinson`, `retina` | A1-A3, A17 | healthy and disease contrasts (SCZ, BD, AD, DLB, ADNC, PD, age) |
| healthy organs | `organ_adipose`, `organ_bone_marrow`, `organ_breast`, `organ_colon`, `organ_eye_anterior`, `organ_ileum`, `organ_muscle`, `organ_optic_nerve`, `tabula_sapiens` | A9 part 2, A17 | one healthy target each |
| immune aging, gut | `ia_jejunum`, `ia_lung`, `ia_mesenteric_ln`, `ia_spleen`, `ia_thoracic_ln`, `ia_crohn`, `kong_ileum`, `kong_colon`, `ibd_ileum`, `ibd_rectum` | A14-A16 | healthy; `ia_crohn`, `ibd_crohn*` |
| lung | `lung_hlca`, `lung_natri`, `lung_natri_control`, `hlca_fibrotic`, `hlca_tstates`, `natri_tstates`, `hlca_tstates_a7`, `natri_tstates_a7` (A7 panel and stress genes: whole T class plus the other lung classes; A6 stress, A7) | A4-A7, A12 | `lung_hlca`, `lung_ipf_hlca`, `ild_natri`, `hlca_fibrotic`; T states (figure 7) |
| heart, kidney, islets | `heart_global`, `heart_dcm`, `heart_dcm_healthy`, `heart_myocarditis`, `kidney_sc`, `islets`, `islets_a96` | A8, A13, A96 | `heart`, `dcm`, `dcm_healthy`, `myocarditis`, `kidney_sc`, `ckd_sc`, `islets`, `islets_t2d` |
| blood, skin, endometrium | `blood_onek1k`, `blood_aida`, `blood_lupus`, `skin_ad`, `skin_ad_nonlesional`, `skin_ssc`, `endometrium` | A6, A9, A11, A20 | `blood_onek1k`, `blood_aida`, `lupus`, `skin_ad*`, `skin_ssc*`, `endometrium*` |
| immune-damping inputs (clock genes + 359 score genes) | `*_scored` (9 datasets), `kong_cd`, `kong_cd_scored`, `masld_pooled`, `masld_pooled_scored`, `masld_pooled_a96`, `periodontitis`, `periodontitis_scored`, `taurus`, `taurus_scored` | A20, A43-A45, A51, A96 | |
| skin single cell (psoriasis, healthy differentiation bins, body sites) | `psoriasis_sc`, `psoriasis_sc_bins`, `psoriasis_sc_kc_states`, `healthy_skin_bins`, `healthy_skin_sites` (shared reader `_skin_sc/skin_sc.py`) | A58, A60, A61, A63-A66, A69 | |
| mouse colon (known ZT) | `colon_chrono`, `colon_chrono_genome`, `colon_chrono_raw` | A21, A31 | |
| mouse, known ZT (V13) | `colon_prepub` (pre-publication colon h5ad, 8 animals x 2 regions, regenerating colon), `skin_dermis_zt` (GSE223109, 18 animals; `build(parts=k)` splits the fibroblasts) | A19 part 4 (V13) | |
| gene-rhythm maps (all genes) | `gene_map_dlpfc`, `gene_map_lung`, `lung_hlca_genome`, `heart_global_genome`, `organ_adipose_genome`, `organ_breast_genome`, `organ_colon_genome`, `organ_muscle_genome`, `liver_xu_genome`, `breast_hbca_genome`, `skin_ganier_genome` | A27, A31b, A39, A50 | |

Skin single-cell builders (`_skin_sc`): they read the processed pseudobulk tables of the legacy raw builders
(`$SCCC_DATA/psoriasis_sc/pseudobulk{,_bins_both,_kc_states}.csv`, `healthy_skin/pseudobulk_bins{,_sites}.csv`: QC, the A58
marker-and-Leiden class rule with random_state 0, keratinocyte differentiation quartile bins, raw count sums per sample x class).
Those annotation steps (GEO / Census / ArrayExpress downloads, scanpy clustering, hours of CPU) are not re-implemented; the
builders only standardise the table (`donor` = unit key, `compartment` = class, `status`, `stratum` = study, `n_cells`, `depth`,
14 clock genes) and `skin_sc.allgenes()` opens the all-gene count sums (`pseudobulk_all.h5ad`, `pseudobulk_kc_states_all.h5ad`)
used by the A60 / A61 class signatures. Check: `python -c` loading all five builders (rows 363, 355, 131, 514, 174; units 37,
36, 36, 57, 22) and the compute scripts' point estimates against the legacy outputs (psoriasis and
psoriasis_details READMEs).

## Cue-following units (A85)
| folder | source | table | used by | notes |
|---|---|---|---|---|
| `cue_units` | `$SCCC_DATA/atlas_wg/<atlas>/pseudobulk.h5ad` (9 genome-wide Census atlases) | unit (donor x tissue) x coarse class (`COARSE` of the legacy `run_atlas`) pseudobulks of any gene list, majority assay, healthy donors, >= 20 cells, classes with >= 20 units | A85 (`cancer/compute_cue_following.py`) | `table(name, genes)`; no stored table; var names of the files are Ensembl ids, genes are matched on `feature_name` |
