# Bulk, paired and non-human datasets: builders

One folder per dataset, `build.py` with `build() -> (pb, donors)` (pseudobulk tables: donor, compartment, [stratum], depth,
14 clock genes) or `(samples, units)` for paired designs (one row per sample: `unit`, `condition`, `cohort`, `value_type`,
`depth`, clock genes in `clocksync` symbols). All build from the processed files under `$SCCC_DATA`; nothing is saved to
`$SCCC_DATA/pseudobulk` (no `save` calls), no `from_legacy` (every table rebuilds in seconds to a minute). The module docstring of each
builder gives source, accession, the legacy builder replaced and the choices. `_bulk_common/bulk_tables.py` holds the shared reading helpers
and `_bulk_common/check_legacy.py` the equivalence checks (`python -I check_legacy.py [skin|tcga|timed_skin|mouse_ageing|blood|breast|baboon_sheep|mouse_skin]`).

Paired tables are used as `cs.paired(samples, unit="unit", condition="condition", levels=LEVELS, cohort="cohort",
counts="value_type", pseudo_means=gene_means())` (mixes array and count cohorts; `gene_means()` is the legacy `match_means`
input of the pseudo-clock control). `depth` is NaN for arrays. Array values are log2 intensities as in the processed files.

| folder | source | units | compartments / conditions | used by | check | notes |
|---|---|---|---|---|---|---|
| `skin_paired_discovery` | GEO, `skin_lesional/processed` (17 series: psoriasis, AD) | patient (1,568 samples, 784 units) | nonlesional / lesional; cohort `series\|disease` | A53, A54 (`psoriasis_discovery`: 9 series), A56-A57 baselines, A78-A84, A88 | PASS: 18 cohorts (pairs, clock genes, depth) vs `a53_skin.series`, 9 series vs `a54_skin.load`, 18 `match_means` | GSE121212 holds two diseases |
| `skin_paired_heldout` | GEO, `skin_lesional/processed_heldout` (7 psoriasis series) | patient (282 samples) | nonlesional / lesional | A54 replication, A78-A84 | PASS: 7 cohorts vs `a54_skin.load("heldout")`, 7 means (GSE13355 and GSE153007 use GSE30999 reference means) | |
| `skin_paired_other` | GEO, `processed_other` (14 series) + `processed_other_linear` (3, column `source`) | patient (500 samples) | nonlesional / lesional; 10 conditions | A59, A79, A82, A84, A88 | PASS: 17 cohorts vs `a59_other.load` (incl. linear), 17 means (GSE207744 reference means) | linear series are outside the A59 rule (descriptive) |
| `skin_paired_treatment` | GEO, `processed_treatment` (GSE117468 w12, GSE201827 w12 and w52) | patient (844 samples) | visit baseline / post x nonlesional / lesional; `arm`, `treatment` | A56 (H56.2, post hoc interaction), A57, A74, A80 | PASS: 3 files, all four samples per patient vs `a56_treatment.visits` | unit = patient; filter `visit` for a lesional-vs-uninvolved pair |
| `skin_layers` | GEO Visium (GSE206391, GSE202011) and GeoMx (GSE310145), `skin_layers/processed` | patient | epidermis / dermis layers x nonlesional / lesional; `layer_pairs()` aggregates | A62, A77 | PASS: epidermis and dermis, 3 series vs `a62_layers.load` (pairs, counts, depth) | `build()` is unaggregated (GeoMx has several segments per patient) |
| `gut_bulk_ibd` | GEO, `ibd_bulk/processed` (6 CD / UC series) | patient (1,200 samples) | non-inflamed (`nonlesional`) / inflamed (`lesional`) | A55, A79, A82, A84, A88 | PASS: 6 series, CD and UC splits, 6 means vs `a55_ibd.load` | GSE179285 log2 ratio uses GSE193677 means |
| `tcga_paired` | recount3 G029, `tcga/recount3` (10 cancers) | patient (BLCA 19 ... BRCA 112 pairs) | normal / tumour; clinical fields; `tumours()` unpaired | A52, A80, A83, A84, A87, A97, A88 | PASS: 10 cancers (counts, depth, `match` means) vs `a52_tcga.pairs` | counts scaled as recount3 (4e7 / AUC), duplicated symbols summed |
| `skin_timed_gse205155` | GEO GSE205155 (Agilent, 11 subjects x 7 times) | subject x time | epidermis / dermis; donors: subject, clock_h | A37, A67, A76, A78, A82 | PASS vs `skin_paired.table` (counts, depth, stratum, clock_h) | array turned into pseudo-counts (30 M reads) |
| `mouse_skin_timed` | GEO, `mouse_skin_timed/timed_*.parquet` (8 series) | sample | stem / whole_epidermis, `zt` | A72, A78 | PASS vs the parquets read by `a72_mouse_timed` | mouse symbols already mapped to the panel |
| `baboon` | GEO GSE98965 (FPKM) + ENA counts + Ensembl 90 | 12 animals (ZT0-22) | 63 of 64 tissues | A36, A40, A41, A42 | PASS: exact vs `baboon/output/clock_counts.csv` (kept tissues) | `kept_only=False` returns all 64 |
| `sheep` | SheepGTEx v0 (counts, Table S1) | 343 animals | 45 tissues / cell types, 38 BioProjects | A34, A41, A42 | PASS: exact vs `farm_gtex/output/sheep_clock.csv` | 124 animals have several strata (see gaps); about 30 s |
| `mouse_ageing_wolff` | GEO GSE201207 | 36 age x CT pools | 6 tissues | A48 | PASS vs `a48_wolff.load` | `tissue_qc`, `drop_units` provided |
| `mouse_ageing_deota` | GEO GSE190389 | 48 mice (arm x ZT) | 22 tissues | A49 | PASS vs `a49_deota.load` | stratum = arm |
| `blood_timed_gse39445` | GEO GSE39445 (array, 438 samples) | sample / block (subject x protocol) | whole blood | A32 (external), A35 | PASS vs `blood_clock.load` (log-scale matrix, clock time, blocks) | `gene_matrix()` for marker scores; ~45 s |
| `breast_timed_gse233242` | GEO GSE233242 (TPM) + SI Table S1 (`resection_times.csv`) | 42 patients (43 with `join_typo`) | normal / tumour | V25.5, E4 | PASS: 4 variants (N = 10 M, 30 M; typo on/off) vs `v25_5_breast_timed.pseudobulk` | times are recorded hospital hours |

## Not built here (and why)
- `bulk_clock` (A29), `cross_species` (A42), `farm_gtex/human_truth`: inputs are stored outputs (A17 couplings, GTEx gains) and the
  baboon / sheep tables above, plus GTEx (`common/gtex.py`); no other raw input.
- `extended_clock`: E4 is `breast_timed_gse233242`; the D10 atlases are single-cell (other builders); GTEx halves are `common/gtex.py`.
- `skin_invitro` (A75): no deposited transcriptome, only transcribed periods (`periods.csv`).
- Single-cell status-split skin (A58-A69), TAURUS (A53 gut), `a85` Census atlases, `a86` CNV pseudobulks, `celltype_truth` colon chrono-atlas:
  single-cell, belong to the cell-level builders.
- Pig and cattle (FarmGTEx): infeasible (no animal key), as in the legacy.

## Whole-genome GTEx tables
`paired_damping` a81-a102 and `extended_clock/screen` use genome-wide GTEx counts (per-gene cosinor, host regression), `a52`-family
validations use GTEx tissue pairs (`v52_gtex`): `common/gtex.load(genes="all")` is needed, nothing built here.
