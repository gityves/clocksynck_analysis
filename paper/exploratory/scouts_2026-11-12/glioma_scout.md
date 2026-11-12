# Glioma single-cell scout (2026-11-12): metadata and design only

No clock statistic was computed (only UMI sums of the 14-gene panel as one number per cell/pseudobulk, and MGMT UMI). Verified = read from the source page, API or file; inferred = flagged.

## Candidates

| Source | Accession / URL | Assay | Patients | Tumour | Primary / recurrent | Pre-treatment | Malignant + cell types | MGMT status | Time of collection | Size | Licence |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Extended GBmap** (downloaded) | CELLxGENE collection 999f2a15-3d7e-440b-96ae-2c806799c08c, dataset 56c4912d (v 230c8701) | 10x 3' v2/v3, 5' v1 (1.06M cells), Drop-seq, microwell, BD Rhapsody, CEL-seq2, STRT, Smart-seq2 (9k cells); 245k nuclei | 240 donors, 26 studies | glioblastoma (CELLxGENE `disease`); IDH-wt stated for core only | `stage`: unknown 143 donors, Primary 58, Recurrent 26, Recurrent neo-aPD1 13 | not recorded (no TMZ/RT/dexamethasone field); only the neo-aPD1 label | yes: `annotation_level_1` neoplastic/non-neoplastic (SCEVAN CNV, paper), levels 2-3 cell types | `MGMT` column, mixed codes (M, Meth, UM, Un, NM, NT), unknown for 112 of 240 donors | none | 11.2 GB h5ad, raw counts in `raw.X` | not stated in API or paper (CELLxGENE default CC BY 4.0, inferred) |
| Core GBmap | same collection, dataset c888b684 | same, subset | 110 donors | IDH-wt GBM | as above | as above | as above | as above | none | 8.1 GB | idem |
| Brain-tumour atlas, adult diffuse glioma (Powell lab, Med Comm) | CELLxGENE 384e373e-4db1-4605-812f-d9da7cfd246f, dataset 8f3cb9c4 | 10x 3' v3 only | 60 donors (API), 69 tumours (text) | GBM and IDH-mutant | not checked | not checked | not checked | not checked | none expected | 7.3 GB | not checked. Not downloaded, metadata not inspected |
| Abdelfattah 2022 | GSE182109 (GEO text record) | 10x droplet | 18 GBM patients, 44 samples, 201,986 cells (GEO) | GBM | new and recurrent (per GEO summary of the SNU page, partly inferred) | not stated | in GBmap as `Abdelfattah2022` (16 donors, 207,377 cells) | partial in GBmap | none | RAW.tar, not downloaded | GEO |
| Neftel 2019 | GSE131928 | Smart-seq2 (7,930 cells) + 10x (16,201 cells), 28 patients (GEO), 38 on 3CA | GBM adult + paediatric | unknown | untreated not stated | in GBmap as `Neftel2019` (23 donors) | in GBmap partly | none | processed only | GEO |
| Johnson (GSE151506) | GSE151506 | Smart-seq2 (+scRRBS), CD45-enriched, GBM n=7, IDH-mut n=7 | 14 | GBM and IDH-mut | primary | not stated | no (not UMI, immune-depleted) | not checked | none | tpm per sample | GEO; raw on DUOS/EGA |
| Wang 2022, Nat Cancer (downloaded) | GSE174554 (+ EGA EGAS00001004909 raw, restricted) | 10x snRNA, 3'; nuclei | 86 specimens, 49 patients (paper); matrices for 79 samples, 66 patients with Tumor/Normal calls here | IDH-wt GBM (IDH1-mutant excluded) | 36 primary-recurrent pairs + unmatched | recurrent after TMZ + radiotherapy (paper) | only `Tumor` / `Normal` per nucleus, no cell types | not in the GEO files | none | RAW.tar 1.7 GB | CC BY 4.0 (paper) |
| 3CA brain (Weizmann) | https://www.weizmann.ac.il/sites/3CA/brain | mostly Smart-seq2; Wang 2019 10x (10 patients, 16,317 cells); Yuan 2018 microwell (8) | 183 patients, 255,713 cells in total over all brain tumours | GBM, H3K27M, ependymoma, MB, IDH-mut | n/a | n/a | yes (cell types, CNA) | no | none | not downloaded | not stated |
| Wang 2019, Yuan 2018, Neftel, Darmanis, Couturier, Richards, etc. | in GBmap as authors | - | - | - | - | - | - | - | - | - | - |

3CA is a strict subset of what is in GBmap. No candidate records collection or surgery time (verified for GBmap: no time field in `obs`).

## Downloaded
- `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_gbmap/` (extended_gbmap_230c8701.h5ad 11.2 GB, PROVENANCE.md, sha256 6102cc05...).
- `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo_gse174554_wang2022/` (RAW.tar 1.7 GB, Tumor_normal metadata, filelist, PROVENANCE.md). Cap check before download: 408 G used of 500 G cap, 684 G free on the disk. After: about 421 G.

## Mapping of labels to broad classes (GBmap)
malignant = `annotation_level_1` Neoplastic (cell_type "malignant cell"); myeloid = TAM-BDM, TAM-MG, Mono, DC, Neutrophil, Mast; oligodendrocyte lineage = Oligodendrocyte, OPC (non-malignant); T/NK = CD4/CD8, NK; endothelial; mural; neuron; astrocyte; other = B cell, Plasma B, radial glia. Study key = author + cell/nucleus. Patient = `donor_id`; no donor id occurs in two studies. Droplet set = 10x and Drop-seq; second set `umi_all` adds microwell, BD Rhapsody, CEL-seq2, STRT (everything but Smart-seq2). Patient counted once across cell and nucleus runs in the ALL rows.

## Design table, droplet UMI (GBmap, all studies; `output/gbmap_design_table.csv`, set = droplet, study = ALL)

| class | patients >=100 cells | of which >=100 malignant too | median cells / patient | median clock UMI / pseudobulk | median clock UMI / cell |
|---|---|---|---|---|---|
| malignant | 129 | 129 | 929 | 4644 | 2.19 |
| myeloid | 165 | 117 | 1068 | 1903 | 1.13 |
| oligodendrocyte lineage | 69 | 59 | 76 | 1607 | 3.15 |
| T/NK | 76 | 47 | 76 | 923 | 1.36 |
| endothelial | 13 | 9 | 8 | 719 | 4.00 |
| mural | 13 | 12 | 7 | 509 | 2.05 |
| neuron | 11 | 10 | 33 | 7373 | 11.0 |
| astrocyte | 2 | 2 | 6 | 2321 | 16.0 |
| other | 13 | 7 | 8 | 142 | 0.47 |

Median cells per patient is over patients with any cell of the class; the clock columns are over patients with >=100 cells. Of the 129 patients with >=100 malignant cells (`output/gbmap_breadth_droplet.csv`): myeloid >=100 in 117; myeloid and oligodendrocyte in 55; myeloid, oligodendrocyte and T/NK in 20; 8 patients have no other class with >=100 cells, 30 one, 52 two, 39 three or more.

MGMT in malignant cells (droplet): median pseudobulk MGMT UMI 76; 82% of patients (of 129) have >=10 MGMT UMI. Per-study values range from 1 (Wang2020, 3 patients) to 600 (LeBlanc2022); Lee2021 (22 patients, median 119 malignant cells) 59%.

Per-study key rows (patients with >=100 malignant cells; and with >=100 cells of the class): Lee2021 cells 22 (myeloid 22, T/NK 17, oligodendrocyte 1); Wang2021 (CPTAC, nuclei) 18 (myeloid 18, oligodendrocyte 15, T/NK 9, neuron 4); Abdelfattah2022 14 (myeloid 13, oligodendrocyte 8, T/NK 12, endothelial 4, mural 5); Couturier2020 11 (myeloid 9, oligodendrocyte 10); Richards2021 nuclei 10 and cells 6; LeBlanc2022 10; Ruiz2022 nuclei 8 (myeloid 8, oligodendrocyte 8, neuron 3). Full table in the CSV. Nuclei studies have higher clock UMI per cell (Wang2021 malignant 8.4, Ruiz2022 12.2) than cell studies (Lee2021 0.88, Couturier2020 1.06).

Metadata of the 129 patients (`output/gbmap_metadata_summary_ge100malignant_droplet.csv`, MGMT codes M/Meth = methylated; UM/Un/NM = unmethylated is inferred): stage unknown 71, primary 39, recurrent 12, recurrent neo-aPD1 8; MGMT unknown 61, unmethylated 42, methylated 27. So MGMT status is usable in about half (69 patients), clock-state vs MGMT-promoter-status as well as expression.

## Wang 2022 (GSE174554), tumour vs normal only (`output/wang2022_design_table.csv`)

| class | patients with cells | >=100 cells | >=100 and >=100 tumour | median cells | median clock UMI / pseudobulk | median clock UMI / cell |
|---|---|---|---|---|---|---|
| Tumor | 66 | 61 | 61 | 610 | 5299 | 4.67 |
| Normal | 66 | 56 | 51 | 578 | 2291 | 1.83 |

Tumour MGMT: median pseudobulk 48 UMI, 79% of patients >=10 UMI. "Normal" is not split into cell types; the supplement offers none, so a class-wise design needs our own marker annotation (not done).

## Cannot verify / problems
- Licence of the GBmap data and metadata of the 60-donor Med Comm atlas not checked; Johnson 2021 identification of GSE151506 is inferred from its title and design; `Johnson2020` in GBmap is a different 10x set (not traced).
- GBmap MGMT codes are heterogeneous; 10 donors have more than one MGMT label across cells (treated as conflict). Abdelfattah: 16 donors in GBmap versus 18 patients on GEO.
- `is_primary_data` is False for Ruiz2022 and Darmanis2017 cells (duplicates of cells in the core dataset on CELLxGENE); they are kept once here as the extended file contains them once.
- Several GBmap studies are sorted/enriched (e.g. Mathewson 2021 CD45+, Xie, Pombo, Ravi: almost no malignant cells; Lee2021 median 119 malignant cells per patient), so cell-class proportions are not tumour composition.
- Wang 2022: SF8963 matrix is malformed (skipped); 13 samples have >50% of nuclei without a Tumor/Normal call (barcode keys do not match, 31,164 of 228,319 nuclei unmatched, mostly mouse/slice-culture samples and these 13); 'v2' samples taken as the same patient (inferred); the 9 Wang2019 patients of GBmap have SF ids that are in the metadata but not in the matrices of the tar, so there is no double counting between the two sources (checked on ids).
- Whether recurrent/treated status is predictable from `stage` is limited: 143 of 240 donors unknown.
- No time of day exists in any source (GBmap obs verified; Wang not in the GEO files).
