# Inventory GSE266919 (Zhang et al., Cancer Cell 2025, TNBC chemo +/- atezolizumab)

Inventory only (no clock or rhythm analysis). Data: `$SCCC_DATA/tnbc_zhang2025/raw/` (1.36 GB), `PROVENANCE.md` there (URL, date, sha256).
Script: `inventory_GSE266919.py` (+ `rds_parse.py`, pure-python R-serialization reader, `python -I`, nothing executed, R not used since
R 4.1.2 is exposed to the readRDS code-execution CVE-2024-27322). Tables: `output/inventory_tables/` (`summary.json`, per-file sample x cluster counts,
`cells_per_sample_by_file.tsv`, gene lists). All numbers below come from these outputs.

## Files on GEO (series suppl, no per-sample files; GSM entries list "NONE")
| file | size | content |
|---|---|---|
| GSE266919_Bcell.rds.gz | 291 MB | SingleCellExperiment, B/plasma |
| GSE266919_CD4Tcell.rds.gz | 447 MB | SCE, CD4 T |
| GSE266919_CD8Tcell.rds.gz | 363 MB | SCE, CD8 T |
| GSE266919_Myeloid.rds.gz | 255 MB | SCE, myeloid |
| GSE266919_NKcell.rds.gz | 51 MB | SCE, NK |
| GSE266919_Mm_bulk.counts.xls.gz | 7 MB | mouse bulk RNA-seq (OLE2 .xls), not used |

**Only immune compartments are deposited.** No cancer epithelial, fibroblast/CAF, endothelial, pericyte/smooth muscle object exists on GEO
(the series' 66 "Total cells" GSM and 51 "Immune cells" GSM entries have no files). The paper's earlier "previously published" TNBC data
(not in this series) may contain the other compartments; not checked. Raw reads: the submitter states they will go to China GSA (not
requested, no account created). Series: 78 human tumour biopsies from 44 patients per GEO overall design (matches the data below); 117 GSM
entries in total including mouse bulk RNA-seq.

## Format
- Each .rds.gz is gzip inside gzip; an SCE with one assay, `counts`, a dense double matrix (genes x cells), **all values integer, min 0**
  (max 41664 B, 19740 CD4, 7674 CD8, 19164 myeloid, 2017 NK): raw UMI counts, no normalised layer. Barcodes in the matrix match `Barcode`.
- Genes (HGNC symbols plus Ensembl-clone names, GRCh38, Cell Ranger 3.1.0): 16,281 to 21,035 per object, 22,815 in the union (genes filtered per
  compartment, so the gene sets differ between objects; B and myeloid share 19,138).

| object | cells | genes | samples | patients |
|---|---|---|---|---|
| Bcell | 87,595 | 20,151 | 77 | 44 |
| CD4Tcell | 157,219 | 20,514 | 76 | 43 |
| CD8Tcell | 127,782 | 19,701 | 76 | 43 |
| Myeloid | 56,180 | 21,035 | 78 | 44 |
| NKcell | 18,718 | 16,281 | 73 | 43 |
| total | 447,494 | | 78 (union) | 44 |

## Samples, patients, timepoints
- 78 distinct samples = 44 patients; no sample is a repeat library. The preprint's "76 samples" equals the CD4/CD8 object sample count (76),
  not the series total (78).
- Samples are **repeated biopsies of the same patients**: 31 patients have 2 samples (Pre + Post), 11 have Pre only, 1 has Post only, 1 patient
  (P040) has 4 (breast and lymph node, each Pre and Post). Group: 44 Pre-treatment, 34 Post-treatment (no baseline / on-treatment / surgery
  split: only Pre vs Post; Post timing relative to treatment is not recorded).
- Arms (samples Pre / Post; patients): Nab-PTX+atezolizumab 16/12 (16 patients), Nab-PTX 13/9 (12), PTX+atezolizumab 9/6 (9), PTX 6/7 (7).
- Tissue (samples): breast 44, lymph node 16, chest wall 11, liver 6, lung 1 (breast Pre 23 / Post 21).
- Response (samples): R 53, NR 25; Efficacy: SD 38, PR 36, PD 4 (per-patient values not listed here).

## Cell-level metadata columns (identical in all five objects)
`Barcode`, `Sample` (e.g. P002.Pre), `Patient` (P0xx), `Group` (Pre-/Post-treatment), `Tissue`, `Treatment`, `Efficacy` (SD/PR/PD),
`Response` (R/NR), `n_genes`, `n_counts`, `Subset` (fine), `SubCluster`, `MajorCluster`.
**No collection time, clock time, date, procedure type, age, sex, batch, or time-of-day field exists** (GEO sample characteristics: patient,
tissue, group, treatment, response, cell type only). The sample-level metadata are therefore not usable to date the biopsy.

## Annotation levels
| object | MajorCluster | SubCluster | Subset (fine) |
|---|---|---|---|
| Bcell | Bcell | Bfoc, CD19B, Plasma | 14 |
| CD4Tcell | Tcell | CD4 | 16 (Tn, Tcm, Tem, Tfh x3, Treg x2, Tex, Teff, Tisg, Tstr, Tact) |
| CD8Tcell | Tcell | CD8 | 12 (Tn, Tcm, Tem, Tsem, Teff, Tex x4, MAIT, Tisg, Tstr) |
| Myeloid | Myeloid | Macro, Mono, Mast, cDC, pDC | 20 (13 macrophage, 2 monocyte, mast, DC) |
| NKcell | NKcell | CD16, CD56 | 11 |

The compartment is the object (CD4 vs CD8 split is by object; the other levels are sub-states).

## Cells per sample per compartment (`output/inventory_tables/cells_per_sample_by_file.tsv`; 78 samples)
| compartment | samples with >0 | median (all 78) | min | median among >0 | samples >=20 | >=50 |
|---|---|---|---|---|---|---|
| CD4 T | 76 | 1694.5 | 0 (min >0: 45) | 1728 | 76 | 75 |
| CD8 T | 76 | 1322 | 0 (min >0: 23) | 1351 | 76 | 74 |
| Myeloid | 78 | 465.5 | 4 | 465 | 74 | 74 |
| B | 77 | 617 | 0 (min >0: 4) | 651 | 74 | 64 |
| NK | 73 | 147.5 | 0 (min >0: 19) | 167 | 72 | 67 |
| Fibroblast, endothelial, pericyte, cancer epithelial | not deposited | | | | 0 | 0 |

Samples with >=20 cells in CD4, CD8 and myeloid simultaneously: 73 (73 of 78); with B also: 71. Restricted to breast biopsies (44): CD4 44, CD8 44,
myeloid 41, B 44, NK 42 samples have >=20 cells. Total immune cells per sample: median 4,884 (min 8, max 27,495).

## Consequences for A140 (facts only)
- A donor-internal comparison across compartments is possible only among immune compartments (CD4, CD8, myeloid, B, NK) within each biopsy; there is no
  non-immune (stromal-vascular/epithelial) compartment in this deposit.
- Counts are raw integers, so pseudobulks can be built from them directly.
- Pre/Post biopsies per patient exist for 32 patients, but there is no time-of-day information at sample or patient level.
