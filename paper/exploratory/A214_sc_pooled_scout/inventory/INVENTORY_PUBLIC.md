# A214 scout: public single-cell tumour sources (3CA, CELLxGENE Census), 2026-11-06

Scout only; no clock statistic computed, nothing committed. Table: `inventory_public.csv` (29 rows). Census version `2025-11-08` (stable LTS, schema 2.4.0), metadata queries only.

## 3CA
- Portal (weizmann.ac.il/sites/3CA, 124 datasets, >40 cancer types, 2,836 samples per Tyler 2025): **unreachable from this machine** (connection reset on every attempt; WebFetch ECONNRESET; web.archive.org blocked). The per-study table, formats and sizes could not be retrieved. Only secondary statements are available (bioRxiv 2024.10.11.617836 and the Methods page): cells < 1,000 genes removed, malignant calls confirmed by inferCNA-style CNA, non-malignant by markers, file-name suffix says UMI (10x) vs TPM (Smart-seq). Unverified by me.
- Zenodo 5834545 "3CA Data" (ORCESTRA export, v1, 2022-01-10, CC0, 89 .rds files, 10.3 GB; PDAC_1 and PDAC_2 identical md5) is reachable. Downloaded and inspected 10 files (249 MB, `$SCCC_DATA/3ca/zenodo5834545_<study>/`, PROVENANCE.md with sha256), parsed with a pure-python RDS reader (`python -I`, nothing executed; R `readRDS` was blocked by the permission classifier).
  - Format: `SummarizedExperiment`, one assay `data` (dgTMatrix), colData = cell IDs only, **no malignant label, no patient, no CNA, no gene symbols** (Dimnames NULL, rowData empty; "ARNTL" string absent from the objects). Median summed clock counts therefore not computable (no gene identity).
  - Matrix: integer counts in PDAC_4, CRC, BRCA_3, Prostate_3, NeuroBlastoma_5_10X, Lung_5_Patients, BRCA_6, Ovarian_1; non-integer TPM in CRC_3 and BRCA_4 (plate IDs). Two integer files have max values > 1e6 (not UMI-like).
  - Verdict: v1 Zenodo files are unusable for this project. Usable 3CA data would have to come from the v2 portal (unverified) or from the original GEO deposits.

## CELLxGENE Census (obs summaries; `cell_type == 'malignant cell'`, primary data)
1.94 M malignant cells in 41 dataset entries, 17 distinct primary datasets. "Normal epithelium" = regex on cell_type (epithelial, colonocyte, hepatocyte, club, AT2, ...). Donors passing >= 100 malignant and >= 100 epithelial (Census labels only; lineage match checked by tissue):

| dataset (collection) | assay | donors passing | lineage-matched | h5ad size | on disk |
|---|---|---|---|---|---|
| 829a3cd1 Moorman 2025 CRC metastasis (doi 10.1038/s41586-024-08560-0 erratum link) | 10x 3' v3 | 14 | yes (colonocytes/secretory, colon tumour samples) | 0.75 GB | no |
| 1e6a6ef9 LuCA (Salcher 2022, Cancer Cell) | 10x v1-v3, BD, GEXSCOPE, inDrop, few Smart-seq2 | 41 | yes (AT2, club, ciliated) | 17.6 GB (SOMA subset smaller) | partial (re-aggregates e.g. GSE131907) |
| b617ee1b multi-tissue TME atlas (Nat Commun 2024) | 10x 3' v2, 5' v2 | 26 | mostly no (uveal melanoma, ovary "epithelial cell") | 2.6 GB | unresolved |
| 34f5307e / 6c87755e HTAPP breast mets | 10x v2/v3 (nuclei / cells) | 10 / 4 | no (hepatocytes, keratinocytes) | n/m | no |
| e3ed2ba4 CRLM perfusion | BD Rhapsody | 2 | no (hepatocytes) | 0.62 GB | no |
| 1e4214ce FFPE Flex lung | 10x Flex (probe-based) | 1 | yes | n/m | no |

Other datasets (GBmap, breast atlases, pleural effusion, lymphomas, ICB cancer-cell-only, COVID) have no non-malignant epithelium or are epithelium-only. Median raw_sum per cell: Moorman malignant 14.8k, epithelium 13.8k; LuCA 3.8k / 3.4k. Census X holds raw counts by design; I did not check integer-ness per dataset. CNA/clone labels are not in Census obs (unknown).

## Overlap with disk
Moorman: absent (GEO accession unverified; A149 text lists "Moorman 2025" as not added: format/pairs unverified). LuCA: overlaps on-disk studies (GSE131907, others unresolved: sub-study column needs the h5ad obs). Census datasets were not mapped to SRP IDs.

## Could not verify
3CA v2 study list/formats/sizes; per-cell CNA in Census; Moorman GEO accession; LuCA constituent accessions and which are already on disk; raw-integer status of Census X; identity of the Zenodo studies.
