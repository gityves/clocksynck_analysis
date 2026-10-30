# Scout: single-cell breast tumours with ER status (2026-10-31)
Metadata only. No clock-gene value, summary or correlation was computed or read (cell-count columns only). Scripts in `scripts/`; tables `wu_patient_x_class_ER.csv`, `wu_er_summary.txt`, `pal_summary.txt`.
Note: the archive cell tables (`cells/cancer/*.csv.gz`) hold clock-gene columns; only `resection/kind/compartment/n_cells` were loaded.

## Ranked recommendation
| rank | study | ER+ pts | ER- pts | malignant call | stroma / endothelium | verdict |
|---|---|---|---|---|---|---|
| 1 | Wu 2021, GSE176078 (downloaded) | 14 (12 ER+, 2 HER2+/ER+) | 12 (9 TNBC, 3 HER2+) | author `Cancer Epithelial` (CNV + scSubtype) | CAFs, Endothelial, PVL per patient | best available; power limited |
| 2 | Pal 2021, GSE161529 | 18 ER+ (+1 PR+) | 8 TNBC (4 BRCA1) + 6 HER2+ | inferCNV cluster blocks, Figshare annotation (not fetched) | unverified until annotation fetched | larger arms, needs Figshare Seurat .rds |
| 3 | Bassez 2021, Nat Med | HR+ | TNBC | cancer cells annotated (unverified) | unverified | 40 patients, biopsies; data via VIB/EGA, access unverified |
| - | ERP130387 (Tietscher 2023) | ~12 luminal | 1 TNBC | none usable | rich | no ER contrast |
| - | SRP273482 (GSE155109) | unknown | unknown | sorted stroma/EC, <=417 epithelial cells | yes | unsuitable |
| - | SRP378451, SRP429995 | HR+ only | none | - | - | no ER- arm |
| - | SRP537773 (GSE279219), SRP514714 | none | TNBC only | - | - | no ER+ arm |
| - | HBCA / Kumar (SRP445031, atlas/breast_kumar) | - | - | - | - | normal breast (Reed 2024), not tumour |

## 1. Wu 2021 (GSE176078, PMID 34493872), downloaded
Series: 26 primary tumours, "11 ER+, 5 HER2+, 10 TNBC" (GEO overall design). Raw counts in `count_matrix_sparse.mtx`; raw reads on EGA (GSE page). Size 533 MB tar.gz; data at `$SCCC_DATA/sc_breast_er/GSE176078/` (see PROVENANCE.md, sha256).
ER status from GEO `clinical_subtype`: ER+ 12 (includes CID3963, TNBC in metadata.csv), HER2+/ER+ 2, HER2+ 3, TNBC 9 (26 total). Cells per class from `metadata.csv` `celltype_major`:

| arm (GEO clinical_subtype) | n pts | median Cancer Epi | median CAFs | median Endothelial | n>=50 Cancer / CAF / Endo | n with all three >=50 |
|---|---|---|---|---|---|---|
| ER+ (incl. HER2+/ER+) | 14 | 241 | 132 | 184 | 11 / 9 / 13 | 6 (7 at >=30) |
| ER- (TNBC + HER2+) | 12 | 853 | 195 | 135 | 9 / 9 / 9 | 5 (6 at >=30) |
By author label `subtype` (ER+/HER2+/TNBC): ER+ 11, HER2+ 5, TNBC 10; n>=50 Cancer 9/3/8, CAF 7/4/7, Endo 10/4/8; all three 5/2/4.
Caveats:
- Discrepancy: CID3963 is TNBC in `metadata.csv` but ER+ in GEO; CID3586 and CID4066 are HER2+ in metadata but HER2+/ER+ in GEO. The ER label must be fixed (source: paper supplementary) before registration. Not verified against the paper table.
- Several patients have 0 author `Cancer Epithelial` cells (e.g. CID3586, CID3838, CID3946, CID4040, CID4398, CID44041; malignant cells labelled Normal Epithelial or absent), so effective n per arm is about 9-11, and only 5-6 per arm have malignant, CAF and endothelial all >=50.
- Malignant cells are one cluster-level author call per patient; CNV use per sample not verified here.
- Per-class counts are in `wu_patient_x_class_ER.csv` (all major classes). Also Normal Epithelial, PVL, immune available.
- Census (CELLxGENE) holds this atlas split by cell type (collection "A single-cell and spatially resolved atlas of human breast cancers"); not needed.

## 2. Pal 2021 (GSE161529, PMID 33950524; EMBO J, doi 10.15252/embj.2020104063)
GEO design: 69 profiles, 52 patients; tumour samples (Total cells): ER+ tumour 18, PR+ tumour 1, HER2+ 6, TNBC 4, TNBC BRCA1 4; plus 7 ER+ lymph-node metastases, 11 normal epithelial, 13 normal total, 4 BRCA1 pre-neoplastic (`pal_summary.txt`; unique patients by label). Counts: GSE161529_RAW.tar 2.2 GB (10x mtx per sample, not downloaded). Raw reads not deposited (privacy).
Annotation: Chen et al. 2022 Sci Data (doi 10.1038/s41597-022-01236-2) deposit Seurat .rds objects on Figshare (per subtype, epithelial and non-epithelial separately, ER+/HER2/TNBC) and inferCNV block annotation marking tumour vs normal clusters; so malignant calls are cluster-level and CAF/endothelial labels exist, but per-patient cell counts could not be tabulated without downloading those .rds files (pure-python reader needed). Figshare URL and sizes: unverified.
Advantage: ER+ 18 vs TNBC 8 (14 with HER2+, the ER-negative arm if HER2+ counted as ER-; per-tumour ER status for HER2+ not verified). Fewer cells per tumour (median ~5-6k per sample) than Wu.

## 3. Bassez 2021 (Nat Med 27:820; PMID 33958794)
Abstract: HR+ and TNBC, paired pre/on-treatment biopsies, n = 29 treatment-naive anti-PD1 + 11 after neoadjuvant chemotherapy (40 patients). Open access PMC full text not found; per-arm n, malignant calling, access route (EGA raw, VIB portal processed) all unverified. Biopsy samples; likely few stromal cells. Pre-treatment samples only would give one sample per patient.

## 4. Archive breast studies
`scbc_full/` contains 20 studies; breast tumour ones identified from `inventory/inventory_samples.csv` (cancer_type, subtype):
- ERP130387 = PRJEB46175 = E-MTAB-10607, Tietscher et al. Nat Commun 2023 (doi 10.1038/s41467-022-35238-w): 14 tumours selected by exhausted vs non-exhausted immune environment; archive has 12. Paper text: "mostly luminal ... a single triple-negative and a single luminal B-Her2+"; ER+/ER- per patient not in ENA/SDRF (checked E-MTAB-10607 SDRF: age, sex, immune environment only). No usable ER contrast. Archive cells (malignant `epithelial` file / host file): 12 patients, malignant 171-5971 cells, fibroblast 162-2506, endothelial 0-1795 (per patient table from `cells/cancer`).
- SRP273482 = GSE155109 (endothelial atlas, healthy and malignant breast): 9 patients, sorted/enriched stroma and EC (labels stroma/TEC/pEC), 14 tumour samples; epithelial-compartment cells per tumour 8-417 (adjacent normal also present). ER status not in inventory; paper not checked. Unsuitable for malignant-cell contrast.
- Other breast tumour studies in `inventory_samples.csv` (not in scbc_full): SRP378451 (GSE205472, 6 pts HR+), SRP429995 (GSE228499, 9 pts HR+/HER2-), SRP537773 (GSE279219, 12 pts TNBC), SRP514714 (6 pts TNBC metastases), SRP572743 (GSE292824, 4 tumours). Wu and Pal are not in the inventory.
- How the archive identified malignant cells (A86/CNV): not re-read; the archive `cnv/` folder has per-study CSVs (no ERP130387 / SRP273482 files listed in the first rows seen). Unverified, to read before registration.

## Blockers and open items
1. ER label reconciliation for Wu (CID3963, HER2+/ER+ cases) needs the paper's supplementary table.
2. Power: about 9-11 informative patients per arm in Wu; coupling vs host needs the host (CAF + endothelium) pseudobulk with enough cells, only 5-6 patients per arm qualify at >=50 cells in all three classes. Pal adds ~26 tumours if annotation is fetched.
3. Pool Wu + Pal (+ Bassez) for ~40-50 tumours is the realistic design; check Wu and Pal are independent patients (separate cohorts, Sydney vs Melbourne; not verified).
4. mtx for Wu not yet read (2.4 GB, counts only after registration).
