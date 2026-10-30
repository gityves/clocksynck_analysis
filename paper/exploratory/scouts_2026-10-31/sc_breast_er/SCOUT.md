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

---
# Part 2 (2026-10-31): ER reconciliation and per-patient cell counts (metadata only)
No expression value, clock gene or clock statistic was read. Cell-class counts and cluster ids only. Scripts: `scripts/wu_reconcile.py`, `rds_parse.py` (pure-python R reader, python -I, copy of A140 with an altrep fix), `pal_meta_cols.py`, `pal_cluster_counts.py`, `pal_patient_table.py`, `pooled_table.py`. Outputs: `wu_reconciled.csv`, `wu_reconcile_summary.txt`, `pal_cluster_counts.csv`, `pal_sample_tumour_microenv.csv`, `pooled_patient_table.csv`, `pooled_summary.txt`.

## 2.1 Wu 2021: ER label reconciled against the paper
Source: Nat Genet Supplementary Table 1 (ESM4 xlsx; ER, PR, HER2 IHC/ISH, subtype by IHC). Paper ids with a suffix map to cell-table ids (CID4404-1 -> CID44041, 4497-1 -> 44971, 4499-1 -> 44991, 4517-1 -> 45171, 4290 -> 4290A, 4530 -> 4530N). All 26 patients matched.
- Paper IHC subtype = GEO `clinical_subtype` for all 26 patients (0 disagreements). Final ER: 14 ER+ (12 ER+, 2 HER2+/ER+), 12 ER- (9 TNBC, 3 HER2+).
- `metadata.csv` `subtype` disagrees for exactly CID3586 and CID4066 (metadata HER2+, paper HER2+/ER+: ER 100% 2-3+, 70% 3+) and CID3963 (metadata TNBC, paper ER+). Resolved to the paper.
- Flag: CID3963 has ER 30% 1+ (weak; treated, probable recurrence), PR 0, HER2 0; ER+ by the paper but near the ER-low range. CID4066 and CID4398 are post neoadjuvant chemotherapy (treated samples: CID3963, 4066, 4398, 4513, 4523; the last two TNBC).
- Final per-patient label and source are in `wu_reconciled.csv` (`final_label`, `final_ER`, `source`).

## 2.2 Pal 2021 (GSE161529): what the annotation provides
- Chen et al. 2022 figshare doi 10.6084/m9.figshare.17058077 (verified through the figshare API; md5 matches). 26 Seurat .rds files (0.03 to 2.0 GB each, all < 3 GB) plus HumanBreast10X.zip (2.42 GB: per-sample 10x counts, R code, `InferCNV-Annotation.txt`).
- meta.data of the objects holds only `orig.ident, nCount_RNA, nFeature_RNA, group` (sample), `integrated_snn_res.*` and `seurat_clusters`. There is no cell-type label: no CAF, endothelial or pericyte column anywhere in the Total objects. Cluster identities are only in the paper figures.
- Malignant call is cluster-level. From the R code (`Tum <- Cluster %in% c(...)`, 1-based, = seurat_clusters 0 and 4 for ER-Total, 0 and 3 for HER2, 0 and 2 for TNBC) and the inferCNV block file (blocks marking clusters 0/4/5/6 ER, 0/3/7 HER2, 0/2 TNBC; the file does not state tumour vs normal per block). I use the `Tum` clusters as malignant. "Micro-environment" = clusters kept in the code's `Sub` set (non-tumour, non-epithelial; mixture of immune, fibroblast, endothelial, pericyte). CAF vs endothelial cannot be split from metadata alone. Splitting needs a marker-gene step on expression (the Sub objects give finer clusters, still unlabelled); not done here by instruction. Objects downloaded: ERTotal, HER2, TNBC (includes BRCA1 TNBC), BRCA1Tum (not parsed, redundant with TNBC), and the zip.
- Objects cover only 13 ER+ tumours (ER-Total; includes ER-0319 "ERlo, PR+") + 6 HER2 + 8 TNBC. The other ER+ samples (ER-0029-7C/9C, 0056, 0064, 0173 in `PairedER`/`TumLN`, and male ER-0068, 0178 in `Male`) have no malignant-cluster definition in the R code or the inferCNV file, so they are not usable without a new CNV or marker step; not downloaded.
- ER status: ER+ from GEO title and Table EV3 (IHC, e.g. ER-0125 60% 2+). ER-0319 is ER-low with PR+ (10% 2+). The six HER2 tumours (Table EV2) report HER2 amplification only; **ER status of the HER2+ arm is not reported in the files read (unverified)**. TNBC from Table EV2 (4 TNBC, 4 BRCA1 TNBC).
- Patient 0114 contributes two tumours: ER-0114-T3 (ER+) and TN-0114-T2 (TNBC), both age 84 (Tables EV2/EV3). It is in both arms; counted once per sample below.

## 2.3 Pooled counts per ER arm (`pooled_summary.txt`, thresholds >= 50, >= 30 cells per class)
| arm | cohort | n | malignant >= 50 | CAF >= 50 | endothelial >= 50 | all three >= 50 (>= 30) | median malignant / CAF / endothelial |
|---|---|---|---|---|---|---|---|
| ER+ | Wu | 14 | 11 | 9 | 13 | 6 (7) | 241 / 132 / 184 |
| ER- (TNBC + HER2+) | Wu | 12 | 9 | 9 | 9 | 5 (6) | 853 / 195 / 135 |
| ER+ (incl. 1 ER-low/PR+) | Pal | 13 | 13 | n/a | n/a | n/a | tumour 3678; micro-environment (unsplit) 1287 |
| ER- (TNBC only) | Pal | 8 | 8 | n/a | n/a | n/a | tumour 1956; micro-environment 1138 |
| ER unknown (HER2+) | Pal | 6 | 6 | n/a | n/a | n/a | tumour 3620; micro-environment 1590 |
Pal low-count cases: TN-0106 has 54 tumour cells (825 micro-environment), ER-0032 417 tumour cells; all others > 175 tumour. All 27 Pal samples have >= 50 micro-environment cells, but CAF/endothelial subsets are unknown. Pooled potential: ER+ 27, ER- 20 (Wu 12 + Pal TNBC 8), +6 Pal HER2 with unknown ER, with Pal stromal classes pending a marker step. Several Wu patients have 0 author "Cancer Epithelial" cells (CID3586, 3838, 3946, 4040, 4398, 44041), so effective Wu n for malignant >= 50 is 11 / 9.

## 2.4 Same patients across Wu and Pal?
Different cohorts: Wu from Sydney (Garvan, St Vincent's, Chris O'Brien Lifehouse; affiliations on the paper page), Pal from Melbourne (Royal Melbourne Hospital tissue bank, Victorian Cancer Biobank, kConFab; Chen 2022 Methods). Different id systems (CID#### vs 4-digit), different institutions and ethics. No overlap evidence; no direct patient-level cross-check is possible from public metadata, so independence is assumed, not proven. Within Pal, patient 0114 appears twice (ER+ and TNBC tumours); 0064 and 0023 also have normal or BRCA1 samples but those are not tumour samples used here.

## 2.5 Downloads (this part) and blockers
Downloaded: Wu ESM (pdf x3 + xlsx, 18.6 MB) and article html (0.6 MB); Pal figshare: 4 .rds + zip = 6.64 GB; Table EV4 xlsx. Total about 6.7 GB (limit 12). GSE161529_RAW.tar not downloaded (zip contains the same matrices). Provenance files updated in the data folders.
Blockers: (1) Pal CAF/endothelial classes need a marker-based labelling of microenvironment cells, which reads expression; it requires registering the labelling rule first. (2) ER status of the 6 HER2+ Pal tumours unreported. (3) Seven Pal ER+ samples (incl. two males) lack a malignant definition. (4) Wu ER-low case CID3963. (5) inferCNV block file does not state tumour vs normal; the `Tum` clusters come from the R code.
