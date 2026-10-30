# METABRIC and TCGA-BRCA feasibility scout (2026-10-30)

Metadata only; no clock statistic. Scripts: `scout_metabric.py` (-> `scout_metabric.out`), `scout_tcga_brca.py` (-> `scout_tcga_brca.out`), both `python -I`.

## METABRIC (cBioPortal `brca_metabric`)
- Data: `$SCCC_DATA/metabric/brca_metabric/` (PROVENANCE.md there: URL, 2026-10-30, sha256). Total ~1.2 GB.
- Licence: ODbL 1.0 (study LICENSE file): share and adapt with attribution, derived datasets stay open, share-alike. Cite Curtis 2012, Pereira 2016, Rueda 2019.
- Expression file: `data_mrna_illumina_microarray.txt` = Illumina HT-12 v3 log2 intensities (values 4.56-14.88; gene-mean SD 1.53, per-sample means nearly constant: mean 6.71, SD 0.011, so already quantile-normalised). A z-scored companion (`..._zscores_ref_diploid_samples.txt`) exists in the tarball; do NOT use it. 20,603 gene rows x 1,980 tumours (Entrez id column dropped); 218 symbols have >1 probe row (none among clock or purity genes; G2M: BIRC5, CDKN2C, KPNA2, average or take max). 16 NaN cells in 16 samples.
- Samples: 2,509 patients/samples in the clinical files, all Primary, one sample per patient, IDs MB-xxxx; only 1,980 have expression (the Curtis 2012 discovery + validation microarray set; `case_lists/cases_mRNA.txt` = 1,980). Female only.
- Gene coverage by symbol (no alias needed): clock 14/14 (single probe each), Hallmark G2M (minus clock) 200/200, purity 11/11.
- ER: `ER_IHC` (patient file; note the typo level "Positve"): Positve 1,498, Negative 439, NA 43. `ER_STATUS` (sample file, integrates expression/IHC per cBioPortal): Positive 1,506, Negative 474, no NA. They disagree on 128 tumours (IHC- / status+ 55; IHC+ / status- 73). Recommend `ER_IHC` primary, `ER_STATUS` as sensitivity.
- HER2: `HER2_SNP6` (NEUTRAL 1,436, GAIN 438, LOSS 101, UNDEF 5); `HER2_STATUS` Positive 247 / Negative 1,733.
- PAM50 + Claudin-low (`CLAUDIN_SUBTYPE`): LumA 700, LumB 475, Her2 224, claudin-low 218, Basal 209, Normal 148, NC 6.
- Cohort/batch: `COHORT` 1-5 (n 521, 288, 763, 238, 170); also `INTCLUST` (IntClust 1-10, 4ER+/4ER-), `THREEGENE`. No explicit array batch/centre variable beyond COHORT (cohorts = recruitment centres/series, array processing in different sets; treat COHORT as batch covariate).
- Outcomes (1,980 with expression): OS_MONTHS/OS_STATUS 1,980 (1,143 deaths; VITAL_STATUS: Died of Disease 646, Other causes 496, Living 837, 1 NA => disease-specific survival derivable); RFS_MONTHS/RFS_STATUS 1,979 evaluable (803 recurred); median follow-up OS 116 months. Treatment: HORMONE_THERAPY yes 1,216 / no 764; CHEMOTHERAPY yes 412 / no 1,568; RADIO_THERAPY yes 1,173. ER+ (IHC) tumours: hormone therapy yes 1,104, of which chemo-treated 133; no hormone therapy 394 (chemo 7). Also grade (1: 169, 2: 771, 3: 952, 88 NA), cellularity, NPI, age, menopausal state, stage, size.
- Caveat: treatment is non-randomised, hormone therapy is almost exclusively ER+, so endocrine-response analyses are confounded by indication.

PAM50 x ER (IHC), tumours with expression:

| PAM50 | ER- | ER+ | NA | all |
|---|---|---|---|---|
| Basal | 172 | 33 | 4 | 209 |
| Her2 | 115 | 103 | 6 | 224 |
| LumA | 13 | 677 | 10 | 700 |
| LumB | 11 | 456 | 8 | 475 |
| Normal | 28 | 117 | 3 | 148 |
| claudin-low | 100 | 106 | 12 | 218 |
| NC | 0 | 6 | 0 | 6 |
| all | 439 | 1,498 | 43 | 1,980 |

The ER- cells in LumA/LumB (13, 11) are too small for subtype-stratified contrasts; the stratified test has power only in Basal, Her2, claudin-low, Normal. PAM50 x COHORT and ER x COHORT tables are in `scout_metabric.out` (every cohort has >= 44 ER- tumours).
- Proliferation tertiles (G2M) not computed here.
- Overlap with TCGA / SCAN-B: none expected. Basis: independent series (from the study descriptions, not computed from files: METABRIC = UK/Canada series, microarray; SCAN-B = Swedish, RNA-seq; TCGA = US-led consortium; IDs MB-xxxx vs TCGA-xx vs SCAN-B GSM), no shared identifiers. Not verified by genotype/expression matching.

## TCGA-BRCA on disk
- Expression: `$SCCC_DATA/tcga/recount3/BRCA/tcga.gene_sums.BRCA.G029.gz` (recount3, GENCODE v29): gene-level base-pair coverage sums, NOT read counts. Per `tcga/recount3/PROVENANCE.md`: AUC-scaled counts = round(raw x 4e7 / `recount_qc.bc_auc.all_reads_all_bases`), or read counts = round(raw / `star.average_mapped_length`). Recommended: AUC-scaled counts then log CPM, as in the existing pipeline.
- Samples: 1,256 (Primary Tumor 1,127; Solid Tissue Normal 112; Metastatic 7); primaries cover 1,093 patients (34 patients with a second primary aliquot; need one-per-patient rule). Histology of primaries: ductal 799, lobular 212, other types 116 (mucinous 17, metaplastic 13, medullary 6, mixed 29, other 49, NOS 1, NA 1).
- ER (IHC) ON DISK: `tcga/recount3/BRCA/tcga.tcga.BRCA.MD.gz`, column `xml_breast_carcinoma_estrogen_receptor_status` (one primary per patient, n=1,093): Positive 806, Negative 237, Indeterminate 2, NA 48. Detection method is heterogeneous (`xml_er_detection_method_text`, mostly missing). PR and HER2 IHC columns also present.
- PAM50: NOT on disk; fetched `$SCCC_DATA/tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt` (column `SUBTYPE`, PanCanAtlas, via cBioPortal; the GDC does not host the TCGASubtype file). Coverage on the 1,093 patients: LumA 499, LumB 197, Basal 171, Her2 78, Normal 36, missing 112. Licence: Broad GDAC TCGA pipeline licence (open, provisional data). Sources recorded in `tcga/pancan_subtypes/PROVENANCE.md`.

| PAM50 | ER- | ER+ | Indet. | NA | all |
|---|---|---|---|---|---|
| Basal | 144 | 20 | 0 | 7 | 171 |
| Her2 | 46 | 26 | 2 | 4 | 78 |
| LumA | 11 | 469 | 0 | 19 | 499 |
| LumB | 3 | 182 | 0 | 12 | 197 |
| Normal | 13 | 20 | 0 | 3 | 36 |
| no label | 20 | 89 | 0 | 3 | 112 |
| all | 237 | 806 | 2 | 48 | 1,093 |

## Recommendations and caveats
- METABRIC: use `data_mrna_illumina_microarray.txt` (log2 intensity). TCGA: recount3 AUC-scaled counts -> log CPM. SCAN-B: as in `a127_lib.py`. Microarray intensities have a minimum of 4.56 (background floor) and probe-based; coherence (correlation-type) measures are scale-free but not floor-free, so check low-expressed clock genes (NPAS2, PER-family, CIART) against the floor before comparing across platforms. Compare direction/effect on rank or correlation scale only, not absolute values.
- METABRIC ER- tumours are mostly Basal/Her2/claudin-low; a within-LumA/LumB ER contrast is not feasible (13 and 11 ER-); the same holds in TCGA (11 and 3).
- Use COHORT as a covariate/stratifier in METABRIC.
- No blocker. Survival and endocrine analyses need PI approval before any scoring.
