# A126 data scout (2026-10-29)

Scope: data and metadata only. No clock statistic, coherence, phase or clock-gene summary was computed. Facts marked [V] were checked in a file on disk or a downloaded metadata file; [U] were not (web page or memory only).
Downloaded (metadata only, parsed with `python -I`, never executed): `$SCCC_DATA/scout_a126/{gse198732,gse111976,gse161529}` (13 MB; Murrow metadata rds parsed with the repo's `stream_rds.py` reader, copied to the scratchpad; no R).
Cell-count and donor-count tables below are metadata tabulations only.

## Headline

1. **Murrow GSE198732 carries no per-donor hormone-phase annotation.** [V] The metadata data frame has `Age, BMI, Obesity, Parity, Parity.Y.N, HC_use, Type, Sort, Batch` and nothing on cycle day, progesterone or luteal state. The authors' "hormone state" is a high/low ER/PR-signalling donor clustering from expression (paper text, [V] via PMC; not stored in GEO). So the luteal score has to be computed by us, as planned, and the "authors' annotation" arm of the pre-registration cannot be used for Murrow.
2. **Murrow alone gives 9 donors per tertile, below the >=10 rule.** [V] 28 donors; 27 have >=50 HR+ luminal cells. Adding the Komen Tissue Bank (KTB) donors in the same GEO series gives 34 (7 more with >=50 HR+ cells), i.e. 11/11/12. Pooling Murrow + Kumar + Reed (study as stratum) is comfortably feasible (below).
3. **Only the endometrium has an authors' phase label per donor** (GSE266265: all 60 donors; GSE111976: 10 donors on 10x; E-MTAB-10287: 5 donors). No breast dataset records cycle phase.
4. **On disk now:** Kumar HBCA, Reed HBCA (both as Census donor x cell-type pseudobulks with the full 61,497-gene panel; PGR, TNFSF11, WNT4, AREG, MKI67 present as var names [V]), GSE266265 (full per-cell counts export). **Not on disk:** Murrow, Pal (GSE161529), Wang (GSE111976), Garcia-Alonso (E-MTAB-10287), Ulrich 2024.
5. **Caveat on the luteal score for the on-disk pseudobulks:** Kumar's pseudobulk is grouped by Census `cell_type`, where HR+ luminal and secretory/progenitor luminal share the term "luminal epithelial cell of mammary gland" [V: 114 donors with n>=50]. The authors' labels (`author_cell_type` e.g. LummHR-SCGB, `broad_cell_type` lumhr) are per-cell columns that the pseudobulk dropped [V `obs_columns.csv`], so a Kumar luteal score needs a rebuild grouped by `broad_cell_type`. Reed's Census terms already separate "luminal hormone-sensing cell" (54 donors n>=50), "luminal adaptive secretory precursor" (51) and "luminal epithelial cell" (38) [V], so Reed works as built.

## Per-dataset tables

### 1. Murrow et al. 2022 (Cell Systems), GSE198732: primary

| item | finding |
|---|---|
| (a) accession, files | GEO GSE198732, two platform series matrices (GPL20301 4.1 KB, GPL24676 3.3 KB) [V]. Supplementary, all `.rds.gz` (R objects, no h5ad) [V listing]: `breast.data.rds.gz` 1.7 GB; `RM.raw.counts.matrix.rds.gz` 441 MB; `KTB.raw.counts.matrix.rds.gz` 56 MB; `RM.metadata.rds.gz` 1.7 MB and `KTB.metadata.rds.gz` 0.17 MB (downloaded, parsed); MULTI-seq barcode count files (small). Not on CELLxGENE Census [V: absent from `census_datasets.csv`]. Raw counts matrix needs an R-serialization reader (precedent: `stream_rds.py`, GSE266265); `breast.data.rds` is probably the Seurat object [U]. FASTQ also on GEO/SRA (not needed). |
| (b) per-donor metadata | RM: 28 donors (RM142 ... RM307), 86,136 cells, ages 19 to 39 (median 23) [V]; all premenopausal by design. BMI present for 16/28 donors only [V]; `Obesity` (obese 10, non-obese 6, NA 12) [V]; Parity (0 / 1 / 2 / 3) for 22/28, `Parity.Y.N` no 11, yes 11, unknown 6 [V]. `HC_use`: none 16, combined 9, progestin 3 [V] (progestin or combined contraception will perturb PR targets: a confounder for the luteal score). No cycle day, no hormone assay, no BRCA (healthy reduction mammoplasty). KTB: 9 donors, 8,485 cells, no age/BMI/parity columns; older (37 to 47), parous, BMI 20.7 to 28.3 per paper text [U]. |
| (c) cell annotation | `Type` [V]: HRpos_Luminal 22,274; Secretory_Luminal 16,788 (this is the luminal progenitor/secretory lineage); Basal 19,021; Fibroblast 14,863; Vascular_Endothelial 7,661; Vascular_Accessory (pericyte/SMC) 2,476; Lymphocyte; Macrophage; Lymphatic_Endothelial. Donors with >=50 cells: HRpos_Luminal 27, Secretory_Luminal 26, Basal 21, Fibroblast 27, Vascular_Endothelial 26, Lymphocyte 13. Sort gates: Live 44,546; Epithelial 19,037; Basal 11,758; Luminal 10,795 (sorted gates distort compartment proportions, not within-type expression, but gate must be a stratum). |
| (d) premenopausal donors; tertiles | 28 premenopausal (RM) + 9 KTB (older, status not stated in the metadata [U]). Luteal tertiles from our score: 27 donors usable (HR+ >=50 cells and a stromal host with >=50 cells: 26) gives 9/9/9 (8/9/9 with host). **Infeasible alone against the >=10 rule; feasible only pooled** (Murrow+KTB 34 gives 11/11/12). The authors' own high/low ER/PR split is a binary clustering (per the text: 2 of 7 nulliparous and 0 of 7 parous non-HC-users were "high" ) [U], i.e. few donors in the high arm. |
| (e) on disk | No. Builder would be new (`breast_murrow`). |
| suitability | **Medium.** Best-designed (premenopausal, parity, obesity, HC use, sorted epithelium) and the only source for H126.3 with real obesity contrast (10 vs 6 with known status; BMI for 16). Too small for a standalone tertile test, and BMI is missing in 12/28, so H126.3 would be n = 16 by BMI or 16 by `Obesity`. |

### 2. Kumar et al. 2023 (Nature), HBCA via CELLxGENE: on disk

| item | finding |
|---|---|
| (a) accession, files | CELLxGENE Census 2025-11-08 dataset `842c6f5d-4a94-4eef-8510-8c792d1124bc` ("scRNA-seq data - all cells", 714,331 cells) [V]. On disk: `atlas_wg/breast_kumar/pseudobulk.h5ad` 99 MB (4,178 groups x 61,497 genes), `atlas/breast_kumar/pseudobulk.h5ad` 0.8 MB (small gene panel), `donors.csv`, `obs_columns.csv`; `cells/a152_breast_kumar/cells.parquet` (T cells only). |
| (b) per-donor metadata | 126 donors [V]; `donor_menopausal_status`: premenopausal 91, postmenopausal 25, unknown 10 [V]; BMI for 117 donors (90 of the premenopausal) [V]; `bmi_group` (5 levels); `age_group` young 96 / old 30; `development_stage` ; `procedure_group` reduction mammoplasty 76, cancer mastectomy 35, prophylactic mastectomy 15 [V]; `reported_diseases` none 86, breast cancer 40 [V] (the normal tissue of cancer patients: a covariate); `breast_density`; sites Baylor 73, MD Anderson 31, UC Irvine 22 [V] (site = stratum). No parity, no BRCA, no cycle phase, exact age not given (decade stage) [V]. |
| (c) cell annotation | Census `cell_type` (39 types), `author_cell_type` (58, e.g. LummHR-SCGB) and `broad_cell_type` (lumhr, lumsec, basal, ...) per cell [V `obs_columns.csv`]. In the pseudobulk only Census `cell_type` survives: luminal epithelial (114 donors), basal (105), mammary gland epithelial (86), fibroblast (126), vein endothelial (112), capillary endothelial (103) at n>=50 [V]. HR+ vs progenitor luminal not separable in the existing pseudobulk. |
| (d) premenopausal; tertiles | 91 premenopausal (90 with BMI). Tertile split of luteal score: ~30 per tertile if HR+ cells are rebuilt (luteal score donor-level). Comfortably >=10. |
| (e) on disk | Yes (pseudobulk); needs a rebuild by `broad_cell_type` (via Census, ~25 min, 4 threads as in the A152 extraction) or a per-cell extraction of the 4 score genes plus clock panel for lumhr cells only. |
| suitability | **High** for H126.1 (n), H126.2 and H126.3 (BMI 117). No cycle label, so the luteal score is the sole stratifier; phase is a latent variable, and with menopausal status the post vs premenopausal contrast is a natural negative control (no luteal signal). |

### 3. Reed et al. 2024 (Nature Genetics), HBCA: on disk

| item | finding |
|---|---|
| (a) accession, files | Census dataset `0ba636a1-4754-4786-a8be-7ab3cf760fd6` ("HBCA - global", 803,283 cells; also HBCA-epithelial 352,496, HBCA-stroma 416,825, HBCA-immune 25,382, integrated 2,122,065) [V]. On disk: `atlas_wg/breast_hbca/pseudobulk.h5ad` 36 MB (599 groups x 61,497 genes), `donors.csv`. Authors' primary archive (ArrayExpress/EGA) not verified [U]. |
| (b) per-donor metadata | 55 donors [V]; exact `donor_age` 19 to 65 (median 38); BMI for 23 donors [V]; `parity` 0:17, 1:7, 2:17, 3:7, 4:3, unknown 4 [V]; `brca_status` BRCA1 17, BRCA2 11, WT 2, assume_WT 22, unknown 3 [V]; condition: prophylactic mastectomy 27, reduction mammoplasty 22, contralateral 6 [V]. **No menopausal-status column**; age is the proxy: age <=45: 38, <=50: 43, >=55: 9 [V]. No cycle phase. |
| (c) cell annotation | Census terms [V, donors with n>=50]: luminal hormone-sensing 54, luminal adaptive secretory precursor 51, luminal epithelial 38, basal-myoepithelial 53, fibroblast of mammary gland 54, blood vessel endothelial 52, perivascular 52, lymphatic endothelial 42. Authors' `level1/level2` labels are in the obs (LummHR etc., [V `obs_columns.csv`]) but only Census terms are in the pseudobulk. |
| (d) premenopausal; tertiles | ~38 donors aged <=45 (menopause inferred, [U]). 12 to 13 per tertile if all 38 are used; fewer with BMI. BRCA1/2 carriers 28 vs non-carriers 24 to 27 with HR+ luminal groups: descriptive contrast feasible. |
| (e) on disk | Yes. Can score luteal genes in the "luminal hormone-sensing" pseudobulk without rebuild. |
| suitability | **High** (BRCA contrast, parity, HR+ separation already built). Age-only menopausal proxy and 23/55 BMI are the weaknesses. BRCA mutation carriers' tissue is confounded with procedure type (prophylactic mastectomy), and age. |

### 4. Pal et al. 2021 (EMBO J), GSE161529

| item | finding |
|---|---|
| (a) accession, files | `GSE161529_RAW.tar` 2.2 GB (10x mtx per sample), `features.tsv.gz`, filelist [V]. Counts only; no raw reads in SRA (privacy) [V]. Not downloaded. The authors' re-aligned matrices exist in GSE276609 (VCU) [V via an existing note]; Chen et al. 2025 bioRxiv "Human Breast Cancer Single Cell Atlas" is on Census (global 621,200 cells) [V]. |
| (b) per-donor metadata | Series matrix [V]: normal tissue = 24 samples from 18 patients (Total or Epithelial fractions, some both); menopause Pre 13, Post 5 (incl. 1 oophorectomy, 1 hysterectomy); parity nulliparous 11 / parous 7; no age, no BMI, no cycle phase in GEO characteristics; 4 BRCA1 pre-neoplastic samples (2 pre, 2 post-oophorectomy) [V]. |
| (c) cell annotation | Not in GEO; Chen et al. 2022 (Sci Data) provides Seurat objects with cell labels [U]. |
| (d) premenopausal; tertiles | 13 premenopausal normals. **Too few for a tertile test**; at most added as a stratum. |
| (e) on disk | No (a note in the scout of earlier sessions lists it as "normals, not paired"; not in `geo/`, `scbc_full`). |
| suitability | **Low** for A126 (13 premenopausal, no BMI/age/phase). Useful for BRCA1 pre-neoplastic descriptives only (n=4). Still needed for A127 H127.3. |

### 5. Endometrium across the cycle

| item | GSE266265 (endometriosis atlas) | Wang et al. 2020, GSE111976 | Garcia-Alonso et al. 2021, E-MTAB-10287 |
|---|---|---|---|
| (a) files | `raw/` Seurat rds 41.7 GB (+ epithelial monocle3 rds 0.2 GB); **already exported**: `export/counts_{i,x}` 4.75 GB each, `meta.csv` 115 MB, 34,727 genes x 466,371 cells [V]. Pseudobulk 148 KB (68 clock-panel genes only) [V]. | `GSE111976_ct_endo_10x.rds.gz` 408 MB (10x, 10 women), `ct.csv.gz` 12 MB (C1 Smart-seq-like, 2,148 cells, 19 women), cell-annotation csvs (0.6 MB), donor-phase csvs [V listing; small csvs downloaded]. | Per-sample Matrix Market (11 samples, `.mtx` 18 to 232 MB each, ~1.6 GB total) + sdrf/idf [V via BioStudies]; no h5ad. Census mirror: dataset `9ddea8d9` "Non-pregnant Uterus (Endometrium - All)", 100,307 cells, listed as 15 "donors" in Census (sample-level), ENA ERP128125 [V]. |
| (b) donor metadata | 60 participants, one sample each; cycle phase per donor: proliferative 28, periovulatory 18, secretory 14 (early 4, mid 7, late 3); `CycleDay` and `ProgesteroneSerum` for all 60 [V]; endometriosis 35 (mild 15, severe 20) vs control 25 [V]; phase x status: controls 13 prolif / 5 periov / 7 secr [V]. | 10x: 10 donors with `phase_canonical` (proliferative 2, secretory early 3, early-mid 3, mid 1, late 1) and `day` (16 to 26) [V]; C1: 19 donors incl. menstrual 4, proliferative early 3, late 3, secretory early 3, mid 3, late 2 [V]. 2 women in both. | 5 donors (3 biopsies + 2 organ donors), phase: proliferative 2 donors (A13, E3), early-secretory 1 (A30), mid-secretory 1 (E2), late-secretory 1 (E1) [V via BioStudies]. No cycle day. |
| (c) cell annotation | `AnnotationMain` [V]: stromal 411,534; epithelial 26,270; endothelial 15,855; myeloid; lymphocyte. `AnnotationRefined` [V]: glandular 11,390; mid-secretory 6,027; early-secretory 1,570; ciliated 4,338; lumenal ~630; TP63+/KRT5+ 50; fibroblast subtypes, EC, pericyte. Epithelial cells per donor: median 310, min 85; 58/60 donors >=100 [V]. A `cycling.fibroblast` flag exists (25,309 TRUE) [V]; no cycling-epithelial flag (G2M score needed [U]). | 10x [V]: stromal fibroblasts 22,102; unciliated epithelia 1 and 2 (18,826; 16,399); ciliated 4,393; lymphocytes; endothelia; SMC; macrophages. | Census terms [V]: stromal 37,111; epithelial 17,898; glandular secretory epithelial 13,100; ciliated 3,189; endothelial; smooth muscle. |
| (d) premenopausal; tertiles | All premenopausal by design (endometrial biopsies; ages not in metadata checked). Phase tertiles feasible by progesterone serum or cycle day (20/20/20). Phase-labelled: 28 prolif vs 14 secretory. | 10 (10x) or 19 (C1) donors. C1 cells per donor ~110: too sparse for pseudobulk. | 5 donors. |
| (e) on disk | **Yes** (`geo/GSE266265`; builders `paper/common/datasets/endometrium` and `endometrium_scored`). Phase is already the `phase` column (proliferative/other/secretory) of `pseudobulk/donors.csv`, where `other` = periovulatory [V]. | No. | No (Census route available via `census_env`). |
| suitability | **High** for H126.0 (largest, full per-cell counts and PGR, TNFSF11, WNT4, AREG, MKI67, TOP2A all in the 34,727-gene list [V]), but a descriptive caveat: epithelial cells per donor are few (median 310) so cycling cells per donor may be a few tens: **donor-level power to be checked before registration**. Endometrial epithelium is a hormone-driven proliferative tissue: a good test of proliferation vs clock independent of breast. | **Low-medium**: 10 donors, secondary replication only. | **Low** (5 donors). Ulrich et al. 2024 PNAS (Census `8a961a54`, 20 donors, 167,910 cells, endometrium and myometrium) is a better external endometrium set than Wang or Garcia-Alonso, with phase labels in the paper [U]. |

### 6. GTEx breast (counts only; controlled access, aggregates)

| item | finding |
|---|---|
| Source | `paper/data/Gtex_phenotypes` (dbGaP phs000424.v8; subject and sample attributes) [V]. |
| Counts | Women with a Breast - Mammary Tissue RNA-seq sample (freeze RNASEQ): 168. Aged <45: 46 (16 with BMI >=30, BMI known for all 46); aged <50: 71. All 46 have a BMI value. Including non-RNA-seq analytes: 174 women, 48 <45. |
| Cycle information | No menopause, pregnancy or contraceptive variable is among the 190 subject variables (only SEX, AGE, BMI plus medical-history items) [V]. A luteal score can only be inferred from expression (breast is whole tissue, a mixture of adipose and stroma with few HR+ luminal cells). |
| Suitability | **Low-medium, descriptive.** 46 women under 45: a tertile split of an expression-derived score gives ~15 per tertile, but the HR+ luminal fraction in whole tissue is low and varies with adiposity. Time of death comes from the controlled file; report aggregates only. |

## Recommended primary design

**Cohort for H126.1/H126.2/H126.3:** pool premenopausal donors of Kumar (91), Reed (about 38 by age <=45; add rebuilt Murrow RM, 27 usable). Study as stratum (Kumar site as sub-stratum; Murrow sort gate and HC use as covariates). Murrow KTB (7 usable donors) is older and parous: include only as a descriptive stratum, since menopausal status is not recorded.

**Luteal score and tertiles.** Per donor, the mean z of PGR, TNFSF11, WNT4, AREG across the HR+ luminal cells (Murrow `HRpos_Luminal`, Kumar `broad_cell_type == lumhr`, Reed "luminal hormone-sensing"). Z-scoring and tertile cutoffs **within study** (library chemistry differs: v2 vs v3), then top vs bottom tertile pooled. Require >=50 HR+ cells and >=50 host (fibroblast/endothelial) cells per donor. Expected tertile sizes: Kumar ~30, Reed ~12, Murrow ~9 (pooled ~50 per top/bottom tertile). Expectation to fix in the registration: progestin or combined contraception (Murrow 12/28 donors) shifts the score; run with and without HC users (HC status exists for Murrow only, so the cross-study analysis should not condition on it; treat as a Murrow sensitivity). Postmenopausal Kumar donors (25 + 10 unknown) can serve as the within-pipeline negative control (score ~ flat).
- Score validity gate before any clock reading: that the score is bimodal or has a spread separating tertiles, that it is not a depth or cell-number artefact (HR+ cell count per donor), and that it correlates with age <45 vs older in the expected direction. This can be checked without any clock gene.
- No dataset offers an authors' phase label in breast. Murrow's authors' high/low ER/PR state (binary) could be recomputed only by their clustering; do not rely on it.

**H126.0 (cycling vs non-cycling epithelium within donor).**
- *Endometrium (GSE266265): feasible.* 60 donors, epithelial cells per donor median 310 (58/60 >=100). Cycling cells flagged by the A97 G2M set with the A24 top-quartile rule leaves maybe ~75 cycling and ~75 non-cycling cells for a median donor and only ~25 for the smallest donors; the number of donors with >=30 cells in each group must be counted before registration (G2M scoring is expression-based and was deliberately not done here). Phase-confounded: report within phase strata (proliferative has the cycling cells), so the paired comparison (same donor, cycling minus non-cycling) is the right design and proliferative-minus-secretory remains descriptive.
- *Breast: feasible with caveats.* Murrow `HRpos_Luminal` + `Secretory_Luminal` + `Basal` per donor yield 17 to 25 donors with >=100 cells of a lineage; cycling cells are rare in normal breast epithelium (<1 to 2% typical [U]), so with ~1,000 epithelial cells per donor ~10 cycling cells per donor is expected: probably underpowered for within-donor cycling vs non-cycling. Kumar/Reed per-cell data have more cells (Kumar 714k, Reed 803k); the epithelial subset with cycling flag needs the counts to be extracted (not in pseudobulk). The registration should state the minimum cycling cells per donor (suggest >=30) and report the number of eligible donors per cohort, with H126.0 for breast declared infeasible if fewer than ~15 donors qualify.

**Secondary reads.** Epithelial phase relative to stromal host by luteal tertile (Kumar, Reed, Murrow have host fibroblast/endothelium for >=50 cells in 27/28, 52 to 54, 112 to 126 donors). BRCA carriers vs non-carriers in Reed (28 vs 24 to 27; age and procedure confounded). Parity: Murrow 11 vs 11 plus Reed (17 nulliparous vs 34 parous with known parity). BMI: Kumar 117 (premenopausal 90) is the main H126.3 cohort; Murrow obese 10 vs non-obese 6 is a replication; Reed 23.

**Not recommended for the pre-registered primary:** Pal, Wang, Garcia-Alonso (too few premenopausal donors or too few cells), and GTEx breast (descriptive only).

## Things to check or decide before registering

1. Kumar rebuild by `broad_cell_type` (HR+ luminal) and per-cell extraction of the G2M set, via `census_env`, 4 threads (the existing `cells/a152_breast_kumar` extractor is the template).
2. Murrow counts need an R-serialization reader (precedent `stream_rds.py`); `RM.raw.counts.matrix.rds.gz` is 441 MB compressed (RAM at parse is bounded by the streaming reader, not by the file). Confirm gene ID type (symbol vs Ensembl) and integer counts.
3. Overlap of donors across Kumar, Reed and Murrow was not checked. Murrow (UCSF, Komen Tissue Bank) and Kumar (Baylor/MD Anderson/UC Irvine) are different cohorts per the metadata [V for Kumar sites]; Reed's sites are not in the metadata [U].
4. Reed menopausal status is missing: choose an age rule (<=45) before scoring; label any other as post hoc.
5. Endometrium per-donor cycling-cell counts (power) and whether `phase` "other" (periovulatory) is merged.
6. Controlled data rule: GTEx numbers above are aggregates; do not add per-donor values.
