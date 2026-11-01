# Scout C4 single-cell: clock locked to the cell cycle in malignant cells? (feasibility, 2026-11-01)

Read-only scout. Only cell counts, UMI depth and detection of the 14 clock genes were computed; no clock statistic across
cell-cycle bins. Scripts and CSVs are in this folder (`scout_*.py`, run with `python -I`, sccc env, 4 threads).

## 1. Data behind the A87 single-cell set
`sh.load_sc` -> `data.load("cancer_sc")` (currently the legacy conversion, `paper/common/datasets/cancer_sc`) = 29 cohorts with
copy-number (CNV) calls, 418 patients; **251 patients in 29 cohorts have `tumour:malignant`** (191,188 cells); the A87 set (79
patients, `paired_epithelium`) has 57,855 malignant cells. Malignant label = `cnv_call == "malignant"` of
`$SCCC_DATA/cnv/<study>.csv.gz` (sample, barcode, patient, cnv_score, clonality, cnv_call). Raw UMI counts: GEO h5ad
`$SCCC_DATA/geo/<GSE>/<GSE>.h5ad` or scBaseCount per-sample h5ad `scbc_full/<SRP>/*.h5ad` (integer counts, per the builder docs).
Per-cell tables already exist: `$SCCC_DATA/cells/cancer/<study>_epithelial.csv.gz` (14 clock counts, S, G2M, no barcodes) and
**CoPhaser per-cell cell-cycle tables of 17 studies** at `/home/maxine/Documents/paychere/scCircadianCancer/outputs/states/cc/<study>.csv.gz`
(theta, cycling_logit, S, G2M, clock counts, patient, kind, compartment, total_umi; `$SCCC_DATA/cells/cancer_cycle/` is NOT populated, the
builder `paper/common/datasets/cells_cancer_cycle/build.py` has `from_legacy`). Nine of the 17 have CNV calls
(ERP149257, ERP153272, GSE131907, GSE132465, GSE178341, GSE181919, GSE183904, SRP273482, SRP327758): 198 patients, **110 with >= 300 malignant cells**.

Per-study counts (`per_study_malignant_counts.csv`, all patients with malignant cells; `a87` rows in the CSV):

| study | cancer | pts | malignant cells median (min) | pts >= 300 | pts >= 500 | median UMI/cell* | mean detection, 14 genes (cells) |
|---|---|---|---|---|---|---|---|
| GSE178341 | colorectal | 60 | 848 (96) | 56 | 46 | 19,326 | 0.148 |
| SRP327758 | oesophageal | 31 | 258 (57) | 12 | 7 | 16,244 | 0.301 |
| GSE132465 | colorectal | 22 | 459 (157) | 15 | 9 | 19,992 | 0.128 |
| GSE183904 | gastric | 22 | 258 (81) | 9 | 4 | 8,881 | 0.118 |
| GSE131907 | **LUAD** (Kim) | 14 | 332 (54) | 8 | 5 | 11,214 | 0.136 |
| GSE166555 | colorectal | 12 | 450 (85) | 8 | 4 | n/a (no CoPhaser table) | |
| GSE181919 | head and neck | 12 | 267 (68) | 4 | 3 | 23,222 | 0.196 |
| ERP149257 | colorectal | 8 | 442 (60) | 6 | 2 | 11,993 | 0.218 |
| GSE294300 / GSE315534 / GSE200997 | colorectal | 10 / 5 / 8 | 183 / 2,385 / 261 | 4 / 5 / 3 | 3 / 4 / 2 | n/a | |
| SRP416180 | bladder | 6 | 1,214 (133) | 4 | 4 | n/a | |
| GSE212966 | pancreatic | 6 | 680 (497) | 6 | 5 | n/a | |
| 11 other cohorts | kidney, prostate, liver, thyroid, ampullary, gastric, breast (SRP572743, 4 pts) | 4-7 each | 142-705 | 0-3 each | | n/a | |
| **Total (29 cohorts)** | | **251** | | **153** | **106** | | |
*median `total_umi` of malignant cells (CoPhaser tables). All cohorts were chosen with sample median UMI >= 1,500; the actual medians are 8-23k.
Detection per cell is low (mean over the 14 genes 0.12-0.35; CLOCK, CRY1, PER2, DBP highest, TEF, CIART, ARNTL lowest; 68-97 % of cells have >= 1 clock UMI, median
1-4 clock UMI per cell): `per_study_cell_level_detection.csv`. Per-cell clock phase is hopeless, pseudobulk bins only.
No per-patient clinical label is stored with the 29-cohort set beyond cancer type (`donor` table: study, patient, cancer, sample types). Proliferation must be
derived from the cells (G2M/S score or CoPhaser cycling fraction), which makes it non-independent of the cell-cycle binning.

**Wu 2021 breast** (GSE176078, `$SCCC_DATA/sc_breast_er/GSE176078/extract/`, 2.4 GB mtx + metadata; integer counts, 100,064 cells, 26 patients): malignant = author
`Cancer Epithelial` (24,489 cells, 20 patients, 6 patients with 0). Median UMI 8,573; mean detection 0.114 (ARNTL 0.069, TEF 0.022, CRY1 0.220); 67 % of cells >= 1 clock UMI.
Per patient (`wu_per_patient.csv`): n >= 50 = 20, **>= 300 = 14 (ER+ 6, ER- 8)**, >= 500 = 13. ER: `scouts_2026-10-31/sc_breast_er/wu_reconciled.csv`
(final_ER from Suppl. Table 1; IHC subtype, treatment). Read in 4 threads, 1 pass of the mtx (a few minutes).
**Pal 2021 breast** (GSE161529, figshare 17058077, `.../sc_breast_er/GSE161529/`; raw 10x only inside `HumanBreast10X.zip`, counts NOT extracted; Seurat rds with meta only): malignant =
authors' `Tum` clusters (A170). A170 units (`pal_per_patient_malignant.csv`): 21 tumours (ER+ 13, TNBC 8; the 6 HER2+ have no ER label and were excluded), median 2,899
malignant cells, min 54, **19 with >= 300**, 18 >= 500; mean UMI per cell (depth / n) median 8,679, range 2.8-16k. Cell-level clock detection not computed (needs the zip).
Neither Wu nor Pal has a CoPhaser table: CoPhaser would have to be run (GPU; archive fits took 1.6-2.6 min per cohort).

## 2. LUAD
On disk: **GSE131907 (Kim 2020 LUAD)**, 14 patients with malignant cells, 8 with >= 300 (5 >= 500): too few alone. Other lung single-cell on disk: ERP153272 (NSCLC atlas, 14 tumour
samples but only 1 patient with CNV malignant calls; the inventory lists 7 patients), ERP156718 (NSCLC immune/stromal, adenocarcinoma + squamous, 7 pts, median UMI 9.8k), GSE117570 (early NSCLC, tiny),
GSE153935 (Drop-seq NSCLC, fibroblast study), many scBaseCount lung studies with no patient labels (SRP347260 = GSE189357, 9 tumour samples, 126k cells, patient unknown). No LUAD atlas
with harmonised malignant labels. Candidates (from knowledge, not verified here, nothing downloaded):
1. **Salcher 2022 LuCA** (Cancer Cell; CELLxGENE collection `1e6a6ef9`, also in the Census scan `census_scan_2026-10-11.md`): about 1.28 M cells, 29 datasets, 318 patients (LUAD,
   LUSC, others), author tumour-cell annotation; h5ad several GB; mixed chemistries, so a median UMI >= 1,500 filter per sample is needed. Best candidate.
2. **Laughney 2020** (Nat Med, GSE123904; ~ 50 LUAD, 10x, ~ 40k cells) and **Wu F 2021** (Cell Res, GSE148071; 42 NSCLC, 10x v2, ~ 90k cells, malignant annotated), 5-20k UMI typical.
3. **Zhu 2022 / Bischoff 2021** (GSE189357 AIS-MIA-IAC, on disk without patient labels; Bischoff EGA, controlled).

## 3. Earlier work on cell cycle x clock (archive `scCircadianCancer`)
- **A24 cancer arm** (`docs/cancer/ANALYSIS_PLAN.md`; `paper/results/figure_5/scripts/cancer/compute_cycling_a24.py`): phase of S, G2M, cycling epithelial cells minus G1 on the stromal host (joint model); pooled
  -0.49 h [-1.02, 0.03] (7 studies); tumour minus normal uninformative.
- **T1-T3, C1** (2026-10-04 amendment; `compute_cycling_scores.py`): G2M-hi vs non-cycling offset -0.29 h (-0.57, 0.00), p 0.053 (12 studies), not without ERBB2-high (p 0.30); T2 RORE minus E-box +0.06 p 0.025;
  T3 G2/M fraction vs host NR1D1 slope -0.02 ns. Score-based (Tirosh), pseudobulk per patient x state.
- **CoPhaser C2-C4, C2b** (`compute_cycling_cophaser.py`; archive `scripts/s3/s4/s5/s6`; outputs `paper/results/figure_5/output/cancer/{c2_c4_*,c2b_*,c3_genes*,decision}.csv`): one CoPhaser model per cohort
  (cell-cycle mode, 98 packaged cell-cycle genes, 2,000 HVG context, RPE1 Fourier warm start, `auto_hyperparameters(cycle="cell_cycle")`, clock and clock-like genes removed, <= 3,000 cells/sample, train
  <= 80k cells; the CoPhaser repo is at commit 93861a8, version used then not recorded in the JSON). Pre-registered decision rule: CoPhaser G2/M vs non-cycling -0.02 h [-0.55, 0.51] (9 studies) = fail;
  C2b (post hoc, landmark window) -0.10 h, ns; C3 half-life gradient rho -0.06 p 0.83; **not confirmed**. Clock genes along the CoPhaser cycle in 17 studies: 59 % of gene x study pairs vs 53-60 % for pseudo-clock genes (ns).
  CoPhaser tumour vs normal epithelial phase -0.22 h but pseudo-clock -0.09 h (not clock-specific).
- **Caveats on the fits** (`outputs/states/cc/*.json`, `cc_diag/`): the cycling (non-cycling branch) call was off or degenerate in 3 studies; GSE183904, SRP327758 and SRP539998 score < 45 % G1; for GSE183904 and SRP327758 the CNV-malignant cycling_logit
  is > 0 in 100 % of cells (branch off: no usable cycling/non-cycling split); GSE183904 phase fails the S -> histone -> G2M ordering and was dropped in C2b. Phase `u` is oriented post hoc from the G2M-score peak.
- **Bulk C4 (A179)**: not explained by linear G1/S and G2M composition (0 of 5, `paper/exploratory/A179_c4_cellcycle`). Literature: Bieler 2014, Feillet 2014, Droin 2019 (`paper/docs/literature_positioning.md` section 2).
- **Consequence for registration**: the archive already ran clock-gene tests across cell-cycle phase (per-cell offsets, profile amplitude vs pseudo-clock sets) on these cells. The new pseudobulk-by-bin test is not blind
  to that history; the registration must state it, and the old results were null / not clock-specific, which sets the prior (small effect) and argues for a pseudo-clock pool matched to cell-cycle sensitivity.

## 4. Power sketch (cell counts only)
- Malignant >= 300 cells per patient: **153 of 251 patients** (all 29 cohorts), 49 of 79 in the A87 set; Wu 14, Pal 19 (ER+ 13 / TNBC 8 incl. 2 small), LUAD (Kim) 8.
- With CoPhaser tables available now (9 CNV cohorts): 110 patients >= 300 malignant, 94 of them with >= 150 cycling cells (cycling_logit > 0; unreliable in GSE183904 / SRP327758 where the branch is degenerate).
  Real 4-6 bins x >= 50 cells need >= 200-300 cycling cells per patient if only cycling cells are phase-binned (non-cycling G1 as an extra bin) (rough estimate, not computed).
- Most patients cluster in CRC (GSE178341 56 of 110); the cohort-level strata are needed for the high-proliferation contrast (two arms of about 30-50 patients each), ER and tissue confounded with cohort.
- Breast: Wu + Pal 33 patients >= 300 (ER+ 19, ER- 14) before the CoPhaser fits; the 2 x 2 (ER x proliferation) arms are about 8 each.
- Pseudobulk depth per bin: 50-cell bins at 10-20k UMI/cell give 0.5-1 M UMI, but a clock gene is detected in 2-30 % of cells, so about 1-15 UMI per gene per 50-cell bin (TEF, CIART, ARNTL lowest): the bins are Poisson-limited, so a clocksync read per bin needs shrinkage across bins or the gene pool; Run an MDE step (permuted cell-cycle labels within patient, pseudo-clock sets) before registration.
- Mitosis / library-size artefacts (cycling cells have different UMI, histone share, G2/M shutdown of short-lived mRNAs: NR1D1) must be controlled with the same pseudo-clock pool and the C3 half-life gradient.

## Files
`per_study_malignant_counts.csv`, `per_study_cell_level_detection.csv`, `per_patient_malignant_cells_in_cycle_tables.csv`, `per_patient_malignant_cycling_counts.csv`, `per_study_cycling_counts.csv`,
`wu_per_patient.csv`, `wu_cell_level_detection.csv`, `pal_per_patient_malignant.csv`, `malignant_per_patient_pseudobulk_meta.csv`.
