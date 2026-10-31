# Gastric adenocarcinoma and lung histology (LUAD vs squamous) bulk cohorts: scout (metadata only, 2026-10-31)

Blinded: no clock-gene expression value was read, printed or used. Expression matrices were opened only for whole-matrix min/median/max (`scale_diagnostics.csv`).
Platform coverage = annotation tables only (`platform_coverage.csv`; GENES from `a127_lib`, Hallmark G2M = 200 genes, purity = 11 genes). All numbers are in
`series_table.csv`, `platform_coverage.csv`, `scale_diagnostics.csv`, `meta/*.csv`. Scripts (run with `python -I`): `scripts/download.sh`, `parse_hdr.py`, `platform_coverage.py`, `meta_tables.py`, `series_table.py`.
Data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_gastric_lung/<SERIES>/` (series matrix + PROVENANCE.md with sha256; `platforms/`). Scout data folder 555 MB; `cancer_circadian` 367 GB (cap 400 GB).
GSE41271 and the GPL570, GPL96, GPL6884 tables were copied from `scout_luad_lihc`, not re-downloaded. GPL8432 exists on ftp only as a 897 MB family SOFT: its probe table was extracted (`platforms/GPL8432.txt`) and the family file deleted (sha256 in PROVENANCE).

## A. Gastric adenocarcinoma
| series | n tumours | n normal | platform | scale (whole matrix) | clock | G2M | purity | downloaded |
|---|---|---|---|---|---|---|---|---|
| GSE62254 (ACRG) | 300 | 0 | GPL570 | min 0.57, med 1.48, max 4.41: NOT log2 (fits log10 MAS5, inference) | 14/14 | 192 | 11/11 | yes |
| GSE66229 (ACRG + normals) | same 300 | 100 | GPL570 | not downloaded | 14/14 | 192 | 11/11 | no: OVERLAPS GSE62254 |
| GSE15459 (Singapore) | 200 (192 adenocarcinoma after authors' exclusions) | 0 | GPL570 | linear MAS5-like (max 384695) | 14/14 | 192 | 11/11 | yes |
| GSE84437 (Yonsei) | 483 (433 with clinical fields) | 0 | GPL6947 HT-12 v3 | linear, raw-like (4.5 to 38106; 1.6% NaN) | 14/14 | 191 | 11/11 | yes |
| GSE26253 (Samsung, FFPE) | 432 | 0 | GPL8432 DASL | log2-like (6.6 to 15.0) | 13/14 (no PER3) | 178 | 11/11 | yes |
| GSE26899 (Korea Univ) | 93 | 12 | GPL6947 | not downloaded | 14/14 | 191 | 11/11 | no (<150) |
| GSE29272 (NCI) | 134 | 134 paired | GPL96 | not downloaded | 13/14 (no CIART) | 190 | 11/11 | no (<150) |
| GSE27342, GSE57303, GSE34942 | 80, 70, 56 | 80, 0, 0 | GPL5175 / GPL570 / GPL570 | not downloaded | | | | no (<150) |

Annotations in GEO fields (the finding that matters): none of the four >=150 series deposits Lauren, ACRG/TCGA molecular subtype, MSI or EBV. GSE62254/66229 carry only tissue and patient ID
(ACRG subtypes, Lauren, stage, EBV, MSI are in the publication supplement / authors, PMID 25894828 and 29725014; not verified, not downloaded). GSE15459: Lauren and outcome are in a supplementary
`GSE15459_outcome.xls` (not downloaded, content unchecked). GSE84437: age, sex, pT, pN, survival only (433 of 483). GSE26253: stage (2: 167, 3A: 111, 1B: 68, 4: 67, 3B: 19), recurrence.
Site, Lauren and stage are deposited only in GSE26899 (93 tumours: Lauren intestinal 59 / diffuse 31 / mixed 2; antrum 55, body 29, cardia 7, fundus 2; stage 1-4: 11/18/27/36) and cardia vs non-cardia in GSE29272 (62 / 72). Paired normals: GSE66229 (100), GSE29272 (134), GSE27342 (80), GSE26899 (12).
Overlaps: GSE62254 = GSE66229 tumours (use one). GSE15459 and GSE34942 are Singapore batches A and B (sample names differ, not provable). GSE26253 (Samsung) may share patients with ACRG (Samsung Medical Center): not verifiable. GSE84437 is Yonsei.

## B. Lung histology (adenocarcinoma vs squamous, same platform)
| series | ADC | squamous | other tumour | normal | platform | scale | clock | G2M | downloaded |
|---|---|---|---|---|---|---|---|---|---|
| GSE37745 | 106 | 66 | large 24 | 0 | GPL570 | log2-like (2.3 to 15.2) | 14/14 | 192 | yes |
| GSE30219 | 85 | 61 (SQC) | LCNE 56, basaloid 39, carcinoid 24, "SCC" 21 (small-cell, inferred), LCC 3, other 4 | 14 | GPL570 | log2-like (2.2 to 14.6) | 14/14 | 192 | yes |
| GSE41271 | 183 | 80 | 12 | 0 | GPL6884 Illumina | log2-like (1.5 to 15.2) | 14/14 | 192 | yes (copied) |
| GSE81089 (RNA-seq) | 108 | 67 | large cell/NOS 24 | 19 | GPL16791 HiSeq | series matrix has no expression table; FPKM only as supplementary tsv (not downloaded) | whole transcriptome | | matrix only |
| GSE50081 | 127 | 42 | 12 | 0 | GPL570 | not downloaded | 14/14 | 192 | no (squamous <60) |
| GSE3141 | 58 | 53 | 0 | 0 | GPL570 | not downloaded | | | no (<60) |
| GSE19188 | 45 | 27 | LCC 19 | 65 | GPL570 | not downloaded | | | no (<60) |
| GSE33532 | 40 | 16 | mixed 24 | 20 | GPL570 | not downloaded | | | no (<60) |

Oesophagus: no series has >=50 adenocarcinoma and >=50 squamous. GSE13898 has 75 EAC and no ESCC; GSE26886 has 21 EAC / 9 ESCC; GSE53624 (119), GSE23400 (53 + 51 on GPL96 + GPL97) and GSE20347 (17) are ESCC only. Not downloaded.
Squamous in GSE30219 is exactly 61 (SQC); the "SCC" label (21) is read as small-cell carcinoma from the CIT publication, not from GEO fields: confirm before pooling. Basaloid (39) is a separate class.
Overlaps: GSE37745 (Uppsala 1995-2005) and GSE81089 (Uppsala 2006-2010) are different patients by operation dates. GSE41271 (MD Anderson), GSE30219 (Grenoble) and GSE37745 are independent institutions. TCGA not expected in any.

## Recommendation
- A (gastric): primary GSE62254 (300 ACRG tumours, single platform, 14/14 clock, 192/200 G2M; subtype/Lauren/MSI/EBV must come from the published supplement, so check availability before registering any subtype contrast; scale must be fixed first because the deposited values are not log2) and GSE84437 (483, Illumina, linear, 14/14; no histology/subtype annotation, a second independent cohort for pan-gastric statements). GSE15459 (192 adenocarcinoma, Affymetrix, same platform as ACRG; Lauren likely in the supplementary xls) as third. GSE26253 is a weaker option (FFPE DASL, PER3 missing, 178/200 G2M). Use GSE62254 only, not GSE66229, to avoid double counting. Lauren-stratified tests would need GSE26899 (93) or the supplementary tables.
- B (lung histology): GSE37745 (106 ADC / 66 squamous, GPL570, RMA-scale, single institution) and GSE41271 (183 / 80, Illumina) are the clean pair. GSE30219 (85 / 61, plus normals, GPL570) is third but borderline on squamous n and has mixed neuroendocrine classes. GSE81089 (108 / 67 RNA-seq) would be a platform-independent check but needs the FPKM supplement (outside the allowed file types here). Minimum squamous n across the three arrays is 61, so power for a squamous-specific contrast is modest.
- Oesophagus: not feasible at >=50 per histology in GEO bulk series.

## Not verified
- ACRG subtype/MSI/EBV/Lauren availability and coding (not in GEO fields); GSE15459_outcome.xls contents; whether "SCC" in GSE30219 is small-cell; GSE62254 scale (log10 inferred from range); possible Samsung/ACRG patient overlap with GSE26253; GSE84437 sub-series structure and the 50 samples lacking clinical fields; GSE81089 histology/stage/smoking codes (taken from the design text); coverage for GPL5175, GPL6102, GPL18109, GPL571, GPL97 (annotation not downloaded) and GPL16791 (RNA-seq).
- Normals in GSE81089 (19) assumed paired (design text), matched status unchecked.
