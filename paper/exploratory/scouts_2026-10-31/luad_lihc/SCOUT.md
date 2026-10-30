# LUAD / LIHC independent bulk cohorts: scout (metadata only, 2026-10-31)

Blinded: no clock-gene expression value was read, printed or used. Expression matrices were opened only for whole-matrix min/median/max
(`scale_diagnostics.csv`) and probe counts. Platform coverage = annotation tables only (`platform_coverage.csv`; GENES from `a127_lib`, Hallmark G2M = 200
genes, purity = 11 genes). All numbers below are in `series_table.csv`, `platform_coverage.csv`, `scale_diagnostics.csv`, `meta/*.csv`.
Scripts: `scripts/parse_hdr.py`, `platform_coverage.py`, `meta_tables.py`, `series_table.py` (run with `python -I`).

Data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_luad_lihc/<SERIES>/` (series matrices + PROVENANCE.md with sha256; `platforms/` annotations).
Scout folder 921 MB; `du -sh` of cancer_circadian: 365 GB before, 366 GB after.

## Per-series table (full text of annotations and overlap in `series_table.csv`)
| series | cancer | n tumours | n normal | platform | scale | clock genes | G2M | purity | downloaded |
|---|---|---|---|---|---|---|---|---|---|
| GSE72094 | LUAD | 442 | 0 | GPL15048 Rosetta/Merck Affy custom | log2-like, IRON + RNA-quality PLS corrected | 14/14 | 189/200 | 11/11 | yes |
| GSE68465 | LUAD | 443 | 19 (Stratagene reference) | GPL96 U133A | linear MAS5 (max 130623) | 13/14 (no CIART) | 190/200 | 11/11 | yes |
| GSE31210 | LUAD | 226 | 20 | GPL570 U133Plus2 | linear MAS5 (max 53814) | 14/14 | 192/200 | 11/11 | yes |
| GSE41271 | LUAD | 183 (of 275 NSCLC) | 0 | GPL6884 Illumina WG6v3 | log2-like | 14/14 | 192/200 | 11/11 | yes |
| GSE30219 | LUAD | 85 (of 293 tumours) | 14 | GPL570 | not checked | 14/14 | 192/200 | 11/11 | no |
| GSE50081 | LUAD | 128 (of 181 NSCLC) | 0 | GPL570 | not checked | 14/14 | 192/200 | 11/11 | no |
| GSE11969 | LUAD | 90 (of 163) | 5 | GPL7015 Agilent 21.6K | not checked | 12/14 (no NR1D2, CIART) | 151/200 | 11/11 | no |
| GSE14520 | LIHC | 247 (225 GPL3921 + 22 GPL571) | 241 (239 paired non-tumour + 2 healthy) | GPL3921 / GPL571 Affy U133A | log2-like | 13/14 (no CIART) | 191/200 | 11/11 | yes |
| GSE25097 | LIHC | 268 | 249 (243 adjacent + 6 healthy; plus 40 cirrhotic) | GPL10687 Rosetta/Merck Affy | NOT log2-like (min 0.005, median 0.11, max 66) | 14/14 | 174/200 | 10/11 | yes |
| GSE36376 | LIHC | 240 | 193 adjacent | GPL10558 Illumina HT-12v4 | log2-like | 14/14 | 192/200 | 11/11 | yes |
| GSE63898 | LIHC | 228 | 168 (cirrhosis, not paired normal) | GPL13667 Affy U219 | log2-like (RMA) | 14/14 | 184/200 | 11/11 | yes |
| GSE76427 | LIHC | 115 | 52 adjacent | GPL10558 | not checked | 14/14 | 192/200 | 11/11 | no |
| GSE54236 | LIHC | 81 | 80 cirrhotic | GPL6480 Agilent | not checked | 14/14 | 191/200 | 11/11 | no |
| GSE10143 (GSE10141 = its 80-tumour subset) | LIHC | 80 | 307 | GPL5474 DASL 6k | not checked | 9/14 | 140/200 | 10/11 | no |
| GSE64041 | LIHC | 60 | 65 | GPL6244 HuGene 1.0 ST | not checked | 14/14 | 190/200 | 11/11 | no |
| GSE39791 | LIHC | 72 | 72 | GPL10558 | not checked | 14/14 | 192/200 | 11/11 | no |

Not retrievable / not pursued (no download; details unverified, from a web search only):
- LIRI-JP (ICGC): the RIKEN liver RNA-seq is listed by NBDC as controlled access (hum0158, JGA); I could not confirm an open tier for expression. Treat as blocked.
- CPTAC-LUAD (Gillette 2020, about 110 tumours: below the 150 target) and CPTAC-HCC / Gao 2019 (159 HBV-HCC, paired): RNA-seq via PDC/CPTAC portal or NODE (registration); not checked further.
- OncoSG (EGA, controlled) and Chen 2020 East Asian LUAD (about 100 patients): not pursued.

## Annotations as natural subgroups (downloaded series)
- GSE72094: KRAS 154 mut / 288 WT, EGFR 47 / 395, STK11 68 / 374, TP53 111 / 331; smoking Ever 335 / Never 33 / missing 74; stage, age, sex. No TRU/PP/PI label deposited (would have to be derived from expression, which is not blinded).
- GSE68465: four sites (MI 178, MSKCC 104, DFCI 82, HLM 79), stage (pN/pT), grade, smoking history (Never 49 of 443), adjuvant chemo/RT, survival. No mutations.
- GSE31210: stage I-II only; EGFR 127, KRAS 20, ALK-fusion 11, triple-negative 68; ever/never smokers 123/123; Brinkman index; 15 of the 20 normals labelled patient-linked.
- GSE41271: histology (ADC 183), stage, tobacco history, sex, race. No mutations.
- GSE14520: HBV throughout (viral status CC 160, AVR-CC 58, N 6, 18 unknown), cirrhosis Y 223 / N 19, BCLC / TNM / CLIP, AFP, tumour size, multinodular, survival; paired non-tumour for 239; two array batches (platform column).
- GSE25097: tissue class only (no clinical fields in GEO).
- GSE36376: AJCC stage (over 433 arrays: 1 184, 2 179, 3a 44, 3b 17, 3c 9); adjacent non-tumour. Aetiology not in GEO fields.
- GSE63898: sex, BCLC (A 178, 0 19, B 17, C 11), comparator is cirrhosis. Aetiology not in GEO fields.
Proliferation (G2M) and purity are computable on all downloaded platforms; coverage is 174-192/200 G2M. Missing clock genes: CIART on GPL96 and GPL3921/GPL571.

## Overlap notes
- None of the series is expected to contain TCGA samples (independent institutional cohorts; no TCGA IDs in the headers). Not provable from metadata alone.
- GSE10141 is the 80-tumour training subset of GSE10143. GSE63898 and GSE10141/10143 cite the same meta-analysis (PMID 31344396): patient overlap not verifiable.
- GSE68465, GSE72094, GSE31210, GSE41271 come from different institutions (multi-site, Moffitt, NCC Japan, MD Anderson). GSE68465 contains four sites: use site as a batch/stratum. GSE30219 and GSE50081 are also different cohorts.
- GSE14520 mixes two platforms (GPL3921, GPL571): analyse the 225-tumour GPL3921 batch, or per platform.

## Recommendation
- LUAD primary: GSE72094 (442 tumours, single platform, 14/14 clock genes, KRAS/EGFR/STK11/TP53 calls, smoking, stage; processed log2) and GSE68465 (443 tumours; linear MAS5, CIART missing, multi-site so site-stratified; the most independent platform). Descriptive: GSE31210 (226, stage I-II, EGFR-heavy, single centre), GSE41271 (183 ADC, Illumina; histology filter needed), GSE30219 / GSE50081 (below 150 ADC).
- LIHC primary: GSE14520 (225 tumours on GPL3921 + 22 on GPL571; HBV, cirrhosis annotated, paired non-tumour; CIART missing) and GSE36376 (240 tumours, Illumina single platform, all 14 clock genes, adjacent non-tumour; limited annotation). Descriptive: GSE63898 (228, HCV-cirrhosis background, Affy U219), GSE25097 (268 but scale unclear and weakest annotation; check scale before use), GSE76427 (115).
- Aetiology split is not available inside a series in GEO fields except GSE14520 (HBV-only); LIHC aetiology contrast would need external cohorts.

## Blockers and caveats
- No open expression for LIRI-JP, CPTAC-HCC or OncoSG found in this scout; counts for them are from a web search and unverified.
- GSE25097 values are not on a recognisable log2 scale: a transformation rule must be registered before scoring.
- GSE68465 and GSE31210 are linear MAS5: log2 with a floor (as A168 `expression()`) before use.
- Illumina (GSE41271, GSE36376) and Affymetrix probe-averaging per symbol follows the A166 convention (`prep_a168.py expression()`); only GPL annotation tables were needed from `platforms/` (GPL96, GPL13667, GPL6244, GPL15048, GPL10687, GPL5474, GPL7015 as acc.cgi text).
