# C6 blind GEO cohort scout (metadata only), 2026-11-01

Scope: HANDOFF_2026-10-31 item 6 (C6). Only GEO accession text pages (series and GSM brief views), FTP directory listings
and platform title pages were read. No series matrix, count table, RAW.tar or CEL was downloaded or opened; no expression
value was seen. Sizes are from the FTP listings (matrix = `<GSE>_series_matrix.txt.gz`, expression values in the file body).
Blind check = `find`/`ls` and `grep -ril` over `$SCCC_DATA` (names, PROVENANCE.md), `paper/`, and the archive
(`*.md *.py *.txt *.tsv`). "0 hits" means none of them contains the accession.
Parser family column: GPL570, GPL6947, GPL6884, GPL6244, GPL96, RNA-seq counts are in the family; everything else needs a new
probe-to-gene map (unchecked clock-gene coverage is stated per row; coverage cannot be verified without the platform table
or the matrix, so "expected" means from the platform type only).

Gaps up front: no qualifying blind cohort (>= 150 primary tumours, one platform, processed values on GEO) was found for
**pancreatic** or **endometrial**; **cervical** has only a DASL cohort (risky). Prostate has one in-budget cohort.

## Colorectal (TCGA value weak or null)

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE39582 | GPL570 (in family) | 566 colon cancers (443 discovery + 123 validation) | 19 non-tumoral mucosa | `GSE39582_series_matrix.txt.gz` (164 MB) | expected 14/14 (U133 Plus 2; CIART absent as in other GPL570 cohorts) | 0 hits | **take** (CIT, French multicentre, frozen) |
| GSE14333 | GPL570 (in family) | 290 primary CRC | none | `GSE14333_series_matrix.txt.gz` (81 MB) | expected as above | 0 hits | **take** (Jorissen, Melbourne; MAS5 quantile-normalised) |
| GSE17536 | GPL570 | 177 (Moffitt) | none | `GSE17536_series_matrix.txt.gz` (50 MB) | expected as above | 0 hits | optional third |
| GSE13294 | GPL570 | 155 | none | `GSE13294_series_matrix.txt.gz` (44 MB) | expected | 0 hits | skip: same group and processing as GSE14333, possible sample overlap not checked |
| GSE41258 | GPL96 | 390 samples, mixed (normal, adenoma, primary, metastasis) | yes | not listed | | 0 hits | skip: tumour count after exclusion unverified, mixed stages and sites |

## Pancreatic (TCGA value weak or null) — no qualifying cohort

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE71729 | GPL20769 (Agilent 4x44K, not in family) | 145 primary (+ 61 metastatic, 17 cell lines, 134 normals) | 46 pancreas + 88 distant | `GSE71729_series_matrix.txt.gz` (15 MB) | expected (Agilent genome-wide) | 0 hits | **near miss (145 < 150)**: take only if the PI accepts 145; metastatic/cell-line/normal samples must be filtered out by `tissue` |
| GSE62452 | GPL6244 (in family) | 69 | 61 adjacent | RAW only | | 0 hits | no (< 150) |
| GSE183795 | GPL6244 | 139 | 102 adjacent + 3 | `GSE183795_normalized_matrix.txt.gz` (19.6 MB) | | NOT BLIND: named in `scouts_2026-11-01/SCOUT.md` (metadata scout only) | dropped, and < 150 |
| GSE21501 | GPL4133 | 132 | | RAW only | | 0 hits | no (< 150) |
| GSE131050 | GPL16791/GPL18573 RNA-seq | 97 primary FF/FFPE/FNA + 66 biopsies, mixed protocols and 2 platforms | | 2 supplementary TSV (5-9 MB) | | 0 hits | no (no >= 150 on one platform, mixed sample types) |

## Endometrial (TCGA value weak or null) — no qualifying cohort

| GSE | GPL | n primary tumours | normals | processed file | blind check | recommendation |
|---|---|---|---|---|---|---|
| GSE17025 | GPL570 | 91 stage I | 12 atrophic | matrix (12 MB) | 0 hits | no (< 150) |
| GSE106191 | GPL570 | 64 carcinoma (+ 33 hyperplasia) | | matrix | 0 hits | no (< 150) |
| GSE2109 (expO) | GPL570 | unknown endometrial subset of 2,158 mixed-tumour samples | | per-batch tars, no single matrix listed | 0 hits | not usable without per-sample tissue annotation (would need the series-matrix header or GSM pages for 2,158 samples); not pursued |

On disk and therefore not blind: `endometrium_gse51981` (GSE51981, GSE6364, GSE141549).

## Prostate (TCGA value weak or null)

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE116918 | GPL25318 (Almac Prostate Disease Specific Array, custom FFPE-biopsy array; not in family) | 248 (diagnostic biopsies, radiotherapy cohort) | none | `GSE116918_series_matrix.txt.gz` (42 MB) | UNCERTAIN: custom array, 14/14 not verifiable from metadata; check on download | 0 hits | **take**, with the coverage caveat |
| GSE62116 | GPL5188 (Affymetrix Human Exon 1.0 ST, not in family) | 235 FFPE radical prostatectomy | none | `GSE62116_series_matrix.txt.gz` (1.6 GB) | expected (core exon-level, needs gene summarisation) | 0 hits | optional second (large file, exon-level matrix) |
| GSE46691 | GPL5188 | 545 FFPE RP | none | matrix has no values (20 KB); values only in `GSE46691_quantile_normalized.txt.gz` (5.7 GB) | | 0 hits | no: file alone exceeds the budget |
| GSE21034 | GPL5188 / GPL10264 (same 131 primaries on both) | 131 primary per platform | 29 adjacent | matrix GPL10264 51 MB, GPL5188 1.0 GB | | 0 hits | no (131 < 150) |
| GSE70768 | GPL10558 | 125 tumour + 74 benign | | 48 MB | | 0 hits | no (< 150) |

## Squamous controls

### Head and neck

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE65858 | GPL10558 (Illumina HT-12 v4; sister of GPL6947/GPL6884, not strictly in family) | 270 tumour arrays (253 labelled Primary; HPV+ and HPV- mix, Heidelberg/Uppsala-type multicentre) | none | `GSE65858_series_matrix.txt.gz` (20 MB) | expected 13-14/14 (HT-12 covers the clock genes in GPL6947) | 0 hits | **take** (filter to `tumor_type: Primary`) |
| GSE39366 | GPL9053 (Illumina GA) | 138 | | 12 MB suppl | | 0 hits | no (< 150, same group as GSE65858) |
| GSE41613, GSE42743, GSE27020 | GPL570/GPL570/GPL96 | 97, 74-103, 109 | | | | 0 hits | no (< 150) |

### Cervical squamous

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE44001 | GPL14951 (Illumina HT-12 WG-DASL v4, FFPE) | 300 early-stage cervical cancers (mostly squamous; histology column not verified) | none | `GSE44001_series_matrix.txt.gz` (43 MB) | RISKY: DASL, partial transcript coverage, clock gene detection unverified; low-expression genes (NPAS2, CRY2, NR1D2) often below detection in DASL | 0 hits | take only as a flagged optional; a failed gate here is uninformative |
| GSE63514 | GPL570 | 28 cancers of 128 samples (CIN, normal) | 24 | 36 MB | | 0 hits | no (< 150) |

### Oesophageal squamous

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE53625 | GPL18109 (Agilent lncRNA + mRNA V2.0, not in family) | 179 ESCC (paired) | 179 paired | `GSE53625_series_matrix.txt.gz` (103 MB) | expected (mRNA probes), unverified | FLAG: accession appears only as the relation line in `scouts_2026-10-31/gastric_lung_histology/hdr/GSE53624_series_matrix.hdr`; its sub-series GSE53624 (119 pairs) was scouted as "not downloaded" and no expression was read | **take, with the flag**; the PI may rule it not blind on the grep criterion, in which case oesophageal has no cohort |
| GSE47404 | GPL6480 | 71 | | | | 0 hits | no (< 150) |
| GSE23400, GSE45670, GSE75241 | | 53, 38, 30 | | | | | no (< 150) |

### Lung squamous (not GSE37745, GSE41271, GSE30219)

| GSE | GPL | n primary tumours | normals | processed file (size) | clock-gene coverage | blind check | recommendation |
|---|---|---|---|---|---|---|---|
| GSE157009 (Alliance cohort I) | GPL570 (in family) | 249 SCC, early stage, multi-institution | none | `GSE157009_series_matrix.txt.gz` (57 MB) | expected as GPL570 | 0 hits | **take** |
| GSE157010 (Alliance cohort II) | GPL570 | 235 SCC | none | `GSE157010_series_matrix.txt.gz` (65 MB) | as above | 0 hits | **take** (independent second cohort; GSE157011 is their SuperSeries, do not download it too) |
| GSE4573 | GPL96 | 130 SCC | | 6.3 MB | | 0 hits | no (< 150) |
| GSE43580 | GPL570 | 150 NSCLC, of which about 73 SCC | | 33 MB | | 0 hits | no (SCC subset < 150) |
| GSE50081 | GPL570 | 181 NSCLC, SCC minority | | | | NOT BLIND (named in `scouts_2026-10-31/luad_lihc`) | dropped |

## Exclusions (not blind or not usable)

- Not blind: GSE81089 (RNA-seq NSCLC, on disk and in `paper/`), GSE19188 (4 hits in `paper/`), GSE50081, GSE183795, GSE51981, GSE6364, GSE141549 (endometrium on disk); GSE53624/GSE53625 flagged above.
- Used by A173 and excluded by the brief: GSE37745, GSE41271, GSE30219.
- Below 150 tumours: GSE13294 (155, kept out for overlap), GSE21034 (131 per platform), GSE70768 (125), GSE71729 (145, near miss), all other series above.
- Mixed-type or no single-platform subset: GSE2109, GSE131050, GSE41258.

## Selected (default block)

| type | GSE | n | platform | download |
|---|---|---|---|---|
| colorectal | GSE39582 | 566 | GPL570 | 164 MB |
| colorectal | GSE14333 | 290 | GPL570 | 81 MB |
| prostate | GSE116918 | 248 | GPL25318 (Almac DSA) | 42 MB |
| head and neck | GSE65858 | 270 | GPL10558 | 20 MB |
| oesophageal squamous | GSE53625 (flag) | 179 | GPL18109 | 103 MB |
| lung squamous | GSE157009 | 249 | GPL570 | 57 MB |
| lung squamous | GSE157010 | 235 | GPL570 | 65 MB |
| pancreatic | none >= 150 | | | |
| endometrial | none >= 150 | | | |
| cervical | GSE44001 (optional, DASL risk) | 300 | GPL14951 | 43 MB |

Default block total: about 532 MB. Optional additions: GSE44001 43 MB, GSE71729 15 MB (145 < 150), GSE17536 50 MB,
GSE62116 1.6 GB; everything together is about 2.3 GB, under the 5 GB cap.

## Download commands

```bash
D=/mnt/max8tb/data/CoPhaser_data/cancer_circadian; G=https://ftp.ncbi.nlm.nih.gov/geo/series

# colorectal
mkdir -p $D/c6_colorectal_gse39582; curl -fL --retry 3 -o $D/c6_colorectal_gse39582/GSE39582_series_matrix.txt.gz $G/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz
mkdir -p $D/c6_colorectal_gse14333; curl -fL --retry 3 -o $D/c6_colorectal_gse14333/GSE14333_series_matrix.txt.gz $G/GSE14nnn/GSE14333/matrix/GSE14333_series_matrix.txt.gz
# prostate
mkdir -p $D/c6_prostate_gse116918; curl -fL --retry 3 -o $D/c6_prostate_gse116918/GSE116918_series_matrix.txt.gz $G/GSE116nnn/GSE116918/matrix/GSE116918_series_matrix.txt.gz
# head and neck squamous
mkdir -p $D/c6_hnsc_gse65858; curl -fL --retry 3 -o $D/c6_hnsc_gse65858/GSE65858_series_matrix.txt.gz $G/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz
# oesophageal squamous (blind flag, see table)
mkdir -p $D/c6_escc_gse53625; curl -fL --retry 3 -o $D/c6_escc_gse53625/GSE53625_series_matrix.txt.gz $G/GSE53nnn/GSE53625/matrix/GSE53625_series_matrix.txt.gz
# lung squamous (Alliance cohorts I and II)
mkdir -p $D/c6_lusc_gse157009; curl -fL --retry 3 -o $D/c6_lusc_gse157009/GSE157009_series_matrix.txt.gz $G/GSE157nnn/GSE157009/matrix/GSE157009_series_matrix.txt.gz
mkdir -p $D/c6_lusc_gse157010; curl -fL --retry 3 -o $D/c6_lusc_gse157010/GSE157010_series_matrix.txt.gz $G/GSE157nnn/GSE157010/matrix/GSE157010_series_matrix.txt.gz

# optional
mkdir -p $D/c6_cervical_gse44001; curl -fL --retry 3 -o $D/c6_cervical_gse44001/GSE44001_series_matrix.txt.gz $G/GSE44nnn/GSE44001/matrix/GSE44001_series_matrix.txt.gz
mkdir -p $D/c6_pancreatic_gse71729; curl -fL --retry 3 -o $D/c6_pancreatic_gse71729/GSE71729_series_matrix.txt.gz $G/GSE71nnn/GSE71729/matrix/GSE71729_series_matrix.txt.gz
mkdir -p $D/c6_colorectal_gse17536; curl -fL --retry 3 -o $D/c6_colorectal_gse17536/GSE17536_series_matrix.txt.gz $G/GSE17nnn/GSE17536/matrix/GSE17536_series_matrix.txt.gz
mkdir -p $D/c6_prostate_gse62116; curl -fL --retry 3 -o $D/c6_prostate_gse62116/GSE62116_series_matrix.txt.gz $G/GSE62nnn/GSE62116/matrix/GSE62116_series_matrix.txt.gz
```

Notes for the registration: every selected series was identified from titles, sample counts and sample annotations only;
tumour counts above come from GSM source-name/characteristics fields and should be re-confirmed on the downloaded header
before the amendment fixes n. Samples that are normals, metastases, cell lines or xenografts (GSE39582, GSE65858,
GSE53625, GSE71729) must be filtered out by annotation.
