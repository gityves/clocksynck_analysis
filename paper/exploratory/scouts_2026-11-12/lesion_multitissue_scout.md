# Scout 2026-11-12: paired lesion / non-lesion tissue (Topic A) and living multi-tissue cohorts (Topic B)

Scope: sample sheets and detection summaries only. No phase, rhythm fit or lesion-vs-normal comparison of any gene was computed. Data are under
`/mnt/max8tb/data/CoPhaser_data/cancer_circadian/paired_lesion/<GSE>/` and `.../living_multitissue/<GSE or E-TABM-1140>/` (PROVENANCE.md with URL, date 2026-11-12, sha256 per file in each). Downloaded: 245 MB (A) + 6.3 GB (B) = 6.6 GB.
Scripts: `scripts/gdssearch.py` (E-utilities search), `gsemeta.py` (header scan), `download.py`, `download_muther.py`, `scoutlib.py`, `build_samples.py`.
Outputs: `output/paired_lesion_<GSE>_samples.csv`, `output/living_multitissue_<GSE>_samples.csv`, `*_clockgenes.csv` (present / probes / median percentile per clock gene), `output/lesion_multitissue_summary.csv`, `output/living_multitissue_GSE24293-4-7_identity_check_sexgenes.csv`.

Definitions. "Complete pairs" = patients with at least one sample of each of the two statuses in the `pair_def` column of the summary. Median percentile = per sample, rank of the gene among all genes (probes collapsed to gene symbols by the mean; RNA-seq as log2 CPM, zero-count genes included in the denominator), then median over samples; the value in the tables is the median over the 14 clock genes of that per-gene number. Denominators differ across datasets (column "genes ranked"), so the percentiles are not comparable between platforms.

Clock time of collection: a keyword scan (time of day / hour / am-pm / morning / fasting / circadian / ZT) of every `!Sample_*` and `!Series_*` field of all 18 downloaded series matrices found nothing recorded (GSE320593 mentions dexamethasone timing only). No dataset has sampling clock time in GEO or in the MuTHER SDRF. Constant-offset cancellation is therefore what a clocksync analysis would have to rely on.

## Topic A: paired lesion vs non-lesion, outside skin-psoriasis/eczema, IBD and tumour-normal

Ranked by usable paired n, same tissue type, verified patient link, data type. "DL" = downloaded.

| # | Accession | Tissue and pair | n paired patients | Data type | DL |
|---|---|---|---|---|---|
| 1 | [GSE16134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16134) | gingiva, periodontitis site vs healthy site, same mouth | 69 (120 patients; 241 diseased, 69 healthy samples) | Affymetrix U133 Plus 2 array | DL |
| 2 | [GSE10334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10334) | gingiva, same design, same group | 64 (90 patients; 183 diseased, 64 healthy samples) | Affymetrix U133 Plus 2 array | DL |
| 3 | [GSE39491](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39491) | Barrett metaplasia vs squamous oesophagus and vs gastric cardia, same patient (premalignant) | 40 (three tissues each) | Affymetrix U133A 2.0 array | DL |
| 4 | [GSE57218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57218) | OA vs preserved cartilage, same joint (RAAK) | 33 (40 patients; 7 healthy donors unpaired) | Illumina HT-12 v3 array | DL |
| 5 | [GSE43292](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43292) | carotid plaque vs adjacent macroscopically intact artery | 32 | Affymetrix HuGene 1.0 ST array | DL |
| 6 | [GSE192352](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192352) (array twin: GSE192354) | leiomyoma (benign tumour) vs own myometrium | 30 | RNA-seq raw counts | DL |
| 7 | [GSE119717](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119717) | AAA: dilated vs non-dilated perivascular fat (not aortic wall); plus omental and subcutaneous fat | 30 (all male, Italy) | Illumina HT-12 v4 array | DL |
| 8 | [GSE26051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26051) | tendinopathy lesional vs non-lesional tendon | 23 donors, but different tendon in 19 of them | Affymetrix U133 Plus 2 array | DL |
| 9 | [GSE151243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151243) | hidradenitis suppurativa lesion vs perilesional skin | 20 | RNA-seq raw counts | DL |
| 10 | [GSE37837](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37837) | endometriosis: autologous eutopic vs ectopic (endometrioma) | 18 (subject id present) | Agilent array | no (below 20) |
| 11 | [GSE68309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68309) / GSE68372 | varicose vs normal vein, same patient | about 7 (14 samples in the expression series; GSE68372 titles it "20 paired" but its 20 samples are methylation) | Illumina array | no |
| 12 | [GSE311770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311770) | vitiligo skin, lesional / perilesional / healthy 1 mm biopsies | 16 patients with at least two sites (13 with three); 105 samples | RNA-seq | no (skin, below 20) |
| 13 | [GSE136825](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136825) | nasal polyp (42), polyp-side inferior turbinate (33), control turbinate (28) | not determinable: titles carry sample numbers only, no patient field | RNA-seq | no |
| 14 | [GSE141549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141549) | endometriosis lesions, endometrium, peritoneum (408 samples, two platforms) | not determinable: "SAMPLE n" is a sample number, no patient id | Illumina array | no |
| 15 | [GSE25628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25628) | endometriosis ectopic (8), eutopic (8), normal (6) | no patient id; below 20 | Affymetrix array | no |
| 16 | [GSE298871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298871) | vitiligo lesional vs non-lesional skin | 5 (17 samples) | RNA-seq | no |
| 17 | [GSE120521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120521) | stable vs unstable segments of the same carotid plaque | 4 | RNA-seq | no |
| 18 | [GSE99621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99621) | IPF scarred vs normal-appearing lung, same patient | 3 IPF patients (26 samples incl. controls) | RNA-seq | no |
| 19 | [GSE163154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163154) | carotid plaque with vs without intraplaque haemorrhage (MaasHPS, 43) | none: different plaques, not segments of one plaque | Agilent array | no |
| 20 | [GSE186476](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186476) | cutaneous lupus lesional / non-lesional (29 samples) | no patient id | RNA-seq | no |

Rejected on inspection: GSE7305 (endometriosis; 10 disease vs 10 normal samples from different women), GSE104140 (32 plaque sections by pathology class, no patient id), GSE44270 (keloid: cultured cells), GSE75152 (Dupuytren: unpaired control tissues), GSE51588 (OA subchondral bone, 40 OA vs 10 normal, unpaired), GSE338426 (nasal polyp vs separate controls).
Not found on GEO with the E-utilities queries tried (several wordings each): Crohn's stricture vs non-strictured ileum, COPD/IPF regions with at least 20 patients, diabetic-foot edge vs own intact skin, RAAK or other OA RNA-seq with patient-matched damaged/preserved cartilage, subchondral bone pairs, keloid vs own normal skin tissue (only fibroblast cultures), tendinopathy other than GSE26051. Other fibroid cohorts (GSE199849 n=136, GSE207350 n=108, GSE244187 n=18) were seen in search results but their pairing was not checked.

### Per-dataset notes (downloaded)

Clock detection (from `output/lesion_multitissue_summary.csv`):

| Dataset | samples | patients | complete pairs | clock genes present | median percentile | genes ranked |
|---|---|---|---|---|---|---|
| GSE16134 | 310 | 120 | 69 | 14 | 65.8 | 22836 |
| GSE10334 | 247 | 90 | 64 | 14 | 64.9 | 22836 |
| GSE39491 | 120 | 40 | 40 | 13 (CIART absent) | 50.5 | 13909 |
| GSE57218 | 73 | 40 | 33 | 14 | 68.7 | 19608 |
| GSE43292 | 64 | 32 | 32 | 14 | 66.4 | 21147 |
| GSE192352 | 62 | 30 | 30 | 14 | 76.1 | 25167 |
| GSE119717 | 119 | 30 | 30 | 12 (NPAS2, NR1D1 absent) | 49.2 | 11812 |
| GSE26051 | 46 | 23 | 23 | 14 | 66.6 | 22836 |
| GSE151243 | 40 | 20 | 20 | 14 | 89.0 | 55863 |

- GSE16134 / GSE10334: patient id parsed from "Patient N" in source_name. The two series have the same design and platform and look like one group's work; numbering is per series and whether patients overlap between them is not documented, so do not pool them without a check. Up to 3 sites per patient and status (replicates). Only 69 and 64 patients have a healthy site; the other patients are diseased-only. Age, sex, smoking are not in the sheet. Healthy vs diseased is a site classification (BoP, probing depth, attachment loss), sampled at one visit.
- GSE39491: sheet title "Case NC/BE/NE k" gives the case; 40 complete triplets in the sheet (the series text says 43 BE patients). Premalignant lesion, flag. Age and sex are not in the sheet.
- GSE57218 (RAAK): patient from "RAAK_n" in the title; 26 F / 14 M. Series-matrix values are normalised; non-normalised table also downloaded. Array chip/batch is not in the sheet, so batch versus OA status cannot be examined.
- GSE43292: all patients have the label "hypertension"; age and sex not recorded. Pairs are from the same endarterectomy. Plaque and intact tissue were both processed in the same series, batch not listed.
- GSE192352: 62 samples, 30 patients; 29 patients contribute one leiomyoma and one myometrium sample, patient 10 contributes three leiomyomas (ME10_S3, ME10_2_S17, ME10_3_S22) and one myometrium. Benign tumour, flag. Age 31 to 49. Two sheet titles do not match count-column names exactly (replicate suffix naming); the count table has 62 sample columns. GSE192354 is the same cohort on array and methylation platforms.
- GSE119717: perivascular adipose, not the aortic wall; "non-dilated" is the neck PVAT of the same aneurysm patient. 30 dilated, 30 non-dilated, 30 subcutaneous (arrays duplicated for subjects 50 and 51), 29 omental samples; all male, Italy. Also qualifies for Topic B (four depots, but an aneurysm cohort, not metabolic).
- GSE26051: donor field links the two samples, but the lesional and the non-lesional sample come from different tendons in 19 of 23 donors (e.g. donor 5: biceps lesional, brachialis non-lesional). That breaks the "same tissue type" requirement; usable only as a within-person, different-site contrast. 15 M / 8 F.
- GSE151243: perilesional (not distant) skin; 17 F / 3 M; counts use Ensembl ids (mapped with the local GENCODE v29). Skin, so lower priority. Percentile denominator contains 55,863 genes including unexpressed ones, hence the high value.

## Topic B: living multi-tissue metabolic cohorts

| # | Accession | Tissues (same person, same visit/operation) | n persons with at least 2 tissues | Metabolic status recorded | Data type | DL |
|---|---|---|---|---|---|---|
| 1 | [E-TABM-1140 (MuTHER)](https://www.ebi.ac.uk/biostudies/studies/E-TABM-1140) | subcutaneous fat, skin (punch biopsies), lymphoblastoid cell line | 856 women; 681 with fat + skin; 655 with fat + skin + LCL | no: sample sheet has sex, zygosity, twin id only | Illumina HT-12 v3, log2 normalised (processed zips) | DL |
| 2 | [GSE24293 / GSE24294 / GSE24297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24335) (union GSE24335) | liver (651), omental fat (848), subcutaneous fat (701) at RYGB, MGH 2000-2007 | cannot be established: person link is not in GEO | not in GEO | Agilent two-colour, log ratio vs pool | DL (all three plus the 2200-sample union) |
| 3 | [GSE73439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73439) | subcutaneous and visceral adipose (pregnancy / labour study) | 51 women (25 no labour, 15 labour, 11 non-pregnant) | no BMI or T2D | Affymetrix HuEx gene-level (exon-level series skipped, 219 MB) | DL |
| 4 | [GSE24883](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24883) | subcutaneous and visceral adipose | 32 women, 8 each lean / overweight / obese / obese with metabolic syndrome | group only | Agilent two-colour, log ratio vs reference pool | DL |
| 5 | [GSE133979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133979) | rectus abdominis muscle, subcutaneous adipose | 35 persons, 33 with both (24 with PDAC diagnosis; 11 without a diagnosis entry) | cancer/cachexia, not metabolic | variance-stabilised counts | DL |
| 6 | [GSE119717](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119717) | 2 perivascular depots, omental, subcutaneous fat (AAA) | 30 men | aneurysm | HT-12 v4 array | DL (stored under paired_lesion) |
| 7 | [GSE320593](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320593) | liver (26 patients) and subcutaneous fat (60 patients) | none: no patient has both tissues | no | RNA-seq raw counts | DL, but not a multi-tissue set (see notes) |
| 8 | [GSE136822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136822) | periaortic, omental, subcutaneous fat (peripheral artery disease) | 11 subjects (below 30) | PAD | Illumina array | no |
| 9 | [GSE145284](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145284) | subcutaneous + omental fat of obese girls at sleeve gastrectomy, plus follow-up after weight loss (24 samples) | below 30 | obesity | RNA-seq | no |
| 10 | GSE15773, GSE29718, GSE73108 | omental + subcutaneous (19, 20 and 12 samples) | below 30 | insulin resistance / obesity | arrays | no |

Same tissue before and after weight loss (flagged separately, not downloaded, only titles and sample counts checked): GSE19494 (88 samples), GSE72158 (84), GSE351786 (38), GSE66921 (8), GSE9007 (10). GSE62199 (201 samples) is overfeeding, not weight loss. GSE113823 is mostly a rodent multi-tissue RYGB set with a human arm of 5 pooled-patient pairs.
Excluded as single-tissue: GSE71416 (omental only), GSE16415 (visceral only), GSE135251 (liver only). GTEx excluded as post-mortem. The stomach series of the Greenawalt cohort (GSE24335 does not contain it) was not located.

### Per-dataset notes (downloaded)

| Dataset | samples | persons | persons with >=2 tissues / pair | clock genes present | median percentile | genes ranked |
|---|---|---|---|---|---|---|
| E-TABM-1140 | 2356 | 856 | 681 (fat + skin) | 14 | 59.6 (fat file) | 19573 |
| GSE24293-4-7 | 2200 | n/a | n/a (index-based 651 is invalid) | 14 (presence only) | not computed | 19523 |
| GSE73439 | 102 | 51 | 51 | 14 | 65.8 | 19932 |
| GSE24883 | 64 | 32 | 32 (unverified link) | 14 (presence only) | not computed | 19297 |
| GSE133979 | 68 | 35 | 33 | 13 (RORA absent) | 58.0 | 12908 |
| GSE320593 | 172 | 86 | 0 multi-tissue; 86 surgical + non-surgical | 14 | 70.9 (adipose file) | 19419 |

- Greenawalt (GSE24293, GSE24294, GSE24297, union GSE24335): the series give only "MGH_<tissue>_<k>" with k running 1..651 / 1..848 / 1..701 per tissue; there is no clinical field, no patient id and no linking file. I tested whether index k is the same person across tissues with sex genes (XIST minus the mean of 5 Y-linked genes, per sample; not a clock gene). Within each tissue the two components are strongly anti-correlated (r about -0.92 / -0.96 / -0.97), so the score is valid, but across tissues at the same index r is 0.053 (liver vs omental, 651 indices), 0.039 (liver vs subcutaneous, 651) and -0.012 (omental vs subcutaneous, 701), against permuted-index means of 0.004, -0.001, 0.004 (`output/living_multitissue_GSE24293-4-7_identity_check_sexgenes.csv`). The indices are therefore not person identifiers, and a same-person multi-tissue design cannot be built from GEO alone. A person key would have to come from the original MGH/Rosetta data (Greenawalt 2011; Schadt 2008) or from the authors. Data are two-colour ratios against a tissue pool, so within-sample percentiles are not meaningful (only gene presence was recorded).
- MuTHER (E-TABM-1140): processed per-individual matrices (replicate arrays merged) for fat, LCL and skin, downloaded as the three processed zips (left zipped, 0.8 GB); raw intensity files (3 GB) not downloaded. Person key = SDRF "Individual" (TWPID); 856 individuals; 826 fat, 705 skin, 825 LCL samples; 655 individuals have all three, 190 have two, 11 have one. All women (twins; zygosity and co-twin id given). No BMI, diabetes or NAFLD field in the sheet, so metabolic status would need the TwinsUK phenotype release. LCL is a cell line, not a tissue. Twins are correlated, so pair-level analyses need a twin-level block. Clock summary uses the fat file mapped via the GPL6947 (HT-12 v3) annotation.
- GSE73439: pregnancy and labour cohort (operation type not stated in the sheet); not a metabolic cohort; no age, no BMI. Subject id is explicit (51 subjects, one subcutaneous and one visceral sample each).
- GSE24883: GEO states paired SAT and VAT in 32 women, but the sheet has no person id; I paired by group plus the running index in the title (for example "Overweight_3" in both depots). That link could not be verified (all women, no sex-gene check possible) and the Greenawalt result shows that indices alone can mislead. Two-colour log ratio against a mixed liver, fat and muscle pool.
- GSE133979: patients undergoing abdominal surgery; muscle and adipose columns share the S/M/C+number prefix. 35 persons by prefix, 33 with both tissues; the diagnosis field is "Pancreatic ductal adenocarcinoma" for 24 and empty for the 11 "C" controls. 18 F / 17 M. Cancer cohort, cachexia focus; flag.
- GSE320593: the two biopsies of one patient are a non-surgical (percutaneous) and a surgical biopsy of the same tissue (title prefix is the patient code); 60 patients for adipose and 26 for liver, with no patient present in both. Perioperative dexamethasone is the study's variable (DEX vs NO_DEX), and sequencing batch is confounded with biopsy type (batch D contains only surgical samples; counts per batch and type are in the sheet, e.g. batch A 44 non-surgical vs 31 surgical, batch B 16 vs 29). Useful for same-tissue, same-person comparison, not for cross-tissue.

## Problems and cautions (summary)

- Pair timing: periodontitis sites, OA joints, plaque segments, myometrium/leiomyoma, Barrett biopsies and AAA fat depots are all taken at one visit or operation. GSE320593 and the MuTHER skin/fat punch biopsies are the only designs where two samples of one person may differ in timing, and no clock time is recorded for any of them.
- Missing or unreliable person keys: Greenawalt (index not a person id), GSE24883 (title index, unverified), GSE141549, GSE136825, GSE25628, GSE186476 (no patient field).
- Different tissue under the "paired" label: GSE26051 (different tendons), GSE119717 (perivascular fat, not wall), GSE39491 (control tissues are two other organs).
- Benign or premalignant lesions: GSE192352 (fibroid), GSE39491 (Barrett).
- Batch confounding that can be seen in the sheets: GSE320593 (biopsy type vs batch A-D). Array batch is absent from GSE57218, GSE43292, GSE10334, GSE16134, GSE26051, GSE39491 sheets, so it cannot be assessed from metadata.
- Clock-gene presence: all 14 present in 12 of the 15 downloaded datasets (the Greenawalt trio counted once); CIART absent in GSE39491 (U133A), NPAS2 and NR1D1 absent in the 11,812-gene ranked set of GSE119717, RORA absent in GSE133979 (VST table). Percentiles are within-dataset only.
- A search-process note: a first full run of the sample builder read the GPL5175 family SOFT file past the platform table and exhausted memory on the shared machine; I killed it, fixed the reader to stop at `!platform_table_end`, and reran everything under an address-space limit.
