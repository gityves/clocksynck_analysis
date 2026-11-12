# Scout 2026-11-12: cancer datasets with clock-time information

Verified = read from the source page / GEO record / file listing today. Inferred = stated as such. Data root `$SCCC_DATA` = `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/` (420 G before, 427 G after; about 7 GB downloaded). No expression statistic was computed.

## 1. Diamantopoulou et al., Nature 2022 (CTC, Smart-seq2)
| field | value |
|---|---|
| accession | [GSE180097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180097), PMID 35732738 (verified) |
| organism / assay | human CTCs from 1 patient (GEO donor "Patient") and two human CTC-derived xenografts (BR16, LM2) in NSG mice; Smart-seq2, 306 samples (verified) |
| samples | patient 36 samples (single CTC / CTC cluster / CTC-WBC cluster); BR16 and LM2 xenograft samples (single, cluster, WBC-attached; some are pools of 10 cells) |
| time info | field `timepoint:` = "resting" / "active" for patient and for xenograft at 10h/22h; sample title carries the clock hour (xenograft "at 10h", "at 22h"; pooled LM2 samples "at 02h/06h/10h/18h/22h"; some timepoint fields are clock numbers 200, 600, 1000, 1800, 2200); patient samples are titled "day_N" / "night_N" (title only). No ZT label in GEO: 10h and 22h match ZT4/ZT16 of the paper only by inference. Patient identity per sample: only one patient in GEO; xenograft identity in `donor:`. Mouse identity is only in the sample title string (e.g. "from mouse 0200_1" for pools). |
| files | `GSE180097_RAW.tar` (36.6 MB): one `*.counts.txt.gz` per sample; no annotation file (metadata via family soft) |
| access / licence | open, no licence stated |
| downloaded | `diamantopoulou2022_gse180097/` 37 MB |

## 2. Wang et al., Cell 2024 (Scheiermann lab; GEO title "...dictate a temporal response to immunotherapy")
| field | value |
|---|---|
| accession | [GSE260641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260641), PMID 38723627 (verified); code: github.com/zqun1/circadian_immune_mouse_melanoma (verified, 10 R scripts, no data) |
| organism / assay | mouse, B16-F10-OVA melanoma in C57BL/6N females 8-9 wk; 10x scRNA-seq of live CD45+ TILs; 8 samples (verified) |
| time info | `time: ZT1/ZT7/ZT13/ZT19`, 2 replicates each (also in title) |
| files | single processed Seurat object `GSE260641_seu_obj_all_merged_filtered_processed_pop_def_man_1_alt_obj.rds.gz` (1.24 GB); no counts matrices, no raw. Seurat metadata includes `timepoint` and `mouseID` (inferred from the GitHub Step0 script) |
| not in GEO | endothelial data, Bmal1 endothelial / T-cell deletion data: no GEO record found (search of GEO by author/topic; the paper's Data availability statement could not be read: Cell 403, PMC blocked). Human data: the code uses TCGA-SKCM and zeitZeiger-predicted times (inferred from code and press), i.e. no recorded clock times |
| access | open, no licence stated |
| downloaded | `wang2024_gse260641/` 1.2 GB |

## 3. Human tumour cohorts with recorded clock time of collection (other than GSE233242)
| candidate | finding |
|---|---|
| REMAGUS02 [GSE26639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26639) | 226 patients, Affymetrix U133 Plus 2 (verified). Sample fields: tissue, ER, PR, HER2, grade, centre only. **No time field** (verified from the family soft). |
| Vasilyev et al. 2026 bioRxiv ([10.64898/2026.03.15.711880](https://doi.org/10.64898/2026.03.15.711880)) | Verified in text: six cohorts (REMAGUS 226, Mainz 200, Rotterdam 286, Transbig 198, Unt 125, Upp 251; 1,286 total); "Only 108 samples had sampling time available, all from the REMAGUS trial". The main text carries no data-availability statement naming a repository for the times; SI Notes (not read) hold dataset details. **Times not found public**; ask the authors (Vasilyev / Rand, Warwick; Giacchetti, Levi for REMAGUS). |
| TimeTeller Vlachou et al. 2024 PLoS Comp Biol ([10.1371/journal.pcbi.1011779](https://doi.org/10.1371/journal.pcbi.1011779)) | Statement quoted: "All data used except that of Bjarnasson et al. is publicly available from the sites given in the Supplementary Information"; Bjarnason oral mucosa data ship with the R package github.com/VadimVasilyev1994/TimeTeller. Oral mucosa of healthy volunteers (also Boyle smokers, Feng normal/dysplasia/OSCC without verified times). No tumour data with clock times. |
| GEO text searches (esearch, 2026-11-12) | "time of surgery", "time of day", "collection time", "time of biopsy", "sampling time", "time of collection", "resection time", "circadian" combined with cancer/tumor/Homo sapiens: scanned all hits for the narrow phrases and the circadian list (97). **No human tumour series with a per-sample clock-time field except GSE233242.** "Circadian" hits are cell lines/organoids (GSE262627 PDAC lines + 6 organoids, dex time courses; GSE119668 U2OS; GSE221103/GSE237608 neuroblastoma lines), GSE11469 (healthy oral mucosa, 8 men, 24 h, not tumour), TCGA analyses. The narrow phrases matched mostly unrelated series (peri-surgical biopsy handling GSE105777 / GSE73235-7 concern sampling method, not clock time). Limit: GEO free-text search is incomplete; sample-level fields are not all indexed. |
| literature | PeerJ 2018 (Ruben et al.): "For all data from TCGA and GEO, the times of day of surgery are not available." (search snippet, not verified in full text). |
| CPTAC / TCGA | **Not verified.** No evidence found that collection clock time is released; I did not inspect the CPTAC clinical dictionaries. |

Negative: nothing new to download for item 3.

## 4. scRNA-seq cohorts with metastasis and normal tissue of the host organ, same patients
| cohort | cells | patients x sample types | counts / annotation | access | downloaded |
|---|---|---|---|---|---|
| Liu et al. 2022 Cancer Cell, [GSE164522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522) | sorted CD45+ only (metadata: 3,940 CD45- cells of 257,843) | metadata has 10 patients (patient08-17), each with PBMC, metastasis tumour (MT), metastasis-adjacent normal liver (MN), primary tumour (PT), primary-adjacent normal colon (PN); 9 also lymph node. Cells: PBMC 46,973; PN 29,795; PT 37,529; MN 58,664; MT 37,359; LN 47,523 (verified from metadata csv). GEO lists 17 patients and 101 samples in total (series says; the released metadata covers the 10). | six per-tissue count csv (0.4-0.8 GB each), metadata csv with celltype_global/major/sub, VDJ | open | yes, `liu2022_gse164522/` 3.8 GB |
| GSE315534 (CRC primary, adjacent normal colon, liver mets; PMID 42082451, released Apr 2026) | all cells (no sorting stated; unverified) | 6 patients: T and N colon for 6, liver metastasis for 3. **No normal liver.** | 10x matrices per sample (15 samples), no annotations | open | yes, `hu2026_gse315534/` 0.84 GB (first author not checked) |
| Che et al. 2021, [GSE178318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178318) | all cells (140,281 barcodes) | 6 patients: primary CRC + matched liver metastasis; 3 PBMC; **no normal liver or colon**; treated (8) / naive (7) | one pooled 10x matrix, no annotation file | open on GEO, "for usage terms contact the PI (Jianming Li)" | yes, `che2021_gse178318/` 0.52 GB (kept as the largest all-cell CRLM set) |
| Wu et al. 2022 Cancer Discov (10.1158/2159-8290.CD-21-0316) | scRNA + spatial | no GEO accession found. A NODE (biosino.org) project OEP001756 lists 12 human CRLM scRNA samples (from a search result; the page itself did not render, so access, content and normal liver are **unverified**) | unknown | likely application-gated (unverified) | no |
| GSE298084 (PMID 40828155) | CD45+ colon and liver + PBMC | 3 synchronous CRLM patients, peritumour colon and peritumour liver (not a tumour sample) | 10x h5, no annotation | open | yes, `gse298084_crlm_macrophage/` 0.15 GB |
| GSE235057 | immune cells, HCC and CRLM | 10 HCC patients with tumour+adjacent liver; CRLM tumours only (no paired normal in GEO) | RAW.tar 0.96 GB; raw in HRA004767 | open | no |

Best option for the same-patient design: GSE164522 (immune compartment only). No all-cell cohort with normal host liver was found.

## 5. Human blood series dissociating circadian phase from sleep-wake / wall time
| accession | design | samples | time fields (verified) | melatonin phase | files |
|---|---|---|---|---|---|
| [GSE48113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48113) Archer 2014 PNAS | forced desynchrony 28 h, leukocyte arrays GPL15331 | 287 samples, 22 subjects (14 each, 7 per day 1 and day 4) | `time sample taken` (clock, 174 values), `day sample was taken in protocol` (1/4), `time point` (1-7), `sleep timing` = in phase / out of phase with respect to melatonin (2 levels), `per3 genotype` (4/4, 5/5), `subject` | **Only the binary in/out-of-phase label; no DLMO-relative time per sample in GEO** | series matrix (73 MB) + soft; RAW.tar 2.6 GB not downloaded |
| [GSE107537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107537) Kervezee 2018 PNAS | simulated night shift (day- vs night-oriented schedule), PBMC arrays GPL24324 | 103 samples, 8 subjects | `subject`, `condition` (Day/NightOrientedSchedule), `relclocktime` (hours relative to habitual wake/schedule, 42 values) | **none in GEO** (DLMO not a field; relclocktime is clock-relative) | series matrix 15 MB + soft 32 MB; RAW.tar not downloaded |
| [GSE122725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122725) Kervezee 2019 | same protocol with bright-light vs control group | 141 samples, 11 subjects (control 103, bright 38; the control group likely overlaps GSE107537, inferred from counts) | as above + `experimentalgroup` | none in GEO | series matrix 16 MB + soft 36 MB |
| GSE39445 Moller-Levet 2013 | on disk `blood_timed/GSE39445` | skipped | | | |

Downloaded: `archer2014_gse48113/` 189 MB, `kervezee2018_gse107537/` 47 MB, `kervezee2019_gse122725/` 52 MB. Per-subject DLMO values are in the papers' supplements (not checked).

## 6. Mouse tumour datasets
| study | what exists (verified unless stated) | access |
|---|---|---|
| Hadadi 2020 Nat Commun ([10.1038/s41467-020-16890-6](https://doi.org/10.1038/s41467-020-16890-6)) | PyMT mice, chronic jet lag (10 wk) vs LD; mRNA-seq of dissociated primary tumour and bone-marrow mononuclear cells, 5 JL + 5 LD mice; ENA [PRJEB33802](https://www.ebi.ac.uk/ena/browser/view/PRJEB33802) (19 runs, FASTQ). Methods: LD mice sampled at ZT3-ZT4, "in the case of JL mice, ZT was not followed": single time of day, no time series. | raw FASTQ only in ENA (not downloaded: rule); FPKM in paper supplement (unverified) |
| Aiello 2020 Sci Adv ([10.1126/sciadv.aaz4530](https://doi.org/10.1126/sciadv.aaz4530)) | murine melanoma, circadian disruption; cytokines, flow, qPCR | Data statement: "present in the paper and/or the Supplementary Materials"; no sequencing deposit |
| Liang et al. bioRxiv 2026 "Circadian misalignment underlies immune escape in breast cancer" ([10.64898/2026.05.26.726543](https://doi.org/10.64898/2026.05.26.726543), posted 2026-05-29; Myles Brown / Kellis labs) | exists. 4T1 orthotopic tumours, BALB/cJ, 12:12 LD, CT4/10/16/22, n = 4 tumours per time point, snMultiome (101,321 nuclei). Human TNBC: CYCLOPS2.0-inferred phases on Zhang et al. scRNA-seq (GSE266919, 76 samples) and Wu 2021 (GSE176078); no recorded human clock times. Data statement: GEO deposit "will be publicly released as of the date of publication"; processed data URL https://compbio.mit.edu/Circadian_BRCA/ returned an empty directory on 2026-11-12. | **not yet public** |
| Masri 2016 Cell ([GSE73222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73222)) | liver (not tumour) of WT and KrasLSL-G12D;p53fl/fl lung-tumour-bearing mice, ZT0/4/8/12/16/20 x 3 replicates, Agilent arrays, 36 samples | open; series matrix 6 MB, RAW 305 MB; not downloaded |
| Kettner 2016 Cancer Cell ([GSE75475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75475)) | WT liver, normal schedule vs chronic jet lag, 12 and 30 wk, ZT2/10/18 (2-3 mice per cell), 32 arrays; no tumour samples in the series | open; series matrix 3 MB; not downloaded |

## Not verified / limits
- Wang 2024 Data availability statement and key resources table (Cell 403, PMC blocked): claims about absent endothelial / Bmal1-deletion deposits rest on GEO search only.
- Vasilyev 2026 SI notes and any repository for the 108 REMAGUS sampling times.
- Wu 2022 NODE project content and access; CPTAC / TCGA collection time fields; time field content of GSE294391 / GSE193525 style blood-cancer series (not opened).
- GEO free-text searches are incomplete for sample-level fields.
- Diamantopoulou: mapping 10h/22h to ZT4/ZT16 is inferred, not read in GEO.
