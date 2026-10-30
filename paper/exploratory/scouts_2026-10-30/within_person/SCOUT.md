# Within-person disease designs: data scout (roadmap 6.6), 2026-10-30

Metadata only: sample tables and cell counts per sample x class. No expression value, clock gene or coupling was read.
Script: `scripts/scout_paired.py` (run `python -I -W ignore scripts/scout_paired.py`); full printout `scout_output.txt`; tables are the CSVs in this folder.
Eligible families: fibroblast_stromal (fib) and endothelial (endo). Threshold: >= 50 cells in the sample x class unit.
Paired patient = same patient with both conditions, each with >= 50 cells in the class.

## Ranking

| rank | candidate | design | n patients with both states | paired patients with fib >= 50 / endo >= 50 / both | eligible families per sample (of 2) | on disk |
|---|---|---|---|---|---|---|
| 1 | TAURUS (Thomas 2024), UC + CD | inflamed vs non-inflamed biopsies | 29 | 24 / 28 / 24 (patient x state pooled over sites); 23 / 28 / 23 (sample level) | 2 in 171 of 216 samples, 1 in 40, 0 in 5 | yes (raw counts, 12.7 GB) |
| 2 | Psoriasis sc (GSE173706, GSE228421, GSE278330) | lesional vs non-lesional | 30 (visits pooled), 31 (same visit) | 12 / 20 / 12 (pooled); 12 / 21 / 12 (same visit) | see section 2 | yes (per-GSM matrices) |
| 3 | Kong 2025 TI stroma (Census) | Crohn vs normal tissue, same donor | 6 | 4 / 3 / 3 | stroma only | pseudobulk (68 genes) only |
| 4 | Kong 2025 colon stroma (Census) | same | 4 | 2 / 2 / 2 | stroma only | pseudobulk only |
| - | GSE147424 AD (He 2020) | AD lesional vs non-lesional | not determinable from GEO metadata | not counted | not counted | yes (RAW tar, 17 samples) |
| - | GSE266546 (endoscopy / surgical IBD) | histology grades | 0 same tissue; 29 patients span tissues | stroma too sparse | see section 4 | annotation only |
| - | IBDverse E-MTAB-16986 / 16999 | one TI sample per individual | none | no cell counts | n/a | metadata only |
| - | bulk IBD (7 series) and bulk skin (about 45 series) | paired mucosa / skin | various | one compartment per sample | 0 (not eligible: no compartments) | yes |

## 1. TAURUS (Thomas et al. 2024, Nat Immunol; Zenodo 14007626)
- Source tables: `ibd/TAURUS/pseudobulk/samples.csv`, `class_counts.csv`; `PROVENANCE.md` in `ibd/TAURUS`.
- 987,743 cells, 216 samples (biopsies), 41 donors: 38 IBD (UC 108 samples, CD 96 samples) and 3 healthy (12 samples). Inflamed 104, non-inflamed 100.
- 29 patients have both states, 12 have one. Both states at the same timepoint: 21 patients. Pre- and post-treatment samples in the same patient: 34. Patient x site units with both states: 37 (inflamed and non-inflamed at the same site, at different timepoints).
- Timepoint x state: Pre 67 inflamed / 37 non-inflamed; Post 37 / 63. Inflamed vs non-inflamed is therefore partly confounded with treatment timepoint (and with site).
- Sites: descending colon 60, rectum 58, ascending colon 38, terminal ileum 36, sigmoid 24.
- Patients with fib and endo >= 50 in both states (pooled over sites): 24 (12 UC, 12 CD). Sample level: 23.
- Median cells per sample, inflamed / non-inflamed: fib 137 / 366, endo 252.5 / 153, mural 119 / 104, epithelial 716.5 / 2128.5, T 1168.5 / 576.5, myeloid 169 / 75.5. Fibroblasts are the thin class in inflamed biopsies.
- Class definitions: fib = LP, Epi, C3pos, THY1pos_FAPpos_PDPNpos and generic Fibroblast; mural = Pericyte; Myofibroblast excluded (as in A45).
- Separate libraries per sample: yes (sample ids CID...; chemistry 3' v3.1 164 samples, 3' V3 52; Batch column present). Raw counts: yes (X integer). Time of day: none in the metadata. Pretreatment flag: 112 False, 104 True.
- Licence: Zenodo record; licence text not read in this scout. Prior use: A45 (28 patients, 165 biopsies at the time; primary not met; stroma coupling 0.80-0.91 in both states per archive RESULTS_preliminary).
- Tables: `taurus_sample_class_cells.csv`, `taurus_patient_condition_cells.csv`.

## 2. Psoriasis single-cell (A58 pseudobulk, `psoriasis_sc/pseudobulk.csv`)
- 3 GEO series, 36 patients; samples are per patient x visit x status (libraries merged per sample in the pseudobulk). Patients per series: GSE173706 10, GSE228421 5, GSE278330 21. 68 sample units.
- Patients with both lesional and non-lesional: 30 (visits pooled), 31 (same visit).
- Samples with fib >= 50 / endo >= 50, lesional vs non-lesional: GSE173706 9 / 9 vs 8 / 9; GSE228421 5 / 5 vs 5 / 5; GSE278330 6 / 13 vs 1 / 8.
- Paired with both families: 12 patients (GSE173706 7, GSE228421 5); none from GSE278330 (non-lesional fibroblasts almost absent there).
- Median fib cells per sample: lesional 104.5, non-lesional 35.5. Mural 0 in the median sample (not eligible).
- Time of day: none known. Raw counts on disk (mtx / csv per GSM). Public GEO. Overlap with the psoriasis data already used in A54/A58 (not blind to the earlier analysis).
- Table: `psoriasis_sample_compartment_cells.csv`.

## 3. Kong et al. 2025 Crohn atlas (CELLxGENE Census pseudobulk, 68 genes)
- Donor-level pairing exists: donors labelled "Crohn disease; normal" (diseased and normal tissue of the same donor).

| dataset | donors | donors with CD and normal | paired fib >= 50 | paired endo >= 50 | both |
|---|---|---|---|---|---|
| TI stroma | 40 | 6 | 4 | 3 | 3 |
| colon stroma | 34 | 4 | 2 | 2 | 2 |
| ileum immune | 42 | 6 | n/a | n/a | n/a |

- Pseudobulks only (clock-panel genes); the raw Census h5ads are not on disk (available from the Census S3). Inflamed vs non-inflamed is not annotated at the pseudobulk level, only "Crohn disease" vs "normal". Library structure not checkable from the pseudobulk. Too few donors. Table: `kong_summary.csv`.

## 4. GSE266546 (IBD endoscopy and surgical TI / AC)
- Only the cell annotation files are on disk (no count matrix). 202,013 cells, 170 samples, 85 patients.
- Histology classes per sample: Normal_CD 58, Control 44, Mild 31, Quiescent 17, Severe 11, Moderate 9. No patient has two histology grades in the same tissue (0 patient x tissue); the 29 multi-grade patients differ by tissue (AC vs TI), so there is no same-site contrast.
- Fibroblast or endothelial cells: median 2 per sample; only 7 of 170 samples have fib >= 50 and 8 have endo >= 50. Epithelial-dominated. Not usable.

## 5. IBDverse (E-MTAB-16986, E-MTAB-16999)
- Metadata only: 421 individuals (CD 125, healthy 296) and 343 individuals (CD 111, healthy 232). One TI biopsy per individual in 16999; 16986 has TI, rectum and blood (311 individuals with 2 samples). `ti_ses_inflamed` is an individual-level label. No same-person inflamed / non-inflamed contrast and no cell counts on disk.

## 6. GSE147424 (AD, He 2020; `healthy_skin/GSE147424`)
- 17 samples: 4 lesional, 5 non-lesional AD, 8 healthy. Sample titles (S1_LS, S3_NL...) do not link lesional and non-lesional samples to patients; the pairing needs the paper (BioProject says 5 AD patients). Cell counts per class not computed (matrices not parsed). Most likely at most about 4-5 pairs.

## 7. Bulk series (not eligible)
- `ibd_bulk`: 7 series (GSE107499, GSE179285, GSE193677, GSE208303, GSE243625, GSE66407, GSE95095), `skin_lesional`: about 45 series (psoriasis, AD and other lesional vs non-lesional). Bulk gives one compartment per sample, so no host-referenced reading is possible. Skin: GSE205155 is epidermis vs dermis every 4 h in 11 healthy subjects (microarray), not a disease design.

## 8. Web search (titles only, nothing downloaded, not verified on GEO)
- Smillie et al. 2019 (Cell, PMC6662628): 18 UC patients and 12 healthy, paired inflamed and adjacent non-inflamed biopsies from the same procedure, epithelial and lamina propria fractions; 366,650 cells; accession to be taken from the paper (Single Cell Portal expected). Closest additional candidate to TAURUS; not on disk.
- GSE95446 (Kinchen 2018, colonic mesenchyme in UC): stroma-focused scRNA-seq; whether pairs exist is unverified.
- GSE166924: sorted bulk-style RNA-seq of IBD cell fractions, inflamed vs non-inflamed; not single-cell.
- Search did not return another human AD scRNA study with paired lesional / non-lesional samples beyond GSE147424.

## Recommendation
1. TAURUS is the only candidate with enough pairs (24 paired patients with both eligible families, 28 for endothelium alone) and raw counts on disk. Its caveats: it is already analysed in A45 (not blind), inflamed vs non-inflamed is confounded with treatment timepoint and partly with site, and inflamed biopsies have few fibroblasts (median 137).
2. Psoriasis sc: 12 patients with both families; on disk but already used in A54/A58.
3. Kong: 3 and 2 paired donors, insufficient.
4. Acquisition to consider (ask the user first): Smillie 2019 UC (independent of TAURUS).

Blockers: none on disk-level for items 1-2. No time-of-day information in any candidate. TAURUS licence text and the GSE147424 patient pairing were not verified.
