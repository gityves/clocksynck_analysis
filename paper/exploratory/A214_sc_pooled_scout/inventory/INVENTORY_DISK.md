# A214 disk inventory: malignant vs normal epithelium, same patient (counts only)

Script `../scripts/inventory_disk.py`; tables `inventory_studies.csv` (61 rows), `inventory_patients.csv`. Eligible patient = >= 100 malignant (copy-number call) and >= 100 normal epithelial cells. (a) `adjacent_normal:epithelial`; (b) `tumour:epithelial_nonmalignant`. No clock statistic computed.

## Totals
| set | studies | patients (all) | eligible (a) | eligible (b) | eligible either |
|---|---|---|---|---|---|
| per-cell tables, 29 with calls + barcodes | 29 (of 38 tables) | 421 | 108 (24 studies) | 115 (20 studies) | 168 |
| A149 new cohorts, pseudobulk only (17) | 17 | 187 | 48 (33 with copy-number-called malignant; 15 with unsplit `tumour:epithelial`) | not available | 33 (called only) |
| Wu GSE176078 (author labels) | 1 | 26 | 0 (tumour samples only) | 4 | 4 |

Cells over eligible patients, per-cell tables: (a) 92,063 malignant / 103,169 normal; (b) 116,654 / 57,433.
Most patients: GSE178341 (a 31, b 26), SRP327758 (a 1, b 25), GSE183904 (a 5, b 13), GSE166555 (a 10, b 8), GSE131907 (a 7, b 1).

## Missing or ambiguous
- 9 per-cell tables cannot be split (no copy-number calls and/or no barcodes in the legacy file; 132 patients): ERP130387, ERP156718, SRP301444, SRP329496, SRP335254, SRP407590, SRP467578, SRP539998, SRP567887. Their (a) and (b) are NA, not 0.
- No per-cell table on disk (cancer_sc pseudobulk only): SRP309698, SRP353200, SRP414286, SRP429995, SRP514981.
- Cohorts with a call file but 0 eligible patients are real zeros (e.g. GSE181919, ERP153272, SRP273482, GSE236696, GSE271690).
- A149 new cohorts have no per-cell tables (only `output/pb/*.parquet`); the 12 `cnv_a149` call files hold calls only, no counts. Clock UMI there is the group mean per patient, not a per-cell median; (b) does not exist in these pseudobulks. Cohorts with unsplit tumour epithelium (GSE167297, GSE176031, GSE193337, GSE200972, GSE310916; GSE160269 and others have no adjacent normal): the "malignant" group is all tumour epithelium.
- Assay: all 43 cancer_sc cohorts are droplet 10x (3' or 5'; SRP327758 includes "other" chemistry). A149 non-10x droplet/microwell: GSE242299 (inDrops), GSE293189 and GSE176031 (Seq-Well), GSE282701 (SeekOne), GSE242889 (CeleScope). No plate-based cohort on disk among these.
- Pal GSE161529 (breast, CoPhaser fit valid, 21 tumours, 43,693 malignant cells by author `Tum` clusters): normals are reduction mammoplasty, not paired; no per-patient counts made, eligible 0 by design. Wu clock UMI not read (author labels, nCount_RNA only).
- CoPhaser fit validity is hand-copied from `C4_single_cell/fits/FITS.md` (y, y_low, suspect, weak, n). Cycling annotation possible = valid fit (y, y_low) or median genome UMI per malignant cell > 1,500 (true for nearly all with a split).
- Clock UMI per cell is mostly 0-5 (median 0-2 in many cohorts), see `*_clock_umi_median_*`.
