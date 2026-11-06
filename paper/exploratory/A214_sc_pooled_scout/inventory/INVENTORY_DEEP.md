# A214 deep-cohort scout (2026-11-06): >= 8 clock UMI per cell in malignant and normal epithelium

Table: `inventory_deep.csv` (98 rows). Counts only; no clock statistic. Scout, nothing committed.

## Verdict
No study found that meets >= 8 clock UMI per cell in both compartments with >= 100 + 100 cells per patient. Nothing verified at >= 16.

## Measured
- **Moorman 2025 CRC metastasis (CELLxGENE 829a3cd1; 10x 3' v3; raw layer; 14 clock genes)**: 14 patients have >= 100 malignant and >= 100 colonocyte/secretory cells. Total UMI per cell ~15k (malignant 15.2k, epithelium 14.5k, Census raw_sum) but clock UMI median 2 (malignant mean 3.3, epithelium mean 2.3), about 1.3 per 10k, well below the 3.5-4 per 10k of our data. It fails. Downloaded 0.75 GB to `cellxgene_829a3cd1_moorman2025/` (PROVENANCE.md written).
- **scRNAseq_database**: 5,144 per-sample panel files (12 genes, no TEF/CIART), 4,556 SRX/ERX plus 588 GSM. The values in `scbasecount_human_studies.csv` match `obs` (n cells, median UMI: exact). Median clock density over all samples is 3.2 per 10k (panel, intron-inclusive GeneFull_Ex50pAS counts, so not directly comparable with our exonic counts). Solid-tumour samples >= 15k UMI, excluding cell lines and blood cancers: 70 samples in 39 studies (all 10x 3'), sample-level median clock UMI 1-20 (all cells, malignant and stroma mixed). Files carry no cell labels, so compartment-level clock UMI cannot be computed from them; full h5ads (GCS) would be needed. Higher-clock samples are mostly non-epithelial (glioblastoma, rhabdomyosarcoma, hepatoblastoma SRP496487 median 19) or cell lines/organoids (SRP417718 PDAC cells 20). Candidates with possible epithelial pairs: SRP376674 (CRC, adjacent normal, 3 samples, median 8), SRP329496 (ovarian, on disk), SRP573418 (Barrett/EAC, 1 sample, 8), SRP396819 (PDAC adjacent, 1 sample, 3), SRP414514 (prostate, 10). Too few samples per study to reach 100 patients.
- Database bugs/caveats: the ERX/NRX/GSM keying is not uniform (96 `NRX*` files, 4 SRX files missing from `map_srx_study.csv`; study IDs for 3,596 samples were absent from the local maps, resolved through ENA for the 70 tumour samples). GSM files include Smart-seq2/ SS-like plates (median "UMI" 100k-360k are reads, e.g. GSE300150 ependymoma, GSE239854, GSE131309) and the GSM tumour set is dominated by CNS tumours, organoids and cell lines.

## Census (obs only, malignant cell labels, primary data)
Median raw_sum per malignant cell by assay: Smart-seq2 (reads, 0.2-1.4M), 10x 3' v3 glioblastoma 28k/22k (1 donor, no epithelium), BD Rhapsody LuCA 18.7k malignant vs 8.9k epithelium (9 patients with >= 100 each), BD CRLM 16.8k vs 7.5k (hepatocyte/cholangiocyte, 2 patients), 10x v3 Moorman 14.8k vs 14.6k (measured above), others <= 7.5k (TruDrop polyps 7.5k, breast 3.5-7k). Smart-seq2 LuCA: 3,187 malignant, 1,251 epithelial, 1 patient with both >= 100 (reads, not UMI). Depth in epithelium is the limit, not malignant.

## Alternative chemistries (user confirmed plate is eligible; non-UMI = reads, Poisson corrections need UMI)
- **Plate**: GSE81861 (CRC, 11 patients, C1 SMARTer) and GSE103322 (HNSCC, 18 patients, Smart-seq2, matrix log2 TPM) are the known paired designs; both have too few cells per patient for 100 + 100 and are read-based. GEO search for Smart-seq3, FLASH-seq, VASA-seq and Parse found no tumour with malignant + normal epithelium.
- **Flex**: GSE338505 (FFPE colon, 1.7 GB), GSE318840 (osteosarcoma), 10x 900k NSCLC tumour/normal-adjacent Flex aggregate; probe-set membership of the 14 clock genes could not be checked (cf.10xgenomics.com returns 403, www.10xgenomics.com 429). Depth reported in Census for Flex is 1.5-2.5k UMI: worse.
- **Imaging**: CosMx GSE348495 (CRC, 3 adjacent), GSE346380 (LUAD vs normal, unpaired), GSE312757 (CRLM), Xenium GSE326226 (MIBC). Panel content (Xenium Prime 5K, CosMx 6K) not verified; spill-over between adjacent malignant and normal cells affects clock genes strongly. Not measured.

## Shortlist
Empty at >= 8. Nearest: none measured; hypotheses only: (1) intron-inclusive reprocessing of Moorman (would lift counts, as scBaseCount GeneFull vs exonic suggests; untested), (2) BD Rhapsody LuCA (9 patients, malignant 18.7k UMI but epithelium 8.9k, clock not measured, 17.6 GB), (3) full scBaseCount h5ad of SRP376674 (3 samples).

## Not verified
10x pages and probe lists, 3CA v2, Xenium/CosMx gene lists, per-patient clock UMI for any dataset except Moorman, Flex/BD/imaging counts, Smart-seq3 tumour studies (none found, absence not proven).
