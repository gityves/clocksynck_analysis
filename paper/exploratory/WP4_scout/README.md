# WP4 scout (A113, A123, A130, A131)

Scout only: no hypothesis statistic was computed (no phase against time, no offsets against truth). Table: `output/scout_wp4.csv`
(git-ignored), per-tissue GTEx counts `output/gtex_window_counts.csv`. Script: `scout_gtex_counts.py` (reads `obs` only, aggregates only).

## Datasets (A113)
| accession | design | genes | age, sex | matrix | builder |
|---|---|---|---|---|---|
| GSE108539 (Perrin 2018; same biopsies as Loizides-Mangold 2017) | 10 volunteers (9 M/1 F), vastus lateralis, 6 biopsies each every 4 h from 12:00 (57 libraries; one subject has 3); hourly isocaloric meals, lights 07-22 h | 13 of 14 (ARNTL absent); log2 edgeR RPKM exon and intron, not counts | cohort aggregate only (22-54 y); GEO subject ids not linkable to Table S1 | public | `common/datasets/muscle_timed_gse108539/build.py` |
| GSE104674 (Stenvers 2019) | 12 men (6 lean, 6 obese T2D), adipose, 4 biopsies each, ZT 15:30 (d2), 0:15, 5:45, 11:15 (d3) | 14 of 14 (CIART = C1orf51); raw counts | per subject, age 40-71 | public | `common/datasets/adipose_timed_gse104674/build.py` |
| Smith 2024 JCEM | same cohort summary as Perrin (9 M/1 F, 30 +- 10 y, BMI 24.1 +- 2.7); TaqMan array, 45 genes | partial | aggregate | not public (on request) | none |

Data: `$SCCC_DATA/muscle_timed/GSE108539/`, `$SCCC_DATA/adipose_timed/GSE104674/` (PROVENANCE.md with sha256). Project total 309 GB of the 400 GB cap.

## A130
V25.5 (GSE233242): 42 primary patients with age (33-93 y, SI Table S1, all female). Pooled with GSE104674: 54 subjects with per-subject age, range 33-93 (60 y); the rule (>= 30, >= 25 y) holds. Muscle enters only as aggregate.

## A123
GTEx v10 fast deaths with time of death: 244 donors, 111 died 08-18 h. Tissues with >= 20 such donors: 40 of 48 (benchmark 8 tissues: 79-98). CHIRAL commit 1f78a9a: only `sample_metadata.txt` (per-sample CT) is local; the per-tissue tables sit on a Drive link.

## A131
Public GTEx dictionaries (dbGaP pht002743, v8 and v10 SampleAttributesDD) give SMCENTER only as "Code for BSS collection site" (B1, C1, D1, with A1 variants). No locations. Two usable centres in fast deaths (B1, C1).
