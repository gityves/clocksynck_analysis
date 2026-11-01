# CPTAC scout for C1 replication (metadata only, 2026-11-01)

Source: GDC REST API (`/files`, `/cases`; program CPTAC; JSON metadata only, no expression file downloaded).
Selection: data_type "Gene Expression Quantification", workflow "STAR - Counts" (RNA-Seq), sample type Primary Tumor or Solid Tissue Normal, pair = same case with both.

## Findings
- CPTAC on GDC = CPTAC-2 (342 cases) and CPTAC-3 (1866 cases). 3308 STAR gene-count files, all **open access** (the controlled files are splice-junction files, not needed). Data type: STAR raw/augmented gene counts (`*.rna_seq.augmented_star_gene_counts.tsv`, ~4.2 MB each), same pipeline as the TCGA STAR counts already used.
- **CPTAC-2 (breast 133, colon 105, ovary 71, ...) has tumour RNA-seq only, zero normal-adjacent files**: no pairs. CPTAC-3 has no breast or colon tumour RNA-seq on GDC.
- **ER status for BRCA: not usable.** No BRCA pairs, and GDC clinical for CPTAC carries no ER/PR/HER2 field (ER only in PDC / paper supplements, for the CPTAC-2 tumours without NAT). H175.2 (ER) cannot be replicated from CPTAC.
- **Overlap with TCGA:** CPTAC-3 cases are C3L-/C3N- IDs, independent patients (no TCGA IDs in the pair sets; all pair cases start with C3). CPTAC-2 breast/colon/ovary were TCGA-analysed tumours, but irrelevant here (no NAT). Independent replication is therefore valid.
- Cancer-type assignment is from GDC primary site and `primary_diagnosis` (LUAD = adenocarcinoma, LSCC = squamous, CCRCC = "renal cell carcinoma, nos"); approximate, to be refined from the CPTAC clinical tables. Some cases have several tumour files (multiple aliquots); pairs counted per case.

## Pairs, size, rough MDE
Rough MDE = 0.194 x sqrt(517 / n_pairs) (A175 H175.1 pooled MDE 0.194 at 517 pairs; `A175_c1_pairing_strata/output/power.csv`). This is only a scaling, ignoring per-cancer arm balance and noise differences; the D needed (0.158) requires n >= ~780 pooled pairs.

| type | tumour cases | NAT cases | pairs | files (T+N) | size | rough MDE (alone) | TCGA counterpart |
|---|---|---|---|---|---|---|---|
| LUAD | 228 | 215 | 215 | 430 | 1.8 GB | 0.30 | LUAD (in A175 pool) |
| LSCC | 109 | 96 | 96 | 193 | 0.8 GB | 0.45 | LUSC (not in pool) |
| CCRCC | 223 | 151 | 151 | 383 | 1.6 GB | 0.36 | KIRC (in pool) |
| UCEC | 240 | 45 | 45 | 134 | 0.6 GB | 0.66 | UCEC (in pool) |
| HNSCC (oral/larynx/pharynx sites) | 110 | 61 | 61 | 162 | 0.7 GB | 0.57 | HNSC |
| PDAC | 161 | 59 | 52 | 128 | 0.5 GB | 0.61 | PAAD |
| SKCM | 162 | 66 | 66 | 132 | 0.6 GB | 0.54 | SKCM |
| STAD | 151 | 27 | 27 | 54 | 0.2 GB | 0.85 | STAD |
| GBM | 200 | 17 | 1 | - | - | - | - |
| BRCA, COAD, OV | CPTAC-2 only | 0 | 0 | - | - | - | - |
| **Pool of the 8 types** | | | **713** | 1616 | **6.7 GB** | **0.165** | |
| LUAD+LSCC+CCRCC | | | 462 | 1006 | 4.3 GB | 0.205 | |

(Another 17 pairs in other kidney histologies, excluded.)

## Recommendation
- **Not feasible as a clearly powered independent C1 replication for H175.1 in the strict sense**: the pooled 713 pairs give a rough MDE ~0.165, just above the target 0.158 (effect D +0.158), and the pool mixes 8 types with very unequal n (LUAD/CCRCC/LSCC carry 65%). Any single type is far underpowered.
- **H175.2 (ER) not replicable** (no BRCA pairs).
- Reasonable use: a pre-registered pooled confirmation at MDE ~0.165 (a directional/replication-of-sign test of D > 0, expected power ~45-50% at the true effect 0.158 if equal to the TCGA point estimate, rough), or a combined TCGA + CPTAC meta-analysis (1230 pairs, rough MDE ~0.126) with CPTAC as an independent-cohort consistency check. Downloads are cheap (6.7 GB), so the option is low cost; the decision is the PI's. Check the A175 pool definition (it used 10 cancers, 517 pairs) before registering, and that CPTAC proliferation arms can be defined the same way.

## Download plan (user runs; nothing downloaded here)
- Disk: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian` = 369 GB used (cap 400 GB); the full pool needs ~6.7 GB, leaving ~24 GB headroom. Volume /mnt/max8tb has 724 GB free.
- Manifests (GDC format, id/filename/md5/size/state) and sample sheets (file_id, case, sample, sample_type) in `cptac_manifests/<TYPE>_manifest.tsv` and `<TYPE>_samples.tsv` for LUAD LSCC CCRCC UCEC HNSCC PDAC SKCM STAD.
- Run (md5-verified, resumable, plain curl on the open GDC data endpoint, no token):
  `paper/exploratory/scouts_2026-11-01/cptac_download.sh LUAD LSCC CCRCC UCEC HNSCC PDAC SKCM STAD`
  into `$SCCC_DATA/cptac_<type>/` (the manifest and sample sheet are copied there too).
- Alternative: `gdc-client download -m cptac_manifests/LUAD_manifest.tsv -d $SCCC_DATA/cptac_LUAD/`.
- Still needed before analysis: CPTAC clinical tables (stage, histology correction) from PDC/GDC clinical for proliferation-independent covariates; gene-ID mapping identical to TCGA STAR (same pipeline).
