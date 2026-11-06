# A226 scout (2026-11-06): new paired tumour / adjacent-normal single-cell cohorts with stroma

Metadata only until the download step; no clock, pseudo-clock or pairing statistic computed. Exclusion lists: A149 / A151 `output/units.csv`.
Table of candidates and decisions: `scout_candidates.csv`; raw screens: `db_screen.csv`, `census_*.csv`, `geo_series.csv`, `geo_filtered.csv`.

## Sources searched
1. **User database** (`scRNAseq_database/output/datasets`, 691 GEO studies without SRA raw reads; `scBaseCount` table) and the curated cancer inventory
   (`$SCCC_DATA/inventory/inventory_samples.csv`, 903 studies): every study with tumour and adjacent-normal samples of the same patients is already a
   `cancer_sc` / A149 cohort (units.csv); the rest have <= 3-5 pairs (SRP431612 = GSE229253, SRP367723 CD45+ only, SRP307611 3 pairs).
2. **On disk** (`$SCCC_DATA`): `cellxgene_829a3cd1_moorman2025` is epithelium only (no fibroblast / endothelial cells); LuCA parts and `cells/` have no matched normal;
   A214 inventories list nothing new.
3. **CELLxGENE Census 2025-11-08** (`scripts/scout_census.py`): primary cells of 203 datasets with non-normal disease; donors with fibroblasts / endothelium (>= 20 cells) in both a
   tumour and a normal sample. Datasets with >= 3 such donors: Pelka (GSE178341, A149), the Thomas gut atlas (19053a82; donors are A149 cohorts), OPSCC (b6b5ea88, kept),
   RCC Li 2022 (5af90777, kept), CRLM atlas (e3ed2ba4; 83 normal fibroblasts in total), the multi-tissue TME atlas (b617ee1b; its Qian 2020 donors are not in A149, kept).
4. **GEO E-utilities** (`scripts/scout_geo.py`; 6,226 human single-cell/nucleus series 2018-2026 matching tumour + adjacent/normal terms, 198 with >= 3 tumour-like and >= 3 normal-like
   sample titles, 138 after dropping A149/A151 cohorts and non-single-cell, triaged by sample titles / summary): kept GSE308103 (LUAD snRNA, 23 patients with normal + LUAD,
   treatment-naive), GSE183219 (NSCLC sorted stromal cells, 12 pairs), GSE155109 (breast EC-enriched + stroma fractions, 9 patients). Others rejected (see csv): treated normals
   (GSE205506), < 4 pairs, bulk / spatial, T-cell sorted, Seurat-only objects, no sample mapping (GSE244267 downloaded, not used).

## Outcome
No large new paired single-cell cohort exists outside A149 / A151: the public paired cohorts with stroma were used there. After building the seven new cohorts with the
A149 rules (`output/eligibility.csv`), patients with both compartments are: fibroblasts 11 (GSE308103 5, Qian 2020 CRC 6), endothelium 13 (GSE308103 5, GSE183219 8).
The loss is mostly the A149 sample depth rule (median UMI >= 1,500; GSE308103 is FFPE snRNA at about 1,400-1,600 UMI per nucleus: 7 of 23 patients pass) and the
compartment presence rule (>= 60 % of patients or >= 6). Without both rules (sensitivity upper bound): fibroblasts 48 patients / 6 cohorts, endothelium 47 / 5.
Downloads: 9.9 GB in `$SCCC_DATA/a226/` (PROVENANCE.md per source, sha256, parsed only).
