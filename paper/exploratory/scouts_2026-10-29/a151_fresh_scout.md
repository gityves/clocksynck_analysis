# A151 fresh scout: paired tumour / adjacent-normal single-cell cohorts not used by A149/A151

Date 2026-10-29. Metadata only: no clock statistic, no expression read, no count matrix downloaded.
[V] verified from a GEO sample-title list, GEO series text or the CELLxGENE index; [U] unverified.
Downloaded (about 20 MB, in `$SCCC_DATA/scout_a151/`, parsed with `python -I` or the sccc env for one xlsx): GEO sample lists
(`fam_*.txt`), GEO esearch/esummary JSON, CELLxGENE collection and dataset index JSON, the GSE268238 metadata xlsx.

## Excluded as already examined
- The 61 units of `A149/output/units.csv` (= `A151/output/units.csv`), including the GSE164789 / GSE229253 lung and GSE150290,
  GSE167297, GSE183904, GSE270680 gastric cohorts. GSE229706 is the NYU superseries of GSE229253 (same nuclei).
- GSE206785, GSE290298 (log-normalised only) and GSE333596 (libraries not separable) were dropped by A149 note 3 rules, so they
  are not fresh. They are on disk under `geo/`. GSE178341 is Pelka (CELLxGENE 4d82bd8e), already used.
- Moorman 2025 and Joanito 2022 CRC: see the table (epithelium only / Synapse login).

## Candidates

| accession | cancer | paired patients | patient mapping | annotation / cells | format, size | Census | on disk | verdict |
|---|---|---|---|---|---|---|---|---|
| GSE268238 (+ GSE201347, GSE231539, GSE254188; Kyushu) | gastric | **16 untreated pairs** [V] (G40-43, 49-52, 55 in GSE201347; G47 T in 201347, N in 268238; G66, 67, 69, 71, 74, 75 T in GSE254188, N in 268238). G54 pair has chemo T (excluded) | GSM titles `gastric cancer [GxxT]` / `non-malignant gastric mucosa [GxxN]` (the Gnn number is the patient) [V] | none in GEO text; Seurat object has author annotation [U]; cells per patient unknown [U] | `GSE268238_GC.PNI.rds.gz` (about 22 GB), `GSE201347_GC.integrated.no.NAC.rds.gz` (18 GB); counts slot unverified [U]; BD Rhapsody (RSEC MolsPerCell) and 10x mix [U]; raw FASTQ in SRA PRJNA1115225 (not used) | no [U] | no | **register** (largest new gastric; heavy parse of an RDS with the pure-python reader used for TNBC; stromal content needs checking [U]) |
| CELLxGENE f7cecffa (Li 2022 Cancer Cell, ccRCC; EGAD00001008030 for raw) | kidney (ccRCC) | up to 12 donors with tumour core / interface / normal kidney [V donors = 12]; pairs not counted [U] | `donor_id` in the h5ad obs | 270,855 cells (dataset "Single-cell transcriptomic datasets of Renal cell carcinoma"); fibroblast, endothelial and lymphatic endothelial types present [V]; disease labels include `normal` | h5ad (CELLxGENE, open download, size not listed [U]); 10x 5' | yes | no | **register** (stroma-rich, open h5ad; includes oncocytoma / benign, filter to ccRCC) |
| CELLxGENE 3c34e6f1 + Zenodo 14284357 (OPSCC) | oropharyngeal SCC | at most 7 [U] (10 donors, `normal` and tumour) | `donor_id` | 82,844 cells; fibroblast, endothelial, lymphatic endothelial [V] | h5ad (open), 10x 3' v3 | yes | no | pool with others (< 10 pairs) |
| GSE318418 | HCC | 4 of 7 patients [V] | titles `patient A-xx, tumor tissue / tissue adjacent to tumor` | live-sorted and CD45-sorted libraries; stroma not guaranteed | RAW (size not checked) | no | no | too small, sorted |
| GSE196964 | colorectal | 2 pairs (Colon0008/0040 vs ColorectalTumor0008/0040) [V]; 9 more tumours only | GSM titles | unknown | RAW.tar; raw also deposited in dbGaP | no | no | too small |
| GSE201348 (HTAN VUMC snRNA, 72 samples) | colorectal | 2 patients have CRC + unaffected [V] | title prefix A001 etc. | polyp-centred | RAW.tar | HTAN a48f5033 | no | no |
| GSE221561 | oesophageal SCC | 2 pairs (S6417, S9478), rest treated [V] | titles | n/a | n/a | no | no | no |
| GSE173950 | oesophageal adenocarcinoma | about 10-12 N/T pairs by title [U] | titles `OAC44N / OAC44T2` | Drop-seq, split Naive/Chemo/CRT files | 3 DGE files | no | no | register only after the naive-only subset is counted [U] |
| CELLxGENE 1cbfb478 (Moorman 2025) | colorectal | 31 trios [U] | donor_id | CELLxGENE sets are epithelial subsets only (`..._epithelial`, `_non-tumor_epith`, `_untreated_epith`); no fibroblast type [V] | h5ad | yes | no | fails stroma criterion |
| Joanito 2022 (Synapse syn26844071) | colorectal | 63 patients, about 39 normal; SG1/SG2 + KUL5 new [U] | n/a | n/a | Synapse login | no | partly (GSE132465, GSE144735) | login, excluded by rule |
| CELLxGENE a18474f4 (Busslinger / Gastric Cancer 2023) | gastric | 70 donors in the cancer integration; normal and cancer may be from different donors [U] | donor_id | 5,754 fibroblasts, endothelial subsets present [V] | h5ad | yes | no | check pairing from obs before use |
| GSE164522 | CRC + liver metastasis | CRLM01-17, tumour/mets design [U] | titles | n/a | n/a | no | no | not a normal pairing |
| GSE156626, GSE156625 | HCC | 3 patients, Nanostring GeoMx (156626) | n/a | not single cell | n/a | no | no | no |
| OMIX006159 (CNCB-NGDC) | breast | 7 patients, tumour + adjacent [per search snippet, U] | n/a | fibroblast paper | NGDC, access by application | no | no | controlled, listed only |
| 3CA, TISCH2 | many | not checked in depth (3CA website reset the connection; TISCH2 holds mostly tumour-only GEO series) [U] | | | | | | |

## Findings
- GEO esearch (single-cell AND tumour AND adjacent/paired normal, 2021-2026, 138 + about 375 series) and the CELLxGENE index
  (402 collections, 2,244 datasets) give very little beyond what A149/A151 already took. Most remaining series are tumour-only,
  immune-sorted, treated, or have 2-4 pairs.
- Verified fresh pairs: gastric 16 (Kyushu). Plausible but unverified: ccRCC up to 12, OPSCC at most 7, oesophageal adenocarcinoma
  about 10. Total reachable about 16 verified, about 45 with the unverified ones.
- Colorectal, liver and breast are not extended: CRC fresh sets are epithelial-only (Moorman), tiny (GSE196964, GSE201348) or
  login-gated (Joanito); no new HCC set has more than 4 paired patients with unsorted cells; breast has no open paired
  adjacent-normal set with stroma (the large breast atlases CELLxGENE 9432ae97 and dea97145 are tumour-only; OMIX006159 is gated).
- Within the already-used set, the limit is annotation, not accession: A151 units.csv shows 0 fibroblast pairs for many large
  cohorts (for example GSE150290 with 23 pairs, GSE242299 with 8, GSE178481 with 9, GSE144240 with 10, GSE199654 with 13 and
  GSE160269). A re-annotation of those (a marker-based fibroblast call) would add more patients than any fresh cohort, but
  they are not fresh by the stated rule. Raise this with the PI.

## Recommendation (3 cohorts; a fourth conditional)
1. Kyushu gastric (GSE268238 with GSE201347 / GSE254188): 16 untreated pairs. Gastric was 17 patients in 2 unflagged cohorts
   in the bounded A149 pairing (MDE 0.785), so this at least doubles it. Check the object for fibroblast / endothelial counts per
   patient before registration.
2. Li 2022 ccRCC (CELLxGENE f7cecffa): up to 12 pairs, stroma-rich, open h5ad. Kidney was not estimable. Count pairs from the obs
   (`donor_id`, tissue, disease) first.
3. OPSCC (CELLxGENE 3c34e6f1): at most 7 pairs, head and neck, adds a new cancer type; pool, do not run alone.
4. Conditional: GSE173950 (oesophageal adenocarcinoma) once a naive-only N/T subset of at least 8 is confirmed.
Not recommended: Moorman, GSE318418, GSE196964, GSE201348, Joanito (login), OMIX006159 (access request).
Registration must state per-patient minimum cell counts for fibroblast and endothelial compartments in both tissues, because
the bounded statistic flags units with 4-11 patients.
