# Scout 2026-11-12: multi-region psychiatric brain (topic A) and transplant grafts (topic B)

Scope: sample sheets and detection summaries only. No clock statistic was computed. Numbers below come from
`output/psych_*_samples.csv`, `output/psych_design_tables.csv`, `output/transplant_GSE43974_*.csv` and
`output/clock_gene_detection.csv`, built by `scripts/scout_build.py` (helper: `scripts/geo_esearch.py`).
Data: `$SCCC_DATA/psych_multiregion/<GSE>/` and `$SCCC_DATA/transplant/GSE43974/`, each with a PROVENANCE.md (sha256).
Total new data on disk: 126 MB (brain) + 120 MB (transplant).

## A. Multi-region post-mortem brain

Ranked by usefulness for the host-referenced design (same donor, >= 3 regions, cases and controls, RNA-seq, donor link).

| rank | accession | design | donors with >= 2 regions (by group) | regions | donor id link | time of death | data | downloaded |
|---|---|---|---|---|---|---|---|---|
| 1 | [GSE80655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) Ramaker 2017 | 96 donors, 24 per group (control, MDD, bipolar, schizophrenia) | all 96 have >= 2 regions; with 3 regions: control 22, MDD 21, bipolar 23, schizophrenia 23 (89 donors); 2 regions: 7 | AnCg 96, DLPFC 94, nAcc 91 samples | yes, title prefix X#### (age, sex, diagnosis, PMI identical across a donor's samples) | no | RNA-seq counts, 56,755 genes | yes |
| 2 | [GSE202537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202537) Ketchesin | 72 donors: matched control 36, schizophrenia 28, bipolar psychosis 8 | 71 donors with all 3 regions (control 36, SZ 28, BP 7), 1 donor with 2 | NAc 71, caudate 72, putamen 72 | yes, `patients_id` field | yes (`corrected.tod`, 215 of 215 samples, range -5.6 to 17.1 h) | RNA-seq counts, 59,134 genes | yes |
| 3 | [GSE102556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102556) Labonte 2017 | 48 donors: MDD 26, control 22 (no schizophrenia, no bipolar) | all 48 have >= 4 regions; 6 regions: control 14, MDD 13; 5: control 6, MDD 11; 4: control 2, MDD 2 | OFC 48, dlPFC 48, aINS 48, NAc 48, Sub 43, Cg25 28 (GEO sheet) | yes, title prefix (donor number, attributes identical across regions) | no | RNA-seq FPKM table (log-scale values), 56,623 genes | yes |
| 4 | [GSE53987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987) Lanz | 75 inferred donors: control 19, bipolar 19, MDD 18, schizophrenia 19 | with 3 regions: control 17, bipolar 15, MDD 14, SZ 12 (58); 2 regions: 14 donors; 1 region: 3 | PFC (BA46) 68, hippocampus 68, associative striatum 69 | inferred only: title index within group; age, sex, PMI, pH identical across the 3 regions | no | Affymetrix U133 Plus 2 (series matrix), 24,442 symbols | yes (series matrix, arrays) |
| 5 | GSE54567/68/70/71/72, GSE54564/65/66 Chang 2014 (Sibille) | matched control/MDD pairs, one GSE per region x sex x cohort | see note 5 | BA9, NY_BA9, BA25 (sgACC), MD1_ACC, amygdala | pair number only, shared across GSEs (inferred, no age/sex to verify) | no | U133 arrays | no |
| 6 | [GSE160521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160521) Ketchesin 2021 | 59 donors x NAc, caudate, putamen (177 samples) | 59 with 3 regions | 3 | yes (`hu`) | yes | RNA-seq read counts | no (donor ids overlap GSE202537, no disease field in sheet) |
| 7 | [GSE92538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92538) Pritzker/Akil-Bunney (Kim/Vawter, cell-type deconvolution paper) | 172 subjects, 363 arrays (replicates and two array platforms), MDD, bipolar, schizophrenia, control | none: DLPFC only | DLPFC | yes (`subject id`) | no | Illumina arrays | no |
| 8 | [GSE35978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35978) | 312 samples, cerebellum and parietal cortex, schizophrenia, bipolar, depression, control | 0 donors with 2 regions (291 distinct sample prefixes, each one tissue) | 2 | no | no | Affymetrix arrays | no |
| 9 | GSE12649 (102 BA46), GSE5388 (61 DLPFC) | single region | 0 | 1 | n/a | no | U133A arrays | no |

Notes
1. GSE80655 is the best candidate overall: balanced 24 donors per group, three regions, RNA-seq counts, age/sex/PMI/pH present. No
   RIN and no time of death. Cingulate, DLPFC and nucleus accumbens are not the classic clock-hub regions but the clock genes are well detected
   (median within-sample percentile 81.6 to 98.1 over the 14 genes, all 14 present, non-zero in 100% of samples).
2. GSE202537 is the only dataset that records time of death, which gives the timed (CHIRAL-style) benchmark for the same donors; cases are psychosis
   (schizophrenia 28, bipolar 8, as in Ketchesin 2022) not MDD. Matched controls 36. Median percentiles 81.8 to 97.8, all genes present.
3. GSE102556: 203 of 263 samples are suicide deaths (cause of death field), 89 on medication, including antidepressants (55); smoking and alcohol
   fields are mostly NA. The FPKM table has 281 columns (BA25 46 columns vs 28 GEO samples): 18 BA25 libraries have no GSM entry in the sheet, so
   their phenotype can only be inferred from the donor's other regions (`output/psych_GSE102556_matrix_columns.csv`). Median percentiles 77.9 to 98.1.
   Only MDD vs control; the table in this file counts GEO-sheet samples only.
4. GSE53987: no donor id in GEO; matching of regions rests on identical age, sex, PMI and pH across the three titled `<group>_<region>_<index>` samples
   (no conflicts found in any index). Probe-to-symbol collapse by highest-mean probe (GPL570 annotation from GEO); all 14 genes present,
   median percentiles 59.4 (CIART) to 88.2; values are normalised log intensities so "non-zero fraction" is not reported.
5. Chang 2014: each GSE is one region/sex/cohort of matched pairs, pair numbers appear to be a global list. Pair numbers of the BA25 series
   (GSE54571 + GSE54572, 25 pairs) are all contained in the BA9 series (GSE54567 + GSE54568, 29 pairs): 25 pairs (25 MDD + 25 control) with DLPFC and
   sgACC; amygdala MD1 (GSE54566, 14 pairs) lies fully within BA9 pair numbers (14 of 14), MD3 amygdala (GSE54564, 21 pairs) 15 of 21; MD1_ACC (GSE54565,
   16 pairs) 16 of 16. This is inference from pair numbering only: GEO gives no age, sex or PMI per sample, so donor identity across regions is
   unverified (a data-owner request or the paper's supplement would settle it). MDD only (no bipolar or schizophrenia), 2 to 3 regions, arrays.
   Not downloaded: fewer donors than the top three and unverifiable linking.
6. Searches by E-utilities for other multi-region MDD deposits (Pritzker, Akil-Bunney) returned only GSE92538 (single region) and GSE45642 (already in
   repo, `output/GSE45642_*`). Li 2013 stays not usable (as stated in the brief). GSE102556 plus GSE80655 are the realistic MDD multi-region sets.

Recommended order of use: GSE80655 (all four diagnoses, power), GSE202537 (timed anchor, psychosis), GSE102556 (MDD, 4 to 6 regions, includes
Cg25 and subiculum), GSE53987 (arrays, cross-check).

Clock-gene detection (all four downloaded sets): 14 of 14 genes present in each. Per-gene numbers in `output/clock_gene_detection.csv`
(columns: present, matrix_name, frac_samples_nonzero, median_pctile). CIART is carried under its old symbol C1orf51 in GPL10558 (kidney set only).

## B. Organ transplants

Question: does a graft carry the recipient's or the donor's clock time?

### B(i) single-cell/nucleus graft data with donor vs recipient assignment: no dataset meets the threshold

Required: >= 8 patients each with >= 50 recipient and >= 200 donor cells. None found on GEO.

| accession | content | patients | genotype donor/recipient split | verdict |
|---|---|---|---|---|
| [GSE145927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145927) Malone 2020 | kidney transplant needle biopsies, 10x, plus WES of donor and recipient PBMC | 7 biopsies (7 patients at most; time points 5 d to 7 y) | yes (WES-based genotyping, MGI metadata) | below 8 patients; not downloaded |
| [GSE109564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109564) | one kidney allograft biopsy, inDrop | 1 | yes in paper | no |
| [GSE220797](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220797) | lung biopsies after cold preservation and 2 h reperfusion, scRNA | 6 patients x 2 | sex mismatch possible (donor sex all M; recipient 4 F, 8 samples M) | 6 patients; not downloaded |
| [GSE274199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274199) | paediatric/young-adult lung allograft transbronchial biopsies, 10x | 8 patients, 69 GSMs (incl. TCR/ATAC) | not annotated | recipient-only origin unknown; not downloaded |
| [GSE189536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189536) | kidney graft biopsies, CITE-seq | 11 biopsies | not in GEO sheet | not downloaded |
| [GSE195718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195718), [GSE228300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228300), [GSE151671](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151671) | kidney allograft snRNA/scRNA | 9, 3, 3 | no | no |
| [GSE317012](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317012) | kidney graft biopsies at transplant, 3 and 12 months (BK virus), scRNA | 26 libraries, patients not enumerated in sheet | not annotated | origin split would need own genotype calling from raw reads; not downloaded |
| [GSE203548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203548) | heart transplant endomyocardial biopsy / explant snRNA (CAV vs control) | 10 samples | no | no |
| [GSE334923](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334923) | chronic AMR: PBMC (41) and matched kidney biopsy (7), scRNA | 29 recipients, 7 with graft biopsy + blood | n/a (recipient blood vs graft) | design (ii) in single cell form, 7 pairs < 20; not downloaded |
| GSE242909 | kidney graft TRM cells, 4 patients, donor and recipient origin | 4 | yes per paper | below 8 |
| Byrne 2020 (lung BAL), Snyder/Zuber (intestine) | not found as processed GEO series by title search | | | not on GEO |

Cell counts per patient could not be checked for any set because none reaches 8 patients; per-cell donor/recipient labels are mostly in
the papers' supplements, not in GEO. Nothing downloaded for (i).

### B(ii) graft biopsy paired with native tissue or blood of the same recipient (>= 20 patients)

No bulk dataset found. GSE334923 has 7 matched PBMC and biopsy pairs (single cell). GSE145927 has PBMC for WES only.

### B(iii) paired pre-implantation and post-reperfusion biopsies, cold ischaemia time

| rank | accession | graft | paired grafts | cold ischaemia time in GEO sheet | donor type | clock times | data | downloaded |
|---|---|---|---|---|---|---|---|---|
| 1 | [GSE126805](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126805) | kidney protocol biopsies, pre, post, 3 months, 1 year | 42 patients; pre and post both present in 41; all four time points in 37 | no (only `patient`, `time point`) | no | none | RNA-seq counts and RPKM | no (not requested) |
| 2 | [GSE151648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151648) (on disk) | liver, pre and post | 40 recipients (all 40 with both), IRI+ 23, IRI- 17 | no (field `iri` only) | no | none | RNA-seq counts | already on disk |
| 3 | [GSE43974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43974) | kidney T1 donation, T2 after cold ischaemia, T3 60 min after reperfusion | T2 and T3 pairs: brain-dead 67, DCD 29 (96 total), flagged by a `pair` code of 2 (T2) and 3 (T3) | no | brain-dead, DCD, living | none | Illumina HT-12 v4 per-sample matrix, 554 biopsies | yes (series matrix, 120 MB) |
| 4 | [GSE15480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15480) | liver, pre-implantation and post-reperfusion, with or without ischaemic preconditioning | 12 grafts (6 per arm) | no | deceased | none | Affymetrix | no (< 20) |
| 5 | [GSE87487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87487) | liver pre and post reperfusion | 10 grafts (IRI- 6, IRI+ 4) | no | no | none | RNA-seq counts | no (< 20) |

GSE43974 problems: the sample sheet has no graft/donor id. The title suffix is a per-group, per-time-point counter (BD: 101, 110, 105 samples at T1, T2, T3)
so it cannot be used to link T2 to T3 of the same kidney; the T1 (before retrieval) biopsies come from a different cohort than T2/T3 according to the
series design text. The pair code identifies that 96 T2 and 96 T3 biopsies are in pairs but not which T2 goes with which T3. The linking table would have
to come from the authors or from a supplementary file (`GSE43974_paired_analyses_DBD_DCD_Living.xls.gz` holds analysis results, not checked). Cold ischaemia
time, circadian time of retrieval or implantation and recipient id are all absent from GEO. Detection: all 14 genes present on the array
(medians 41.0 to 59.9 percentile, i.e. low detection on this platform and kidney).
Verdict for the "graft carries recipient or donor time" question: none of the bulk sets records clock time of retrieval or reperfusion, so only the
untimed donor-referenced design (does the post-reperfusion graft move toward the recipient's clock state) could be asked, and only with the authors' pair table.
Cold ischaemia time was not found in the GEO metadata of any of the checked series; it may exist in the papers (not checked against full text).

## Files

- Samples: `output/psych_GSE80655_samples.csv`, `psych_GSE102556_samples.csv`, `psych_GSE202537_samples.csv`, `psych_GSE53987_samples.csv`, `transplant_GSE43974_samples.csv`
- Design tables: `output/psych_design_tables.csv`, `output/transplant_GSE43974_design.csv`
- Clock gene presence/percentile: `output/clock_gene_detection.csv`; GSE102556 matrix column map: `output/psych_GSE102556_matrix_columns.csv`
- Scripts: `scripts/scout_build.py` (needs GPL570.txt and GPL10558.txt platform tables from GEO in a scratch dir, passed as argv[1]), `scripts/geo_esearch.py`.

## Caveats

- Counts of donors rely on GEO sample sheets (GSE80655, GSE102556) or inferred links (GSE53987, Chang): see donor_id_derived column in each samples file.
- GSE102556 FPKM scale is log-like (values around 3 for expressed genes), so percentiles are rank-based and unaffected.
- No full-text check of papers was done; statements about papers (Malone WES genotyping, cold ischaemia in papers) are from memory or GEO text.
- The Chang, GSE35978, GSE92538, GSE12649, GSE5388 and the B(i) candidates were examined from sample sheets only; their expression files were not downloaded.
