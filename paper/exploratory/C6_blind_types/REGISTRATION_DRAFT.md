# A182 registration DRAFT (not committed, not in registered_experiments.md). Items marked [OPEN] need the PI / user before the commit.

## A182 (2026-11-01, written after a metadata-only scout (`exploratory/scouts_2026-11-01/C6_blind_cohorts.md`) and a header-only summary (`exploratory/C6_blind_types/output/header_summary.csv`: header lines, the ID_REF column and platform annotation; no expression value, no clock-gene value, no score) and before any expression value of the series below is read; the TCGA values of the same types are seen (A171 `tcga_all.csv`, A172); HANDOFF_2026-10-31 C6): is the proliferation clock-coherence loss specific to the cancer types where it replicated (breast, lung adenocarcinoma, gastric), or does it appear in any tumour type, including types whose TCGA value was weak or null and squamous carcinomas?

**Why:** the claim names types (HANDOFF stop rule: no new type enters the claim unless scored blind). Breast (A127, A166), lung adenocarcinoma (A171) and gastric (A173: 2 of 3) were chosen or confirmed after TCGA values; the lung squamous contrast was in the predicted direction but failed (A173 H173.4 p 0.058). C6 scores types whose TCGA value was weak or null, plus squamous carcinomas, in cohorts never read. Readout level **L1** (single compartment per person: coherence only; no host, no L2/L3 claim).

**Data** (`$SCCC_DATA/c6_*`, GEO series matrices, parsed with the A171 parser (`A171_luad_lihc/scripts/prep_a171.py`); n and fields from the headers, `header_summary.csv`):

| family | cohort | type | platform | samples | tumours kept (field) | clock genes of 14 |
|---|---|---|---|---|---|---|
| 1 weak/null | GSE39582 (CIT, frozen, RMA + ComBat) | colon | GPL570 | 585 | **566** (`dataset` discovery / validation; 19 non-tumoral mucosa dropped) | 14 |
| 1 weak/null | GSE17536 [OPEN] (Moffitt) | colon | GPL570 | 177 | **177** (all; no tumour-type field) | 14 |
| 1 weak/null | GSE116918 [OPEN] (diagnostic biopsies, radiotherapy cohort; Almac DSA) | prostate | GPL25318 | 248 | **248** (`tissue` Primary Prostate Cancer) | **13 (no NR1D1)**; purity proxy 10 of 11 (no CD3E) |
| 2 squamous | GSE65858 (Heidelberg; RSN + ComBat) | head and neck | GPL10558 | 270 | **253** (`tumor_type` Primary; 14 Secondary, 3 Relapse dropped) | 14 |
| 2 squamous | GSE157009 (Alliance I; MAS5) | lung squamous | GPL570 | 249 | **249** (all SCC) | 14 |
| 2 squamous | GSE157010 (Alliance II; MAS5) | lung squamous | GPL570 | 235 | **235** (all SCC; 2 with Stage NA kept) | 14 |
| 2 squamous | GSE53625 [OPEN] | oesophageal squamous | GPL18109 | 358 | 179 (source_name `cancer tissue`; 179 paired normals dropped) | **not mappable** |

Pipeline's 14 clock genes = ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART (`a168_lib.GENES`); BHLHE40 and BHLHE41 (in the HANDOFF list) are not in it and are not used [OPEN]. G2M set (Hallmark G2M without clock genes, 200 genes): 192 present on GPL570, 187 on GPL10558, 186 on GPL25318. Purity proxy: 11 of 11 except GSE116918 (10; no CD3E). Pseudo-gene pool (probes with a symbol, minus CLOCK_LIKE, G2M, purity and clock genes): 22,596 (GPL570), 16,718 (GPL10558), 19,396 (GPL25318). No two cohorts share a GEO sample accession (GSE39582, GSE17536, GSE14333, others; `header_overlap.csv`); patient-level overlap cannot be checked from headers.

**Definitions** (= A171 / A173): linear matrices (max > 100) → log2 with the smallest positive value as floor; probes averaged per symbol; Y = ln 2 × log2 expression of the clock genes present; proliferation = Hallmark G2M mean z, purity proxy = A97 11-gene mean z, both z-scored over the tumours of the cohort. Coherence = `a127_lib.coh`.
**Tests (one-sided for a loss, 2,000 permutations):** the H171 statistic per cohort: strata = purity-proxy tertiles (≥ 30 tumours); T = mean over strata of coh(high proliferation tertile) − coh(low); proliferation permuted within stratum. Holm within each family.
- **Family 1 (weak/null TCGA value): H182.1 GSE39582, H182.2 GSE17536 [OPEN], H182.3 GSE116918 [OPEN].**
- **Family 2 (squamous controls): H182.4 GSE65858, H182.5 GSE157009, H182.6 GSE157010, H182.7 GSE53625 [OPEN, only if mappable and ruled blind].**

**Seen TCGA value of each type (same statistic, one-sided p for a loss; `A171_luad_lihc/output/tcga_all.csv`, `A172_tcga_unseen/output/per_project.csv`):**

| type | TCGA | n | T | p | pseudo percentile | source |
|---|---|---|---|---|---|---|
| colon | COAD | 458 | −0.108 | 0.146 | not run | A171 tcga_all.csv |
| rectum (no GEO cohort) | READ | 166 | −0.029 | 0.416 | 45.0 | A172 per_project.csv |
| prostate | PRAD | 497 | −0.101 | 0.159 | not run | A171 tcga_all.csv |
| head and neck | HNSC | 502 | −0.170 | 0.079 | 12.5 | A172 per_project.csv |
| oesophagus, any histology | ESCA | 184 | −0.149 | 0.095 | 16.5 | A172 per_project.csv (not split by histology; squamous subset not computed) |
| lung squamous | LUSC | 501 | +0.212 | 0.95 | 96.0 | A172 per_project.csv |
| lung squamous, GEO (seen, descriptive) | GSE37745 / GSE41271 / GSE30219 | 66 / 80 / 61 | −0.061 / −0.087 / −0.072 (mean −0.073, p 0.24) | | | A173 descriptive.csv |

For scale, the replicated types: BRCA −0.118 (TCGA, A171 tcga_all.csv), LUAD −0.361 (TCGA) and −0.242 / −0.293 (GSE72094, GSE68465; A171 tests.csv), gastric −0.187 / −0.533 (GSE62254, GSE15459; A173 tests.csv). Weak/null TCGA values here (−0.03 to −0.11) are below every MDE of A171-A173 (0.17-0.36), so a non-significant result in Family 1 is an absence of evidence of a loss of the replicated size, not equivalence. HNSC and ESCA TCGA values are not null (−0.15 to −0.17, 12th-17th pseudo percentile, Holm-uncorrected p 0.08-0.09): "squamous control" is a prediction, not a known null.

**Predictions (HANDOFF C6):** the types where the loss replicated lose coherence with proliferation; these types do not, unless the TCGA value was a false negative.
- **P1 (Family 1):** no H182.1-3 passes clock-specifically; if one does, the TCGA value for that type was a false negative and the type joins the claim.
- **P2 (Family 2):** no H182.4-7 passes clock-specifically; a pass means the loss follows proliferating carcinomas beyond adenocarcinoma, and the "type-specific, not histology-specific" statement is retracted.
- Direction check (descriptive): the mean T of Family 1 and of Family 2 are compared with the replicated-type cohort values (−0.19 to −0.53, passing A171 / A173 cohorts).

**Specificity companion:** 200 pseudo-clock sets per cohort (A168 construction; pool without clock-like, G2M, purity-proxy or clock genes; number of clock genes present, 13 in prostate); same permutations. "Clock-specific" only if T is below the 5th pseudo percentile. Evaluated for every cohort, since a pass-but-not-specific result is the "general co-regulated gene sets" reading; also a pseudo-set companion for the null reading: a fail is read only with the pseudo-set SD for the cohort.

**Power (after the commit, before the tests; `output/power.csv`):** 2,000 permutations; MDE = 2.49 × null SD per cohort, listed next to the TCGA T of the type and the replicated-type reference values above. Cohorts with MDE > 0.30 are reported as "not estimable at this n" and not read (HANDOFF stop rule); a pooled-per-type descriptive (unweighted mean of the type's cohorts) is added with its own null.

**Descriptive (no Holm, no claim):** each test without PER1 and NR1D1; T per stratum; GSE39582 by `dataset` (discovery / validation), MMR status and tumour location; GSE65858 by HPV status (hpv_dna) and tumour site; GSE116918 with 13 genes versus without the other poorly measured clock genes; pseudo-percentile only for Family 2 pooled.

**Reading:**
- All cohorts fail (or are clock-nonspecific) at powered designs → the loss is specific to the types that replicated (breast, lung adenocarcinoma, gastric): type-specific, not histology-specific, not pan-adenocarcinoma. Bounded at each cohort's MDE.
- A Family 1 cohort passes clock-specifically → its TCGA value was a false negative; add the type with its blind cohort as the evidence.
- A Family 2 cohort passes clock-specifically → proliferating squamous carcinomas also lose coherence; the claim drops the histology wording (cf. A173 H173.4).
- A cohort is underpowered (MDE > TCGA value of its type and > 0.30) → "not estimable"; this cannot support type specificity.
- A failed gate-like check (e.g. fewer than 13 clock genes detected, platform map failure) is uninformative, never "no loss".

**Code:** `paper/exploratory/C6_blind_types/` (`scripts/header_c6.py` is header-only; `prep_c6.py` and `run_c6.py` to be written after the commit, from `A171_luad_lihc` / `A173_gastric_lung_histology`).

## Open design choices (for the user)

1. **GSE17536 inclusion.** Moffitt, 177 colon cancers, GPL570. Its only shared GSM overlap with GSE39582 and GSE14333 is zero (header accessions), but patient-level overlap cannot be seen from headers; GSE14333 (Jorissen, 290) could not be downloaded, so colorectal rests on GSE39582 plus GSE17536. GSE17536 is a SubSeries of GSE17538 (with the VMC training set, which would be a different cohort). Include as H182.2 (separate test, Holm family 1) or as a descriptive replication only? n = 177 is above the ≥ 150 rule.
2. **GSE39582 is ComBat-corrected** (batch correction by the authors: RMA + ComBat). Same for GSE65858 (RSN + ComBat). Registered as is; ComBat may remove batch structure shared across genes (unknown effect on coherence). Keep, or flag as a caveat?
3. **Prostate custom array.** GSE116918 (Almac DSA, FFPE biopsies; radiotherapy cohort, not surgery) maps 13 of the 14 clock genes (NR1D1 has no annotated probe in the GPL25318 table) and 10 of 11 purity genes. Include with 13 genes (as GSE26253 and GSE68465 were handled), or drop? There is no second prostate cohort ≥ 150 (scout).
4. **ESCC (GSE53625).** (a) Blindness: it appears as a relation line (SuperSeries of GSE53624) in `scouts_2026-10-31/gastric_lung_histology/hdr/GSE53624_series_matrix.hdr`; no value was read. (b) Probe mapping: the GEO platform GPL18109 uses Agilent feature numbers and has no gene symbols (only sequences), GPL19748 also none, so no clock-gene coverage can be verified from the allowed annotation. Mapping would need a vendor annotation file (not permitted by the brief) or sequence alignment. Recommendation: drop ESCC (oesophageal has no cohort), or ask for an annotation source. TCGA ESCA (184, mixed histology) is not a squamous-only value.
5. **Gene definition.** The pipeline's 14 (with RORA and CIART) differs from the HANDOFF list (BHLHE40 and BHLHE41 instead of RORA and CIART). The draft uses the pipeline's 14 for comparability with A171-A173; confirm.
6. **Missing types.** Pancreatic, endometrial and cervical have no clean blind cohort (scout), so the C6 weak/null family is colon and prostate only (TCGA READ not testable; PAAD −0.097, UCEC −0.143, CESC −0.175 stay unscored blind). Accept the narrower claim, or admit GSE71729 (145 pancreatic, near miss) and GSE44001 (cervical DASL, risky)?
7. **Squamous prediction vs TCGA.** HNSC and ESCA TCGA values are negative (−0.17, −0.15), so a loss in HNSC is not implausible; the "no loss" prediction is weak for HNSC. Option: state HNSC as a two-sided descriptive control rather than a Family 2 test.
8. **Null reading.** Weak/null TCGA T (−0.03 to −0.11) is below the MDE (expected ~0.18-0.25 in 250-566 tumours). Pre-specify an equivalence bound (e.g. T > −0.20 with a 90 % permutation interval) or accept "bounded at MDE"? Without one a non-significant cohort cannot support "does not lose".
9. **Two lung cohorts.** GSE157009 and GSE157010 are two arms of one Alliance study (SuperSeries GSE157011): separate tests (current) or pool?
10. **Scale check.** Whether each matrix is already log2 is not visible from headers; the registered A171 rule (max > 100 → log2) is applied at prep time. GSE157009 / GSE157010 are MAS5 (linear expected).
11. **Per-cohort Holm.** Holm within each family as in A171 / A173; no overall family across all seven, since the two families make opposite predictions.
