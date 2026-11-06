# Literature tally: clinical cancer studies claiming that clock readouts predict outcome

Question (user, 2026-11-06): the cancer-clock literature is dominated by breast cancer and almost all clock-outcome
claims are positive. Is that consistent with publication / survivor bias, and were claims adjusted for proliferation
and molecular subtype? Not pre-registered; exploratory, descriptive, post hoc. Nothing here is a hypothesis test.

Files: `tally.csv` (79 rows, one per study or per-paper cohort group), `build_tally.py` (the curated table),
`tally_stats.py` -> `stats_output.txt` (all counts below). Per-row confidence is in `read_level` (F full-text card of the
project; A abstract / search snippet; M study known only from a meta-analysis characteristics table; I inferred from
title or memory) and `claim_basis`. **Fields marked "n" for proliferation / random-set control on A, M and I rows are
inferred (default when the snippet does not mention them), not read.** Only 6 rows are full-text reads.

## 1. Method and search log (reproducible, all 2026-11-06)

Sources: project cards (`paper/manuscript/literature/`, 148 cards; prognosis-relevant: cadenas2014loss, li2024tumor,
vasilyev2026robust, vlachou2024timeteller, ye2018genomic, wu2019pan, hammarlund2025cyclops), `literature_positioning.md`
section 10, `NOVELTY_ASSESSMENT_2026-11-06.md` section 1, archive `scCircadianCancer/docs/literature/literature.md`.
The PubMed and bioRxiv MCP tools were listed but need authentication; I did not authenticate. Web searches used
WebSearch (standard mode, summaries of result snippets only, US-only), WebFetch for two meta-analyses. Search strings:

1. circadian clock gene expression prognosis survival breast cancer PER1 PER2 PER3 CRY BMAL1 cohort
2. circadian clock gene signature prognostic risk score hepatocellular carcinoma TCGA validation
3. circadian rhythm-related genes prognostic signature lung adenocarcinoma survival
4. clock gene expression colorectal cancer prognosis PER3 CRY1 BMAL1 immunohistochemistry survival
5. circadian clock genes pan-cancer prognosis TCGA survival analysis clock disruption
6. circadian clock gene prognostic gastric OR pancreatic OR glioma OR renal OR ovarian cancer signature survival
7. PER2 expression breast cancer prognosis overall survival clock gene negative result OR no association
8. circadian clock genes glioma prognosis signature CGGA TCGA
9. circadian rhythm gene signature colon cancer OR colorectal prognostic model immune
10. circadian clock gene expression prostate OR bladder OR kidney cancer survival prognosis bioinformatics
11. BMAL1 ARNTL expression prognostic marker cancer patients survival meta-analysis clock genes
12. NPAS2 OR CLOCK OR CRY1 expression prognosis breast cancer OR non-small cell lung cancer survival
13. breast cancer circadian rhythm gene signature prognostic model TCGA METABRIC risk score
14. clock gene expression head and neck OR esophageal OR cervical OR endometrial cancer prognosis PER1 PER3 survival
15. circadian clock acute myeloid leukemia OR lymphoma OR myeloma clock gene expression survival prognosis
16. tumor circadian clock disruption score inferred phase predicts survival immunotherapy response ordering CYCLOPS OR TimeTeller OR CHIRAL
17. clock gene expression not associated with survival cancer cohort null result PER BMAL1 prognosis "no significant"
18. circadian rhythm gene signature melanoma OR pancreatic OR gastric OR bladder OR esophageal prognostic model immunotherapy 2023 2024
19. PER2 breast cancer ER-positive immunohistochemistry survival cohort circadian clock protein expression outcome
(plus one lookup of PMID 30791227). WebFetch: PMC7241715 (Zhang et al. meta-analysis of low clock-gene expression, 36
studies, search to 1 Jan 2020) and PMC7436287 (PER-family meta-analysis, search to 2 Apr 2019). The 37 `M` rows come from the
PMC7241715 Table 1 (n, cancer, gene, method) and PMC7436287; their per-study HRs sit in forest-plot images I could not read.

Inclusion: clinical tumour material or cohort data; a clock gene, clock score or inferred phase; an outcome (survival,
recurrence, response). Excluded: cell/animal-only work, SNP/epidemiology studies, expression-only (tumour vs normal),
the meta-analyses themselves (used as sources), trials of dosing time (separate question; Landre 2024 and Qian 2021 are
carded elsewhere). Not done: full PubMed query, dual screening, second-reviewer extraction, Embase, Chinese databases.

## 2. Counts (79 rows, 2007-2026; year bins 2005-09: 4, 2010-14: 17, 2015-19: 25, 2020-24: 22, 2025-26: 11)

Claim coding: pos = headline association / prediction claimed; mixed = positive headline with explicitly non-significant
genes, cohorts or models; null = main clock readout not associated; unresolved = per-study result not obtainable.

| cancer group | rows | pos | mixed | null | unresolved |
|---|---|---|---|---|---|
| Colorectal | 18 | 14 | 0 | 2 | 2 |
| Breast | 17 | 13 | 1 | 0 | 3 |
| Lung | 6 | 5 | 1 | 0 | 0 |
| Liver | 6 | 6 | 0 | 0 | 0 |
| Pan-cancer | 5 | 4 | 1 | 0 | 0 |
| Glioma | 5 | 5 | 0 | 0 | 0 |
| Gastric | 5 | 4 | 1 | 0 | 0 |
| Head and neck | 5 | 4 | 0 | 0 | 1 |
| Gynaecologic | 4 | 2 | 1 | 0 | 1 |
| Urologic, melanoma, haematological | 2 each | 6 | 0 | 0 | 0 |
| Pancreatic | 2 | 1 | 0 | 0 | 1 |

- Fraction positive: 69 of 71 resolved rows (97%) have a positive or mixed headline; strictly positive without a
  stated null element 64/71 (90%); mixed 5; null 2/71 (3%; both colorectal, small series); 8 unresolved.
- Breast is the largest single group but not dominant in this set: 17/79 (22%), colorectal 18/79 (23%). Breast
  headline positives 14/14 resolved (100%); non-breast 55/57 (96%). The user's premise "dominated by breast" is only
  partly supported by this search; breast dominates the inferred-phase / ordering-based end (TimeTeller, CYCLOPS,
  COFE), where all four carded outcome papers are breast.
- Readout: 52 rows single gene/protein (42/44 resolved positive or mixed), 27 signature/score/index rows (27/27).
  Signature papers are positive by construction (LASSO-Cox selection on the outcome, then a median split).
- Ordering-based phase / disruption readouts (CYCLOPS magnitude, TimeTeller x2): 3 rows, all breast, all positive, all
  from 2019-2026, directions disagree (CYCLOPS: strong rhythm worse; TimeTeller REMAGUS: good function better in grade
  1-2; Vasilyev 2026: good function antagonistic to survival).
- Independent validation cohort: yes 14, no 49, unclear 16 (breast: 3 yes, 9 no, 5 unclear). Mostly signature papers.
- Adjusted for proliferation: none fully (0/79). Partial (grade or MKI67 as a term, no proliferation metagene or G2M):
  3 (Cadenas grade; Li MKI67; Vlachou grade). Unclear 2. Rest "n" and mostly inferred. Several signature papers report
  cell-cycle enrichment in the high-risk group rather than adjusting for it (snippets).
- Adjusted for molecular subtype: yes 1 (Vlachou: ER, PR, HER2, TN), partial 8 (stratified, restricted or molecular-profile-aware: Cadenas, Li luminal A, Ye, Chen, three glioma, one breast subtype paper), no 24, unclear 46. Among 17 breast rows: yes 1, partial 3, unclear 13.
- Random-gene-set or matched control: 0/79 (read for the 6 full-text rows; inferred for the rest).
- Paper reports its own null genes/cohorts: yes 8, partial 4, unclear 1, no 66 (no/unclear mostly inferred for A, M rows).

### Funnel-style table (rows with n known; 28 rows lack n, mostly TCGA signature papers)

| n | rows | pos | mixed | null | unresolved | headline positive of resolved |
|---|---|---|---|---|---|---|
| < 100 | 19 | 11 | 1 | 1 | 6 | 12/13 (92%) |
| 100-249 | 18 | 14 | 1 | 1 | 2 | 15/16 (94%) |
| 250-599 | 5 | 5 | 0 | 0 | 0 | 5/5 |
| >= 600 | 9 | 7 | 2 | 0 | 0 | 9/9 |

No small-study gradient is visible: 92% small vs 100% large. This is not evidence against small-study bias; the large
rows are TCGA/pooled signature papers whose positivity is built in, and 6 of 19 small rows are unresolved.

### What the two meta-analyses add (read, not mine)

- Zhang 2020 (PMC7241715; 36 studies, OS): pooled HR (low expression) Per1 1.35 (1.06-1.72), Per2 1.43 (1.10-1.85),
  Per3 1.32 (0.99-1.76); Cry1, Cry2, Npas2, Bmal1, CLOCK all n.s. (pooled HR 0.79-1.10); I2 60-86%. Begg p 0.23-0.76;
  Egger p 0.12-0.91 except Per3 p 0.033. Authors read the funnels as symmetric; they list unpublished and non-English/
  non-Chinese studies as not captured. So at gene level, 2 of 8 genes show a pooled effect even though almost every
  individual paper has a positive headline: headline positivity is multi-gene, multi-subgroup selection, not a
  consistent gene-level effect. Their introduction cites at least six null findings (Hu 2014, Oshima 2011, Tokunaga 2008,
  Qiu 2019 lung, Bmal1/CLOCK in lung and CRC).
- PER-family meta-analysis (PMC7436287; 20 gene-level studies, 1,369 patients): protein PER1/2/3 HR 0.67/0.55/0.47,
  PER2 mRNA 0.92 (0.80-1.06) and 0.63 after dropping one study; Begg p all > 0.7, Egger not reported. Protein (IHC)
  results are stronger than mRNA ones, and IHC series are small (29-246).

## 3. Reading for the user's question

Is this consistent with publication or survivor bias? It is consistent, but it does not show it. What fits: 97% positive or mixed
headlines among resolved rows; 2 nulls found only through other papers' introductions or snippets; zero random-set or
proliferation-matched controls; small IHC/qPCR series with flexible cut-points (SP scores, median splits); signatures
selected on outcome; gene-level pooled effects far weaker than paper-level positivity (a selective-reporting pattern:
authors report the genes and subgroups that reach significance). What does not fit or cannot be judged: no small-study
gradient in this sample; meta-analytic Egger tests were non-significant except Per3 (low power at I2 60-86%);
and positivity has a second, non-bias explanation. Clock genes (TIMELESS, CRY2, PER2/3, NPAS2) track the proliferation
metagene and ER status, and proliferation predicts outcome in breast cancer (project card Cadenas; positioning section 10;
A192: the Cadenas score sits at the 7th percentile of proliferation-matched random sets; HR 0.81 falls to 0.88 with G2M).
So "positive" is expected without any bias, which is why the missing proliferation adjustment and missing random-set
controls matter more than the positive fraction. "Survivor bias" in the strict sense (tumour samples from patients
alive at resection, left truncation) was not assessable from snippets; cohorts such as Cadenas are node-negative,
untreated, resected, so conditional on surgery.

Adjustment: proliferation 0/79 fully adjusted (3 partial); subtype 1 yes, 8 partial, 46 unclear. The strongest
breast papers (Cadenas, Vlachou, Li) adjust for grade/ER/HER2 or restrict to luminal A, and none uses a proliferation
metagene, G2M or matched gene sets. Our A191-A194 reruns (positioning section 10) found no clock-specific information
beyond proliferation in METABRIC/SCAN-B/TCGA, direction reproduced.

## 4. Limits of this tally (honest)

- Selection: the search strings targeted claims of prediction, found by a snippet-summarising web search; no PubMed query
  with counts, no screening log, no second reviewer. Studies found are the findable and mostly open-access (PMC, Frontiers,
  MDPI, PeerJ) ones, which favours bioinformatic signature papers (27 rows; the 2021-2026 wave). The project's
  carded set is not a sample of this literature (it was built for methods and biology; only 6 rows are outcome papers).
- 37 of 79 rows come from meta-analysis tables; for 24 rows the claim direction is inferred from title or memory
  (`claim_basis` = I) and for 8 unresolved. Adjustment flags for A, M and I rows are mostly inferred defaults. Rerun with
  full texts before any number goes into the manuscript (project rule: numbers in text come from output tables).
- Null results are structurally hard to find by search ("not associated" is rarely in a title); the 2-null count is a
  floor on retrieval, not an estimate of the null rate. A funnel plot needs per-study effect sizes (not extracted here).
- Overlap: several signature papers reuse TCGA, so rows are not independent; the 4 Qiu cohorts are one paper.
- Ordering-based outcome papers (3 rows) are too few for any rate.
- Next steps if wanted: full PubMed export (authenticated), extract HR/CI per study for a funnel plot / Egger, two-reviewer
  adjustment coding on full texts, and a registered amendment before scoring anything.
