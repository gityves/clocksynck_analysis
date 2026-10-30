# Scout: neoadjuvant endocrine therapy, paired pre / on-treatment biopsies (metadata only)

Date 2026-10-30. No clock statistic, clock-gene correlation or clock-gene expression value was computed or read. Only series-matrix header
lines, sample tables, platform annotations (probe existence) and the probe-ID column of the matrices were parsed (`python -I`).
Scripts: `scripts/read_meta.py`, `summarise_cohorts.py`, `probe_check.py`, `design_checks.py`, `batch_dates.py` (outputs `*.out` beside them).
Data: `/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_endocrine/<GSE>/` (221 MB, `PROVENANCE.md` with URL, date, sha256).
Licence: GEO carries no per-series licence statement; all are public GEO series (no controlled access). Verified = read from the downloaded
metadata or the GEO page; "unverified" = from memory or inference.

## Ranked recommendation (paired patients x platform quality)

| rank | cohort (GEO) | drug | n paired patients | time points | platform | clock genes (of 14 + BHLHE40/41) | on disk |
|---|---|---|---|---|---|---|---|
| 1 | POETIC, GSE105777 | aromatase inhibitor (agent per patient not in metadata) | 157 AI pairs (B, S); plus 56 untreated control pairs (D, S) | baseline core biopsy vs surgery after 2 wk pre-surgical AI | Illumina HT-12 v4 (GPL10558), RSN-normalised | 16/16, all in matrix | yes (matrix only) |
| 2 | Edinburgh/Novartis letrozole, GSE20181 (contains GSE5462) | letrozole 2.5 mg | 58 pre+2wk; 56 pre+3m; 56 with all three (62 patients) | pre, 10-14 d, 90 d | Affymetrix U133A (GPL96) | 15/16, CIART has no probe | yes |
| 3 | Z1031B, GSE87411 | AI (exemestane/letrozole/anastrozole in Z1031; arm not labelled in GEO, unverified) | 109 patients (92 labelled BL + post; 17 "re-analyzed" + post, see notes) | BL vs 2-4 wk | Agilent 4x44K two-colour vs universal reference (GPL6480), log2 ratios | 16/16 | yes |
| 4 | NeoPalAna, GSE93204 | anastrozole alone for 4 wk (C1D1), then + palbociclib (C1D15) | 22 BL+C1D1 (anastrozole only); 22 BL+C1D15; 18 with all three | BL, C1D1, C1D15, surgery | Agilent two-colour (GPL6480) | 16/16 | yes |
| 5 | Edinburgh long-term, GSE111563 | letrozole, extended | 37 pre (tp1) + long-term (tp4); 24 with tp1, tp2, tp4 | tp1 <0 d, tp2 0-120 d, tp4 >120 d | Illumina HT-12 v4 (GPL10558), probe-filtered | 14/16 in matrix (NR1D1, CIART absent) | yes |
| 6 | Edinburgh Illumina letrozole, GSE59515 | letrozole | 25 (all three time points) | pre, 2 wk, 3 mo | HT-12 v4 | 16/16 | yes |
| 7 | Edinburgh ILC/IDC letrozole, GSE55374 | letrozole | 13 (10 with 2 wk) | pre, 2 wk, 3 mo | HT-12 v4 | 16/16 | yes |
| 8 | Baylor AF +/- gefitinib, GSE33658 | anastrozole + fulvestrant (+ gefitinib in 6) | 11 | pre vs 3 wk | Affymetrix U133 Plus 2 (GPL570) | 16/16 | yes |

n values are from `scripts/summarise_cohorts.out` and `design_checks.out`. Only ranks 1, 2, 3 reach the >= ~40 pairs wanted; rank 2 and 5-7 are the same
patient pool type (Edinburgh) on different platforms and cannot simply be pooled (see overlaps).

## Per-cohort notes

### GSE105777 (POETIC, Institute of Cancer Research / Royal Marsden)
- Verified: 426 samples = 157 treated patients x (baseline B, surgery S) + 56 untreated patients x (diagnosis D, surgery S, titled "[reanalysis]").
  Title: "Major impact of sampling methodology on gene expression in ER+ breast cancer" (no pubmed id in the GEO record). Sample-to-patient link: title
  index (`AI.Treated.<n>B/S`, `<n>D/S [reanalysis]`). Matrix 79 MB (82.4 MB on disk). Public Oct 2017.
- Only ER+ tumours (disease field "primary ER+ breast carcinoma"); the design tests ER-signalling blockade inside ER+, it does not contain ER-negative tumours.
- The on-treatment sample is the surgical specimen, the baseline a core-cut biopsy. The authors' point (summary) is that sampling method changes expression
  (FOS, DUSP1, RGS1, FOSB) equally in treated and untreated. The 56 untreated D/S pairs are therefore the required specificity companion for any
  treated change (same permutations: difference-in-difference). Immediate-early / clock-adjacent genes are at risk; no clock value was inspected.
- Treated arm: summary says "2-week's presurgical aromatase inhibitor (AI, n=157)"; which AI, and Ki67 / response labels, are not in the GEO metadata (unverified: POETIC randomised AI vs none, outcome data held by the trial).
- Matrix: 47,323 probe rows, 16/16 clock-gene probes present (RORA 8 probes). Control cohort overlap: see below.

### GSE20181 and GSE5462 (Edinburgh/Novartis letrozole, Miller et al.)
- Verified: GSE5462 (116 samples, 58 patients, pre and 10-14 d) is fully contained in GSE20181 (identical 116 GSMs; GSE20181 adds 60 day-90 samples).
  GSE20181: 62 patients, 176 samples; 56 patients have all three. Responder labels in the titles (responder 110, nonresponder 44, not assessable 22 samples; per-patient counts not computed).
- Patient link: leading number of the title, suffix A pre / B 10-14 d / C 90 d.
- Platform GPL96: no CIART probe; matrix 11 MB. Public Jul 2010.
- Batch confound (verified from submission dates): A and B arrays submitted 4 Aug 2006, C (90 d) arrays 3 Feb 2010. Pre vs 3 months is plausibly confounded with array batch; pre vs 2 wk is within one submission. Hybridisation dates are in the CEL files (not downloaded): unverified.

### GSE59515 / GSE55374 / GSE111563 (Edinburgh, Illumina HT-12 v4; Turnbull, Dunbier)
- GSE59515: 25 patients x 3 (pre, 2 wk, 3 mo), 17 responders / 8 non-responders (verified, from patient characteristics). All 75 samples submitted 17 Jul 2014.
- GSE55374: 13 patients (ILC 7, IDC 6 by pairs; 18 ILC + 18 IDC samples), 13 pre, 10 two-week, 13 three-month; histology and response in metadata.
- GSE111563: 39 patients, 101 samples: 37 pre, 25 early, 39 long-term (>120 d). Timepoints are windows (<0, 0-120, >120 d), not fixed days; drug "letrozole" per the summary. The deposited matrix keeps 22,119 probes (NR1D1 and CIART probes absent); a non-normalized file exists in the supplement (not downloaded).
- Patient overlap between them: not proven either way. Titles and GSMs are disjoint and the ID formats differ (`patient150`, `X109`); the GSE111563 design text says its data were "analysed with" 14 patients (GSE59515) and 9 patients (GSE55374), which suggests the three series are one extended Edinburgh pool whose new patients are in GSE111563, but the numeric ID overlap (3 with GSE59515, 1 with GSE55374) is weak evidence and could be coincidence. Needs the Turnbull/Dunbier supplementary patient list.
- GSE55374 design text states 7 ILC + 8 IDC patients were profiled on Affymetrix and "see GSE20181": those 15 would be distinct from the 13 Illumina patients; GSE20181 IDs (2-69) share no number with GSE55374 or GSE59515 (150-374).

### GSE87411 (ACOSOG Z1031B, Washington University)
- 109 patients, 218 arrays; text says it re-analyses GSE29442, GSE35186, GSE35191. Patient ID in the title; time labels BL / w2 / W2 / M / Redo variants, and 17 patients carry a baseline-like array labelled "re-analyzed" (13 pair with w2, 3 with W2, 1 with M). Reading it as the baseline is an inference (unverified).
- Two-colour log2 ratios against a Stratagene universal reference spiked with MCF7 and ME16C RNA; the reference contains clock genes (cell-line pool) so ratios are pre vs reference rather than absolute, and are loess-normalised (not an obstacle to correlations, but the reference noise is shared across arrays).
- No arm, response or Ki67 field in GEO (the trial randomised three AIs: unverified). Post sample at 2-4 wk. Matrix 9.6 MB, 29,032 probe rows.

### GSE93204 (NeoPalAna, NCT01723774)
- Verified: baseline n=32, C1D1 n=33, C1D15 n=29, surgery n=24; summary says C1D1 follows anastrozole and C1D15 follows anastrozole + palbociclib. 22 patients have BL + C1D1 (anastrozole alone; endocrine-only effect), 22 BL + C1D15. 46 patients in total. Cell-cycle inhibition after C1D1 makes C1D15 and surgery unsuitable for an endocrine-only contrast. Two-colour as GSE87411. The 4-week C1D1 time is from the trial design (unverified in GEO text).

### GSE33658 (Baylor, Anastrozole + Fulvestrant +/- gefitinib)
- 11 patients (AF 5, AFG 6), pre and 3 wk, all ER+; response PD/SD/PR/CR labelled; baseline age and tumour size present. The only series with fulvestrant. Too small alone.

## Overlaps between GEO series (verified unless noted)
- GSE5462 subset of GSE20181 (116 shared GSMs). Never count both.
- GSE105777 untreated control (112 samples, 56 D/S pairs) = GSE73237 GPL10558 subseries re-deposited: all 112 control titles match GSE73237 titles. GSE73237 (superseries of GSE73235 and GSE73236; 158 samples, GPL10558 124 + GPL6884 34 samples) also has 12 extra samples (A/B pairs "after excision" / "after x-ray", sampling-method experiment, 6 patients) and the GPL6884 earlier-chip part; both parts are untreated (GPL6884 part not inspected).
- GSE87411 re-analyses GSE29442, GSE35186, GSE35191 (stated in the record; not checked). 
- GSE59515, GSE55374, GSE111563: see above (unresolved).
- Edinburgh letrozole across platforms: GSE20181 (Affymetrix) vs the three Illumina series: no shared GSM, no shared patient IDs among 20181 / 59515 / 55374.

## Design notes
- Natural primary: the 2-week contrast (pre vs 2 wk). It exists in GSE105777 (157 pairs), GSE20181 (58), GSE59515 (25), GSE55374 (10), is within one array submission for GSE20181 and is the time at which Ki67 / ER-target changes are largest (unverified claim from the literature, not from the files). A 3-month contrast is secondary (batch confound in GSE20181; Illumina series are single-batch by submission date, hybridisation dates unverified).
- Primary dataset: POETIC (GSE105777) with its untreated 56-pair control as the specificity companion; replication in GSE20181 pre vs 2 wk (Affymetrix, independent platform and centre). Z1031B as a third, with the caveat of the reference design. Edinburgh Illumina series as an Illumina-only 3-month secondary (38 patients in GSE59515 + GSE55374, up to 75 with GSE111563 tp1+tp4 if disjoint: unresolved, do not count them yet).
- Coherence within a patient pair needs a cohort-level coherence per time point (tumours are the units), so the paired design is: coherence of the on-treatment set vs coherence of the pre set, with permutation of time labels within patient. Power is set by n (157 in POETIC, 56-58 in GSE20181).
- Confounds: (i) sampling method differs between baseline (core cut) and on-treatment (surgical specimen) in POETIC, hence the control; (ii) ER-signalling blockade reduces proliferation and changes tumour purity and stromal fraction, so coherence changes must be checked against purity/proliferation scores on the same permutations; (iii) responder status is partly defined by the on-treatment sample (Ki67 in Edinburgh), do not use it as a baseline stratifier; (iv) GSE20181 3-month batch; (v) two-colour series are ratios to a shared reference.
- ER-positive only: all cohorts are ER+ (except possibly a handful of HER2+ in POETIC, 44 sample entries labelled HER2-positive). The contrast is "ER blockade within ER+", which is a different question from the ER-negative vs ER-positive comparison in A127/A166; the prediction (if ER signalling supports clock coherence, coherence drops on-treatment) should be stated before scoring, with the untreated D/S control as the null.

## In vitro (listed, not ranked above patients)
- GSE104872: MCF7, 100 nM 4-OHT for 2, 6 and 24 h vs vehicle 24 h, 6 replicates (24 samples, RNA-seq GPL20301). Not dense.
- GSE76369: serum-shock entrained MCF10A (non-cancerous), 0-28 h every 4 h, 8 samples (GPL21250). Dense in time but no endocrine arm, and not MCF7 in this record.
- GSE22593: estradiol time course 1-24 h in ER-engineered MDA-MB-231 (24 samples, GPL96). Not dense.
- PRJCA000628 (BIG Data Center, not GEO): MCF7 +/- tamoxifen time series (search result only, unverified).
- No dense-time endocrine perturbation in GEO was found.

## Not found / not verified
- No tamoxifen or exemestane window-of-opportunity series with paired expression turned up in the searches done (GEO eutils search failed with the tool; two web searches and a list of candidate accessions only: not exhaustive). GSE32072 (FNA pre vs surgical post, 50 samples, GPL96) is not clearly paired (title-level read only) and GSE41994 (first-line AI, metastatic, reference-pool design) is not paired.
- ACOSOG Z1031 arm labels, POETIC drug identity, Ki67 data: not in the GEO metadata.
- Raw CEL/IDAT scan dates: not downloaded.
