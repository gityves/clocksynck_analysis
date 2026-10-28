# Data scout 2026-10-28 (roadmap Phases 6.1, 6.2, 6.6, 6.9 and 7.3)

Metadata only. Nothing was downloaded except small directory listings (GEO FTP `suppl/` and `matrix/` indexes). No code run, nothing scored.

**How this was checked**
- GEO `acc.cgi` pages return a verification wall. Verification used OmicsDI mirrors of the GEO records, GEO FTP directory listings (file names and sizes), PMC full text, and web-search abstracts.
- **[V]** = read from a record, listing or full text this session. **[S]** = from a search snippet or abstract only. **[U]** = unverified, or inferred. **[prior]** = taken from an earlier scout file (`scCircadianCancer/legacy/desync/validation/inflammation_scout_2026-10-16.md`), not re-checked.
- Already in the project (excluded as "new" below; on disk under `$SCCC_DATA`): GSE30999, GSE13355, GSE121212, GSE201827, GSE268958, GSE66511, GSE283313, GSE137218, GSE136757, GSE69967, GSE249027, GSE173706, GSE228421, GSE278330, GSE147424 (raw tar), GSE39445 (`blood_timed`), GSE104674 (`adipose_timed`), GSE108539 (`muscle_timed`), TAURUS and the IBD bulk series.
- Caveat: GSE106800 and GSE107537 are named in `provisory_plan_paper.md` (A124) but are not on disk.

## 1. Human circadian perturbations (6.1, 6.2)

| accession | citation | design | n paired subjects | tissues / compartments | times recorded | platform | processed data public | suitability for clocksync |
|---|---|---|---|---|---|---|---|---|
| **GSE106800** | Wefers et al. 2018 PNAS (PMC6065021; NCT02580513) | Randomised crossover: 3-day control protocol vs 3.5-day misalignment (12 h behavioural shift); vastus lateralis biopsies | 14 men [V] | skeletal muscle only | **Yes: 2 biopsies per condition, 8 AM and 8 PM** (Methods; the Results text says 7 AM/7 PM, an unexplained inconsistency) [V]. 4 biopsies per man, 47 arrays after QC (not 56) [V] | Affymetrix Human Gene 2.1 ST [V] | y: series matrix 13 MB, RAW.tar 226 MB (CEL) [V listing] | **High** for A124. Within-subject paired design with two timed points per condition. Single tissue, so the readout is the core-vs-cue offset within one tissue (A115 estimator), not host-referenced coupling. Only 2 time points per condition limits a per-condition cosinor. Which of the 47 arrays are missing is unknown [U]. |
| **GSE107537** | Kervezee et al. 2018 PNAS (simulated night shift, PBMC) | Two 24 h periods in the same people: baseline vs after 3 nights on a 10 h delayed sleep schedule; constant posture | 8 [V] | PBMC | **Yes: every 4 h, 7 per 24 h, per period** (up to 112 arrays) [V] | Microarray (platform not named) [U] | y: series matrix 15 MB [V listing] | **Medium**. Densely timed, within-person; blood readout is weak (A35). Kervezee's own result is that rhythmic probe sets fall from 3.8% to 0.5% [S]. Good as a descriptive secondary. |
| **GSE48113** | Archer et al. 2014 PNAS (forced desynchrony, blood) | 28 h forced-desynchrony; blood on day 1 and day 4 of the protocol | 22 [V] | leukocytes (whole blood) | **Yes: 7 samples on each of D1 and D4** [V]; constant-routine style timing [U] | Human whole-genome microarray (model not named) [U]; the accession is not stated by the paper in what I read, but the OmicsDI title matches the PNAS title [V] | y: series matrix 70 MB [V listing] | **Medium**. Largest timed within-person misalignment series (22 x 2 x 7). Blood only. Whether the PNAS design is forced desynchrony on day 4 or a later cycle is unclear in the record [U]. Already named in the 2026-10-15 scout as a blood series. |
| GSE39445 (on disk) | Möller-Levet et al. 2013 PNAS (sleep restriction) | Cross-balanced: 1 week sufficient (8.5 h) vs 1 week insufficient (5.7 h), then total sleep deprivation, 10 samples per period | 26 [S] | leukocytes | Yes (10 per condition) | Human whole-genome microarray [U] | y (on disk) | **Medium**. On disk. Within-person, timed. Blood only. |
| GSE87761 | Christou et al. 2019 Sci Rep (PMC6389935), human white adipose clock | 5 subcutaneous adipose biopsies per man, 6 h apart (constant routine, 5 h after start) | 7 men, no perturbation [V] | subcutaneous adipose | Relative times (clock times not listed) [V] | Agilent 4x44k custom (GPL15331) [V] | y: series matrix 7.5 MB [V listing] | **Low** for 6.1 (no perturbation). Useful only as a timed human adipose truth. Not paired with another tissue. |
| GSE122541 | Nurses on day vs night shift, PBMC (reanalysed in a 2026 Front Public Health paper) | Between-person day vs night shift | not paired [S] | PBMC | Not stated [U] | not checked | not checked | **Low**. Between-person, so no within-donor control. |
| Cedernaes 2015 JCEM (no accession found) | Sleep vs overnight wakefulness, 15 men | Crossover; muscle and adipose biopsies from the same men | 15 men [S] | skeletal muscle and subcutaneous adipose | Single morning biopsy | targeted qPCR of clock genes, DNA methylation [S] | **n** (no transcriptome matrix found) | **Low**: the only same-person two-tissue sleep-loss design found, but no public genome-wide matrix. A request to the authors would be outward-facing: ask the user first. |
| NCT05180279 | Colonic biopsies at aligned and misaligned circadian states in UC (ClinicalTrials.gov) | Registered, single tissue | not stated | colon | not stated | RNA-seq planned [S] | **n** (no results posted) | **Low** (not available). |

**Not found.**
- No public multi-tissue human misalignment, shift-work or jet-lag series (same donor, at least 2 tissues) after two extended searches. Human multi-tissue data stay at GTEx level.
- No public timed glucocorticoid-administration transcriptome in human tissue or blood. The nearest are the DREAM trial ancillary study (modified-release hydrocortisone, clock genes in immune cells; qPCR, no accession found) and a PBMC 8 AM/8 PM study of CLOCK and the glucocorticoid receptor (Charmandari-type design; qPCR) [S]. Both are low.
- The Wefers adipose and the Skene/Kervezee metabolomics are not transcriptomic.

## 2. Within-person disease designs (6.6)

IBD was scouted in depth on 2026-10-16 (full table in the prior file). Only the paired-design essentials are repeated; none of it was re-verified here. Only new findings are added below the line.

| accession | citation | design | n paired subjects | tissues / compartments | times recorded | platform | processed data public | suitability for clocksync |
|---|---|---|---|---|---|---|---|---|
| TAURUS: Zenodo 10.5281/zenodo.13768607; GSE282122 | Thomas et al. 2024 Nat Immunol | UC/CD, multi-site, pre/post adalimumab, 3 healthy | **28** patients with inflamed + non-inflamed (15 CD, 13 UC); 18 at baseline; 37 same-site pairs in 20 patients [prior V] | whole biopsy: endothelial, fibroblast, mural (median ~190/191/120 per sample), epithelium, immune [prior V] | no | 10x 3' v3 [prior] | y: raw counts h5ad 12.7 GB, per-compartment files 0.04-6 GB [prior V]; on disk as `ibd/TAURUS` | **High**. Already on disk. Site confounded with inflammation in some pairs (use same-site pairs or site covariate). |
| SCP259 | Smillie et al. 2019 Cell | 18 UC, 12 healthy, Epi + LP fractions of the same biopsy | **~16-18** [prior P/U] | epithelium, immune, stroma (~400 stromal cells per biopsy) | no | 10x v1/v2 | y with free Single Cell Portal login (matrices 0.5-2.4 GB) [prior V listing] | **High**. Best matched design (adjacent biopsies). Not on disk. Needs SCP login: ask the user. |
| SCP3755 | Eshghi et al. 2026 Front Immunol (10.3389/fimmu.2026.1705328; bioRxiv 2024.07.02.601267) | 20 UC, no healthy; inflamed vs uninflamed segments | **18 matched pairs** [prior P] | Epi + LP fractions; fibroblast, endothelial, pericyte, SMC | no | BD Rhapsody WTA (not 10x) [prior P] | y with SCP login (Seurat RDS 1.86 GB) [prior V] | **Medium-High**. Separate chemistry stratum. The accession could not be re-confirmed in this session [U]. |
| GSE266616 | CD, TI + ascending colon (paper not identified) | inflamed, adjacent, non-inflamed, control | **14 CD** [prior V] | TI, colon | no | 10x | y: RAW.tar 3.7 GB, no annotation deposited [prior V] | **Medium**: needs own annotation. |
| SCP1884; CELLxGENE | Kong et al. 2023 Immunity | 46 CD, 25 non-IBD | few: 5 TI, 3 colon in the CELLxGENE subset [prior V] | Epi + LP | no | 10x | y (SCP login) | **Low-Medium** (mostly between-patient). Already examined as A16 (not replicated). |
| GSE134809 | Martin et al. 2019 Cell | 11 ileal CD resections, involved vs uninvolved | **11** [prior V] | ileum LP cells, stroma sparse | no | 10x v1/v2 | y | **Low-Medium**: stroma marginal. |
| E-MTAB-16999/16986, GSE266546 | IBDverse; Vanderbilt LND | large | no pairs | stroma depleted (about 1 anchor cell per sample) | no | 10x; inDrop | y | **Low**: not usable for the anchor design [prior V]. |
| **GSE173706** (on disk) | Reynolds-type psoriasis atlas, Nat Commun 2023 (keratinocytes and fibroblasts amplify inflammation) | Healthy (NS), clinically normal skin of patients (PN), lesional (PP) | **11 matched PN-PP pairs**; 33 libraries: 8 NS, 11 PN, 14 PP; 14 patients, 11 healthy donors [S] | whole skin: keratinocytes, fibroblasts, immune | no | 10x | y (on disk) | **High** for 6.9 and 6.6. Already used. Keratinocyte-only pseudobulk is feasible (A61/A63 used this series). |
| GSE228421 (on disk) | risankizumab single-cell response | lesional + non-lesional, day 3 and 14 | 5 individuals [S] | skin | day 3 / 14 | 10x | y (on disk) | **Low**: n = 5. |
| **GSE147424** (raw tar on disk, 305 MB) | Reynolds/He-type AD atlas, JACI 2020 | AD lesional + non-lesional vs healthy | 5 AD patients (4 lesional, 5 non-lesional libraries), 7 healthy; 39,042 cells [S] | skin: fibroblast, keratinocyte, immune | no | 10x | y | **Low-Medium**: at most 4 pairs. Whether the patients overlap between lesional and non-lesional is not stated [U]. Keratinocyte counts per library unknown. |
| **Nat Commun 2026, AD atlas (s41467-026-69587-7, "multi-modal skin atlas"; possibly also s41467-026-75407-9)** | Not fetched (paywall redirect) | 5' scRNA-seq, 116 biopsies: 54 lesional, 39 non-lesional from 42 AD patients, 23 healthy [S] | **37 of 42 patients have paired NAD + LAD** [S] | skin | no | 10x 5' | accession **unverified [U]** | **High if public** (largest paired skin design found; healthy controls too). Which of the two 2026 papers carries these numbers is unclear [U]. Find the accession in the Data availability section. |
| GSE249027 (on disk) | HS lesional vs non-lesional punch biopsies | bulk RNA-seq | not stated [S] | skin | no | RNA-seq | y | **Medium** (bulk). Already on disk. |
| HS-OmicsDB (bioRxiv 2023.01.26.524600) | meta-resource | 71 L, 19 perilesional, 29 NL, 11 healthy (bulk and sc, mixed studies) [S] | not stated | skin | no | mixed | y (curated) | **Low**: curated across studies, not one paired cohort. |
| Cutaneous / SLE skin scRNA-seq | not found | | | | | | | **Not found** in two searches. |
| GSE38713, E-GEOD-6731 | UC bulk paired mucosa | 46-48 paired inflamed / non-inflamed sigmoid (E-GEOD-6731) [S] | ~46 | colon (bulk, one compartment) | no | arrays | y | **Low**: bulk, one tissue, no host. (The project already uses ibd_bulk series for A55.) |

## 3. New psoriasis bulk cohorts and placebo-controlled treatment series (6.9, future A116)

Excluded as already used: GSE30999, GSE13355, GSE66511, GSE121212, GSE201827, GSE268958, GSE283313 (roflumilast PSORRO sub-study), GSE137218, GSE136757, GSE69967.

| accession | citation | design | n paired subjects | tissues / compartments | times recorded | platform | processed data public | suitability for clocksync |
|---|---|---|---|---|---|---|---|---|
| **GSE31652** | Krueger et al. 2012 JACI (anti-IL-17A ixekizumab precursor LY2439821) | Phase 1, randomised, double-blind, **placebo-controlled**, dose-escalation (5/15/50/150 mg + placebo, 8 per arm) | 40 [V] | lesional skin, same site | Biopsies at weeks 0, 2, 6 [V] | Whole-genome microarray (Affymetrix, model not named) [U] | y (GEO series) [U size] | **High** for A116: true placebo arm, repeated lesional biopsies, drug that acts through IL-17 (the psoriatic delay, A54). Only 8 per arm; lesional only (no non-lesional paired). Without sampling clock time the delay change must be read as an offset against the host or core-cue sets. |
| **GSE51440** | Sofen et al. 2014 JACI (guselkumab first-in-human) | Phase 1, randomised, double-blind, placebo-controlled | 59 biopsies (lesional and non-lesional; baseline, week 1, week 12) [S] | skin | Weeks 0, 1, 12 | microarray | y [U size] | **High-Medium**. Placebo arm and non-lesional skin from the same patient. Platform and arm sizes unverified [U]. |
| **GSE53552** | Russell et al. 2014 (brodalumab, lesional and non-lesional) | Dose-ranging brodalumab; lesional and non-lesional skin over time | 25 patients, 99 samples [S] | skin | over time | Affymetrix [U] | y | **Medium**. Placebo arm not confirmed in the record [U]. |
| **GSE287022** | ObePso-S adipose study (2026): cutaneous adipose tissue in psoriasis, partly IL-17 dependent | Randomised double-blind placebo-controlled phase 4 (secukinumab vs placebo, 12 wk) | 82 patients + 15 healthy [V] | **cutaneous adipose tissue** (lesional and non-lesional) | 12 wk; biopsy timing [U] | RNA-seq (platform not named) | y: raw counts per gene, 10 MB [V listing] | **Medium-High**. Adds a second compartment (adipose) beside the skin readout. The bulk skin of the same trial is under PRJNA833143 (ENA, raw only [S]). No sampling times. Pairing of adipose with skin per patient is unconfirmed [U]. |
| **GSE328285** | SIK2/3 inhibitor GLPG3970 trial; psoriasis (2026) | Bulk RNA-seq: ex vivo treated lesional biopsies; trial responders vs placebo | not stated; count matrix 1.1 MB [V listing] | skin | not stated | RNA-seq | y | **Low-Medium**: small, design not described in the record [U]. Recent (2026). |
| PRJNA1193084 | Roflumilast vs placebo, skin RNA-seq | Phase-2 style PSORRO study | not stated | skin | | RNA-seq | ENA project only; likely the same samples as GSE283313 (already used) [U] | **Low** (probable duplicate). |
| phs001688 | Secukinumab early skin biopsies | baseline to wk 12; 36 patients (24 treated, 12 placebo) [S] | 36 | skin | | | dbGaP controlled access | **Medium** but controlled access: outward-facing, ask the user. |
| GSE47751 / GSE11903 / E-MTAB-6556 | Etanercept early response (20 pts, days 1-14); etanercept 15 pts (wk 0,1,2,4,12) | no placebo | 10-20 | lesional (+ blood in E-MTAB-6556) | days / weeks | microarray, RNA-seq | y | **Low-Medium**: no placebo; useful as time-course controls. |
| E-MTAB-5262; GSE50614; E-GEOD-57376 | Adalimumab early (10 pts, days 0, 4, 14); SRT2104 (11 pts pre/post); infliximab (5 pts) | no placebo documented [U] | 5-11 | skin | | | y | **Low**: small. |
| GSE41662 | Lesional vs non-lesional, 24 pts | no treatment | 24 | skin | | Affymetrix | y | **Low-Medium**: another paired psoriasis cohort for the frozen pipeline. No healthy skin. |
| GSE14905 | Yao et al.; lesional, non-lesional, normal | | 54 subjects, 82 arrays [S] | skin | | | y | **Medium**: has healthy controls and pairs. Status in project unclear: not on disk [U]. |
| GSE54456 | Li et al. 2014 RNA-seq (92 psoriatic, 82 normal) | case-control | not stated | skin | | RNA-seq | y | **Medium**: large, healthy controls (are the 92 lesional/uninvolved paired? [U]). |

## 4. Tumour transcriptomes with recorded surgery time of day (7.3)

| accession | citation | design | n paired subjects | tissues / compartments | times recorded | platform | processed data public | suitability for clocksync |
|---|---|---|---|---|---|---|---|---|
| GSE233242 (the A146 cohort) | Li et al. 2024 PNAS (tumour circadian clock strength, luminal A) | 43 tumour + matched normal pairs (29 LumA, 3 LumB, 2 HER2, 9 TNBC) [V]; the project holds 42 | 43 pairs [V] | breast tumour and adjacent normal, bulk | resection times recorded for all (08-18 h in the project's join) [V] | bulk RNA-seq (platform not named) | y: processed TPM 4.3 MB [V listing] | **This is the first cohort, not a replication.** |
| none | | | | | | | | **No second public breast (or other) tumour cohort with recorded surgery time was found** after three searches (breast, colorectal, lung, liver, glioma, head and neck). Most "circadian signature" papers use TCGA/GEO with no sampling time. |

Leads that are not cohorts, for the PI to judge:
- The PNAS paper's TCGA/GTEx use infers order without times, so it cannot replicate A146.
- Oral mucosa (11 volunteers, every 4 h, RNA-seq) is a timed normal tissue, not a tumour [S].
- Mouse 4T1 at CT4/10/16/22 (bioRxiv 2026.05.26.726543) is the source of the A146 predictions, not a test.
- A replication would need either a request to a clinical group that recorded theatre times (outward-facing: ask the user) or surgical-list timestamps in a biobank. TCGA clinical tables do not carry them [U].

## Reading notes
- Wefers has only 2 timed points per condition per man, in a single tissue, so the clocksync readout is an offset change, not coupling. That matches the A124 design (core-cue offset, paired).
- GSE48113, GSE107537 and GSE39445 are blood-only, but their dense timing allows a within-subject cosinor readout alongside the untimed one. This could validate that clocksync reads misalignment on blood in a design where truth is known.
- For 6.9, the prior single-cell series (GSE173706, 11 PN-PP pairs, 8 healthy) already supports a keratinocyte-only pseudobulk. The 2026 AD atlas, if public, would be the largest paired skin design and is the one accession to locate.

## Status and not yet scouted (scout paused 2026-10-28)

Complete enough to act on: section 1 (Wefers, Kervezee, Archer verified), section 4 (negative result).
Partial: section 2 (IBD relies on the 2026-10-16 scout and was not re-verified; the 2026 AD atlas accession is open; HS and lupus single-cell not found), section 3 (arm sizes, placebo arms and platforms unverified for GSE51440, GSE53552, GSE328285).

Not yet scouted:
- Accession and data availability of the Nat Commun 2026 AD atlas (s41467-026-69587-7 / s41467-026-75407-9).
- Whether GSE287022 adipose samples pair with skin per patient; biopsy times in ObePso.
- Sample tables (patient IDs, array exclusions) of GSE106800, GSE48113, GSE107537 (series-matrix headers not read).
- Placebo arms of GSE53552 and GSE51440; platforms of GSE31652.
- Single-cell hidradenitis and cutaneous lupus; further human perturbation series (jet lag, glucocorticoid with public matrices).
- GSE14905, GSE54456 and GSE41662 pairing and status in the project.
- Eshghi SCP3755 accession; a second timed tumour cohort in non-breast cancers.
