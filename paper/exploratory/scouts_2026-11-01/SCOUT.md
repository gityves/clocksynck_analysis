# Metadata-only scout, 2026-11-01 (HANDOFF Priority 2: H1, H3, H4, H5; Priority 3: D1, D2, D3)

**Method.** NCBI E-utilities (`db=gds`) searches, GEO accession text pages (series and sample metadata) and series-matrix
headers; file sizes from HTTP HEAD requests. No expression value was read, nothing was downloaded except these metadata
pages, nothing was committed, `paper/docs/` was not edited. Helper scripts and raw metadata are in the session scratchpad
(`.../scratchpad/tools`, `.../scratchpad/scouts`). Persons with >= 2 compartments were counted from sample metadata
(patient or subject id x tissue or state); where the id had to be parsed from titles this is stated.
**Hard rule applied:** a design qualifies only if each person contributes >= 2 compartments from the same procedure or
time. Same tissue at two time points (pre/post ischaemia, serial biopsies) is NOT a second compartment; such series are
listed as excluded.
**Already used (grep of `paper/docs/*.md` and the archive `docs/`, plus `$SCCC_DATA` directory names).** "y" below means
the accession appears in the docs or has a data directory under `$SCCC_DATA`.
**Collection time:** "no" means no clock time, date or time-of-day field in the series or sample metadata (checked the
`!Sample_characteristics` keys); the paper may still contain it.

---

## H1. Paired within-person biopsies (IBD, psoriasis, tumour-normal)

### Already used by the project (do not re-scout)
- **IBD (A55, `legacy/desync/validation/paired_damping/ibd_scout.md`):** 55 series screened; 6 built, 563 gated pairs:
  GSE193677 (454 pairs, MSCCR), GSE179285, GSE66407, GSE95095 and GSE208303 (both failed the normal gate), GSE107499.
  About 25 more were excluded (no patient id, one biopsy per donor, different patients per group, < 10 pairs). Same-segment
  pairs are the limit (31 of 454 in GSE193677). **Ileum vs colon from the same patient is not a recorded design in any of
  them**: not found in GEO bulk. Single-cell IBD atlases (Smillie SCP259, Garrido-Trigo GSE214695) were handled by
  A153/A17 and are underpowered.
- **Psoriasis and skin (A115/A116/A148/A160, `$SCCC_DATA/skin_lesional/`, 43 series on disk):** GSE13355, GSE30999, GSE121212,
  GSE201827, GSE117468, GSE69967, GSE137218, GSE193128, GSE66511, GSE153007, GSE75343, GSE283313, GSE268958, GSE197056,
  GSE278330, GSE252029 (psoriasis, lesional/non-lesional); GSE173706 and GSE228421 (single-cell, 12 patients with both
  families); atopic dermatitis sets GSE133477, GSE130588, GSE157194, GSE140684 and others on disk. About 130 gated pairs.
- **Tumour-normal:** TCGA 517 pairs, 10 gated cancers (H138.4); A171 GSE14520 (213 pairs) and GSE36376; GSE31210
  (descriptive); A177/A175 stopped GSE14520 and GSE66229 (MDE above TCGA D); GSE62254, GSE15459, GSE84437 and others on
  disk (scout_gastric_lung); GSE236696, GSE183904, GSE271690 (single-cell, 6 to 65 pairs).

### New candidates (none appears in any doc or under `$SCCC_DATA`)
| accession | design | n persons >= 2 comp. | compartments | platform | collection time recorded | processed download | used |
|---|---|---|---|---|---|---|---|
| **GSE154200** (superset GSE308055: 381 patients, 21 skin diseases, `.h5` / `.rds`, size not exposed) | one centre (Munich), skin punch biopsy from lesional and non-lesional skin of the same patient, same sitting | **186** (psoriasis 88, nummular eczema 52, atopic eczema 46) | lesional, non-lesional skin | bulk RNA-seq, count matrix | no | 17.1 MB `countMat.txt.gz` | n |
| GSE117239 | psoriasis, ustekinumab vs etanercept RCT; lesional and non-lesional at baseline, weeks 1 and 12 | 83 (baseline pairs <= 83) | lesional, non-lesional | Affymetrix GPL570 | visit schedule only | 58.9 MB series matrix | n |
| GSE193309 | AD, 1.5 mm punches, two lesional and two non-lesional sites per patient plus 30 matched healthy | 28-30 AD (+30 HC) | lesional, non-lesional | bulk RNA-seq | visit date recorded | not measured | n |
| **GSE93601 + GSE115577** (Nurses' Health Study I/II, breast) | invasive breast cancer, FFPE tumour and tumour-adjacent normal, ER, PAM50 and proliferation score given | 445 (GSE93601) and 328 (GSE115577); possible overlap between the two deposits (unresolved), so union about 450-770 | tumour, adjacent normal | Illumina DASL array (FFPE), clock-gene coverage unverified | no (diagnosis year only) | 92 MB (GSE93601 matrix); 118 MB data csv (GSE115577) | n |
| GSE87211 | rectal cancer, tumour and matched mucosa control (pre-operative radiochemotherapy given) | 160 | tumour, mucosa | Agilent | no | 35 MB series matrix | n |
| GSE29272 | gastric cancer, cardia and non-cardia, tumour and adjacent normal glands | 134 | tumour, adjacent normal | Affymetrix U133A | no | 27.7 MB series matrix | n |
| GSE183795 | pancreatic ductal adenocarcinoma, tumour and adjacent non-tumour | about 100 (139 T, 102 adjacent; pairing key not verified) | tumour, adjacent | Affymetrix GPL6244 | no | 19.6 MB normalised matrix | n |
| CPTAC-3 (LUAD, CCRCC, LSCC, UCEC, HNSCC, COAD, BRCA) | tumour + normal adjacent, about 100 per cohort | about 100 each (from memory, **not verified here**: outside GEO, GDC/PDC) | tumour, NAT | RNA-seq | procurement times recorded in some clinical tables (unverified) | not measured | n |

**Notes.** (i) GSE154200 uses three diseases in one pipeline, so lesional-vs-non-lesional can be compared across diseases at
a common spread d*. The 88 psoriasis pairs add to about 130 already used. (ii) NHS breast is the only new tumour-normal series
at >= 150 pairs with ER and a proliferation score in the metadata, which are the C1 stratifiers (ER-, proliferation tertile);
its risk is FFPE DASL data (probe set and clock-gene coverage must be checked on the platform annotation before any
expression is opened). (iii) Excluded: GSE229705/GSE229706 (tumour-adjacent lung only), GSE276755/GSE276760 (breast tissue,
pairing not verified).

**Rough MDE on d* (scaling only).** A153 anchor: single-cell atlases with n about 20-40 had MDE >= 0.214 (smallest in
`A153_corrected_contrasts/output/mde_table.csv`), scaled as 1/sqrt(n). Taking 0.214 at n = 40 (optimistic end of the range):

| candidate | n pairs | rough MDE on d* |
|---|---|---|
| GSE154200, all three diseases pooled | 186 | 0.10 |
| GSE154200, psoriasis only | 88 | 0.14 |
| GSE154200 + the 130 psoriasis pairs already used | about 218 | 0.09 |
| NHS breast, all pairs | about 450 | 0.06 |
| NHS breast, ER- stratum (about 25 percent of tumours) | about 110 | 0.13 |

This is a rough scaling, not a power calculation: the anchor comes from single-cell compartments with many cells and a
different noise structure, the paired bulk design has only two compartments (one L2 pair), and MDE also depends on the
baseline coupling c0 and the host spread. The registered MDE step (`mde_table.py` logic) must be run per series before any
contrast. On this scaling the HANDOFF target (MDE <= 0.15) is reached from about 80 pairs.

**Recommendation H1: feasible.** Top candidate GSE154200 (186 pairs, clean counts, 17 MB, all new). NHS breast
(GSE93601 + GSE115577, about 450 pairs) is feasible as a second, larger C1 replication of the ER and proliferation
strata, conditional on the DASL clock-gene check. IBD and psoriasis-in-GEO are otherwise exhausted by A55 and A115.

---

## H3. Ischaemia time, surgery duration or anaesthesia with >= 2 compartments per person

| accession | design | n persons >= 2 comp. | compartments | platform | collection time recorded | download | used | verdict |
|---|---|---|---|---|---|---|---|---|
| GSE274774 | ccRCC, normal kidney and tumour sampled in vivo and after resection (warm ischaemia), 1-3 regions | about 8 patients (42 samples; id parsed from titles) | kidney, tumour, each at 2 ischaemia states | bulk RNA-seq (GPL18573) | ischaemia state (in vivo / ischaemic) only | 8.8 MB RAW.tar | n | qualifies on design, far too small |
| GSE320593 | liver (26 pts) and subcutaneous adipose (60 pts) biopsies, non-surgical vs surgical, with or without perioperative dexamethasone | 0: each patient has one tissue (ids A1..., L1... disjoint), two states | one tissue per person | bulk RNA-seq counts | DEX yes/no, surgical vs not | 5 MB | n | **excluded (single compartment)**; closest design to the dexamethasone question, usable only for L1 |
| GSE260643 | aortic valve surgery, left-ventricle biopsy pre- and post-ischaemia, **AM vs PM surgery** (116 AM, 39 PM samples) | 0: one tissue, two time points | left ventricle only | bulk RNA-seq | surgery AM/PM recorded | only DEG tables (1.6 MB), counts not on GEO | n | **excluded**; note the recorded AM/PM time of day, a possible truth set for L1 only |
| GSE273225 | lung transplantation, biopsies at start and end of implantation (rewarming ischaemia time recorded), left lingula and right middle lobe | up to 34 recipients, left and right lung | two lungs of one donor | NanoString Immunology V2 (579 genes) | rewarming time recorded | RAW.tar | n | **not feasible**: clock genes absent from the panel |
| GSE268041 | lung, cold ischaemia vs 2 h after reperfusion | single tissue | lung | Affymetrix | ischaemia stage | 12 MB | n | excluded |
| GTEx | ischaemic time, many tissues per donor | already used | | | | | y | the existing H3 substitute |
| CPTAC-3 (LUAD etc.) | tumour + NAT with procurement and ischaemia times in clinical tables | about 100 per cohort | tumour, NAT | RNA-seq | partly | not measured | n | **unverified**: outside GEO; needs PDC/GDC clinical table check before any decision |

**Recommendation H3: marginal / not feasible from GEO.** No GEO cohort combines ischaemia, anaesthesia or surgery-duration
metadata with >= 2 compartments per person at n >= 20. GSE320593 and GSE260643 are good single-compartment designs
(dexamethasone; AM vs PM surgery) but fail the rule. The only route is CPTAC (tumour + NAT with procurement times); verify the
clinical table fields first (an unmeasured scout step).

---

## H4. Transplant serial biopsies with a recipient reference compartment

| accession | design | n persons >= 2 comp. | compartments | platform | collection time recorded | download | used | note |
|---|---|---|---|---|---|---|---|---|
| **GSE295373** | deceased-donor kidney: zero-hour procurement biopsy of left and right kidney (two grafts of one donor), 3-month protocol biopsy of recipients | 52 donors with both kidneys at time zero; 11 donors with both recipients at 3 months; 17 more with one 3-month biopsy | left kidney, right kidney (same donor, two recipients) | Affymetrix PrimeView | no (3 months, zero hour) | 48 MB series matrix | n | no recipient reference compartment, but the two kidneys of one donor go to two different recipients: shared donor state at time zero, then two hosts. Same-organ compartments, so weak for a clock gate |
| GSE249451 (CTOT-08) | kidney biopsy and matched blood per visit, for-cause or surveillance | 129 patients, 283 blood-biopsy pairs | graft kidney, recipient blood | array (GPL33988) | no | 51.6 MB series matrix | n | blood reference only (archive A35: weak blood clock, rho 0.31); no time-zero biopsy, so no donor-state arm |
| GSE242302 | acute or chronic rejection: biopsy, blood and urine sediment | 20 samples, 7 patients per tissue (pairing not given) | biopsy, blood, urine | RNA-seq | no | 1.6 MB counts | n | too small |
| GSE43974 | kidney biopsies at donation, after cold ischaemia and 60 min after reperfusion (554 biopsies) | 0 | one tissue, 3 time points | Illumina HT-12 | time point | 554 samples | n | excluded (single compartment); no recipient tissue |
| GSE263614 / GSE276531 / GSE151648 | liver normothermic perfusion or ischaemia-reperfusion | single tissue | liver | RNA-seq | perfusion hours | 4.6 MB (GSE263614) | n | excluded |
| GSE293039, GSE256141, GSE277334, GSE193135 | liver transplant biopsies, rejection phenotypes | cross-sectional, one tissue | liver | various | no | | n | excluded |

**Recommendation H4: not feasible as specified; marginal for one variant.** No series has time-zero graft, later biopsy and
a second recipient solid tissue. GSE295373 is the only series with time-zero biopsy and later protocol biopsy, and its
two-kidney, two-recipient structure supports a different prediction (graft pairs decouple after transplant), but only 11 donors
have both recipients and both compartments are the same organ. GSE249451 has the n (129) but a blood reference.

---

## H5. ICU and sepsis cohorts with > 1 compartment per person

| accession | design | n persons >= 2 comp. | compartments | platform | collection time recorded | download | used | note |
|---|---|---|---|---|---|---|---|---|
| GSE237861 | autopsy: 7 septic-shock deaths and 7 non-infectious critically ill deaths, 6 tissues each | 14 (82 samples) | cortex, hippocampus, heart, lung, kidney, colon | RNA-seq | no | **no count matrix on GEO** (per-patient DE tables only, 0.4-0.6 MB each; counts only via SRA) | n | right design, n = 14, and processed counts missing: would require FASTQ mapping, which the project excludes |
| GSE198103 | CRC patients who died of sepsis vs non-sepsis death: FFPE tumour and adjacent normal | about 25 pairs of 30 | tumour, normal | array (GPL23126) | no | small | n | sepsis is a mixed label in cancer patients; marginal |
| GSE290679 | CD4 T-cell subsets from ICU patients with and without sepsis, 48-72 h after admission | PBMC only (single compartment) | PBMC | RNA-seq | 48-72 h after ICU admission | | n | raw data under dbGaP; excluded |
| GSE302306, GSE300129, GSE196117 and about 20 further sepsis or COVID-ICU series | whole blood or neutrophils | 1 | blood | | | | n | L1 only |
| COVID autopsy sets (Melms, Delorey, others) | multi-organ nuclei, fatal COVID | already analysed | | | | | y (A105, A105b) | existing ICU-like result |

**Recommendation H5: not feasible.** The only multi-tissue ICU series (GSE237861, n = 14) has no processed counts on GEO;
other sepsis series are single-compartment blood. Keep A105/A105b (fatal COVID, ventilated) as the ICU evidence.

---

## D1. Descartes fetal atlas (GSE156793)

**Already used and on disk (A107b, `$SCCC_DATA/natexp/cao2020/`).** 28 fetuses, 72-129 days post-conception, 15 organs;
`pseudobulk_organ/` (112 organ x fetus rows) and `pseudobulk_class/` (775 rows) with the 68-gene panel exist (h5ad 0.09 and
0.34 MB), built 2026-10-19 from the loom; the decompressed loom (21.4 GB) was deleted, the raw gz (13.0 GB) and cell
metadata (0.35 GB) remain. A107b scored 4 gates (no Holm survivor). So D1's "scout first" is done; the new work is the depth
control and the registered negative prediction, not a download.

| GEO supplementary file | size | use |
|---|---|---|
| `GSE156793_S3_gene_count.loom.gz` | 13,007,903,046 B (13.0 GB) | only per-cell counts; required for any pseudobulk with a new gene set |
| `GSE156793_S1_metadata_cells.txt.gz` | 349 MB | cell to fetus, organ, cell type |
| `GSE156793_S4_gene_expression_tissue.txt.gz` | 5.7 MB | tissue-level means, **not per fetus**: not usable |
| `GSE156793_S2_Metadata_genes.txt.gz` | 0.58 MB | gene table |

**Fetuses with >= 3 organs** (from local `samples.csv`, organs with n_cells >= threshold; metadata only):

| minimum cells per fetus x organ | fetuses | >= 2 organs | >= 3 organs | >= 4 organs |
|---|---|---|---|---|
| any | 28 | 23 | 18 | 15 |
| 200 | 28 | 22 | 18 | 12 |
| 1,000 | 28 | 22 | 17 | 10 |
| 5,000 | 26 | 17 | 9 | 5 |

Mean depth about 1,057 counts per cell (sci-RNA-seq3, 93.5 percent nuclei), so a 1,000-cell pseudobulk has about 1 M counts.
Platform: single-nucleus combinatorial indexing; collection time: none (gestational day only).

**Recommendation D1: feasible, no new download.** Use the existing pseudobulks (17 fetuses with >= 3 organs at 1,000
cells); re-streaming the 13 GB gz is needed only if genes outside the 68-gene panel are wanted. The negative control still
needs the depth-matched adult downsampling before the gate is read.

---

## D2. Mother and fetus at term (same delivery)

Archive A106 used first-trimester decidua/placenta from the Arutyunyan single-cell atlas (16 pregnancies), A107/A107b fetal
organs. **No term bulk series was used.**

| accession | design | n persons >= 2 comp. | compartments | platform | labour vs caesarean | collection time | download | used |
|---|---|---|---|---|---|---|---|---|
| **GSE27272** | 91 women, Ceske Budejovice, term delivery; smoker status | **67** (47 with all three) | maternal peripheral blood, cord blood, term placenta | Illumina HumanRef-8 v3 (GPL6883) | mode of delivery: all vaginal (no planned caesarean arm) | no | 24 MB (3 txt.gz: 8.3, 8.5, 6.9) | n |
| GSE73685 (NICHD) | 35 women, 8 tissues, outcomes: term and preterm, with and without labour, pPROM | 28 (>= 2 tissues); term: no labour 7, labour 6 | amnion, chorion, decidua, lower segment, fundus, cord blood, maternal blood | Affymetrix GPL6244 | **labour vs no labour recorded** (six outcome classes) | no | 20 MB series matrix | n |
| GSE30032 | 60 women (Prague), passive smoke | about 52 | cord blood, placenta (both fetal) | Illumina GPL6883 | vaginal 85, caesarean 19 samples | no | 13.5 MB | n |
| GSE185557 | SARS-CoV-2 vs control, peripheral and cord blood | about 19 pairs | maternal blood, cord blood | DNBseq RNA-seq | labour yes (26) / no (12), route recorded | no | 1.1 MB counts | n |
| GSE334075 | maternal BMI, term | at most 14 trios (17 PBMC, 14 placenta, 18 cord) | maternal PBMC, placental villi, cord blood mononuclear cells | RNA-seq | term, no labour field | no | 1.5 MB counts | n |
| GSE181750 | RNA-seq maternal blood, cord blood, placenta | 60 mothers: maternal blood was drawn in the **first trimester**, not at delivery | | RNA-seq | no | no | 7.8 MB | n |

**Recommendation D2: marginal.** GSE27272 is the only adequately sized series (67 persons, 47 trios, fetal placenta and cord
blood plus maternal blood), all vaginal deliveries, so labour is not a factor, and blood and placenta carry weak or untested
clock signal (gate must pass first; archive A35: blood rho 0.31). GSE73685 has the labour factor but 13 term women. Both are
small downloads (about 24 and 20 MB).

---

## D3. Childhood maturation from routine surgical tissue

| accession | design | n persons >= 2 comp. | compartments | platform | age span | collection time | download | used |
|---|---|---|---|---|---|---|---|---|
| GSE126124 | children recruited at diagnostic colonoscopy (CD 39, UC 18, IBDU 2, non-IBD 39); colon biopsy and whole blood | up to 78 colon samples vs 98 blood: **pairing key not in the metadata** (titles carry only disease, age, sex) | colon biopsy, peripheral blood | Affymetrix GPL6244 | 8-18 y | no | 19.7 MB series matrix | n |
| E-MTAB-8901 (Elmentaite, gut atlas, paediatric arm) | healthy children 4-12 y (n = 8) and Crohn's, terminal ileum biopsy, plus fetal gut 6-10 pcw | 8 healthy | cell types within one ileal biopsy (single-cell) | scRNA-seq | 4-12 y plus fetal | no | not measured (ArrayExpress) | n |
| E-MTAB-9543 / E-MTAB-9536 | the same atlas, adult and fetal gut | adult, fetal | up to 10 gut sites | scRNA-seq | fetal, adult | no | not measured | n |
| GSE181022 | newborn vs adult skin biopsies | 1 site per person | skin | array | newborn, adult | no | | n |
| tonsil atlas (Massoni-Badosa 2024, Zenodo / HCA) | tonsillectomy, about 17 donors, children | cell types within one tonsil | single-cell | children | none | | not measured | n (already noted as "donor n unverified" in the docs) |
| GSE328314 | tonsil eQTL, children, 418 samples | single tissue | tonsil | RNA-seq | children | no | | n |
| GSE263248, GSE240858 | gamma-delta T cells across tissues and the lifespan | sorted T cells from organ donors | sorted cells | | lifespan | no | | n |

**Recommendation D3: not feasible.** No series spans infancy to adolescence with >= 2 compartments per child. The nearest
(GSE126124, 8-18 y, colon + blood) lacks a pairing key and an infancy end; the paediatric gut and tonsil atlases are
single-site, single-cell and n <= 17. The gap between the Cao fetal data (to 18 weeks) and GTEx (from 20 y) cannot be bridged
with public data found here; it should be stated as a limit.

---

## Summary

| item | verdict | top candidate | reason |
|---|---|---|---|
| H1 | feasible | GSE154200 (186 lesional / non-lesional pairs, 88 psoriasis, 17 MB counts); second: NHS breast GSE93601 + GSE115577 (about 450 pairs) | new, large, MDE about 0.10 by rough scaling; IBD and old psoriasis already exhausted |
| H3 | marginal | CPTAC tumour + NAT (unverified, outside GEO) | no GEO cohort with ischaemia metadata and >= 2 compartments; GSE320593 and GSE260643 are single-compartment |
| H4 | not feasible (marginal variant) | GSE295373 (11 donors with two grafts in two recipients) | no recipient solid-tissue reference; GSE249451 has blood only |
| H5 | not feasible | GSE237861 (14 patients x 6 tissues) | no processed counts on GEO; rest are blood |
| D1 | feasible, no download | `$SCCC_DATA/natexp/cao2020/` (A107b pseudobulks; 17 fetuses with >= 3 organs at 1,000 cells) | already built; needs depth control |
| D2 | marginal | GSE27272 (67 persons, 47 trios) | all vaginal, blood and placenta weak clock; labour factor only in GSE73685 (13 term women) |
| D3 | not feasible | GSE126124 (best of a poor set) | no cohort spans infancy to adolescence with >= 2 compartments |
