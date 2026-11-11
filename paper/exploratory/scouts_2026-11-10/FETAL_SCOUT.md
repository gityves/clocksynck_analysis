# Scout 2026-11-10: published human embryonic and fetal data in which a clock could be tested

Metadata only: no clock statistic was computed on any dataset listed here. Purpose: the human developmental part of
the paper rests on one organ in one cohort (A189, fetal adrenal, Cao 2020, 7 fetuses) and the mouse counterpart was
downgraded (A236). This scout looks for independent human designs, in any organ.

Counts marked **(v)** were recomputed from the metadata files (sdrf, GEO series-matrix header, CELLxGENE curation
API); the others are as reported by the source and must be checked before use.

## What the method needs

| readout | design | lesson from earlier work |
|---|---|---|
| L1 coherence (clock genes keep the adult phase order across individuals) | one compartment, many individuals | D1b: the number of individuals is the limit, not depth or the sampling window (GTEx at the Cao design: liver n=7 0.805, lung n=11 0.715, heart n=9 0.51, cerebrum n=8 0.07) |
| L2 shared state (compartments of one individual share a clock state) | >= 2 compartments per individual | A181: not estimable in Cao across organs (design-limited) |
| specificity | expression-matched and responsive pseudo-clock sets on the same permutations | A236: within a developmental series, age and any other shared programme pass a clock test unless companions and an age adjustment are part of the design |

A positive needs clock phases that differ between individuals (different collection times, or free-running clocks).
A negative at a design with adult power >= 0.8 is informative; below, it is not.

## Shortlist

| rank | dataset | design (v = verified) | what it can test | access, size |
|---|---|---|---|---|
| 1 | Cardoso-Moreira 2019, E-MTAB-6814, bulk RNA-seq | 207 prenatal libraries in the matrix, 4-20 wpc **(v)**: heart 34, liver 37, hindbrain 33, forebrain 32, kidney 27, testis 26, ovary 18 **(v)**; the individual is not released for 55 of them and was resolved here from the expression data (see below): 151 linked libraries, 66 individuals, 41 with >= 2 organs, 25 with >= 3 **(v)**; same-individual pairs: heart-liver 17, forebrain-hindbrain 14, forebrain-liver 13, forebrain-heart 13, heart-hindbrain 12, heart-testis 12 **(v)**; 90 postnatal libraries (neonate to elderly; 9-25 per organ, no ovary) **(v)** | L2 across organs of one fetus (the paper's question, asked in development); L1 per organ with the postnatal samples of the same organ and platform as positive control; age trend 4-20 wpc | open; CPM (counts recoverable as CPM / smallest non-zero CPM of the library), about 110 MB |
| 2 | Garcia-Alonso 2022, gonads, CELLxGENE 661a402a, scRNA 10x 5' | 51 gonad donors (female 30, male 21), 6-21 wpf **(v)**; donors with >= 200 cells: mesenchymal-stromal 47, supporting 46, epithelial 33, endothelial 24 **(v)**; 47 / 36 / 17 donors with >= 2 / 3 / 4 of these four classes (47 / 37 / 22 over all classes) **(v)**; 39 donors 5' v1, 12 donors 5' v2 **(v)** | L1 per class at n = 24-47 (three to seven times the Cao adrenal n); L2 across cell classes within the gonad | open; h5ad 1.3 + 2.7 GB |
| 2b | Fetal gonad and reproductive tract (2024), CELLxGENE 222e800b, 10x 3' | gonad: 6 female and 5 male donors **(v)** | too few donors for a test of its own: descriptive replication of rank 2 | open |
| 3 | Kang 2011, GSE25219, Affymetrix exon arrays | 57 brains, 1,340 samples **(v)**; prenatal 24 (5.7-37 PCW: 7 at <= 12, 14 at 13-22, 3 at 25-37) **(v)**, 20 of them with >= 6 of the 16 regions sampled at every later stage (the earliest brains carry other regions) **(v)**; postnatal 33 (6 infants 4-12 months, 7 children 2-11 y, 20 from 13 to 82 y) **(v)**; PMI, pH, RIN per sample; the gene-level table lacks PER3 and NR1D1 (12 of the 14 waveform genes) **(v)** | L2 across regions of one brain, prenatal against infant and adult brains on the same platform; only dataset reaching the third trimester and the first postnatal year | open; 73 MB gene-level matrix |
| 4 | Velmeshev 2023, cortex, CELLxGENE bacccb91, snRNA | 106 donors, second trimester to adult **(v, total)**; coarse cell types | trajectory within cortex cell classes; prenatal n unknown until the obs is read | open; 6 GB |
| 5 | RPE / choroid atlas, CELLxGENE 48c15b0c (scRNA dataset) | 60 donors, 10 wpf to adult **(v, total)**; fibroblasts, endothelium, melanocytes | stromal-vascular classes before and after birth in one dataset | open |
| 6 | Fetal retina, CELLxGENE c11009b8 (22 donors, 9-27 wpf) and 5900dda8 (14 donors, multiome) **(v, totals)** | retina has an autonomous clock in the adult (project ED panel) | L1 / L2 across retinal cell classes, discovery and replication | open; 0.7 GB |
| 7 | He 2022, fetal lung, CELLxGENE 2d2e2acd | 12-13 donors, 9-22 wpf **(v, total)** | fibroblast, endothelium, epithelium within lung; at the D1b limit (lung fibroblast 0.815 at n=11) | open; 0.9 GB |
| 8 | Fetal cortex in Down syndrome, CELLxGENE 0e9fd1d3, multiome | 30 donors, 10 wpf onwards **(v, total)** | a disease contrast in fetal cortex; only after a fetal brain positive | open |

Not retained: limb (Zhang 2024; 25 donors but <= 10 wpf), first-trimester brain (Braun 2023; 26 donors, <= 15 wpf,
12 GB; covered by rank 1 and 3), HDBR brain E-MTAB-4840 (190 individuals but FASTQ only: no mapping in this project),
GSE66302 (17 fetuses, DeepSAGE, processed data not located), BrainSpan RNA-seq (same donors as Kang, RPKM only).
Controlled access, would need a request (ask the user first): raw fetal adrenal single-cell data of Kildisiute 2021
and del Valle 2023 (both in EGAD00001008345, Sanger DAC; correction 2026-11-11: the processed Kildisiute object is
open, see below), Jansky 2021 (EGAD00001006624), fetal brain eQTL cohorts (O'Brien 2018, 120 brains, EGA; Walker 2019,
201 donors; BrainVar, 176 donors, Synapse terms not checked). No public fetal or neonatal human dataset with a
recorded time of collection was found.

## Who is who in Cardoso-Moreira 2019 (prep only; clock genes never read)

The matrix names its columns `<Organ>.<stage>.<k>` and gives no library; Supplementary Table 2 lists the libraries
with stage, sex and, for prenatal ones, the HDBR individual. `prep_fetal_bulk/scripts/link_cardoso_individuals.py`
assigns each prenatal column to a library in four tiers (`link_cardoso_summary.csv`): one candidate of that organ,
stage and sex (exact, 82 columns); three deletion-like on / off genes, GSTM1, NPIPB15 and COX20P1 (marker, 29); a
panel of 13 genotype-like genes whose level is concordant between the organs of one individual (panel, 21); the same
for two organ groups taken together (pair, 20); 55 columns stay open and are used only in the `organs_all`
sensitivity. Checks: among same-individual pairs of exact columns the three markers agree in 58 / 58, 14 / 14 and
28 / 28 pairs (chance 0.50, 0.57, 0.61); after all tiers COX20P1 disagrees in 1 of 57 pairs
(`link_cardoso_marker_validation.csv`). Leaving each individual out of the panel selection, 1,183 decisions: right
with the margin 0.953, undecided 0.047, wrong 0 (`link_cardoso_panel_loo.csv`); pair rule, 76 decisions, all right.
A wrong link exchanges two individuals of one organ, stage and sex, which can only weaken an across-organ statistic.
Sex check (Y-linked genes of the column against the sex of the table, 297 libraries): one disagreement among the
linked libraries. Both liver columns at 20 wpc are male, the table lists one male and one female library, and the
panel gave one column to the female individual (female forebrain and kidney columns): the source disagrees with
itself there. A linked library whose sex differs from the majority sex of its individual is treated as unlinked (that
liver library; 151 linked libraries remain).

## Expectations from the literature (scout summary of abstracts; to be verified before citing)

- Non-human primate: clock genes oscillate in the fetal capuchin SCN and adrenal in late gestation (Torres-Farfan 2006,
  doi 10.1210/en.2006-0628). Rat fetal adrenal clock (Torres-Farfan 2011, doi 10.1210/en.2010-1260).
- Rodent: the molecular clock appears with differentiation, in mid to late gestation (Yagita 2010, doi
  10.1073/pnas.0913256107; Umemura 2017, doi 10.1073/pnas.1703170114; Umemura 2022, doi 10.1073/pnas.2114083119).
- Human infants: day-night rhythms of cortisol and melatonin appear in the first 2-4 postnatal months (review:
  Joseph 2015, doi 10.1136/archdischild-2014-306104).
- Hence: no clock expected in the first trimester; a positive in the second trimester would agree with the capuchin
  and A189 adrenal results; the first postnatal year is the transition to adult organisation.

## Feasibility (adult material cut to each design; `feasibility_fetal/`; no fetal expression value read)

Each fetal individual is replaced by one adult donor with the same compartments, who receives its age, sex and
stratum; the registered permutation test is run at alpha 0.05 without the companions; feasible if the pass rate is
>= 0.8 (D1b rule). 200 subsamples (Kang: 100).

| dataset | adult material | readout (n) | adult pass rate | retained as |
|---|---|---|---|---|
| Cardoso-Moreira | GTEx, deaths at 08-18 h | L2 across organs, 9-20 wpc (26) | 0.825 (CI 0.765-0.875) | H237.1 |
| | | L2 across organs, 4-20 (41), 13-20 (11), 4-8 wpc (15) | 0.630, 0.165, 0.140 | descriptive |
| | | L1 4-20 wpc: heart (34), liver (37), hindbrain (33), kidney (27), forebrain (32), ovary (18) | 1.000, 1.000, 1.000, 1.000, 0.965, 0.940 | H237.2 |
| | | L1 4-8 wpc: heart (14), liver (16), hindbrain (16), kidney (9) | 0.990, 0.985, 0.985, 0.875 | H237.3 |
| | | L1 13-20 wpc: heart (9), liver (10), kidney (10) | 0.870, 0.960, 0.870 | H237.3 |
| | | L1 postnatal: forebrain (21), hindbrain (25), heart (10), kidney (9), liver (12) | 0.880, 1.000, 0.955, 0.895, 0.975 | positive controls |
| | | L1 testis (7-26); forebrain 4-8 (12), 13-20 (11); ovary 4-8 (8), 9-20 (10) | 0.215-0.760; 0.640, 0.640; 0.680, 0.775 | descriptive |
| Garcia-Alonso | breast nuclei (10x 3' v3), thinned to the fetal depths | L1 stromal (47), supporting (46), epithelial (33), endothelial (24) | 1.000, 0.995, 0.985, 0.930 | H238.1 |
| | | L2 four classes (47); three classes (44); stromal-supporting pair (44) | 1.000; 0.990; 1.000 | H238.2; descriptive |
| Kang | GTEx brain donors with >= 6 of 13 regions, any hour, 12 genes | L2 across regions: prenatal (20), postnatal (32) | 0.96, 0.98 | H239.1, control |
| | | L1 brain mean: prenatal (20), postnatal (32) | 1.00, 1.00 | H239.2, control |
| | | 13-22 pcw (14), infant and child (13), from 13 y (19): L2; L1 | 0.91, 0.95, 0.95; 0.97, 0.89, 0.98 | descriptive |

Sources: `passrates_cm_day_sex.csv`, `passrates_gonad_sextrend.csv`, `passrates_kang_all.csv` (`output/`). Limits: the
companions are not in the calibration; the Kang calibration is RNA-seq without the sex covariate (the postnatal brains
of the same arrays are its real power check); the gonad stand-ins are another tissue and nuclei.

**The L1 statistic had to change.** The unsigned coherence R (A189) was not feasible in several organs at n > 25
(adult pass rates at 4-20 wpc: forebrain 0.055, kidney 0.495, heart 0.730; `passrates_cm_day.csv`). Cause: R reads a
shift common to the clock genes of a sample as coherence, and such a shift also hides a clock. Simulation
(`shift_artefact.py`, 200 data sets per scenario, `shift_artefact.csv`): without a clock, a common shift of sd 0.2
passes the permutation test of R in 0.625 (n 10) and 1.000 (n 40) of the data sets; a signed coherence with a
per-sample offset fitted out (`clocksync.decompose.signed_coherence`) passes in 0.050 and 0.035 and detects a clock in
1.000 at every shift (R: 0.690 at n 10 with the shift). The fetal amendments use the signed statistic; whether it
replaces R in the earlier results is PI question Q31. The A189 adrenal positive is unchanged by it (post hoc: C 0.926,
p 0.0001; `A189_fetal_within/output/h189_1_signed_posthoc.csv`).

## Plan

1. Download ranks 1-3 (parse only, `PROVENANCE.md` per source) and build the design tables; no clock statistic. Done.
2. Feasibility on adult material cut to each design (D1b rule: pass rate >= 0.8). Done (above).
3. Register one amendment per dataset, with an age adjustment, age-stratified permutations and both companions:
   A237 (Cardoso-Moreira), A238 (Garcia-Alonso), A239 (Kang), `paper/docs/registered_experiments.md`.

**Addendum 2026-11-11: fetal adrenal data for a replication of A189 (links checked; nothing requested, no amendment
registered).** Open: Kildisiute 2021 processed object, https://www.neuroblastomacellatlas.org/
(`adr_all_cellxgene.h5ad`, downloaded for a design check, `natexp/kildisiute2021/PROVENANCE.md`: 57,972 fetal adrenal
cells, 7 samples at 5 gestational ages, w8 to w21, log-normalised values with the total count per cell); Kameneva
2021, GEO GSE147821 (13 samples, weeks 6-14, 10x matrices per sample); Dong 2020, GEO GSE137804 (6 fetal adrenal
samples, ~4 to 14 PCW). Against the age of the Cao adrenals (89-122 days, 12.7-17.4 weeks): Kildisiute's samples
lie on either side (<= w11, w21), Dong has two in range (13+2, 14+4 PCW), Kameneva's oldest reach week 14 (ages per
sample and their convention not checked). The accession
printed in Kildisiute 2021 (EGAD0000100637) is one digit short; ENA PRJEB41516 holds the tumour and organoid reads
only.
