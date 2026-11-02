# C8 draft: tumour internal time and clock state as patient variables (DRAFT for the PI, not registered)

Drafted 2026-11-02. HANDOFF C8 said: "outcomes wait until C1 and C3 are in, then ask again"; the PI added (2026-11-01)
"wait for C6 and the CPTAC scout". C1 (A175, A185), C3 (A176), C6 (A182) and CPTAC (A186) are now scored. Nothing in this
file has been computed: the counts below are metadata only (outcome events, group sizes), read without any clock value.

## Question and framing

Does a tumour's clock state, or the patient's internal time at surgery, carry information about outcome beyond the known
prognostic variables? Exploratory in spirit: SCAN-B already showed no survival association for the off-plane score (A127,
HR 0.94 per SD [0.83, 1.07], p 0.35, 2,969 tumours, 272 events). The question matters for the chronotherapy thread: if a
proliferating tumour's clock is swamped (A176) and follows its host less (A175), timing treatment by the patient's clock
may not reach it.

## What the data allow (the binding constraint)

| readout | needs | cohorts on disk with outcomes |
|---|---|---|
| L1, per-tumour clock signal (in-plane share of the tumour's clock-gene profile; the per-sample counterpart of A176) | one tumour sample | METABRIC (cBioPortal; 2,509 patients, RFS 1,002 events, OS 1,144 deaths; ER+ with hormone therapy 1,104, 438 recurrences), SCAN-B (2,969, 272 deaths; off-plane already seen), NHS GSE93601 (445; 67 recurrences, 46 distant) |
| L3, tumour-to-own-normal phase offset per patient | tumour + adjacent normal, pair gate passing | GSE14520 HCC (about 213 pairs; all patients: 96 deaths, 136 recurrences); NHS (445 pairs, but the pair gate **fails**, p_acrophase 0.098, so L3 is not readable there); TCGA pairs (seen, secondary only) |
| patient internal time at surgery (phase of the normal tissue) | adjacent normal | GSE14520, NHS (normal reference gate passes, p 0.032) |
| endocrine response | pre/post therapy | POETIC GSE105777 (Ki67 response to 2-week AI; seen in A168) |
| recorded resection time (A146 replication) | timed surgery | **none**: GSE233242 is the only cohort; blocked (biobank request, PI) |

## Proposed hypotheses

Primary family (Holm over two), two-sided:
- **H-C8.1 (L1, METABRIC ER+ endocrine-treated, RFS):** per-tumour clock signal (in-plane share, 14 genes, GTEx waveform)
  in a Cox model adjusted for age, size, grade, nodes, a proliferation score (G2M) and PAM50. Prediction: no independent
  association once proliferation is in the model (A176: the signal loss is swamping by proliferation). A significant
  independent HR would be the surprise. MDE about HR 1.14 per SD unadjusted (Schoenfeld, 438 events, 80 %); larger after
  adjustment (collinearity with G2M, to be computed before scoring).
- **H-C8.2 (L3, GSE14520 HCC, recurrence):** per-patient folded distance of the tumour phase from its own liver's phase,
  relative to the cohort's median offset, in a Cox model adjusted for BCLC stage, cirrhosis, AFP, size. Read only if the
  tumour-normal pair gate passes (checked first, before any outcome is read; if it fails, H-C8.2 is not estimable and the
  family is H-C8.1 alone). MDE about HR 1.27 per SD on 136 events (paired subset smaller).

Specificity companion (both): the same Cox z-statistic with the readout computed from 200 expression-matched pseudo-clock
sets; "clock-specific" only if the clock |z| exceeds the pseudo 95th percentile. This matters: any per-tumour "signal"
score can carry proliferation or purity.

Secondary and exploratory (not counted):
- H-C8.1 in all METABRIC (RFS, OS), and in SCAN-B (OS; secondary because the off-plane score was seen) as a direction
  check.
- NHS (L1 only): recurrence and distant recurrence, ER+ (357) and all; MDE about HR 1.41-1.51 per SD, so a fail is
  uninformative.
- **Internal time at surgery** (normal-tissue phase, GSE14520 and NHS): cos/sin terms in Cox (2 df). Read as a proxy of
  the hour of surgery plus individual phase; scheduling (case order, urgency) confounds it, so no causal reading.
- **Endocrine response:** baseline per-tumour clock signal vs Ki67 change after 2 weeks of AI (POETIC, 157 treated),
  adjusted for baseline Ki67; secondary (POETIC seen in A168, though not this association).

## Stop rules and reporting

- MDE above HR 1.5 per SD: "not estimable at this n", not run as a test.
- Failed pair gate: L3 not readable (uninformative, not "decoupled").
- No outcome is read before the gate checks and the adjusted-model MDEs are written into the registration.
- One amendment per question. Report hazard ratios with intervals, never only p.

## Code needed

- clocksync: a per-sample in-plane share (projection of a sample's centred clock-gene profile on the waveform plane over
  its total clock-gene variance) next to `decompose.components`, with a test; per-sample states exist in `paired.states`.
- Cox models: lifelines or statsmodels (check which is in the sccc env).

## Decisions

1. **Family (user delegated the choice, 2026-11-02): both, H-C8.1 and H-C8.2, Holm over two.** H-C8.2 is conditional on
   its pair gate; if the gate fails, the family is H-C8.1 alone (stated before scoring, so no multiplicity is gained).
2. **Internal time at surgery (user, 2026-11-02): kept as exploratory.** Hospital hours remove night and emergency
   surgery, but within 08-18 h the hour still tracks case order (first cases are often scheduled differently), fasting
   and anaesthesia duration and team, and in a 10 h window a patient's internal time is mostly the hour of surgery. Its
   individual part (internal minus recorded time), the part free of scheduling, needs recorded times. So the readout
   stays a proxy of surgery hour plus individual phase, with no causal reading.
3. **Biobank request: avoid if possible (user).** First a targeted search for public cohorts with recorded biopsy or
   resection clock times (T5/H2 scout, 2026-11-02); the request is drafted only if that search finds nothing usable.

## Timed-cohort scout (2026-11-02, Sonnet web and GEO metadata search; nothing downloaded)

Public clinical material with per-patient collection clock time is nearly absent beyond GSE233242.
- GSE111351 (Montaigne et al., Lancet 2018): right atrial appendage at aortic valve replacement, 15 samples, one
  compartment, time only "morning / afternoon"; patient labels inconsistent. Too small and binary for C8 or T5.
- Johanns et al. (bioRxiv 10.1101/2023.04.05.535334; JHEP Reports 2023): 290 bariatric liver biopsies with the exact
  biopsy time (08-16 h); one compartment, histology outcomes, data "upon reasonable request" (no accession found;
  contacting the authors is outward-facing: PI).
- Rejected: GSE81089 (NSCLC tumour + normal + survival, surgery date without hour), GSE112660 (skin, volunteers or
  untimed), cardiac POAF cohorts (no time field in the summaries; sample fields unchecked), TCGA and glioma (no times).
- Not done: a full screen of GEO sample characteristics for fields such as "time of surgery" (E-utilities).
Reading: the A146 replication cannot be done from public data found so far. Options: the Johanns liver data (single
compartment, so internal time from one tissue only), a systematic GEO characteristics screen, or the biobank request.
- **Decision (user, 2026-11-02): the A146 replication is not pursued for now**; A146 is reported as a single-cohort
  result. The morning/afternoon heart validation already in hand is V3 (A174, GSE270258, Jia 2024, 16 left-ventricle
  snRNA-seq patients, 4 compartments), not the Montaigne atrial microarray GSE111351 (bulk, one compartment), which
  would add little.
