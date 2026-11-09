# Method edge: which results need clocksync (2026-11-09)

For each result, could an available method (untimed orderers CYCLOPS / CYCLOPS 2.0 / CHIRAL / CoFE / DCPR; supervised
predictors TimeTeller / ZeitZeiger / TimeSignature / BodyTime; CCD; the Grabe procedure; single-cell orderers Tempo /
CoPhaser / scRitmo) have produced it from the same data? Rows follow `CLAIMS_LEDGER_2026-11-03.md`, plus results scored
after it. Numbers are in `results_preliminary.md` and their `<!-- src: -->` tables; the per-method review of how
published tools handle hospital hours, with our tests of each, is the last section of `paper/presentation/introduction.md`.

**Evidence type:** *demonstrated* = an alternative was run on the data (by us or by its authors); *argued* = the
alternative cannot be applied (sample size, sampling window, no timed training data, platform), not run.

## Facts the arguments rest on

- Untimed orderers lose the GTEx organ map in 08-18 h deaths and recover it with whole-day sampling (A123b, H123.4;
  post hoc: not a reflection or stretch artefact). In a clock-gene-only simulation two-step pipelines are close to
  clocksync even in 08-18 h (A145, H145.1 FAIL after Holm): the edge is shown on real data, not in the simulation.
- CYCLOPS asks for > 250 samples (`anafi2017cyclops`).
- TimeTeller is limited to microarrays for human data: "comprehensive healthy human RNA-seq circadian time-series
  datasets are not yet available" (`vasilyev2026robust`). It works in hospital hours on GTEx (A123b) but is biased when
  the per-gene means of training and test tissue differ (A204 gate failure, post hoc reading; A204b).
- CCD is group-level and "cannot detect a phase difference" (`shilts2018evidence`).
- Single-cell orderers order cells within a sample, where cells share the sampling time.

## A. Method-exclusive positive results

| result | data | why alternatives fail | evidence |
|---|---|---|---|
| 1.1-1.2 organ map from 08-18 h deaths (A123, A123b) | GTEx window donors | orderers fail on the same data | demonstrated |
| 6.2 tumours follow their patient's normal tissue (A138 H138.4) | TCGA, 517 pairs, surgical | CoFE (`ananthasubramaniam2025rhythm`) reported uncoupling; CYCLOPS 2.0 (`li2024tumor`) agreement only in strong-clock luminal A | demonstrated (published methods applied to the question) |
| 1.3 breast resection internal time follows the surgery hour (A200) | GSE233242, 42 patients, 08-18 h | far below orderer sample needs, narrow window; TimeTeller excluded if RNA-seq (platform not checked) | argued |
| 6.1 psoriatic lesion delay vs the same person's skin (A54, A183) | GSE154200, untimed clinic biopsies, paired, about 20 persons per batch | orderers: window and size; a healthy-skin predictor (`wu2018populationlevel`, `wu2020population`; microarray, about 3 h per-sample error) would read lesional mean shifts as timing; clocksync fits mean and phase per compartment separately | argued |
| tumour immune and stromal compartments pair with the patient's normal tissue (A227-A229) | 34 single-cell cohorts, about 7 patients per cohort | sample size, hospital hours, single-cell RNA-seq; CCD phase-blind; TimeTeller microarray-only and no cell-type training data | argued |
| sheep brain partly apart from the body (A42, archive) | SheepGTEx, 343 animals, no slaughter times | window sampling assumed, not recorded; no sheep training data | argued, weak |

Caveats carried by these rows: the psoriasis delay sits at the matched-set specificity line and is not specific
against eczema (A183); microenvironment pairing is lung-weighted, and immune pairing may be cells following body time
through the circulation (A227, A228); A200 is a rerun on data seen before.

## B. Method-exclusive design, answer is a bound or "not detected"

| result | data | note |
|---|---|---|
| 5.2 disease atlases: phase and amplitude equivalence (A180), no effect detected (A165) | single-cell disease atlases, hospital hours | within-person cell types no other method reads; most rows not testable |
| A230 epithelium vs stroma in inflamed gut | TAURUS | not detected (MDE 1.7 h) |

## C. New questions, not method-exclusive

| result | data | alternative that could reach it |
|---|---|---|
| 2.1 shared clock state beyond matched genes (A150, A144) | GTEx, deaths over the whole day | CHIRAL-style two-step (A123b whole day: CHIRAL and CoFE clock genes recover the map; A145 arm F: two-step coupling equivalent) |
| 2.2 human brain apart (A142, A143); 2.3 digestive group | GTEx | two-step with a reliability correction (brain amplitude is lowest); or time of death directly |
| 4.1-4.3 loss of cues, brain death, temperature (A109, A132b, A136, A110) | GTEx | two-step plus time of death; untested in ventilated donors |
| 5.1 T2D and obesity (A164, A202); A233 five chronic conditions; A235 (registered, not scored) | GTEx | two-step coupling |
| 2.4 ageing (A3, A19); A231 coupling vs age; A234 internal time vs age | GTEx, timed brain | CHIRAL donor phase vs time of death; timed cosinor |
| A232 offsets vs proliferation | GTEx | two-step (failed anyway) |
| 3.1 mouse embryo (A108, A121) | collection time recorded | timed cosinor; CCD for the earlier, untimed windows |
| 3.2 human fetal adrenal (A189) | 7 untimed fetuses | CCD-type co-expression (low power); A190 withdrawn |
| 6.3-6.4 tumour coherence loss and reversibility (A127, A166, A171, A173, A176, A168, A196); A212-A225 | bulk tumour cohorts | CCD, CYCLOPS magnitude, TimeTeller Θ; luminal A advantage replicates `li2024tumor`; ΔCCD replicated (A207) |
| 6.5 survival (A191-A194, A198); 6.6 microenvironment vs recorded hour (A146) | breast and pan-cancer cohorts | uses the published readouts or the recorded hour |
| 7.2 aorta after a 6 h advance (A201) | timed mouse | timed cosinor |

## D. Validation and controls (exclusivity not applicable)

1.4 heart AM/PM (A174), adipose (A113), baboon (A36, A40); 1.5 hepatocyte Bmal1 KO (A111, A129) and feeding shift;
7.1 blood limit (A124, A156); 8.1-8.2 method comparisons (A199, A140).

## Reading

- The exclusive positive results are the hospital-hours organ map, tumour-normal pairing, breast resection time,
  psoriasis and the tumour microenvironment: clinical, within-person, hospital-hours material. Only the first two are
  demonstrated; the others rest on the alternatives being inapplicable.
- Most health, cue-loss and ageing results rest on GTEx, where a two-step pipeline on clock genes could have reached
  them; their novelty is the question and the matched-null design.
- Enough for a methods argument (capability on the data class most human disease tissue belongs to). For a biology
  venue the method is the means, and the exclusive biology is the cancer and skin results.

## Possible future tests (not registered)

1. **Psoriasis head-to-head (priority):** CYCLOPS and, if the platform allows, the Wu skin predictor on GSE154200
   lesional vs non-lesional pairs. Turns row 6.1 from argued to demonstrated.
2. Platform of GSE233242 (row 1.3): if microarray, TimeTeller on the 42 patients against the recorded hour.
3. Two-step coupling (orderer then regression, as A145) on GTEx 08-18 h deaths: extends A123b from phase to coupling
   (user, 2026-11-09: phase is enough; optional).
4. CCD and CYCLOPS magnitude on the A215 (TP53) and A219 (subtype) strata: states whether the tumour-intrinsic
   contrasts are method-independent (user, 2026-11-09: not needed now).
5. GTEx brain split with a CHIRAL two-step and a reliability correction: would show row 2.2 is reachable without
   clocksync, or that it is not.
