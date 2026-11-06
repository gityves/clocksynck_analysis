# Novelty assessment: cancer and circadian literature vs our results (2026-11-06, draft for PI review)

Question (user, 2026-11-06): what novelty do we bring, what does the literature do wrong, can we show we do better, and
does this meet the Cell criterion? Sources: `literature_positioning.md`, `manuscript/literature/synthesis.md` (148
carded full texts), `CLAIMS_LEDGER_2026-11-03.md`, its audit, `FINDINGS_MAP_2026-11-03.md`, `VERDICT_TABLE.md`, and
the output tables named per row. Wording follows ledger rule 6: convergence, not refutation.

## 1. What the cancer-clock literature does, and where it fails

| # | practice | examples (card) | problem | our test | outcome |
|---|---|---|---|---|---|
| 1 | Untimed orderers applied to surgical cohorts, then "tumour clocks are uncoupled from the host" | COFE (ananthasubramaniam2025rhythm), TimeTeller REMAGUS (vasilyev2026robust), CYCLOPS 2.0 TNBC immune-escape preprint | Orderers need the whole cycle (anafi2017cyclops, verbatim); surgery is in working hours (grabe2024, gonzalezaponte2025). Random tumour-normal differences arise from the method | A123b, A204, A204b, A140, A138 H138.4 | Orderers lose the GTEx organ map at 08-18 h (CHIRAL −0.165 to CoFE clock −0.417 vs clocksync 0.945); COFE labels collapse (largest gap 0.83 h vs 11 h true) and tumour-normal differences look uniform in 18 of 30 coupled replicates; TimeTeller loses time when training and test tissues differ in mean levels (A204) but not otherwise (A204b); CYCLOPS 2.0 can order along a non-clock factor (A140). In TCGA pairs the tumour clock follows its own normal, clock-specifically: 0.365 vs null 0.122, pseudo-clock q95 0.261, responsive q95 0.310 (517 pairs) |
| 2 | Loss of co-expression or order read as "disruption" / "dysregulation" | shilts2018evidence (ΔCCD), grabe2024quantification, Ye 2018 | Organisation scores; cannot tell a swamped clock from an uncoupled one; no matched random-gene companion | A207, A176, A179 | ΔCCD 0.890 (p 0.001) and Grabe ΔR² −0.413 present in the same tumours that pair with their normal; the coherence loss is a swamped clock plane (7 of 7 below pseudo q05), not linear cell-cycle composition (0 of 5) |
| 3 | Clock genes or clock "strength" as prognostic, without proliferation adjustment or random-set control | cadenas2014loss, li2024tumor (CYCLOPS magnitude), TimeTeller Θ, LASSO risk scores | Clock genes track the proliferation metagene (Cadenas: TIMELESS +0.674, CRY2 −0.612); random signatures are prognostic (Venet 2011) | A191-A194, A198 | Direction reproduces; information beyond proliferation not detected: Cadenas score at the 7th percentile of proliferation-matched sets; CYCLOPS magnitude HR 0.920 [0.761, 1.113]; TimeTeller Θ HR 0.966; no clock × treatment interaction |
| 4 | Lower mean clock-gene expression read as arrhythmia | wu2019pan | No rhythm, phase or co-expression test | none needed | Cite as copy-number down-regulation only |
| 5 | Single-compartment readouts asked to answer a coupling question | all of 1-4 | Coupling and offset to the body are not identifiable from one compartment; per-gene timing not refittable under hospital hours | A139 (V139.1/2 FAIL), A205 | Structural point for the introduction |
| 6 | Absence claims without power or equivalence | "no rhythm", "no relationship" statements across 1-3 | A failed test reported as a negative | equivalence and MDE throughout (A164/A202, A180, A208) | Methodological standard we apply, not a finding |

## 2. What is new (as far as the 2026-11-01 searches and the 148 carded texts show)

1. **Estimand.** Per-compartment coupling c, relative amplitude |z| and offset arg z, read against the clock state the
   same person's other compartments share, from samples taken once in hospital hours. No carded method estimates this
   (synthesis §2; closest precedent: Grabe 2024 tumour vs own normal, without coupling, gate or interval). The
   literature search was limited (positioning §9); repeat before submission.
2. **Specificity standard.** Every contrast against expression-matched pseudo-clock and responsive gene sets on the same
   permutations; gates; equivalence bounds. None of the cancer-clock papers carded uses a matched random-set control.
3. **Biology in cancer.** Tumours pair with their own host clock in bulk; where coherence is lost (proliferating,
   ER− adenocarcinomas) the clock is swamped, not reset or uncoupled; published prognostic clock signals do not add to
   proliferation.
4. **Biology beyond cancer** (the paper's frame, memory "paper framing"): organs share internal time; brain (and partly
   the gut) apart; relations kept when daily cues are lost; metabolic disease within ±0.15 coupling; psoriasis a
   local delay.

## 3. Can we show we do better? Where yes, where no

| claim | evidence | verdict |
|---|---|---|
| Works where orderers cannot (hospital hours) | A123b: Holm 0.0015 on the real organ map; A204: AUC for coupled vs uncoupled 1.0 vs COFE 0.666 [0.518, 0.804] | **yes**, the defensible "better" |
| Better estimator in general | A145 H145.1: MAE 0.060 vs CYCLOPS two-step 0.074, Holm 0.064 (FAIL); with whole days two-step pipelines recover the map (CHIRAL 0.877) | **no**: advantage is partial-day coverage, not exclusivity |
| Fewer false clocks | A145 H145.3: no-clock simulations, gate passes 4 of 30, false-claim 0.127 vs CoFE 0.06 (FAIL) | **no**: a weakness to state and ideally fix before submission |
| Calibrated uncertainty | coupling bootstrap intervals under-cover (A154, A162, A163) | **no**: permutation nulls only |
| Recovers internal time in living patients | breast n 42, Spearman 0.558 above all 200 matched sets (A200); heart AM/PM (A174) | **yes, one cohort**: fresh timed cohort needed |
| Independent truth | A123c donor-disjoint offsets 0.679 (bar 0.7, FAIL; post hoc ceiling about 0.84) | **qualified** |
| Published claims reduce to artefact or proliferation | section 1 rows 1-3 | **yes**, convergent framing; several used our refits of their pipelines (A193 deviations listed) |

## 4. Against the Cell criterion

Cell asks for a conceptual advance or a provocative question on an important biological problem. The cancer results
are mostly corrective: published "uncoupled tumour clocks" and "clock predicts outcome" do not survive a host reference
or a proliferation control. That is valuable, and provocative for a field with many such papers, but on its own it is
a critique plus nulls, with one positive finding (host pairing, bulk). Assessment:

- **Cancer alone: not Cell.** Missing a positive mechanism (what the swamped clock does or does not do), single-cell
  malignant coupling is not estimable (A138 H138.1 c 0.076 [0.014, 0.139], below the 0.083 rule; A149 not
  clock-specific; A187/A188/A195), and the corrective claims partly rest on our refits.
- **Whole paper (health organisation + local breaks, cancer as one chapter): plausible**, if the gaps below close.
  The headline "organs keep a shared internal time; breaks are local and structured" is provocative and supported by
  A142, A109/A136, A164, A138, A176, A54/A183.
- **Nature Methods: reachable now** (hospital-hours validation, orderer benchmark, specificity standard, corrections),
  conditional on the gate false-positive rate and public API, tests and tutorials.

## 5. What is missing for Cell (ranked)

1. **A deliberate multi-organ desynchrony positive** with different baseline phases (in-house jet-lag data; owner and
   PI agreement). Without it the headline rests on nulls plus the Bmal1 KO.
2. **A fresh timed-resection cohort** to replicate A200 (ledger row 1.3) and A146.
3. **Gate false positives** (A145 H145.3): reduce or calibrate, since several rows rest on gate passes (audit item 8).
4. **Audit corrections** to the ledger: headline vs brain death (A132b), "chronic disease" resting on A164, development
   on one fetal organ, row 6.2 wording given H138.1.
5. **A positive biological finding with a mechanism or prediction** in disease, e.g. why brain and gut depart, or
   whether swamped tumour clocks change dosing-time response (no dosing-time data yet; A198 null).
6. **Literature search repeated** (PubMed, bioRxiv 2026) for inter-organ, host-referenced readouts.

## 6. Discussion with the user (2026-11-06)

- **Item 2 (fresh timed cohort):** the user considers heart AM/PM (A174), breast (V25.5/A200) and GTEx restricted to
  hospital hours (A123) sufficient. Breast (user, 2026-11-06): the gate passed at the first analysis (V25.5) and
  nothing was modified after it, so A200 is a registered rerun with an added control, not a second look that shaped
  the analysis. Remaining limits: heart was analysed in the archive (V3) before A174; A123 truth shares donors (A123c
  missed its bar). Downgraded to "would strengthen".
- **Item 3 (gate false positives):** claims rest on gate plus matched pseudo-clock (and responsive) companions, not on
  the gate alone. A145 H145.3 counted gate pass plus c > 0.2 only; the false-claim rate of gate plus companion was
  never simulated. A small registered check (A145 arm N with a pseudo gene set receiving the shared factor) would close it.
- **Item 4 (brain death):** read as loss of internal organisation when the central drive is lost. Precedents: SCN
  lesion destroys inter-tissue phase synchrony within animals (yoo2004periodluciferase, "a form of 'internal
  desynchronization'"); the SCN is required for synchronisation among tissues (sinturel2021circadian); contrary:
  tahara2012in (phases similar in rhythmic SCNx organs). Ventilated, non-brain-dead donors share ICU care and
  continuous feeding and keep the order (A132b 0.678 vs 0.161), and glucocorticoids do not explain it (A136). Caveats:
  post hoc split; A136 on the same donors; other donor-management hormones (thyroid, vasopressin) not tested. The
  headline wording changes to "organs keep a shared internal time unless the central drive is lost".
- **Item 5 (prediction, mechanism):** positives exist (psoriasis, tumour pairing, brain). Proliferation is a
  correlate, not a shown mechanism: physiological proliferation does not reproduce it (A126), no cell-cycle locking
  (A187b, exploratory), single-cell and arrest tests not estimable (A195, A196).

## Limits of this assessment

Numbers are from the verdict table, the ledger and the CSVs re-read here (A145 `h145.csv`, A138 rows, A206, A207); the
audit's number fixes are not yet applied to the ledger. "Not found in the literature" means not in 148 carded texts and
two searches.
