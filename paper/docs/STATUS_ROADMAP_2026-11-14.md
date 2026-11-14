# Status and roadmap (2026-11-14)

Independent assessment requested by the user on 2026-11-14 ("as independent in your opinion as possible"), followed by
a roadmap. Scope read: `results_preliminary.md` (RP), `VERDICT_TABLE.md` (VT), `METHOD.md`, the `clocksync` core and
the output tables named below. `REVIEW_2026-11-10.md`, `provisory_plan_paper.md` and `PI_QUESTIONS.md` were skimmed only
after the assessment was formed. Every number in this file was re-read in the table cited next to it (paths relative to
`paper/`); results cited without a number were read in RP or VT only. Literature references are from memory and are to
be verified before they enter any other document. No code, registered analysis or result entry was changed.

Supersedes `ROADMAP_2026-10-28.md` as the working plan. Phase 8 (naming, scorecard, figures, hardening, tutorials)
stays a PI decision and is not part of this roadmap.

---

## 1. Status

### 1.1 Bottom line

1. The method is sound and its limits are honestly recorded. Its reliable readouts are the **phase offset** and the
   **relative amplitude** of a compartment against the person's shared clock state, plus **pairing statistics with
   companions**. Absolute coupling and coupling contrasts between groups are not reliable (section 1.4).
2. One result is strong and replicated: **the tissues of a living person share one clock state, and that state is
   time**. Several further results are provocative but each rests on one post hoc contrast, one resource or an
   incomplete control (section 1.3).
3. The project is not yet a safe Cell submission. It becomes one with a decisive cross-organ perturbation result (the
   in-house multi-organ jet-lag data are the obvious candidate) and with the open checks of section 2, phase A, closed.
   The Nature Methods alternative is nearly complete in content.
4. Difference from `REVIEW_2026-11-10.md`: lead with the living-person positives and state robustness as equivalence
   of **offsets in hours**, not as coupling equivalences in GTEx followed by brain death.

### 1.2 What stands

| claim | evidence | source table |
|---|---|---|
| Organ offsets and amplitudes are recovered without times, from hospital-hours samples | GTEx, 40 tissues: offsets circular r 0.945, median error 0.35 h; amplitudes r 0.976. CHIRAL on the same samples: -0.165 | `exploratory/A123_working_hours/output/tests.csv` |
| The host state is internal time in a timed primate | baboon, 12 animals: host phase vs ZT R 0.987, median error 0.39 h; 63 tissue offsets vs timed r 0.912, error 0.45 h | `results/figure_1/output/species/a36_readouts.csv` |
| In living women, skin sites and layers share one clock state, clock-specifically | same layer, different sites r 0.628 (stress-responsive sets q95 0.420); layers, same site 0.813 (q95 0.770) | `exploratory/A251_skin_sites/output/score_tests.csv` |
| Fat and skin of a living person share a state that the cultured cell line does not carry | fat-skin r 0.533 (680 women; responsive q95 0.325); fat-LCL 0.042; skin-LCL 0.066; two tissues at one operation 0.686 (83 persons) | `exploratory/A255_living_tissues/output/score_tests.csv` |
| The state advances with the hours between two biopsies | rho 1.26 (exact times) and 2.97 (seven datasets), clock-specific; replicated in fresh studies, rho 2.53 | `exploratory/A252_elapsed_time/output/score_tests.csv`, `exploratory/A256_elapsed_replication/output/score_tests.csv` |
| Brain couples less to the body state than peripheral tissues, in three species | GTEx brain block r* 0.639 [0.530, 0.728] vs reference median 0.978; sheep median 0.205 vs 0.416 (p 0.009); baboon 0.420 vs 0.712 | `exploratory/A142_organ_blocks/output/tests.csv`, `results/figure_2/output/organs/r2_readout.csv` |
| A 12 h shift of behaviour moves the human muscle clock only partly | shift 2.68 h (90 % 1.81 to 3.59), excursion ratio 0.233; 11 subjects | `exploratory/A257_nightshift_muscle/output/score_tests.csv` |
| Mouse aorta follows a shifted light cycle within days | mean shift -2.86 h (90 % -3.95 to -1.77) | `exploratory/A201_aorta_jetlag/output/tests.csv` |
| Tumour and adjacent normal tissue pair on the clock state; immune cells of the tumour do too | TCGA pairs 0.365 vs null mean 0.122 (pseudo-clock q95 0.260; 517 pairs); myeloid 0.581 (p_clock 0.020), T/NK 0.539 (0.030) | `exploratory/A138_cancer_estimability/output/score/tests.csv`, `exploratory/A227_tme_compartment_pairing/output/tests.csv` |

Also standing, read in RP / VT only: the per-person timed validations in breast and adipose, the regenerating colon
shift (A21), the liver clock-knockout and feeding controls (WP3), and rhythmic genes following internal time (A91).
The last three are archive results and must be rerun on current code before they enter a main figure.

### 1.3 Provocative but fragile

| result | why fragile | source table |
|---|---|---|
| Between-organ order is lost in brain-dead organ donors (rho 0.678 in ventilated post-mortem donors vs 0.161 in organ donors, difference 0.517, p 0.008, 79 donors each) | one post hoc contrast in one resource; cause open. Stays descriptive (user decision): replication only, no mechanism amendments | `exploratory/A132_brain_death/output/a132b_rho_difference.csv` |
| The fetal adrenal cortex carries clock-shaped covariation (cortex C 0.608, n 34; medulla -0.048; replication 0.552 late, 0.651, 0.732 early) | not shown to be a running clock: a static between-fetus axis (gestational age, cohort, steroidogenic activity) would pass too; strongest before 12 weeks | `exploratory/A244_fetal_adrenal_classes/output/verdicts_a244.csv`, `exploratory/A242_fetal_adrenal_replication/output/verdicts.csv` |
| Internal time advances with age (+0.64 h per decade, p 0.010) | not clock-specific in its own table (`clock_specific` False); opposite sign in living cohorts (RP) | `exploratory/A234_internal_time_age/output/score_tests.csv` |
| Brain block agreement | needs a stress-responsive companion since A253 showed region pairing in controls is matched by such sets | `exploratory/A142_organ_blocks/output/tests.csv` |
| Organ-map benchmark | passes on overlapping truth only; on disjoint truth offsets r 0.679 (registered threshold 0.7: fail), amplitudes 0.965 | `exploratory/A123_working_hours/output/a123c/tests.csv` |
| Tissue order conserved across species | sheep rank conservation rho 0.331, p 0.067 (fail); baboon vs human offsets not significant | `results/figure_2/output/organs/r3_readout.csv`, `results/figure_1/output/species/a36_readouts.csv` |

### 1.4 Method limits that constrain the claims

- **Absolute coupling c** is attenuated by host reliability and its intervals undercover (A154, A162, A163). Report it
  as a rank or against companions, never as a calibrated fraction.
- **Group contrasts (d, d\*)** are biased by unequal host spread up to the size of the equivalence bound, and label
  permutation is anti-conservative for the intended null (A258; `exploratory/A258_spread_bias/output/summary_bias.csv`).
  The GTEx coupling equivalences therefore move to supplementary with this flag.
- **Gate.** The acrophase null is conservative (two-block gene phases, `clocksync/to_add.txt`), and the gate still
  passes in the clock-knockout aorta (A201). The paper rule should be pairing plus matched and stress-responsive
  companions, not the gate alone.
- **Blood** has no coherent core clock (diurnal composition): an expected exception, not a method failure.

### 1.5 Not informative, closed

Tabula Sapiens pairs (too few), blood forced desynchrony and A124 (detection limits of several hours), inflammatory
bowel disease, mother-fetus, single-cell tumour questions (depth), psychiatric brain regions (A253), glioma (A247),
chronotherapy window (A249, A250), trafficking and adhesion genes (A245), tumour-microenvironment programs (A246),
T-cell night peak (A248), paired lesions other than psoriasis (A254), menstrual cycle, lock ratio, waveform refit.
These go to a supplementary table of negative and uninformative designs.

### 1.6 Housekeeping

- Registered, not scored: **A243** (fetal heart, kidney, lung, liver, pooled open cohorts).
- No background job or subagent running; no uncommitted change to tracked files apart from this document.
- Data folder: 434 GB of the 500 GB cap.

---

## 2. Paper

**Thesis.** A person's tissues share one clock state held together by a systemic signal rather than by each organ's
own cues; in humans, misalignment is between the body and the day, not between organs.

**Target.** Cell if phase C1 delivers a cross-organ perturbation result; otherwise Science (as a sequel to Talamanca
2023) or the Cell Metabolism / Nature Communications tier. Nature Methods remains the method-led fallback.

| figure | content | status |
|---|---|---|
| 1 | Clock relations are readable without times: principle, GTEx organ map, baboon, per-person timed validations | complete; benchmarks to supplementary |
| 2 | One clock state per living person: skin sites, fat-skin vs cell line, elapsed-time rotation; the map; brain and blood as exceptions (three species for brain) | complete once A4 is done |
| 3 | What moves it and what does not: night-shift muscle, aorta jet lag, liver knockout and feeding, offsets equivalent in hours across sex, age and disease; brain-dead donors as a post hoc observation | needs A1, A5; strengthened by B1, B4, C1 |
| 4 | Before birth: adrenal cortex first, as a hypothesis | needs A2 (and C5 if possible); main or supplementary is a PI decision |
| 5 | Tumours: the microenvironment is on body time, malignant cells are not readable; corrections to published signatures | needs A3 |

---

## 3. Roadmap

Tags: **[reg]** amendment written and committed before scoring; **[sim]** simulation only; **[PI]** decision needed;
**[ask]** outward-facing, the user's decision (draft only). Each phase-A item protects a claim already in a figure;
phase B adds the one new line; phase C waits on a decision.

### Phase A. Close the open checks (data on disk)

| # | task | protects | tag |
|---|---|---|---|
| A1 | Equivalence of **offsets in hours** across sex, age, metabolic and chronic disease in GTEx, with a spread-matched null; replaces the coupling equivalences | Fig 3 | [reg] |
| A2 | Fetal adrenal covariate check: gestational age, cohort, a steroidogenic score (STAR, CYP17A1, NR4A1), and the statistic without PER1 and NR1D1 | Fig 4 | [reg] |
| A3 | COFE head-to-head on the TCGA tumour-normal pairs (COFE reports the opposite reading) | Fig 5 | [reg] |
| A4 | Matched and stress-responsive companions for the atlases that lack them, and a stress-responsive companion for the brain block (A142) | Fig 2 | [reg] |
| A5 | Rerun on current code the archive results planned for main figures: A91, sheep brain block (A42 / A103), regenerating colon (A21) | Fig 1-3 | rerun, verdicts unchanged |
| A6 | Gate: acrophase null permuted within the two gene blocks; fix the paper's claim rule as pairing plus companions | all | [sim] then [PI] |
| A7 | Score A243 as registered | Fig 4 | registered |

### Phase B. What synchronises the clocks (new line, in order of priority)

| # | task | data | main risk | tag |
|---|---|---|---|---|
| B1 | **Name the systemic signal.** Two per-sample readers built from genes that stay rhythmic in mouse liver without a local clock (Kornmann 2007): glucocorticoid targets (FKBP5, TSC22D3, KLF9) and temperature targets (CIRBP, RBM3, heat-shock genes). Which one keeps a fixed phase relation to the clock state across tissues and people, and does it loosen in brain? Extension: GTEx long-term steroid users vs other donors (aggregates only, after counting the users) | GTEx; MuTHER fat and skin as confirmation | agonal stress drives both readers in GTEx | [reg] |
| B2 | **Cross-twin pairing.** Twin A's fat against twin B's skin, compared with the within-person 0.533: equal means the shared state is the hour of the visit, lower means a person-specific internal time | MuTHER (E-TABM-1140) | pair identifiers not yet confirmed in the public sample sheet | [reg] |
| B3 | **Organ clock or cell-type clock.** Partition of offsets into organ and cell-type terms over the single-cell atlases; power scout before any residual-sharing test across organs | atlases on disk; cross-tissue immune atlas | few donors with several organs | [reg] |
| B4 | **What moves a human tissue clock, in hours.** Meal-timing and feeding-window interventions with biopsies (Wehrens 2017, Lundell 2020) next to the night-shift result; paired shifts only | check `repeat_biopsy/` first | one tissue per person | [reg] |
| B5 | Lag-amplitude curve on the organ map: a passive follower of a shared driver has relative amplitude cos(lag); tissues off the curve are candidates for an own oscillator | GTEx map, atlases | amplitude also reflects composition | [reg] |
| B6 | Local vs body-wide rhythms, new design: per-gene regression on the host state vs the tissue residual, calibrated on mouse liver clock-knockout data (the gene-set versions failed: A112, A122, A124) | GTEx; Manella 2021, Jouffe 2022 | post-mortem covariance | [reg] |
| B7 | Cell composition by internal time: one registered test | lung atlas | composition confounded with sampling site | [reg] |

Brain death: replication only (kidney transplant biopsies from living, brain-dead and circulatory-death donors,
GSE43974, coherence only; donor-type multi-organ atlases if found). No mechanism amendment.

### Phase C. Decisions and outside data

| # | item | recommendation | tag |
|---|---|---|---|
| C1 | In-house multi-organ jet-lag data | request; this is the result that makes the Cell case and the only route to dynamical modelling | [PI] |
| C2 | STARNET (dbGaP; about 600 living patients, seven tissues, Franzen 2016) | apply: a living organ map to set against the post-mortem one | [ask] |
| C3 | TwinsUK visit times | ask: a timed validation in 680 women | [ask] |
| C4 | Heart cohort GSE260643 | do not map now (one tissue per person, arrest and reperfusion act on PER1 and NR1D1, processed-only rule, 66 GB left); ask the authors for counts | [ask] |
| C5 | Unspliced counts for the fetal adrenal: pre-mRNA leading mRNA by a consistent angle is the one time-free test of a running oscillator | only if released by the authors; mapping breaks the processed-only rule | [PI] |
| C6 | Fetal figure in main text or supplementary; journal fallback order | decide after A2 | [PI] |
| C7 | New atlases: HeartMap (2026) as a scout for within-heart regions, noting that non-failing controls are mostly organ donors; basal ganglia atlas (2026) only if the human part has about 20 donors with several nuclei each | scout, not headline | [PI] |

### Not pursued

Single-cell circadian time mapping on public data (separate methods project), chronotherapy mechanisms (kept as
discussion: immune and endothelial cells on patient time favour host-mediated effects of infusion time), the age
advance in the main text, cross-species conservation of offsets.

### Order of work

A1, A2 and A3 first (two lanes), with B1 and B2 registered in parallel since both use data on disk and test the
thesis directly. A4 to A7 next. B3 to B7 follow in the order listed. C1 to C4 can start at any time and do not block
phases A and B.
