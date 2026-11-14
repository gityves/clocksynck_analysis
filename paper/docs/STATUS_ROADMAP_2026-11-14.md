# Status and roadmap (2026-11-14)

Independent assessment requested by the user on 2026-11-14 ("as independent in your opinion as possible"), followed by
a roadmap. Scope read: `results_preliminary.md` (RP), `VERDICT_TABLE.md` (VT), `METHOD.md`, the `clocksync` core and
the output tables named below. `REVIEW_2026-11-10.md`, `provisory_plan_paper.md` and `PI_QUESTIONS.md` were skimmed only
after the assessment was formed. Every number in this file was re-read in the table cited next to it (paths relative to
`paper/`); results cited without a number were read in RP or VT only. Literature references are from memory and are to
be verified before they enter any other document. No code, registered analysis or result entry was changed.

**Updated 2026-11-15:** sections 2 and 3 restructured around the user's target (plan 1) and its two fallbacks.

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
3. Target and fallbacks (user decision, 2026-11-15; section 2). **Plan 1:** answer what internal time controls in a
   human tissue (the body's time, the tissue's own clock, behaviour, or none); if answered, a tier-one paper.
   **Plan 2:** the biology as it stands, one internal time per person and where it is absent; with only the open
   checks closed its realistic level is Nature Communications or Cell Reports, higher with the offsets bound and the
   systemic signal. **Plan 3:** the method, in Nature Methods, nearly complete in content.
4. Difference from `REVIEW_2026-11-10.md`: lead with the living-person positives and state robustness as equivalence
   of **offsets in hours**, not as coupling equivalences in GTEx followed by brain death. "Organ clocks are aligned" is
   known from rodents and the 2023 population map: on its own it is a confirmation, not the headline.

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

## 2. Paper: target and fallbacks (user decision, 2026-11-15)

User, 2026-11-15: "if we are able to answer this it would be a tier one paper [...] Let's modify the roadmap to include
this as the target, and if it fails a fallback on the biology as it stands now, and if this [is] not interesting enough
the method itself".

| plan | paper | journals | taken when |
|---|---|---|---|
| 1, target | **What internal time controls**: which daily processes of a human tissue follow the body's time, which follow the tissue's own clock, which follow behaviour or the external day, and which follow none | Cell, Science, Nature | gate G2 passes |
| 2, fallback | **The biology as it stands**: one internal time per person, and where it is absent | by tier (section 2.2) | G2 fails, G3 passes |
| 3, fallback | **The method** | Nature Methods | G3 fails |

All three plans share phase A, so no work is lost by aiming at plan 1 first.

### 2.1 Plan 1: what internal time controls

**Why it is open.** The default view is that each organ's clock times that organ's physiology. Mouse work with tissue
clock knockouts and reconstitutions (Kornmann 2007, Koronowski 2019, Greenwell 2019; from memory) shows that only a
minority of liver rhythms need the liver clock. In humans the question could not be asked: there is no knockout, and
"alignment of organ clocks" itself is known from rodents and from the 2023 population map, so confirming it is not a
tier-one result.

**The handle.** Natural variation dissociates several times within one person, and `clocksync` reads them without
time stamps:

| time | read from | dissociated from the others by |
|---|---|---|
| body's internal time | host state shared by the person's compartments | chronotype, shift of the day (GTEx: host state vs time of death) |
| the tissue's own clock | the tissue's deviation from the host state | person-to-person offsets; clock-poor compartments (tumour, blood, cell line) |
| behaviour (sleep, meals) | protocol | night-shift protocol: central time about 0 h, muscle clock 2.7 h, behaviour 12 h |
| external day | recorded hour | elapsed hours between two biopsies; time of death |

**Thesis to test.** Time in a human tissue is largely imposed by the body; the local oscillator is one of several
readers of it. Either outcome (mostly body-driven or mostly locally driven) is an answer; an underpowered partition is
not.

**What already points this way** (none sufficient): rhythmic genes follow internal rather than external time (A91,
archive, to be rerun); muscle moves 2.7 h of a 12 h behaviour shift (A257); blood shows the day through cell
composition without a coherent core clock; a cultured cell line carries no detectable trace of its donor's state (A255). In
tumours the evidence is weak: clock-controlled genes followed the body time of the adjacent normal sample in TCGA
(A249, descriptive), but not clock-specifically in CPTAC, where NFIL3 alone replicated (A250). It is not to be worded
as a finding.

**Figure outline.**

| figure | content | needs |
|---|---|---|
| 1 | Clock relations read without times: principle, GTEx organ offsets, baboon | complete |
| 2 | One internal time per living person, bounded in hours: skin sites, fat-skin vs cell line, rotation between biopsies | A1, A4 |
| 3 | The times of a tissue come apart: body vs external day, tissue vs body, behaviour vs clock | A5 (A91), T1 |
| 4 | What follows which time: genome-wide partition by tissue and by program; the systemic signal | T2, T3, T4 |
| 5 | Where the local clock is weak or absent (brain, fetus, tumour, blood) and what still keeps time there; the cancer chapter | A2, A3, T5 |

### 2.2 Plan 2: the biology as it stands

**Thesis.** A person's tissues share one clock state held together by a systemic signal rather than by each organ's
own cues; in humans, misalignment is between the body and the day, not between organs. Its tier depends on what is
added:

| evidence | realistic target | roadmap items |
|---|---|---|
| open checks closed only | Nature Communications, Cell Reports | phase A |
| offsets bounded in hours, plus the systemic signal named or a living multi-organ map | Cell Metabolism, possibly Science | A1, and T4 or C2 |
| a cross-organ perturbation result on top (optional data, not a requirement) | Cell or Science credible | C1 |

**What justifies aiming high.** A method that removes the need for timed, repeated sampling and is validated against
truth (section 1.2); a lead result replicated in living people with specificity companions; findings that depart from
the rodent picture (the partial shift of muscle, brain as the exception in three species) or are new (fetal adrenal
cortex, brain-dead donors, the tumour microenvironment on body time); pre-registration with negatives reported.

**What could pull it lower, and the item that answers it.**

| risk | answer |
|---|---|
| The headline reads as expected: tissues sampled at the same hour share a clock state. The new content is how tightly they agree, in hours | A1, B1 |
| No mechanism and no perturbation across organs; the systemic signal is a hypothesis | T4, B3, C1 |
| Living evidence covers three superficial tissues (fat, skin, muscle); liver, gut and heart are post-mortem only | C2 |
| Each provocative result is fragile on its own (section 1.3) | A2, A7, brain-death replication |
| No functional consequence: 14 clock genes, no link to physiology | this is plan 1 (phase T) |
| Visible method limits (section 1.4; disjoint-truth benchmark 0.679 against a registered 0.7) | A4, A6; report offsets, not coupling |
| Novelty partly taken: the 2023 GTEx organ map; COFE claims the tumour side with the opposite reading | A3; lead with living people, not the map |
| Scope reads as a collection | frame every chapter as a test of the thesis |

**Outlines.** One paper including the cancer chapter (user, 2026-11-06), in either of two formats:

| figure | long format (Cell) | short format (Science): one rule and its exceptions |
|---|---|---|
| 1 | readable without times: principle, GTEx map, baboon, per-person timed validations | same, condensed |
| 2 | one clock state per living person; the map; brain and blood as exceptions | one time per living person, with the bound in hours across sex, age and disease |
| 3 | what moves it and what does not: night-shift muscle, aorta jet lag, liver knockout and feeding, offsets in hours, brain-dead donors (post hoc) | what moves it: night-shift muscle, aorta jet lag, systemic signal |
| 4 | before birth: adrenal cortex first, as a hypothesis | where it is absent, one panel each: brain, fetus except adrenal cortex, malignant cells vs microenvironment, brain-dead donors |
| 5 | tumours: microenvironment on body time, corrections to published signatures | - (tumour detail in the supplement) |

### 2.3 Plan 3: the method

Untimed, host-referenced inference that works on hospital-hours samples, single-cell atlases and two-tissue cohorts,
where time-assignment methods fail; validated in GTEx, baboon, timed human cohorts and mouse controls; calibration
limits stated (section 1.4). Content is nearly complete after A4 and A6; the public API, tests and tutorials are
Phase 8 (PI).

### 2.4 Gates

Thresholds are fixed in the amendments, before scoring.

- **G1, after the power scouts (T0).** A design of phase T is registered only if its scout shows it can separate the
  times at the program level. If none can, go to G3 directly.
- **G2, after phase T is scored.** Plan 1 is taken if the partition is (i) specific: the body and tissue readers beat
  matched pseudo-clock and stress-responsive readers; (ii) calibrated: the same procedure recovers the known classes in
  mouse liver with and without a local clock; (iii) replicated in at least two independent designs or tissues;
  (iv) quantitative: a share of genes or programs per time, with intervals.
- **G3, after A1.** Plan 2 is taken if the between-tissue offset bound is tight enough to be a headline (bound set in
  the A1 amendment). Otherwise plan 3.

---

## 3. Roadmap

Tags: **[reg]** amendment written and committed before scoring; **[sim]** simulation or scout, no hypothesis scored;
**[PI]** decision needed; **[ask]** outward-facing, the user's decision (draft only).

### Phase A. Open checks, shared by the three plans (data on disk)

| # | task | protects | tag |
|---|---|---|---|
| A1 | **Headline bound.** Between-tissue offset spread within a person, in hours, with its detection limit (GTEx, and the living cohorts where two tissues exist); equivalence of offsets in hours across sex, age, metabolic and chronic disease, with a spread-matched null; replaces the coupling equivalences | plans 1 and 2; gate G3 | [reg] |
| A2 | Fetal adrenal covariate check: gestational age, cohort, a steroidogenic score (STAR, CYP17A1, NR4A1), and the statistic without PER1 and NR1D1 | fetal result | [reg] |
| A3 | COFE head-to-head on the TCGA tumour-normal pairs (COFE reports the opposite reading) | cancer chapter | [reg] |
| A4 | Matched and stress-responsive companions for the atlases that lack them, and a stress-responsive companion for the brain block (A142) | organ map, brain exception; plan 3 | [reg] |
| A5 | Rerun on current code the archive results planned for main figures: **A91 first** (the first piece of plan 1), then the sheep brain block (A42 / A103) and the regenerating colon (A21) | plans 1 and 2 | rerun, verdicts unchanged |
| A6 | Gate: acrophase null permuted within the two gene blocks; fix the paper's claim rule as pairing plus companions | all; plan 3 | [sim] then [PI] |
| A7 | Score A243 as registered | fetal result | registered |

### Phase T. Target line: what internal time controls

| # | task | data | main risk | tag |
|---|---|---|---|---|
| T0 | **Scouts before any registration.** (a) Why the earlier gene-set attempts failed (A112, A122, A124) and what design avoids it. (b) Power of T1 to T3 at the program level, by simulation on the real designs. (c) Calibration set: mouse liver with and without a local clock, fed in and out of phase | Manella 2021, Jouffe 2022 and the WP3 mouse data on disk | - | [sim] |
| T1 | **Behaviour vs tissue clock vs central time.** Genome-wide in the night-shift muscle data: assign programs to the time they follow (about 0 h, 2.7 h, 12 h) | GSE106800 (11 subjects, 44 arrays; RP A124) | few arrays; A124 had a 7.3 h detection limit with gene sets, so genes must be pooled into programs | [reg] |
| T2 | **Body vs tissue vs external day, per gene.** Regress each gene on the host state, on the tissue's deviation from it, and on time of death given the host state (aggregates only). A gene following the deviation is under local clock control | GTEx | sample quality drives both the deviation and the gene: matched and stress-responsive companions required | [reg] |
| T3 | **Rotation in living people.** Which genes rotate with the clock state between two biopsies, and which do not | repeat-biopsy cohorts of A252 / A256 | interventions between biopsies | [reg] |
| T4 | **Name the systemic signal.** Per-sample readers of glucocorticoid targets (FKBP5, TSC22D3, KLF9) and temperature targets (CIRBP, RBM3, heat-shock genes), built from genes that stay rhythmic in mouse liver without a local clock: which keeps a fixed phase relation to the host state across tissues and people, and does it loosen in brain? Extension: GTEx long-term steroid users (aggregates only, after counting the users) | GTEx; MuTHER fat and skin as confirmation | agonal stress drives both readers in GTEx | [reg] |
| T5 | **Compartments with a weak or absent local clock** as natural knockouts: what still follows body time there (tumour, blood, cultured cell line). In tumours the reader must be shown clock-specific first (A250 was not) | TCGA and CPTAC pairs, blood cohorts, MuTHER | low odds in tumours | [reg] |
| T6 | Lag-amplitude curve on the organ map: a passive follower of a shared driver has relative amplitude cos(lag); tissues off the curve are candidates for an own oscillator | GTEx map, atlases | amplitude also reflects composition | [reg] |

### Phase B. Support for plan 2

| # | task | data | main risk | tag |
|---|---|---|---|---|
| B1 | **Cross-twin pairing.** Twin A's fat against twin B's skin, compared with the within-person 0.533: equal means the shared state is the hour of the visit, lower means a person-specific internal time | MuTHER (E-TABM-1140) | pair identifiers not yet confirmed in the public sample sheet | [reg] |
| B2 | **Organ clock or cell-type clock.** Partition of offsets into organ and cell-type terms over the single-cell atlases; power scout before any residual-sharing test across organs | atlases on disk; cross-tissue immune atlas | few donors with several organs | [reg] |
| B3 | **What moves a human tissue clock, in hours.** Meal-timing and feeding-window interventions with biopsies (Wehrens 2017, Lundell 2020) next to the night-shift result; paired shifts only | check `repeat_biopsy/` first | one tissue per person | [reg] |
| B4 | Cell composition by internal time: one registered test | lung atlas | composition confounded with sampling site | [reg] |

Brain death: replication only (kidney transplant biopsies from living, brain-dead and circulatory-death donors,
GSE43974, coherence only; donor-type multi-organ atlases if found). No mechanism amendment.

### Phase C. Decisions and outside data

| # | item | recommendation | tag |
|---|---|---|---|
| C1 | In-house multi-organ jet-lag data | optional, not a requirement (user, 2026-11-02); it would raise plan 2 by a tier and is the only route to dynamical modelling | [PI] |
| C2 | STARNET (dbGaP; about 600 living patients, seven tissues, Franzen 2016) | apply early (access takes weeks): a living organ map with liver and arteries, and a second resource for T2 | [ask] |
| C3 | TwinsUK visit times | ask: a timed validation in 680 women, and the external-day term of T2 in living people | [ask] |
| C4 | Heart cohort GSE260643 | do not map now (one tissue per person, arrest and reperfusion act on PER1 and NR1D1, processed-only rule, 66 GB left); ask the authors for counts | [ask] |
| C5 | Unspliced counts for the fetal adrenal: pre-mRNA leading mRNA by a consistent angle is the one time-free test of a running oscillator | only if released by the authors; mapping breaks the processed-only rule | [PI] |
| C6 | Format and journal once G2 and G3 are known; fetal result in main text or supplementary | decide after A2 and the gates | [PI] |
| C7 | New atlases: HeartMap (2026) as a scout for within-heart regions, noting that non-failing controls are mostly organ donors; basal ganglia atlas (2026) only if the human part has about 20 donors with several nuclei each | scout, not headline | [PI] |

### Not pursued

Single-cell circadian time mapping on public data (separate methods project), chronotherapy mechanisms (kept as
discussion), the age advance in the main text, cross-species conservation of offsets.

### Order of work

1. **Now:** T0 (scouts) and A1, with the A91 rerun of A5. T0 decides which of T1 to T3 are worth registering (gate
   G1); A1 decides whether plan 2 has a headline (gate G3).
2. **Next:** register and score the designs that pass G1, and T4, which serves plans 1 and 2. A2 and A3 run beside
   them.
3. **Then:** gate G2. If it passes, T5 and T6 and the remaining phase-A items complete plan 1. If it fails, phase B
   and the rest of phase A complete plan 2; if G3 failed as well, A4 and A6 complete plan 3.
4. **At any time:** C2 has the longest lead time and should be opened first; C3 and C4 do not block anything.
