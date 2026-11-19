# Questions for the PI (open first, answered below)

Each entry: date, question, what is on hold until it is answered, the evidence (with `src` pointers in
`results_preliminary.md`). Work that does not depend on an answer continues. When answered, move the entry to
"Answered" with the answer and its date.

## Open

### Still with the PI (updated 2026-11-13)
- **New 2026-11-13:** results A245-A256 and what they change are summarised in `REVIEW_2026-11-10.md` §6; decisions asked there:
  the luminal A row after A224, whether to request ages (Wu 2020) and TwinsUK phenotypes (MuTHER), whether A243 is run (Q32),
  and the wording of the cancer chapter (Q24) after the synthesis in the ledger addendum of 2026-11-12.
- **Q42** (below, new 2026-11-19): section 6 two-axis checks scored or closed (A328-A334b); decide S1 timed control, A333 power run, D3 remainder, cell-mechanism follow-up.
- **Q41** (below, new 2026-11-19): Cell Article with brain death as the climax, or Cell Resource? Brain-death checks done (A325, A327, Task 2-3, B1, B4), lymphocyte bound (A326).
- **Q39** (below, new 2026-11-17): tally of the cell-type program (blocks A-C, A276-A283): what the cell-type line now claims, and which of B2, C1, C2 to run next?
- **Q38** (below, new 2026-11-16): two side projects scouted (clock and cell cycle; liver zonation): start either, and with which data?
- **Q37** (below, new 2026-11-16): add the ROSMAP single-nucleus cortex (427 donors) and RADC time of death to the Q12 request, for a human calibration of cell-type offsets?
- **Q35** (below, new 2026-11-16): tally of the biology-paper checks: two exceptions fell (fetal adrenal first; brain as a clock statement), the cancer and robustness lines stand, no systemic signal named. Biology paper or the Nature Methods ordering?
- **Q34** (below, new 2026-11-15): gate G2 (iii) after A262: no second design reproduces A261; blood forced desynchrony is neither central nor behaviour (damped 0.36, shifted 4 h); night-shift PBMC uninformative. How to use the muscle / blood contrast?
- **Q33** (below, new 2026-11-15): gate G3 fails as registered (A259): the bound in hours is not a headline. Which plan does the paper take while phase T continues? **Answered (c)**; update: no muscle headline (clock-core design not estimable).
- **Q32** (below): three human fetal datasets give no adult-like clock organisation before birth (A237-A239): how the development section reads. Updated with the size of the clock the negatives exclude (A240), A241 and **the replication of the fetal adrenal clock on open data (A242, PASS): option (iii) is done, the recommendation changes**.
- **Q31** (below): the unsigned gene coherence R reads a shift common to the clock genes as a clock; replace it by the signed coherence paper-wide?
- **Q30** (below): the mouse embryo timing positive (A108) is not clock-specific (A236): how the paper words it.
- **Q29** (below): which confirmation first for the Q25 headline: timed brain cohorts, living skin cohort, jet-lag data.
- **Q28** (below): brain death: one descriptive table of single-organ clock integrity; full-n number and caveats.
- **Q27** (below): contrasts "at a common spread" (d*) are biased at unequal host spread; simulation and reporting.
- **Q26** (below): A130 sign label looks inverted; A234 finds an age advance in GTEx.
- **Q25** (below): the headline: two kinds of desynchrony; robustness audit; human local (IBD) test.
- **Q24** (below): the cancer chapter after A214-A226: structure and claims for the paper.
- **Q12** (below): brain time-of-death request. The PI reviews and sends the draft themselves; nothing is sent by the agent.
- **A119:** lung T-cell procurement scout (outward-facing / procurement).
- **Outcome analyses (C8):** scored as A191 on the user's instruction (2026-11-02): both primaries fail (HR 0.98 and 0.96 per SD); draft in
  `paper/exploratory/C8_internal_time_outcomes/DRAFT_C8.md`; results in `results_preliminary.md` (A191). PI to review.
- **Replication of published clock-survival claims (C8b, A192-A194; user's instruction 2026-11-02):**
  - Cadenas genes: carried by proliferation, with RORC as a partial exception.
  - TimeTeller Theta: not prognostic.
  - CYCLOPS magnitude (Li et al. 2024): not prognostic (HR 0.92 top tertile), not reproduced in TCGA with our refit.
  - Clock x treatment (A198): no interaction with chemotherapy or endocrine therapy beyond proliferation (MDE 1.19 / 1.27).
  - Outward-facing use (sharing with these groups, wording in the paper) is the PI's call.
- **Phase 8** (naming, scorecard, figure regeneration, hardening, tutorials, METHOD/CHANGELOG update): PI.
- **Closed 2026-11-01:** Q21, Q22, A136 (registered). **Closed 2026-10-31:** Q8, Q13, Q16, Q17, Q18, Q19, Q20, and the
  A121 and A135 placements.

### Q42 (2026-11-19): the two-axis checks of PAPER_PLAN section 6 are run; four decisions remain
**Question.** The section 6 checks (prompt `prompts/two_axis_checks_2026-11-19.md`) and the user's follow-up on the cell
mechanism (A334, A334b) are scored, closed at power, or registered and waiting for data. The claim set stays frozen.
Decide:
- **(a) S1 timed control:** add the Duan 2023 skin set to the server job? It needs FASTQ mapping, outside the approved
  exception.
- **(b) N1:** answered by the user (run): not decidable at this design (median MDE 0.726 against 0.357).
- **(c) D3:** answered by the user (run, A335 part 5): not estimable (no tissue pair with 20 shared donors).
- **(d) Cell mechanism:** search for deep plate-based data with a same-tissue clock knockout, or keep "quieter" without a
  cell mechanism?

**On hold:** Figure 3e (S1), Figure 5 supplement (D3).

| check | amendment | verdict | key numbers (source CSV under `paper/exploratory/`) | figure / wording it changes |
|---|---|---|---|---|
| T1-T5 assembly | - | done | 21 conditions; Delta muscle 0.79, tumours 0.37, ICU 0.18, lesions 0.19 (`two_axis_map/output/two_axis_table.csv`); registry of 219 amendments (`docs/triage/registry_2026-11-19.tsv`) | Figure 5a map; "amplitude to a quarter" corrected to "excursion to about a third" (A257) |
| S1 fetal running loop | A328 | registered blind, waits for server counts | none yet | Figure 3e stays "the open test" |
| G1 organ factor across species | A330 | not decidable at this design | power 0.604 at true rho 0.6, 1.000 at rho 1; mouse 0.455 at 8 organs (`G1_organ_gain_species/output/power.csv`) | Figure 4g: "organ" stays open; no cross-species claim |
| P1 person loudness, living | A331 | PASS (blind, powered); not beyond responsive genes | leverage 0.177 (twin-pair 90 % 0.119-0.236), p 0.0005; matched q95 0.117 (98.5th pct), responsive q95 0.179 (93.5th) (`P1_person_loudness/output/tests.csv`) | Figure 4f holds in living women, with the responsive caveat |
| N2 departure on the core loop | A329 | partly leakage | core excess 0.250 (90 % 0.143-0.361) vs full 0.969; coupling -0.076 (matched -0.001); signal-state R^2 excess -0.009 (-0.013 to 0.022) (`A329_core_departure/output/tests.csv`) | model part 6: "departure carried mainly through signal-inducible clock genes, cause open"; no core-loop coupling readout |
| N1 output retention | A333 | not decidable at this design (power stage run 2026-11-19 on the user's instruction) | median MDE log Phi 0.726 (loss to 0.48 detectable; gate 0.357) (`N1_output_retention/output/power.csv`); archive A91 pooled lock, not blind | Figures 1 and 7: N1 stays an open question |
| D1 skin loudness | A332 | descriptive | psoriasis r +0.499 (94th / 90.5th pct of companions) (`A332_descriptive/output/d1_series.csv`) | section 3, Figure 6d: loudness change not clock-specific |
| D2 brain death per organ | A332 | descriptive | organ amplitudes 0.92-1.13; lung -1.58 h (90 % -3.11 to -0.14), tibial artery +1.09 h, tibial nerve +0.99 h (`A317_braindeath_confirm/output/d2_per_organ.csv`) | Figure 5b arrows; "rearranged, not quieter" |
| D3 ventilation time | A332 | organ-donor part only | loudness falls with ventilation time while the order is already lost (RP A332) | Figure 5 supplement; remainder (c) |
| D4 developmental curve | A332 | descriptive | hindbrain birth-12 y logV 1.800 (99.5th pct), C 0.772; kidney >= 12 wpc logV 1.865 but C -0.115 (`two_axis_map/output/d4_development.csv`) | Figure 3c; V read with C |
| Cell mechanism (user) | A334 | calibration gate FAIL | 8 h mixing: rho ratio 0.489 (expected 0.500), share_R 0.591 (90 % 0.384-0.799; needed >= 0.75); ZT heterogeneity p < 1e-6 (`A334_dephase_vs_damp/output/calibration.csv`, `calib_mix.csv`) | section 3: desynchrony below the resolution stays possible; not decidable with public mouse 10x data |
| Cell mechanism, aorta D sign (user) | A334b | clock check FAIL: not decidable, not scored | smooth-muscle D aligned 3.14, Bmal1 KO 8.42 (2.68x); MDE dlog D 0.92 vs expected damping -0.47 (`A334_dephase_vs_damp/output/b_power.csv`) | as above; decisive design: same-tissue clock knockout with D near 0, plate-based depth (d) |

**Recommendation.**
- (a) Yes, if the server budget allows: without the timed control, S1's lead can only be read against the postnatal
  reference.
- (c) Only if Figure 5 keeps a supplement.
- (d) Keep "quieter", and scout plate-based datasets before any new registration.

### Q41 (2026-11-19): Cell Article with brain death as the climax, or Cell Resource? (after checks B1-B4, L1)
**Update 2026-11-19 (A335, after the independent review `reviews/plan_review_2026-11-19.md`).** PAPER_PLAN_2026-11-19 already moved brain death to one panel of Figure 5. The review's condition for the Article (night-shift loss and psoriatic shift surviving the core loop) is met: night shift s 4.06 h, R 0.239 on the core loop; psoriasis core q 0.049 (at the line). Also: ICU loudness loss robust (core, age / sex matched), time-of-death reading holds with a held-out waveform; living person loudness (A331) does not survive the core loop. Proposed answer: (a') Cell Article on the two-axis frame, brain death one panel.

**Question.** The public-data checks of the brain-death claim are done (prompt `prompts/cell_route_2026-11-19.md`,
REVIEW_2026-11-19 sections 2.1, 5, 7). Which route: (a) Cell Article with brain death as the climax, (b) Cell Resource
(method, specificity standard, atlas of clinical datasets), (c) Article only after new procurement data? **On hold:**
the dataset sweep and Figure 1 of the Resource route (not started; inclusion rule to be registered first), figures.

**Brain-death evidence (GTEx only; every test NOT BLIND, the contrast was found post hoc in A132b).**
| check | result | reading |
|---|---|---|
| A317 H317.1 | D_rho 0.621, relabelling p 0.0025, above all 200 matched and 200 responsive sets | organ donors' arrangement differs from post-mortem ventilated donors, clock-specifically |
| A317 H317.2 | excess C 0.269 vs 0.317 (FAIL) | organs of brain-dead donors still share a clock state: rearranged, not dissolved |
| A318 H318.1 | D_rho 0.584 after linear ischaemic time + RIN adjustment (PASS) | not a linear ischaemia / RNA-quality effect (groups barely overlap: 258 vs 561 min) |
| A318 H318.2 | organ donors with longer ischaemia keep more of the arrangement, +0.410, p 0.040 | gradient in the collection-mode direction |
| Task 2 (A318 post hoc) | common ischaemic range: 0.239 on 36 pairs; noise-matched 0.244-0.473 (q05-q95), R36 at the 5th percentile; organ donors' rho rises beyond noise (0.184 vs <= 0.066); dropped pairs = every pair with heart, lung, blood | borderline; the fall is organ-donor samples with post-mortem-range ischaemia looking post-mortem (H318.2 again) |
| Task 3 (living surgical reference) | precondition: no shared pair's organ-donor minus post-mortem offset beyond the site bound (1.22, -0.01, -0.03, 1.68 h; all intervals contain 0) | cannot decide; stopped |
| B1 (baboon referee) | fast deaths vs baboon over 6 tissues: absolute rho 0.31, p 0.24 | the baboon cannot referee; stopped |
| B2 (A325) | post-mortem ischaemic tertiles (555-585 min apart, hour-stratified) keep the arrangement: H325.2 rho 0.862, p 0.62, power 0.99 | the post-mortem side carries no collection component of the organ-donor size |
| B3 (A327) | D_rho leg-only pairs 0.244, possibly perfused 0.265, mixed 0.456; no perfusion dependence (p 0.66); MDE 1.0 | consistent with brain death in unperfused tissue, underpowered |
| B4 | no GTEx field for donation after circulatory death; all 298 organ donors brain dead with a cross-clamp time | the direct control does not exist in GTEx |
| Literature | McTiernan 2021: brain-dead donor hearts lose PER2 / NFIL3 rhythm and shift ARNTL / NR1D1 against failing hearts (single organ, pseudo-time) | first on the order between organs, with a single-organ precedent |

**Verdict on the Article route (recommendation).** The effect is large, clock-specific and survives everything that can
be tested on the post-mortem side (A318 H318.1, A325). What remains is on the organ-donor side: within organ donors the
arrangement moves toward the post-mortem one with longer ischaemia (H318.2, Task 2), and procurement itself (perfusion,
recovery order, anaesthesia) has no counterpart in GTEx (B4); the checks that could separate it are uninformative or
underpowered (Task 3, B1, B3). Public data cannot settle it. As a climax of a Cell Article it would be the least
secure claim carrying the most weight; reviewers will raise the procurement alias first. **Recommended: not the
climax without new data**; keep it as an open, provocative question ("In GTEx, brain-dead organ donors carry a different
organ arrangement from every post-mortem group; brain death and collection at procurement cannot be separated in this
resource") in either route. With new data (below) the Article route reopens.

**Lymphocytes (L1, A326).** The human designs read a lymphocyte clock only above about 0.7 of the fibroblast amplitude;
a mouse-colon-sized T-cell clock (0.25) would not have been seen. "No readable clock" stays, as a bound, not a species
difference: it is not a provocative result either.

**Resource route (one-page sketch; not started).**
- *Message.* What untimed hospital-hours tissue can and cannot say about human clocks, with a method (clocksync), a
  specificity standard (expression-matched and stress-responsive companions, calibrated: A283, A285, A320) and an atlas
  of public clinical datasets with at least two compartments per person.
- *Figure 1, "what hospital-hours tissue can and cannot say".* Can: one internal time per person across tissues in
  living surgical cohorts (A251-A256); cell classes carrying it (A272-A283); organ relations in GTEx deaths (A293).
  Cannot: per-sample time; lymphocyte clocks weaker than 0.7 of fibroblasts' (A326); brain death vs procurement (above).
  Published claims our checks contradict, and only these: COFE's tumour uncoupling (A265, A267) and the published
  prognostic clock signatures (A191-A194, A220). A123b is a method failing in one sampling regime (untimed orderers in
  08-18 h deaths), not a contradicted claim; A207 shows disruption scores coexisting with tumour-normal pairing. The
  edge is real hospital-hours data (A123b); in simulation A145 found none after Holm.
- *Atlas.* A sweep of public datasets with >= 2 compartments per person (tumour-normal, multi-site surgery, multi-organ
  donors, single-cell atlases): inclusion rule and denominator registered before any run; per dataset the gate, the
  couplings / offsets with companions, and "uninformative" where the gate fails. Outward-facing: the public atlas holds
  GTEx aggregates only.
- *Biology in the Resource.* Cell-class participation, conserved class gains, the gain-over-departure warning
  (coupling mostly measures loudness), chronic-disease nulls, psoriasis, tumours keep the patient's time.
- *Cost.* API, tests, tutorials (planned anyway); the sweep (weeks); no new biology claims beyond the freeze.

**Outward-facing options (drafts only on request; ask first).** (1) Multi-organ procurement biopsies of donors after
brain death and after circulatory death (Swiss transplant centres; UK QUOD biobank, scope to check): the direct test.
(2) Archived multi-organ RNA of rat brain-death models with ventilated sham controls (no public dataset found).
(3) The in-house jet-lag data (Q29): if they become available, register before opening an ordinal predictor: organs
with a higher feeding weight b (A322, A324) lag the SCN longer after a light shift and balanced-weight organs dip most
in amplitude.

**Also for the PI.** The i-TIMES pooled analysis (eight Roche trials, infusion time not randomised; overall survival HR
1.039, 0.925-1.168, matched cohort) is carded from a secondary report only; the LungTIME-C01 retraction (Nat Med,
2026-06-24) is for design and data integrity: cite as fragile clinical evidence, never as a hospital-hours false
positive.

### Q40 (2026-11-18): where the paper stands after the organisation framework, and how the cancer chapter reads in it
**Before / after, 2026-11-15 to 2026-11-18** (STATUS_2026-11-18 sections 1-3, updated today).
- *Before (11-15/16):* no replicated headline for "what internal time controls" (Q33, Q34); the biology thesis had lost
  two exceptions (brain as a clock statement, fetal adrenal first; Q35 recommended the Nature Methods ordering); the
  amplitude-coupling relation was unexplained; atlas couplings failed the responsive control; death-type results were
  archive results without companions.
- *After (11-18):* a coherent account of how the clock is organised in the body. Tissues express one shared time with
  9-fold different gain (shared variance 85-fold) over a clock-specific departure of nearly constant size (A291, A294); within organs
  stromal-vascular, epithelial and parenchymal cells and macrophages carry it, lymphocytes do not (A272-A288);
  composition changes gain, not phase (A295); the organ arrangement holds without alignment to the day or behavioural
  cues and is lost in brain death (A293); peripheral clocks are self-sustained (mouse re-analysis, A297).
- *Fell or stayed open:* adrenal +4.4 h offset (absent in mouse and sheep, A303: an effect of death or of species),
  endothelium-first order (A293), meal-driven shared departures (A302), cortex cell-class coherence (A280).
**Claim tiers.** Tier 1 (living, replicated, clock-specific): shared internal time advancing with elapsed hours;
cell-type participation; one phase per organ; conserved relative gains. Tier 2 (GTEx deaths, clock-specific): gain over
a constant departure; organ arrangement independent of daily alignment. Supported with care: gain not composition
alone; offsets equal across sex, age and metabolic disease; brain-body two domains (descriptive).
**Today's scoring (A304-A306).**
- *Cancer, A304 (TCGA paired, 517 pairs, powered):* proliferating tumours express less of the patient's time (gain ratio
  0.63, p 0.0055, Holm 0.011, below the matched and responsive companions) and do not drift further from it (departure
  rise not detected, power 1.000 at a doubling). TP53 and luminal A: no departure rise (powered, FAIL); ER: not
  powered. Purity does not explain it (the proliferative arm has more stroma, which raises gain).
  Tumour coupling to the body's time, disattenuated by GTEx reliability (descriptive): lung adenocarcinoma 0.54,
  liver 0.38, kidney clear-cell 0.27, prostate 0.26, breast 0.18, papillary kidney 0.12, thyroid 0.11, colon 0.06.
- *Cancer single cell, A305:* not powered (2 cohorts, 3 malignant-stroma pairs; MDE 0.156 > 0.10); closed.
- *Disease, A306 (GTEx, powered for 5 readouts):* T2D, hypertension and smoking change neither the gain nor the
  departure of their organ relative to the rest of the body (MDE 0.26-0.35 on log z2 or log psi; smoking departure
  0.42). With d* and offsets already equal (A164, A233, A235), chronic disease leaves the organ-level organisation
  intact. Lung gain in smokers is lower (-0.51, p 0.068), not powered.
**Recommendation for the cancer chapter (main text).** One framework, three statements:
(1) the tumour keeps the patient's time (A138; COFE's opposite is a sampling artefact, A265, A267);
(2) within the tumour, the participation rule of healthy organs holds: fibroblasts and endothelium carry it clock-
specifically, immune cells share only a responsive state (A227, A271);
(3) as tumours proliferate they express less of that time, without drifting from it (A304): gain falls, the departure
does not rise. This replaces "loss of coherence" (A127-A173, tumour-only, which could not tell damping from rewiring,
A212) with a statement on the body's time, and it is the organ-level analogue of lymphocyte dilution in healthy tissue.
Whether malignant cells carry a weaker clock or simply carry less of it is not resolvable with public single-cell data
(A214, A305). Per-cancer ranking (A3) as an Extended Data panel, descriptive, with the GTEx-reliability caveat.
The disease section reads: clock organisation is robust to chronic disease (couplings, offsets, gains, departures)
and lost in brain death; cancer is the one condition where the organ's gain falls, with proliferation.
**Open decisions (user / PI).** (a) Venue: Cell (biology case now real) vs Nature Methods fallback. (b) Spliced /
unspliced fetal counts: FASTQ mapping is outside the processed-data rule; rule change or an author request (draft only).
(c) Q37 (time of death for brain cohorts / human known-time truth at cell-type level). (d) A290 (the other session's calibration of the amplitude-coupling relation against responsive genes):
rescheduled by its owner as a weekend run (deviation registered, commit 9a6479e); nothing needed from the PI. (e) Is the smoking lung gain worth a powered design (larger smoking cohort)?
**Evidence.** RP A304, A305 power note, A306; STATUS_2026-11-18; VERDICT_TABLE rows A304-A306.

### Q39 (2026-11-17): tally of the cell-type program (blocks A-C)
**Scored (all on the calibrated test of A283 where a companion test is used).**
- Test calibration: q95 rule anti-conservative (A276); percentile against a 5000-draw placebo null calibrated (V283
  PASS); every A272/A274 verdict unchanged (H272.1 still fails against responsive genes; H272.2, H274.1-3 pass).
- Stromal-vascular over lymphocyte (D): passes in three independent body sources (A272 0.277, A274 0.152, IBDverse
  A280 0.413) and in living donors alone (A282: 0.225, p_cal 0.011 / 0.011; post-mortem 0.261).
- Lymphocytes: no readable clock of their own (A281: coherent in 2 of 9 atlases; fail where depth-matched fibroblasts
  pass in 5 of 7), and in living tissue they share responsive genes more than clock genes (A282).
- Untimed offsets unbiased against known time in two mouse tissues (A273 colon, A279 aorta); aorta macrophages not
  less coupled than stromal-vascular cells (H279.1 FAIL).
- Cortex: pooled sharing passes over HBCC, DLPFC, MTG (H283.5) but not in MSSM, the largest cohort (H280.2 FAIL); no
  cohort beats the matched family alone; offset sign unresolved (H280.3 not read).
**Not scored.** B2 residency (design ready from A280/A281 tables); C1 disease contrasts (inventory only); C2 tumour
pairwise (same patients as A149 S3, power not computed); D1-D3 scouts (Q37, Q38).
**Question.** (1) Is the claim for the paper "stromal-vascular cells carry the shared tissue clock, lymphocytes neither
share it nor keep their own; cortex unresolved"? (2) Which next: B2 (residency: do tissue-resident lymphocytes differ),
C1 (disease contrast, with power first), or C2? (3) Move the calibrated pairwise readout into `clocksync` (the program
says: if A1 confirms calibration)?
**Evidence.** RP A276-A283; VERDICT_TABLE rows A272-A283.

### Q38 (2026-11-16): clock and cell cycle (D2), liver zonation (D3): scouts only
**D2. Is the distribution of cell-cycle phases of a donor's cells related to the donor's clock state?** CoPhaser is
used for cell-cycle phase only, at median UMI per cell above 1.5k (memory rule). Mean UMI per cell (pseudobulk depth /
cells; the median per cell is not in the pseudobulk tables) and cells per donor:
- Kong gut epithelium (Crohn's study, healthy and disease donors): colon 34 donors, 5.8k UMI per cell, median 1,686
  epithelial cells per donor; terminal ileum 42 donors, 5.5k, 1,524. Best candidate.
- Colorectal (colon_crc, 62 donors, epithelial 22.8k UMI per cell, 1,870 cells per donor; tumour and normal), KPMP
  kidney (77, 5.0k, 1,114), muscle (32), breast (126, 933 cells per donor).
- T cells: 1.3k to 5.8k UMI per cell; the Kong epithelial fractions carry T cells at 5.3k to 5.8k.
- Power (rough): with 34 to 42 donors a circular-linear relation of about |r| 0.45 is detectable at 80 % (alpha 0.05).
  Needs per-cell counts (Census download per atlas, size not checked) and the CoPhaser runs with their validation
  plots; uses the lab's unpublished tool.
**D3. Do pericentral and periportal hepatocytes of a person share one clock phase?** Healthy liver
(`liver_xu`, Census harmonised liver, 36 donors): donors with both a centrilobular and a periportal hepatocyte
compartment of >= 30 cells: 11 (5 / 3 / 3 by assay), so 8 degrees of freedom after strata, below the project's minimum
of 15 for a pair; pericentral and periportal sinusoidal endothelium: 14 donors (>= 30 cells). Larger single-nucleus
cohorts are mostly MASLD: GSE202379 (47 patients, 2 healthy; per-cell counts exported), GSE289173 (48 donors),
GSE256398 (6 healthy of 17). These have no zone labels: a per-cell zonation score from marker genes fixed in advance
(e.g. CYP2E1, GLUL, CYP1A2 pericentral; HAL, SDS, CPS1, ASS1 periportal, Droin 2021 for the mouse) would split
hepatocytes by zone; in MASLD zonation itself changes, which confounds a disease cohort.
**Recommendation.** Neither is needed for the cell-type section; both are projects of their own. D2 with the Kong
epithelium is the more feasible; D3 is underpowered in healthy liver and confounded in MASLD.
**On hold.** Nothing.

### Q37 (2026-11-16): human truth with known time at cell-type level (cell-type program, D1)
**Question.** Extend the unsent Q12 request (`exploratory/A141_brain_stability/drafts/tod_request_draft.md`, updated
2026-11-16) to the ROSMAP single-nucleus prefrontal cortex (Mathys et al. 2023, 427 participants, Synapse syn52293417
under a Data Use Certificate) with the hour of death from the RADC hub (recorded for ROS and MAP donors, used by Lim et
al. 2017)? Send, edit, or not?
**What it would give.** The first human check of untimed cell-type offsets against known time (A273 is mouse colon
only); the tie-break of the cortex offset sign (PsychAD +1.3 h against SEA-AD -1.5 h for the neuron-vascular pair,
A274); whether the shared cortical cell-class state reads time of death. Brain cohorts on disk (`$SCCC_DATA/brain_tod`,
bulk) cannot give cell classes.
**On hold.** Nothing; outward-facing, so drafts only. Not verified: whether the time-of-death field is released through
the hub for this use (RADC to confirm), and the overlap with `psychad_radc` donors.

### Q35 (2026-11-16): tally of the biology-paper checks (stop rule: two exceptions fell)

**Thesis tested** (`prompts/biology_paper_checks_2026-11-15.md`): the human clock network is robust, and the places
where it breaks are specific: brain, blood, the fetus except the adrenal cortex, malignant cells against their
microenvironment, brain-dead donors. Rule of the prompt: if two or more exceptions fall, the same material is reordered
into a Nature Methods paper.

| exception / line | check | outcome | status |
|---|---|---|---|
| Fetus absent except the adrenal cortex | 1a, 1b: A264, A266 (adrenal); A243, A269 (heart, kidney, lung, liver) | adrenal not clock-specific without PER1 / NR1D1, core loop not significant (C 0.417, p 0.057); lung (C 0.530) and liver (0.583) carry the same covariation, heart excluded at 0.5 x adult; no core loop clock-specific in any organ; not explained by steroid or glucocorticoid activity | **fell** ("adrenal first" and "absent elsewhere" both gone; a weaker multi-organ observation remains, supplement) |
| Brain | 3: A263 Part A (human, baboon, sheep) | brain's lower coupling to the body (human L -0.514, baboon -0.528, sheep -0.707) lies at the 30th-88th percentile of matched and responsive gene sets; companions drop more; clock r* 0.803 vs companion medians 0.42 / 0.38 | **fell as a clock statement** (erratum A142: "brain is less coupled to the body's shared state") |
| Malignant cells vs microenvironment / tumour-normal pairing | 2: A265, A267 | COFE's "uncoupled" reading not reproduced (uniform in 3 of 10 cancers) and uninformative in a simulation on the real design (7 of 10 coupled replicates read "uncoupled"); COFE's own pairing is patient similarity (not above matched sets, not above generic axes); clocksync 0.361 (p 0.0001), 10/10 coupled and 0/10 uncoupled replicates | **stands** (the head-to-head favours clocksync) |
| The systemic signal (mechanism) | 4: A268 | glucocorticoid, cold and heat readers consistent across GTEx tissues (K 0.70-0.79) but not above matched / responsive readers, and at a phase 4-7 h away in living MuTHER tissue | **not named** (clean negative at the stated power) |
| Blood; brain-dead donors | not checked (blood an expected exception; brain death descriptive by decision) | - | unchanged |

**Update 2026-11-16 (run after the stop rule at the user's request: "Still go for Part B"; A270 and A271 registered on
the user's go-ahead).**
- *A263 Part B (H263.4 FAIL):* atlas couplings 0.393, at the 87.5th percentile of the responsive sets (q95 0.423); only
  fibroblast / stromal (and "other") clear both families. Sheep (H263.2) fails too.
- *A270:* in GTEx the responsive sets do not read time of death (17 % of sets; clock R 0.631 above all 400 sets); in
  baboon most sets of either family read ZT. The responsive failures are not explained away as "daily hormone program"
  in human deaths.
- *A271 Part A (H271.1 FAIL):* the tumour-normal myeloid, T/NK and CD8-negative T pairings (A227, A228) are not above
  the responsive sets; fibroblasts and endothelium are.
- *A271 Part B (H271.2 PASS):* "stronger clocks are more coupled" (A104, periphery rho 0.734) is clock-specific: the
  companion sets give rho near 0 (q95 0.258 matched, 0.453 responsive). In GTEx bulk tissues the clock coupling (median
  0.675) is far above both families (0.107, 0.267), unlike the atlas cell types (descriptive).
- *Consequence (post hoc, user and agent):* coupling size shows a shared state, not the clock. Couplings are reported
  as the excess over each cell type's own companions; clock identity rests on phase-based readouts (host phase vs
  recorded time, offsets vs timed truth, r*), which the clock passes far above every companion. The clock-specific
  backbone narrows to fibroblast / stromal cells (both families) in health and in tumours. Figure 3 (H109.1-R) and the
  TCGA pairing (H138.4-R) already passed the responsive family (2026-10-28). This strengthens the Nature Methods
  ordering: the companion discipline is itself a result, and the biology claims that survive are the stromal backbone
  and the timed readouts.

**What stands.** One clock state per living person, the state is time (living tissues, elapsed hours), offsets
recovered from hospital-hours samples, group offsets within an hour, tumour-normal pairing; all with companions where
registered. **What is gone as a clock claim.** The brain exception and the fetal adrenal "first". **Method results
added by these checks.** The companion discipline itself removed two exceptions that a coupling statistic alone would
have kept; in the head-to-head, the only competing untimed method's published conclusion fails on the same data.

**Recommendation: the Nature Methods ordering**, as the prompt's rule says (two fell), with the biology kept as the
application chapters: (1) the method and its calibration (timed recovery, hospital hours, the COFE head-to-head and the
simulation on the real design); (2) one clock state per living person as the flagship application; (3) where the
shared state is weaker and what is not clock (brain, fetus: shared-state statements with companions, as cautionary
results); (4) cancer pairing. The biology paper would rest on robustness alone, with no specific breaks and no named
signal, which does not meet Cell's "conceptual advance" bar. If the PI prefers to keep the biology framing, Part B
(H263.4) should be run first: a failure there would also remove the organ-map couplings.

**On hold:** the plan and journal choice (user and PI); any rewording of roadmap section 2 and of Figure 2.
**Evidence:** RP A243, A263 (Parts A and B), A264-A271; VERDICT_TABLE rows of the same; roadmap progress note 2026-11-16.

### Q34 (2026-11-15): gate G2 (iii) after A262: no second design reproduces A261, and blood differs from muscle

**Question.** G2 (iii) asked for a second misalignment design. No muscle or adipose design with shifted meals or sleep
has public matrices (scout T5). The blood designs give the following:
- **Forced desynchrony (GSE48113, 21 subjects):** the genome-wide rhythmic program does not follow the displaced sleep.
  It is damped to 0.355 (0.203 to 0.621) and displaced 3.96 h (2.23 to 5.69). Central time and behaviour are both
  excluded (p < 0.0001); "lost" is not excluded at p 0.0101, against the level 0.01.
- **Simulated night shift (GSE107537, 8 subjects):** uninformative (rhythm gate not passed, 0.245 vs 0.254).

The registered condition for "which time a tissue's outputs follow is tissue-specific" is not met. Which of these does
the paper take?
- (a) Report muscle (A261) and blood (A262 D1) side by side as a contrast: muscle follows behaviour, blood is neither
  on central time nor on behaviour and loses its excursion. State that neither is replicated and that blood is
  composition.
- (b) Treat G2 (iii) as failed, so plan 1 has no replicated headline, and fall back to plan 2 or plan 3 for the main
  text.
- (c) Seek a replication of A261 in muscle or adipose. That needs data that are not public (Wilms 2019 adipose has no
  accession; contacting the authors is outward-facing, a draft only, your call).

**On hold:** the plan choice (with Q33), any wording of a muscle versus blood contrast, and contacting data owners.
**Evidence.** RP A261 and A262; scouts T1b (the clock-core design is not estimable) and T5;
`literature_positioning.md` section 11.
**Recommendation: (a) for a results section, without headline status.** Combined with the literature, it is
consistent: muscle outputs follow behaviour (Wefers 2018), blood rhythms are lost or stay put (Archer 2014, Kervezee
2018). What is new is the same in-hours program machinery applied to both designs against explicit times. It does not
meet Cell's "conceptual advance" bar alone, so plan 1 stays below headline until T4 is in. Plan 3 stays the floor.

### Q33 (2026-11-15): gate G3 after A259, and where plan 1 stands

**Update 2026-11-15 (after review; item 2 of the follow-up): plan 1 has no headline in muscle from GSE106800.**
Decision recorded: option (c); the plan choice waits for G2, plan 3 stays the floor, roadmap 2.2 and the Figure 2
headline are not reworded. The clock-autonomous core design (do direct clock targets stay with the muscle clock while
the rest follows behaviour?) is **not estimable** on these arrays (scout `exploratory/scouts_2026-11-15/T1b_clock_core/`,
matched stand-ins only, no core gene read): a program clears the A261 rhythm gate in 0.157 of draws at 400
reference-rhythmic genes and 0.827 at 600, against 65 genes for the A260-derived core; and a damped clock-following
program (ratio 0.233) is never told from central time (correct call at most 0.010 at any size). So the three layers
cannot be shown to come apart in human muscle with these data. A261 stands as a program-level, in-hours confirmation of
Wefers 2018 with a masking limit ("follows behaviour, shifted or evoked"). What remains for plan 1: whether blood goes
the other way on the same machinery (A262, registered next; Kervezee 2018 and Archer 2014 suggest it does), which
would support "which time a tissue's outputs follow is tissue-specific" (not written until both are scored), and T4.

**Question.** The registered G3 component fails: a tissue's departure from its person's body time is bounded at 1 to
3 h, not below 1 h, and the departure is not shown to be clock. Do we (a) keep plan 2's "one internal time per
person" with the bound restated as "group offsets agree within 0.2 to 0.6 h; within-person tissue scatter below 1 to
3 h, mostly sample-level", (b) move to plan 3 as the roadmap's rule says when G3 fails, or (c) defer the plan choice
until phase T is scored (G2), since plan 1 has a new positive?
**On hold:** the plan, format and journal (C6); rewording of roadmap section 2.2 and of the Figure 2 headline. Phase T
continues (it does not depend on the answer).
**Evidence.**
- A259 (RP; `exploratory/A259_tissue_deviation/output/score/tests.csv`): median upper bound of sigma over 11 GTEx tissues
  2.20 h (0 of 11 below 1 h); per-tissue sigma 1.03 h (thyroid) to 3.22 h (lung); no tissue above both companion
  families (matched sets carry 0.30 to 0.74 of kappa); MuTHER fat-skin 2.20 h (upper 2.85 h). Power: the estimator
  detects 1.0 h with 0.85 at strong clocks.
- A1 (b) (RP): the group offset equivalences of A208, A231 and A235 stand against the spread bias (expected bias <= about
  0.1 h at the observed spreads); A233 has no interval.
- Plan 1 so far: A260 calibrates the untimed readers on mouse liver (AUC 0.997; shift reader r 0.869); A91 reproduces
  on current code (body vs day per gene); **A261: the genome-wide rhythmic program of human muscle follows behaviour
  (shifted or evoked) after a 12 h shift (shift 11.65 h, ratio 0.854; matched sets 0 / 200) while its clock moves
  2.68 h (A257)**. Limit (2026-11-15): with two clock times, a 12 h shift and an acute fed / awake vs fasted / just-woken
  response predict the same data, and the protocol inverts the behavioural state at each biopsy exactly; Wefers 2018
  already reported on these arrays that the muscle clock did not follow the inverted schedule; T3 fails
  G1; the T2 tissue term is not registrable (A259); Hallmark programs are too weak to partition (A261 H261.2).
**Recommendation: (c).** G3 as registered fails, but the reason is informative rather than a failure of the method:
sample-level variation, not clock, sets the within-person bound, and the between-group offsets remain tight. Plan 1 now
has a headline candidate (in muscle the daily transcriptome follows behaviour, the clock does not) that needs G2 (iii),
a second design. Next, independent of the answer: T4 (systemic signal) and a replication of A261 in another
behaviour-shift or meal-timing design if one exists on disk or public, plus companions for A91. Plan 3 stays the floor.

### Q32 (2026-11-11): the human developmental section after A237-A239
- **Update 2026-11-16 (A243 and A269; tally Q35):** pooled open cohorts give a clock-shaped covariation in the fetal
  lung (C 0.530, Holm p 0.026) and liver (0.583, 0.025), between half and three quarters of the adult size as the
  adrenal; heart negative at an adult pass rate of 0.800 at half the adult size, kidney negative (three quarters
  excluded). Both survive a glucocorticoid-target covariate (A269 H269.2); neither passes on the core loop (lung C
  0.473, p 0.060; liver 0.686, p 0.009, below the responsive q95 0.712). Reading: "the clock genes covary before birth
  in the adrenal, lung and liver, beyond glucocorticoid activity; a core loop is not shown in any". Recommendation:
  supplement with one main-text sentence; the one time-free test of a running oscillator is spliced / unspliced counts
  (roadmap C5: author release, outward-facing, or FASTQ mapping, outside the processed-data rule).
- **Update 2026-11-15 (A264, scored; not blind):** the cortex covariation survives a cubic age trend with cohort-specific
  age slopes (C 0.659, p 0.0021) and the ACTH steroidogenic score (C 0.526, p 0.0121), both clock-specific; without PER1
  and NR1D1 it is C 0.489, p 0.0290, above matched sets but below the responsive q95 (0.634): **not clock-specific**. By
  the registered reading the fetal result moves to the supplement as a hypothesis; the verdict turns on the responsive
  family, whose q95 rises with 12 genes (RP A264, limits). Counted as a fallen exception in the tally (biology-paper
  checks, question to follow after Q34). A243 (other organs) is still to be scored.
- **Update 2026-11-15 (A266, user-requested follow-up; not blind):** on the core loop without PER1, PER2 and NR1D1 the
  cortex covariation is not significant (C 0.417, p 0.057; adult pass rate 1.000, 0.995 at 0.75); the relation of the
  responsive genes to the clock state is not estimable (adult power 0.50). The fetal cortex expresses MC2R / MRAP highly
  and no detectable melatonin receptor. Reading: a hormone-driven (ACTH / cAMP, own cortisol) response of the
  signal-induced clock genes is the simpler explanation; a core loop >= 0.75 x adult is excluded. Recommendation: keep
  the fetal adrenal in the supplement as "clock genes respond to the fetal adrenal's own hormonal axis", not as a
  running clock.
- **Update 2026-11-13 (A244, scored):** the fetal adrenal signal is cortical (C 0.608, Holm 0.0076, clock-specific, also
  after foreign-RNA correction); the medulla shows none at adult power 1.000. Below adult size; strongest before 12 weeks.
  Suggested wording: "a clock-shaped covariation between fetuses is present in the steroidogenic adrenal cortex from the
  first trimester, below adult size; none is detected in the medulla". A243 (other organs, pooled) is still registered,
  not run; it needs a scout and up to 60 GB.

- **Question:** three open human datasets were scored on a common registered design (age and sex adjustment, signed
  coherence, both companions): organs of one fetus (Cardoso-Moreira, 4-20 wpc), gonad cell classes (Garcia-Alonso,
  6-21 wpf), brain regions (Kang, 8-37 pcw). No prenatal primary readout passes; the postnatal L1 control passes in
  all six sets (five organs and the brain, n 9-32, unrecorded times); in the brain the coherence across individuals
  is absent before birth (C -0.166, p 0.58; 13-22 pcw 0.047) and present in infancy and childhood (C 0.735, p 0.003)
  and from 13 years (0.663, p 0.007). Does the paper (i) state this as its developmental result ("the clock
  organisation read by the method is not detected before birth in seven organs, four gonad cell classes and the
  brain, and is present from infancy; the fetal adrenal of A189 is the one exception, in 7 fetuses"), (ii) move
  development to the supplement, or (iii) first seek a replication of the adrenal positive, which needs controlled
  access (Kildisiute 2021 and del Valle 2023, EGA: an access request, your decision)?
- **On hold:** any further fetal amendment: the adrenal replication, a follow-up of the two leads (gonad endothelium,
  Holm 0.056 and clock-specific; fetal forebrain with the RNA quality covariate, descriptive C 0.605, p 0.004), an L2
  companion that works when every gene set is shared (Kang: the postnatal gate passes but is not above the
  responsive sets, so H239.1 is uninformative as registered).
- **Evidence:** `results_preliminary.md` (A237-A239); `A237_fetal_organs`, `A238_fetal_gonad`,
  `A239_fetal_brain/output/verdicts.csv`. Main limit: a negative needs clock phases that differ between fetuses as
  they do between adults sampled in hospital hours; clocks entrained tightly to maternal time and collected in a
  narrow window, or of low amplitude, would be missed.
- **My recommendation:** (i), worded as "not detected" with that limit, with the Kang trajectory as the figure panel
  (same platform on both sides of birth) and the Cardoso-Moreira postnatal controls beside the prenatal organs; (iii)
  only if you want the adrenal exception to carry weight.
- **Update 2026-11-11 (A240, A241; nothing new on hold).** (a) *What the negatives exclude* (A240, adult material
  at the fetal designs with the clock scaled): a clock organised as in adults of at least 0.75 (9 of the 18
  negatives) or 1.0 (9) times the adult between-individual variation would have been seen in 8 of 10 such cohorts; a
  clock of half the adult size is detected in 0.000-0.745 of them (median 0.200) and is not excluded. The bound is
  conservative by construction. Adult pass rates are unchanged when the deaths are restricted to 08-14 h or 08-12 h
  (0.865-1.000): a narrow collection window does not explain the negatives. Proposed wording for (i): "no clock of
  adult size and organisation is detected before birth; a clock below about three quarters of the adult one is not
  excluded". (b) *Adrenal as host inside Cao 2020* (A241): fails and is uninformative (adult pass rate 0.140 at 7
  fetuses); it says nothing on which organs follow the adrenal. (c) *Option (iii) needs less than stated above:*
  the processed fetal adrenal object of Kildisiute 2021 is open (7 samples, w8-w21), as are Kameneva 2021
  (GSE147821, 13 samples, weeks 6-14) and Dong 2020 (GSE137804, 6 samples); only the raw reads and del Valle 2023
  are behind EGA (EGAD00001008345). Few samples lie in the age range of the Cao adrenals (12.7-17.4 weeks); a
  replication would be a new amendment with its own power calculation, not registered. Links:
  `scouts_2026-11-10/FETAL_SCOUT.md` (addendum).
- **Update 2026-11-11 (A242: option (iii) done on open data, no access request).** Nine open cohorts, 59 fetal
  adrenals, none from Cao 2020, pooled with a cohort indicator and within-cohort permutations; adult power 0.99-1.00.
  All three registered readouts pass and are clock-specific: from 12 weeks (19 fetuses, 7 cohorts) C 0.552, p 0.022;
  13 bulk adrenals of a ninth cohort C 0.651, p 0.0024 (Holm 0.0048); before 12 weeks (25 fetuses, 5 cohorts) C
  0.732, p 0.0007. No cohort carries the result alone (leave one out: p 0.0099-0.046). Against adult stand-ins the
  fetal coherence lies between half and three quarters of the adult between-individual clock variation; the two
  cohorts scored alone do not pass (Kameneva p 0.246, Wang p 0.071). *What it changes:* the adrenal positive is no
  longer one cohort of 7 fetuses but two registered positives over ten cohorts. *What it does not show:* that the
  adrenal runs first; the negatives of the other organs exclude a clock of 0.75-1 times the adult one (A240) and the
  adrenal signal is itself below that. *My recommendation now:* (i) with the adrenal as a replicated positive: "a
  clock-shaped covariation between fetuses is present in the adrenal from the first trimester (ten cohorts, 66
  fetuses), below the adult size; no clock of adult size is detected before birth in seven other organs, the gonad
  and the brain"; the Kang trajectory stays the figure panel. *For your decision (each a new amendment, nothing
  started):* (1) the same pooled design for another organ with several open fetal cohorts (heart, kidney, lung,
  liver), which is the only way to compare organs at equal power and ask which starts first; (2) cortex against
  medulla within the pooled adrenals, with the ambient-RNA limit of A190; (3) leave-one-cohort-out of the
  before-12-weeks readout (not registered in A242, so post hoc).
- **Update 2026-11-11 (2): (1) and (2) are registered on the user's instruction, as A243 and A244** (staged: a design
  addendum with the cohorts and the adult power is committed before any fetal statistic; not scored). (3) stays
  open. Their outcome will be added here; the wording of the developmental section remains your decision.

### Q31 (2026-11-10): the unsigned gene coherence R and a shift common to the clock genes

- **Question:** the L1 statistic of the paper, the gene coherence R of `clocksync.decompose.components`, is unsigned
  and fits no per-sample offset. A shift common to the clock genes of a sample (depth, composition, RNA quality,
  normalisation) is then read as coherence: in simulation without any clock, a common shift of sd 0.2 (log units)
  passes the per-gene permutation test of R in 0.625 of the data sets at n 10 and in 1.000 at n 40, and the shift
  hides a real clock at small n (pass rate 0.690 at n 10 against 1.000 without shift). A signed coherence with the
  offset fitted out (`clocksync.decompose.signed_coherence`, added with tests) stays at its nominal level (0.050,
  0.035) and keeps full sensitivity (1.000 at every shift). Do we (i) replace R by the signed coherence in every L1
  result of the paper (re-scoring, labelled as an audit, not as new tests), (ii) keep R and add the signed value as a
  companion column where a claim rests on R, or (iii) keep R as it is and state the limit?
- **On hold:** nothing already scored is changed. The new fetal amendments (A237-A239) are registered on the signed
  statistic. An audit amendment (D1b feasibility and the L1 / `d_coherence` results: A109-A111, A126, A127, A134,
  A135, A140, A161, A176, A179, A184, A189, A190, A212, A213, A215, A219, A225, C4 A195, the cancer coherence
  chapter and the figure 4 / 5 decompositions) is not registered before the answer.
- **Evidence:** `scouts_2026-11-10/feasibility_fetal/output/shift_artefact.csv` (200 simulated data sets per
  scenario). On adult GTEx material cut to the Cardoso-Moreira design the sensitivity of the two statistics differs
  strongly in some organs (forebrain, n 32: pass rate 0.055 unsigned, 0.970 signed; kidney, n 27: 0.495 and 1.000;
  `feasibility_fetal/output/passrates_cm_day.csv` and `passrates_cm_day_L1.csv`). The one positive that rests on R
  alone, the fetal adrenal of A189, is unchanged by the signed statistic (post hoc: C 0.926, p 0.0001, Holm 0.0009,
  above both companions, leave-one-gene-out minimum 0.904; `A189_fetal_within/output/h189_1_signed_posthoc.csv`).
- **My recommendation:** (i), as one registered audit amendment with the list of affected results fixed in advance;
  pass / fail verdicts that change are reported as errata. The gate (L2) and the coupling estimates do not use R.

### Q30 (2026-11-10): the mouse embryo "time-of-day state" after A236

- **Question:** A236 tested whether the late-embryo timing positive (A108 H108.1, R 0.52, p 0.014) reflects a running
  clock. All three registered tests fail: the time-locked variation of the clock genes is not clock-shaped (S1 0.069,
  null mean 0.368; power 0.98), no rotation is seen (rho 0.015; power only 0.33) and the alignment is not specific
  (18.5 % of expression-matched and 41 % of responsive pseudo-clock sets reach p < 0.05; R at the 93rd and 81.5th
  percentile). Does the paper (i) drop the embryo timing panel, (ii) keep it with the companions and the wording
  "time-of-day information, not clock-specific; age and clock time confounded within litter age", or (iii) keep only
  the untimed gate of the late window (p_pairing 0.001) as the developmental mouse statement?
- **On hold:** any use of H108.1 as evidence of a functioning embryonic clock, and the idea of placing embryo organs
  on external time from the inferred phase. Nothing else depends on it; the human fetal adrenal positive (A189) passed
  both companions and stands.
- **Evidence:** `results_preliminary.md` (A236); `A236_embryo_running_clock/output/h236_primary.csv`. The design is
  the limit (one embryo per clock time and litter age, 18 h in 6-h steps): a data set crossing age and clock time
  would be needed to ask the question properly.

### Q29 (2026-11-10): which confirmation first for the Q25 headline

- **Update 2026-11-13 (user's instruction of 2026-11-12: (a) and (b) both).** (b) done: A251, PASS, clock-specific (154
  living women, skin sites and layers), and extended to two organs with a genetic control (A255, MuTHER, PASS); no age
  contrast (ages not public). (a) age: the archive already holds it on the same data (A3); (a) depression: not
  scorable (A253: the patients of Li 2013 are not public, and in four other cohorts responsive genes pair across
  regions as strongly as clock genes). (c) the in-house jet-lag data remain the open route to a desynchrony test. See
  `REVIEW_2026-11-10.md` §6 for the full update and the list of data to request.

- **Question:** the fresh review (`REVIEW_2026-11-10.md` §3-§4) keeps the Q25 headline but finds that its break arm
  rests on one post hoc contrast (A132b) and its robustness arm on GTEx alone. Which confirmation do we pursue first:
  (a) timed post-mortem brain cohorts (GSE71620, Chen 2016: age 16-96, two areas; GSE45642, Li 2013: six regions,
  controls vs MDD), (b) a living skin cohort (Wu 2020: 154 women, three sites x two layers, untimed), or (c) the
  in-house jet-lag data?
- **On hold:** registration of D2-D4 of the review (the planned A117 was never registered); nothing already scored.
- **Evidence:** review §3.1 ("Missing for Cell") and §4 (D2, D4, with power arguments). (a) would also replicate or
  refute the A234 brain advance on other donors (31 subjects under 40 y against 13 in GTEx), which bears on Q26.
  Without (a), (b) or (c) the review's frank reading is a Nature Methods paper.

### Q28 (2026-11-10): brain death, one descriptive companion and the number shown

- **Question:** within the decision "descriptive, no mechanism amendments", (i) is one registered descriptive table in
  scope: within-tissue clock-gene coherence and amplitude in organ donors vs ventilated post-mortem donors (are the
  single-organ clocks intact while their mutual order is lost)? (ii) Should the main number be the full-n readout
  (rho 0.093, 298 organ donors) rather than one equal-n subsample (0.161), with donor management (steroid and hormone
  protocols, time to procurement) named as an open alternative?
- **On hold:** the wording "first human observation consistent with SCN-lesion work"; no analysis is started.
- **Evidence:** review §2.2 and §4 (D5). The registered H132.1 / H132.2 failed; A132b is post hoc; A136 H136.3 uses the
  same donors; A136 H136.1 (glucocorticoid exposure) did not pass its gate (p 0.098). Order lost with rhythms kept is
  the SCN-lesion picture (Yoo 2004); damped clocks would be a different result, and no readout separates them yet.

### Q27 (2026-11-10): contrasts at a common host spread (d*) are biased when the arms differ in spread

- **Update 2026-11-13 (A258, on the user's instruction):** bias measured (none at equal spread; -0.11 to -0.27 at log V ratio
  -0.8 to -1.0 with a noisy host); the arm-specific prior (`contrast(host_prior="arm")`) removes most of it but fails the
  registered acceptance, so no corrected column; the label-permutation test is anti-conservative at unequal spread.
  Flagged rows (expected bias > 0.05): A231 0-3 of 54, A233 6-33 of 302, A235 5-55 of 262. Recommendation: keep the
  registered numbers, mark the flagged rows, drop "d* removes the spread", and register a spread-matched null.

- **Question:** may a calibration amendment be registered (simulation at the real designs of A164, A231, A233 and the
  fibrotic-lung contrasts; candidate arm-specific correction), and if d* moves, do we keep the registered numbers with a
  corrected sensitivity column, or re-register the equivalence family on the corrected estimator?
- **On hold:** any sentence stating that d* "removes the sampling spread"; MDEs of contrasts with unequal spread. No
  published number is changed.
- **Evidence:** review §1, M1 (`_core/host.py:216, 226`; `dataset.py:481-505`). The pooled-prior host posterior means
  are over-shrunk in the arm with the narrower spread, so its |z| and d* read low although the coupling is the same.
  Post hoc feasibility simulation (not a result; `exploratory/scouts_2026-11-10/`): truth d* = 0, bias -0.11 to -0.19 at a spread ratio
  of e^-0.8 to e^-1.0 with a noisy host, ≤ 0.07 at e^-0.3 to e^-0.45, none at equal spread. Exposure of the equivalence
  claims is small (median log spread ratio between -0.32 and +0.24 in A231 / A233); the bias runs against the
  fibrotic-lung T-cell loss; 31 % of 4,214 rows in the exploratory contrast tables have |log V ratio| > 0.4.

### Q26 (2026-11-08): A130 direction label looks inverted; A234 finds the age advance in GTEx

- **Question:** confirm the A130 sign reading below, and whether the ageing clause can now say "internal time advances
  with age (GTEx), organisation unchanged (A231)".
- **On hold:** the erratum wording for A130 and any ageing sentence in the paper; nothing else.
- **Evidence:** A130 computes internal − wall time; with internal phase rising with clock time (H113.1), its -0.339 h per
  decade is a delay (not significant), although registered and summarised as "earlier with age". A234 (registered,
  same convention, direction checked in the data): body +0.64 h per decade (Holm 0.010), brain +0.99 (Holm 0.002), brain
  − body -0.01 (post hoc); estimate probably inflated at low power. Living cohorts (A130: 42 breast-surgery women in
  hospital hours + 12 men) and GTEx fast deaths then disagree in sign. See `results_preliminary.md` A234 and the A130
  sign note.

### Q25 (2026-11-07): the headline — two kinds of desynchrony, and how much "robust" is tested

**Proposal (user and agent, 2026-11-07).** Separate two kinds of desynchrony:
- **Body-wide (between organs):** in mice it follows conflicting cues (day feeding, abrupt light shifts) or loss of the
  central drive (SCN lesion, Yoo 2004). We detect an imposed one with known truth (Xin day feeding, H111.1 PASS,
  Fisher-Lee 0.988). In humans only brain death breaks the inter-organ order (A132b, post hoc: 0.161 vs 0.678).
- **Local (inside a tissue):** follows injury, regeneration, inflammation or transformation, with no cue conflict.
  Known-truth exemplar: regenerating mouse distal colon after DSS (Carmona-Alcocer et al., Nat Commun 2026), epithelium
  shifted against the stroma, recovered blind (−8.95 h untimed vs −9.03 h timed; A21 archive, qualified). Human
  candidates: psoriasis (A115 descriptive; A183 not clock-specific), tumours (stroma keeps host time, A227-A229).
- **Candidate headline:** "Human organs keep a shared internal time through ageing, chronic disease, cancer and the loss
  of daily cues; it breaks body-wide only when the central drive is lost, and otherwise only locally." Provocative
  question: body-wide desynchrony may require conflicting cues or loss of central drive, which mouse models impose and
  most human chronic disease does not. The in-house jet-lag data (conflicting cues) would then be the test of a
  prediction (transient desynchrony expected), not a missing control.

**Audit of the "robust" clauses** (`exploratory/scouts_2026-11-07/robustness_audit.md`, subagent, 22 rows; key numbers
re-checked against `results_preliminary.md`, some line references may be off by a line or two):
- **Equivalence shown (registered):** T2D d* +0.019, 90 % [−0.037, +0.076], obesity +0.004 [−0.083, +0.091], both
  inside ±0.15 (A164 re-read after A202); sex d* +0.020 [−0.003, +0.043] (A208, not yet in the ledger).
- **Not detected only:** ageing (A180 equivalence fails: γ 90 % [−0.147, 0.089] vs ±0.0375 per decade); disease atlases
  (A153 no primaries; "coupling preserved in disease" cannot be stated as tested); loss of daily cues: coupling above
  the pseudo null passes in ventilated donors (H109.1, excess 0.310) but the registered order test fails (H109.2 rho
  0.261, p 0.165; reading "a shared donor state rather than clocks with preserved phase relations"); "keeps the order"
  rests on Hardy 4 (secondary, rho 0.784) and the post hoc A132b split. ICU continuous feeding is not a measured
  variable.
- **Untested:** night shift and forced desynchrony (MDE 7.3 h and 13.1 h), human jet lag, season, other GTEx diseases.
- **Counter-cases:** brain death (contradicts the current ledger headline, CL:30); psoriasis uninvolved skin already
  late (A115), so not purely lesion-local.
- **Ledger fixes found (to apply after your decision):** FINDINGS_MAP still quotes pre-A202 A164 intervals; A190 cited
  at status B though withdrawn; A121 attribution; A124 is muscle, not blood; A54 −1.43 h is the post hoc value (discovery
  −1.50 h); SCN-specific Bmal1 KO (no coupling loss, VIP-Cre, 6 vs 6) absent from the ledger.

**Consequence.** As written, "keeps a shared time through ageing ... and the loss of daily cues" overstates the
evidence: only metabolic disease and sex are equivalence-backed. Honest version: "Inter-organ coupling is unchanged
within ±0.15 in metabolic disease; no body-wide desynchrony was detected with ageing, chronic disease or cancer at the
stated MDEs; the inter-organ order breaks when the central drive is lost."

**Human local test (scout, `exploratory/scouts_2026-11-07/ibd_scout.md`).** Epithelium-vs-stroma offset, inflamed vs
non-inflamed, same patient:
- TAURUS (on disk, 29 patients): stroma read in A45/A53 (coupling 0.80-0.91 in both states); the epithelial offset was
  never read, so a readout is non-blind, secondary only; confounded with treatment timing.
- GSE266616 (Crohn's, open, 3.7 GB): 14 patients with inflamed + adjacent-to-inflamed (1 strictly non-inflamed);
  epithelium and lamina propria in one library; depth unchecked.
- Smillie 2019 (~18 UC) and Eshghi (~18 UC): Single Cell Portal login (PI step, Q19).

**Questions.**
1. Adopt the two-kind framing and the honest headline wording above?
2. Make the jet-lag data a registered prediction test of the body-wide clause (owner's agreement)?
3. Human local test: (a) GSE266616 blind, inflamed vs adjacent, plus TAURUS as non-blind secondary; (b) acquire Smillie
   / Eshghi (login) for a cleaner blind design; (c) leave the local clause to mouse colon and tumours.
4. Apply the ledger fixes now?
- **Agent recommendation:** 1 yes; 2 yes; 3 (a) now, (b) if you can log in; 4 yes.
- **Done (user, 2026-11-07):** 3 (a) registered and scored as A230. GSE266616 not scored (host gate fails, 8 eligible
  patients, adjacent references only). TAURUS (non-blind): epithelial offset change +0.58 h, p 0.39, MDE 1.71 h: no
  mouse-sized shift in inflamed human gut. The fibroblast companion was replaced by the pseudo-clock check before scoring
  (user decision; two anchors give no leave-one-out companion). Option (b) (Smillie / Eshghi, login) remains.

### Q24 (2026-11-06): the cancer chapter after A214-A226 — structure and claims for the paper
The user wants **one paper** that includes the cancer chapter, held together by "how the body keeps its time and when it
loses it". This needs your review (Phase 8: framing and figures). All numbers: `results_preliminary.md` (sections named).
- **Positive, model-based:**
  - The tumour clock state matches the same patient's normal tissue (H138.4: 0.365 vs pseudo q95 0.261; 517 pairs, 7 of 10
    cancers p < 0.01).
  - The match does not depend on purity (H138.8: low − high +0.025, MDE 0.16), so malignant cells contribute. The shared
    state partly reflects the shared surgery hour, so the claim is that both tissues keep body time, not that a coupling
    strength is measured.
  - High-proliferation tumours pair less (A175 H175.1, +0.158, clock-specific).
  - Tumour fibroblasts and endothelium follow the patient's normal stroma (A149 note 5, post hoc, about 100 patients each).
    A confirmatory test in independent cohorts (A226) could not be scored: all sizeable public paired studies are already
    A149 cohorts (11 / 13 new patients). Proposed wording (user agreed, 2026-11-07): post hoc origin, reproduced in
    pre-specified analyses on overlapping patients (A227 anchors: fibroblasts 0.686 vs q95 0.539, endothelium 0.677 vs
    0.512; A228 control: tumour endothelium vs normal fibroblasts 0.642 vs 0.544), independent confirmation not possible
    (A226). Stroma keeps the signal in colorectal, where immune compartments do not.
  - HCC is less coherent than its own adjacent liver (GSE14520, −0.33).
  - Tumour myeloid and T/NK cells keep the patient's normal-tissue clock state, clock-specifically (A227, registered: myeloid
    0.581 vs pseudo q95 0.548, 169 patients; T/NK 0.539 vs 0.510, 227 patients); B / plasma and mast pair without clock
    specificity. Gradient stroma > myeloid ≈ T/NK > B follows residency; immune passes may reflect shared circulation.
  - A228: the clock-specific T/NK pairing holds within CD8-negative T cells (0.532 vs q95 0.501, 214 patients), not CD8 T
    or NK. Body time in tumour immune cells is not detected against normal fibroblasts (T/NK 0.441 vs q95 0.458; myeloid
    0.467 vs 0.487), although the endothelial reference (descriptive) and lung are positive; colorectal is at the pseudo
    level throughout.
  - A229 (Wang et al. Cell 2024 correlate): tumour endothelial VCAM1 (0.426 vs q95 0.402) and, borderline, the tumour T/NK
    share (0.378 vs 0.371) follow the endothelial clock state clock-specifically (secondaries); ICAM1 and the CD8 share
    (primaries) do not. Lung-weighted.
- **Corrective (the method's specificity and design discipline applied to the literature):**
  - Clock "prognosis" is proliferation (A191-A198).
  - Circadian LASSO signatures are no better than random gene lists; they "validate" exactly where random lists do (A220;
    KIRC, LGG, METABRIC 100 % of random lists).
  - Clock-gene mutations and losses occur at background rates; PER1 / PER3 losses ride on 17p (TP53, 98.9 % co-loss) and
    1p36 (A222).
  - "Clock genes are downregulated" is partly real in paired samples (PER1-3, CRY2, TEF, RORA; L −0.33) but not for ARNTL,
    NPAS2 or NR1D1, which GTEx comparators turn into apparent losses, partly through sampling time (A221).
  - Genomic-lesion effects on coherence are subtype composition (A215-A217).
  - Literature tally: 0 of 79 clinical claims adjusted for proliferation, 1 for subtype, 0 used random-set controls (draft,
    abstract level).
- **Open:**
  - Breast LumA tumours are more coherent at matched proliferation (A219, 3 cohorts, clock-specific).
  - Cis copy-number adjustment removes 63 % of that gap (A224), but a simulation shows additive copy-number noise cannot
    produce it (A225). The copy-number entanglement is real but not resolvable in bulk (biology vs tumour state).
  - TP53-mutant tumours without 17p loss have lower PER1 / PER2 (A224 H224.3, −0.16 log2), not adjusted for tumour state.
- **Not estimable:**
  - Host desynchrony (needs >= 2 normal tissues per cancer patient).
  - Single-cell per-cell questions (A214: 5-10x more depth needed).
  - Level-based damping readout (A223).
- **Proposed shape:**
  1. Tumours and their stroma keep body time.
  2. What the field read as clock disruption or prognostic clock information is proliferation, subtype composition,
     genomic position, design or sampling time; one main figure, the rest in the supplement.
  3. Open: tumour-state differences in coherence.
- **Questions:**
  - Is this structure acceptable?
  - How strongly should the literature correction be worded? It names published signatures, CYCLOPS magnitude and
    TimeTeller.
  - Should the A220-A222 tests appear as main or Extended Data?
- **On hold:** the cancer figure and its text.

### Q12 (2026-10-27): request brain time-of-death fields? (A141)
A141 (registered after your Q6 answer):
- 8 of 12 brain contrasts are noise-driven, so not estimable.
- 4 are anchor-dependent beyond noise: PD_mssm, VaD_mssm, ASCVD_hbcc and the PD atlas.
- Your Q6 answer foresaw a time-of-death request in this case. A draft (not sent) covers PsychAD on Synapse and the PD
  atlas on AMP PD, plus optionally ROSMAP: `paper/exploratory/A141_brain_stability/drafts/tod_request_draft.md`.
- **Question:** send the request(s), edit them, or not request?
- **On hold:** nothing else. The brain panels stay in Extended Data with the result stated.
- Evidence: `results_preliminary.md` (A141).

## Answered
### Q36 (2026-11-16): calibration of the pairwise cell-type test (A276-A278)
**Question (for the user first).** Under no pairing, the clock's statistic exceeds the q95 of its companions in 6 to
15 % of 500 placebo draws (A276, band 7.5 %): the clock-plane states are more isotropic and lighter-tailed than the
companion states, which raises null |r| by a few percent; the pooled statistics (means over atlases) turn that into the
largest rates. A per-set null correction (A277) leaves two responsive combinations at 0.080 and 0.078 (the responsive
sets overlap, so their reference is too narrow). The percentile test calibrated against its own placebo distribution
(A278) scatters around 0.05 (0.018 to 0.082, mean 0.043) but misses its own band in one combination, a band that was
too narrow for a split-sample check. **Options:** (a) adopt A278 with the 1000 placebo draws as the null and read the
restated tests H278.1-5 (`calibrated_a278.py --read`); (b) enlarge the null (e.g. 5000 draws) and re-run the A278
calibration with a band set for two-sample error before reading; (c) replace the statistic by one whose null does not
depend on state shape (mean |r|^2 over pairs, E = 1/df), a new amendment and a new calibration. **Recommendation:** (b)
then read: it keeps the registered statistics, and the cost is minutes. **Update 2026-11-17 (null only, not a
registered test):** option (c) is not calibrated either: with the mean |r|^2 over pairs and the q95 rule, 500 placebo
draws give 0.042 to 0.154 (H272.1 responsive 0.154, H272.2 responsive 0.116; mean 0.076), since the overlap of the
responsive sets is untouched (`A276_pairwise_calibration/output/option_c/tests.csv`). Only the percentile read against
its own placebo distribution (a, b) is calibrated by construction. A280 has a ready script for it
(`calibrated_a280.py --q36`), with the real pairing locked until the answer.
**What is at stake.** The scored margins are large for H272.2 (D 0.277 against q95 0.177), H274.1 and H274.2 (clock
above all 200 companions); a calibrated test changes them only if the null rate of "above all companions" exceeds 0.05,
and it is 0.005 to 0.028 (A278). The closest is H274.3 against the matched family (cortex, T_CX 0.497 against q95 0.466,
p 0.020), whose placebo rate at the q95 was 0.102. Every placebo check removes all pairing; the registered null has
companions that share (|r| 0.3 to 0.4), where shape effects should be smaller (not shown).
**On hold.** Every further test of the cell-type line (A2 pooled table, B1, B2, C1-C3 scoring). Scouts, inventories and
genome-table builds continue.
**Evidence.** RP A276, A277 calibration, A278 calibration; `exploratory/A276_pairwise_calibration/output/`.
**Answer (user, 2026-11-17):** "(b) 5000 draws, then read". Registered as A283 (V283 band 0.0685 per combination, both
split directions; restated tests H283.1-5 only if V283 passes); A280 addendum adopts the same form.

### Q21 (2026-11-01): cancer L2 after A175/A176, and the Q8 wording
**Evidence** (`results_preliminary.md`, A175, A176, C1 replication scout):
- A175 H175.1 (TCGA, secondary): high-proliferation tumours pair less with their own normal tissue (D +0.158, p 0.021,
  clock-specific). ER split not estimable (MDE 0.41).
- A176: the coherence loss is "swamped", not "disordered", in 7 of 7 cohorts. In proliferating tumours the clock-plane
  variance across tumours falls clock-specifically. In ER− tumours it is mostly a generic rise in off-plane variance.
- So the lower pairing in A175 may be less clock signal rather than less host following. The two cannot be separated
  without a spread-matched pairing readout.
- Primary replication: GSE14520 + GSE66229 pooled MDE 0.214 > TCGA D 0.158. Stopped by the stop rule.
**Questions:**
1. Q8 wording: replace "clock loss" by "less clock signal relative to non-clock variation (swamped)" in Figure 5?
2. Register a spread-matched pairing readout (pairing at common tumour clock-plane variance, d* logic) on TCGA as a
   secondary? This would be a second amendment on the C1 question, which the one-amendment rule forbids without your OK.
3. Look for CPTAC (LUAD/STAD/BRCA with normal adjacent tissue) for a powered replication? This means a PDC/GDC download
   (open access), size to be checked.
**On hold:** any L2 wording for cancer beyond "secondary, TCGA". Nothing else is paused.
**Update (scout 2026-11-01, `exploratory/scouts_2026-11-01/SCOUT.md`):**
- A powered C1 replication may exist in GEO: the Nurses' Health Study breast tumour/adjacent pairs (623 unique patients with
  tumour + adjacent normal: 445 on GPL22920, 178 on GPL17585 (GSE115577 re-analyses the GSE93601 arrays); ER− 117 /
  ER+ 503; Affymetrix HTA arrays, all 14 clock genes annotated; `scouts_2026-11-01/prep_nhs/output/nhs_counts.csv`).
- For H1, GSE154200 has 186 persons with lesional + non-lesional skin (psoriasis 88, nummular 52, atopic eczema 46).
- Proposal: download both (open GEO, small), check clock-gene coverage and gates, run the MDE step, and bring the
  primary registrations to you before scoring. Go ahead?

**Update (2026-11-01, after the downloads):** the NHS replication is stopped: 445 usable pairs, tumour–normal pair
gate fails (p_acrophase 0.098) and H175.1 MDE 0.198 > TCGA D 0.158 (`results_preliminary.md`, pre-registration prep).
CPTAC (point 3) is the only route left to a powered C1 replication.

### Q22 (2026-11-01): two primary registrations to approve before scoring (C6 blind types, H1 skin)
**A183 (H1, GSE154200 skin)**, draft `exploratory/scouts_2026-11-01/prep_gse154200/REGISTRATION_DRAFT.md`: blind
replication of the psoriatic lesion delay (A54 held-out −0.93 h; MDE 0.81 h, 87 persons) and its specificity against
nummular and atopic eczema in the same pipeline (MDE 0.69 h). d* is not estimable (two compartments), so the HANDOFF
target "MDE on d* ≤ 0.15" is replaced by the paired offset. **Question:** approve as is?

**A182 (C6)**, draft `exploratory/C6_blind_types/REGISTRATION_DRAFT.md` (headers only read). My defaults in brackets:
1. GSE17536 (177 colon, Moffitt; GSE14333 could not be downloaded): [separate test in the weak/null family].
2. GSE39582 and GSE65858 are ComBat-corrected by the authors: [keep, caveat in the text].
3. Prostate GSE116918 (custom FFPE array, 13 clock genes, no NR1D1; the only prostate cohort ≥ 150): [include with 13].
4. ESCC GSE53625 (no symbol annotation on GEO; relation line in an earlier scout header): [drop].
5. Clock genes: [the pipeline's 14 (with RORA, CIART), as in A171-A173], not the HANDOFF list (BHLHE40/41).
6. Pancreatic, endometrial and cervical have no blind cohort ≥ 150: [accept the narrower family (colon, prostate);
   no near-miss cohorts].
7. HNSC and ESCA lean negative in TCGA (−0.17, −0.15): [keep HNSC as a test; state the TCGA value next to it].
8. Null reading: TCGA values of the weak types (−0.03 to −0.11) are below the expected MDE (about 0.2), so a
   non-significant cohort cannot say "no loss". [Add an equivalence test per cohort: T > −0.19 (the smallest
   replicated cohort effect, GSE62254) with a 90 % interval.]
9. GSE157009 / GSE157010 (two arms of one Alliance study): [separate tests].
**D2 (GSE27272, mother–fetus at term), optional third registration:** 60 mother + cord blood pairs; the cord blood ×
mother gate passes (p 0.022 / 0.013), placenta × mother fails. The only estimable question is one-group: is cord blood
coherent with its own mother's blood beyond shuffled pairs, clock-specifically and beyond the responsive-gene
companion (delivery stress)? [Register it as a secondary, marginal design (n 60); smoking not estimable.] Draft on
your yes.
**On hold:** scoring of A182 and A183. Their power steps are done (A183) or run right after the commit (A182).
- **Answer (PI, 2026-11-01):** Q21.1 yes, Figure 5 uses "less clock signal relative to non-clock variation
  (swamped)". Q21.2 yes, the spread-matched pairing readout is registered as A185 (second C1 amendment, allowed).
  Q21.3 scout CPTAC (metadata, commands for the user). Q22: A183 approved as drafted; A182 with the equivalence bound
  −0.19, registered tests GSE39582, GSE17536, GSE116918, GSE157009, GSE157010, GSE65858 descriptive, GSE53625 dropped,
  all other defaults; D2 registered as a secondary (A184). Also: Q12 the PI reviews and sends the draft; C8 not yet;
  A136 proceeds as drafted.

### Q13 (2026-10-28): roadmap phase 0 decisions (ROADMAP_2026-10-28.md)
I am going ahead with phases 1, 2 and 4 (code, data integrity, simulations) and the phase 6 scouts. None of them
depends on the answers below. The defaults in brackets apply when phase 3 starts, unless you say otherwise.
- **0.1 Reported coupling readout:** excess over the pseudo-clock median (Δc, with pseudo q95), or raw c?
  [default: Δc as main, raw c in the supplement]
- **0.2 Eligibility rule for disease primaries:** a compartment enters only if its control Δc exceeds the pseudo q95 and
  the MDE on d* is <= 0.2, registered as A153 before the contrasts are rerun. [default: adopt]
- **0.4 Anchor self-inclusion (C7): A155 meets the decision rule.** Anchor amplitudes are inflated by 9-20 % with 3-4
  anchors (4-6 % with 6), while offsets (<= 0.06 h), coupling (already leave-one-out) and the A123 map (r unchanged) are
  unaffected. Proposal: leave-one-out becomes clocksync's anchor amplitude and offset readout (opt-in now, default
  after your OK), and Figure 2 amplitude panels are rerun in phase 3. [default: adopt]
- **0.6 Parked items:**
  - A136: register now? [default: I draft the registration for your review, not scored]
  - A147 reference-size sensitivity: [default: run, labelled post hoc, 2 min]
  - A139 grid: [default: skip]
  - A138 Q2 simulator: [default: option (b), fix the clock signal-to-noise with a dated note]
  - Q12 and A119: yours (outward-facing / procurement).
- **On hold:** phase 3 onwards (anchors, backbone, A150 extension, disease reruns).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** defaults accepted: Δc is the main readout (raw c in the supplement); family-level eligibility (as A153 note, Q15 a); leave-one-out anchor amplitude and offset become the default readout; 0.6 items as applied (A136 stays a draft for the PI).

### Q14 (2026-10-28): COARSE ordering choices (bone marrow, muscle)
The data audit found three ordering choices in the shared COARSE rules. They are not bugs of the "T cell" kind:
- "granulocyte monocyte progenitor cell" (21,500 cells) goes to `monocyte`;
- "mesenchymal stem cell" (28,275) and "skeletal muscle satellite stem cell" (13,531) go to `progenitor`, because
  progenitor comes before fibroblast_stromal and muscle.
- **Question:** keep them, or move GMP to progenitor and the two stem labels to fibroblast_stromal / muscle?
- **Default:** keep (the archive A9 choice). Changing them would rebuild organ_bone_marrow and organ_muscle and
  relearn their anchors.
- **On hold:** nothing. If you change them, I rebuild before phase 3.
- **Answer (PI, 2026-10-30):** keep.

### Q11 (2026-10-27): A140 Q1 simulator gate fails; rerun with the published CYCLOPS settings and a cleaner factor calibration?
The registered gate (arm D, whole-day sampling, CYCLOPS circular correlation with truth >= 0.7, seed 1,500) fails: median
0.067, 0 of 50 at >= 0.7. In the same replicates the clocksync gate passes 43 of 50, and the three CYCLOPS compartment phases
agree with each other (median |r| about 0.97) while ignoring the clock: CYCLOPS follows the shared non-clock factor. Per the
registration Q1 stops; H140.1-3 are not scored.
Two choices of the run may cause the failure, so the gate does not yet show that the published pipeline cannot work:
- CYCLOPS was trained with 10 models per fit instead of the template's 80, for cost: 80 models take about 80
  worker-minutes per fit, about 4-5 h on 16 workers for the 50 gate replicates.
- The non-clock factor was calibrated on the GSE266919 Pre samples, because the A17 tables hold clock genes only. These
  samples fail the gate, and their leading eigenvalues include Poisson noise, so the factor may be overstated.
  Post hoc diagnostic (8 replicates): factor scaled to 0.1 gives median truth r 0.30 (1 of 8 >= 0.7); without the
  factor, CYCLOPS crashed in its Julia code.
- **Question:**
  - (a) Accept the stop: report "the published pipeline does not recover a coupled clock at this factor size", as
    descriptive.
  - (b) Or register a revision before any rerun: 80 models, the factor calibrated on gate-passing tumour immune
    pseudobulks with all genes (needs those pseudobulks built), and the factor size as a grid. The gate is then rerun,
    and H140.1-3 only if it passes.
- **On hold:** A140 Q1 scoring and Q2 (ii) (the published pipeline on GSE266919).
- Evidence: `paper/exploratory/A140_dip_cyclops/README.md` (Q1), `output/q1/{simulator_gate,diag_summary}.csv`, `results_preliminary.md` (A140 Q1).

- **Answer (PI, 2026-10-27):** (b) revise and rerun: a dated revision registered before any rerun (80 CYCLOPS models; factor calibrated on genome-wide pseudobulks of gate-passing tumour immune compartments; factor size as a grid); gate, then H140.1-3 only if it passes.
### Q10 (2026-10-26): the A19 age model fails its null calibration in the first full count-level run
The archive validated `Dataset.covariate` (A19 part 2) with a model-based simulation (null rejection 0.025-0.08). The
migrated acceptance (M3 full run) keeps the real count tables and permutes ages within stratum: rejection at alpha 0.05
exceeds the registered 0.07 in 4 of 6 compartments for locked amplitude and for gamma (lung AT1 locked amplitude 0.220;
PsychAD excitatory neurons gamma 0.125, OPC 0.145). Lung AT1 and NK/ILC fits mostly sit at coupling 1 on the age subset.
- **Question:** (a) may I investigate (boundary fits on age subsets; the package's simulated null against the permutation
  null), and (b) if the count-level failure holds, should A19 part 2 (H19.1-2 nulls, the neuron gamma secondary, p <= 0.01
  per compartment against the simulated null) be re-read with a permutation p (ages permuted within stratum), as an
  erratum to the archive?
- **On hold:** any use of the A19 age readouts (neuron gamma) in figures or text; the erratum. M3 of resilience runs as is.
- Evidence: `results/extended_data/scripts/calibration/README.md` (M3 full run), `output/calibration/age_acceptance_*.csv`.
- **Answer (PI, 2026-10-26), part (a):** investigate. Diagnostic `paper/exploratory/A19_acceptance_check/check_null.py` (post hoc, real ages not read). Part (b) open.
- **Diagnostic result (2026-10-26):** the simulated null is too narrow in brain (SD ratio permutation / simulation
  1.10-1.28; 8-13 % outside its 95 %) and off-centre in the lung atlas (AT1 locked amplitude t median 1.84 vs 0.17; 34 %
  outside); the age fit is degenerate (own noise at floor) in 78-84 % of AT1 and NK/ILC permutations
  (`paper/exploratory/A19_acceptance_check/README.md`). Proposal for (b): add a within-stratum permutation null to
  `Dataset.covariate` (`null="permutation"`, default; refit per permutation; fits with own noise at the floor reported as
  not estimable for gamma), rerun the acceptance with it, then re-read A19 part 2 on the real ages as an erratum. Needs
  your go (clocksync change).
- **Answer (PI, 2026-10-26), part (b):** add the within-stratum covariate shuffle. Done: clocksync 613e224
  (`Dataset.covariate(null="permutation")`, default; gamma not estimable when own noise is at its bound). Running: the
  acceptance with it, then A19 part 2 on the real ages (erratum to follow). Note: the M3 port of A19 part 2 has no design
  check, so even under the simulated null the archive's neuron gamma secondary is not reproduced (-0.104, SE 0.077,
  p 0.173 vs -0.124, p 0.013; six compartments the legacy check excluded enter; `resilience/README.md`).
- **Result (2026-10-26):** with the permutation null the acceptance holds within Monte Carlo error (null rejection
  0.015-0.080; one cell of 24 at 0.080). A19 primaries are unchanged (null). Neuron gamma -0.116 (SE 0.054, p 0.031),
  with 2 of the 6 archive compartments at p <= 0.01. Erratum written (`results_preliminary.md`). The hold on the neuron
  gamma readout is lifted; it is a nominal secondary. Closed unless you object.

- **Answer (PI, 2026-10-27):** closed; the neuron gamma result goes in the main text as a nominal, uncorrected secondary (with the erratum).
### Q2 (2026-10-23): figure-table numbers used against their archive verdict (WP0 item 1)
- A104 rho +0.73 quoted as support; H104.1 predicted rho < 0 (archive: not supported; H104.2 not replicated).
- A54 psoriasis "replicated in discovery and held-out": the delay replicates (H54.2) but H54.1 and H54.3 fail.
- Mouse embryo "from E15": H108.1 passes in the late window but the window gate does not (p_acrophase 0.08).
- Quoted values that differ from their table: benchmark offset r 0.934 / 0.900 (quoted 0.91 / 0.88); baboon offsets
  0.455 h (0.5 h); healthy atlases 22 targets / 17 studies ("about 30").
- **Question:** remove, restate or keep each claim (with the verdict stated)?
- **On hold:** edits of the Figure 1, 2 and 5 tables.

- **Answer (PI, 2026-10-27):** go through each claim with the proposed text before editing.
  Item by item (PI, 2026-10-27): A104 restated as not supported (post hoc positive rho labelled); embryo restated with the
  borderline gate; psoriasis 'the delay replicates' with H54.1/H54.3 not met; numbers from the tables (0.934 / 0.900; 0.455 h;
  22 atlases, 17 studies, 14 tissues). Applied to the figure tables of `provisory_plan_paper.md`.
### Q3 (2026-10-23): panels not on the anchored A17+ estimator (WP0 item 4)
13 of 30 figure-table rows flagged (`paper/exploratory/WP0_integrity/output/estimator_consistency.csv`). By design
(the text must name the readout): r* (A103, A42), A104, A85 index, A92/A93 coherence C, A33 tissue host. Needing a
decision: Fig 1 A21 colon (an anchor decouples, V17.2 fails); Fig 2 GTEx organ departures (pre-A17 timed) and
neuron/retina panels (pre-A17 free fit); Fig 4 A20 "damping" (withdrawn in A22) and COVID A105 (2 anchors, no A22
confirmation); Fig 5 fibrotic-lung T cells (A5 free-fit 0.45).
- **Question:** drop, recompute on the anchored estimator, or keep with the readout named, per panel?
- **On hold:** recomputation of those panels.

- **Answer (PI, 2026-10-27):** lead's default per panel: by-design rows kept with the readout named; GTEx departures, neuron/retina and fibrotic-lung T cells recomputed on the anchored estimator; A20 'damping' and A21 colon dropped; COVID A105 to Extended Data with '2 anchors' stated.
- **Done (2026-10-27):** neuron/retina and fibrotic-lung T cells recomputed (results_preliminary.md, "WP0 / PI Q3"). PI (2026-10-27): GTEx departures kept with the timed readout named (no anchored version); neuron panel restated in ED (retina under both rules; cortex anchor-defined), message 'the brain's backbone is neuronal' dropped.
### Q6 (2026-10-24): brain disease contrasts are anchor-dependent (WP0 item 8)
With an all-LOO host instead of the A17 anchors, brain contrasts (Figure 4: PD, AD, ASCVD, VaD, SCZ ...) change sign in
1-7 of 9 compartments (block agreement r* 0.42-0.67); the peripheral headline contrasts (ILD, fibrotic lung, T2D, A44,
AD skin) keep their sign. Per-cohort cancer tumour-minus-normal coupling changes sign in 3 of 9 cohorts (near zero).
- **Question:** keep the brain disease panels with the anchor dependence stated, move them to Extended Data, or drop?
- **PI (2026-10-24):** decide evidence- and literature-driven whether one compartment has a good reason to be the brain reference (covariance, clock strength, literature). -> A135 (registered before any brain contrast is recomputed).
- **A135 result (2026-10-24):** no compartment family qualifies (excitatory neurons first in 2 of 5 control sets; rule 4 of 5). By the registered rule the brain disease panels go to Extended Data with the anchor dependence stated. Confirm, or keep open?
- **Time of death for the brain cohorts (checked 2026-10-24, column names only):** absent from every local donor table
  (PsychAD MSSM/HBCC/RADC, SEA-AD DLPFC/MTG, PD atlas) and from the public SEA-AD donor metadata (PMI only). Controlled
  portals not checked (no access from here): PsychAD on Synapse syn52160016; ROSMAP clinical (`ROSMAP_clinical.csv`, AD
  Knowledge Portal); PD atlas (Prashant et al. 2024, AMP PD Knowledge Platform). PI to decide whether to request it.
- **On hold:** Figure 4 brain panels.
- Evidence: results_preliminary.md, WP0 item 8.

- **Answer (PI, 2026-10-27):** brain donors die at any hour, so the brain cohorts likely cover the full cycle and the
  inferred amplitudes and phases are close to absolute. Compare them with the GTEx and simulation amplitudes to check
  whether the amplitude and phase spread is large enough for the model's estimates to be stable. If not, nothing can be
  said. If yes, describe first which compartments agree, and which should agree or are noise-driven. If they are not
  noise-driven, the PI may request the time-of-death fields to learn what follows it. -> to be registered (descriptive
  amendment) before any brain readout.
### Q9 (2026-10-24): direction after WP0-WP4 (step back)
Proposals (see the chat summary of 2026-10-24; each would be registered or written into the plan once agreed):
- (a) Figure 3 claim: from "a self-sustained peripheral network" to "organs stay coupled clocks when cues are disrupted
  (Hardy 4); perturbation acts in two separable ways: a shared displacement that keeps the relations (ventilated, not
  brain dead) and a loss of the relations (brain dead)". Drop the mechanism row of WP2 from the main text.
- (b) Figure 1: add the genetic positive control (A129: knockout liver decouples, specific) and the living cohorts
  (A113, with the post hoc estimator stated); drop "cue/locked split validated in mouse" (A112 failed).
- (c) Method: make "pull vs order" (signed pole-deviation Spearman and shrink ratio) a registered clocksync readout,
  validated on the A133 simulation grid before use.
- (d) Process: every registration states the attainable minimum p and the null-free ceiling of each statistic
  (H111.1 and A133 set unreachable thresholds).
- (e) Next tests, in order: A136 drug vs brain death within organ donors (glucocorticoid-response score), an external
  brain-dead / multi-organ replication, A134-A135 (running), then WP5-WP8 with WP5 reframed by A134.
- (f) Journal: the disease positive for Cell rests on cancer (being reframed) and the anchor-robust peripheral
  contrasts; Nature Methods is the stronger fit for the validated method.
  Correction (2026-10-24, PI comment): the psoriatic local-clock delay is the strongest pre-specified disease positive
  (archive: A54 H54.2 held-out -0.93 h, p 0.004; A78 core minus systemic cues -1.97 h, p 0.0002; A81b GTEx cue set
  -1.69 h, p 0.0003; other skin diseases -0.75 h, p 0.22). It is a delay (lesional clock later), not an advance. Open:
  trait vs state and treatment specificity (WP6) -> proposed to run WP6 before WP5.
- **Question:** which of (a)-(f) do you agree with?
- **Answer so far (PI, 2026-10-25):** (a) agreed; (b) agreed; (c) asked what "pull vs order" means (explained in chat,
  decision pending); (d) agreed (in the spirit of stating power); (e) cancer is a priority for the literature: CoFE
  (Ananthasubramaniam 2025) may have been misled by its method, and many papers used CYCLOPS 2.0 although it works only
  when samples cover the whole day. Determine whether anything can be estimated in tumours, or whether the variance
  between malignant cells and patients is too large. -> WP5 to be redesigned around this (proposal in chat, to register).
  Follow-up (PI, 2026-10-25): (c) agreed: pull vs order becomes a clocksync readout, validated on the A133 grid first.
  Both WP5 (redesigned) and WP6 will run; cancer first (A138). Literature agents may start now.
  A138 choices (PI, 2026-10-25): thresholds left to the lead (20 % / 50 %); CYCLOPS 2.0 tested both plain and with the
  healthy whole-day anchoring of Li et al. 2024 (assumes tumour = healthy Fourier coefficients); the CoFE authors were
  already told their tumour clock is wrong (no email needed). A138 registered (0162cc2).
- **On hold:** plan rewrite; A136 registration.


- **Answer (PI, 2026-10-27):** (f) journal decided after WP6. Plan rewrite and A136 registration: write both now as drafts for review.
### Q4 (2026-10-23): clocksync change
`clocksync.genes.gene_columns` returns -1 silently for a missing symbol (indexing the last column; the A85 bug).
Proposed: `missing="ignore"|"warn"|"raise"` and an all-absent guard
(`paper/exploratory/WP0_integrity/output/gene_audit/proposed_fix.diff`, with a test). clocksync holds your uncommitted
work, so it was not applied.
- **Question:** apply it (own commit, only these files)?
- **Answer (PI, 2026-10-24):** apply it. Done: clocksync 0a288ac (`genes.py` and a new test only; default behaviour kept, `missing="raise"` used in `paper/common/gtex.py` and `natexp_common`).
- **On hold:** nothing (no current result affected).

### Q7 (2026-10-24): what A110 can still say about the SCN
GTEx organ donors have no brain samples (0 of 301), so the registered SCN test (hypothalamus coupling, Hardy 3-4 vs
brain-dead donors) is impossible; only fast vs Hardy 3-4 remains (descriptive, being scored). The illness-duration test
of memory vs maintained coupling (H110.1) is infeasible: the interval is recorded for 16 of 103 Hardy 4 donors.
- **Question:** accept that GTEx cannot separate SCN loss from donor drugs, or look for another dataset (e.g. a
  brain-dead vs non-brain-dead ICU cohort with multi-organ or blood time series)?
- **Answer (PI, 2026-10-24):** accept that GTEx gives no statistically robust SCN result; the SCN reading stays a hypothesis.
- **On hold:** nothing now; the SCN sentence of Figure 3 stays a hypothesis.

### Q8 (2026-10-24): what the cancer "coherence loss" means (WP0 item 2)
A86/A87's lower share, coherence and alignment in malignant cells are not a depth effect (malignant pseudobulks are
deeper; a 4-fold shallower normal epithelium keeps coherence 0.68 vs malignant 0.23). But per-cell D is equal in
malignant and normal cells (3.89 vs 3.90), which rejects scrambled gene phases within cells (simulated 0.42). A87 also
shows 3-fold higher off-plane (non-clock) variance of the clock genes across tumours, with clock-plane variance kept.
- **Proposal:** restate the main cancer claim as "malignant cells keep a coherent clock in each cell; across tumours,
  clock-gene expression is dominated by non-clock variation that breaks the shared-phase structure", after a test:
  (a) A87 coherence within patients across cells (or clones), (b) coherence after removing the copy-number-explained
  part of each clock gene (WP5 H114 dosage). Register as an amendment before scoring.
- **Question:** agree with testing this before WP5, and with the possible restatement of the Figure 5 claim?
- **Answer (PI, 2026-10-24):** agreed -> A134.
- **A134 result (2026-10-24):** inferred copy number explains a small part (malignant coherence 0.23 -> 0.29, p 0.017 vs normal; about 1/10 of the gap); oncogenic programmes nothing; per-cell test uninformative. Figure 5 wording: keep "coherence lost" with the per-cell caveat and the copy-number share, or restate? (open)
- **On hold:** the Figure 5 cancer narrative; WP5 (A114) framing. The A86/A87 numbers themselves stand.
- Evidence: results_preliminary.md, WP0 item 2 (second and third parts).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** the Figure 5 wording sub-question waits for C3 (in-plane vs off-plane variance: "swamped" vs "disordered").

### Q5 (2026-10-23): A132 outcome and the Figure 3 wording
Within ventilated (Hardy 0) donors, the push of pair phases toward 0 / 12 h is seen in both brain-dead organ donors and
ventilated post-mortem cases (not brain dead), with no separation (H132.1, MDE 0.44-0.54 h). Descriptively, only the
non-brain-dead group keeps the fast-death phase relations (rho 0.69, p 0.0006 vs 0.16, p 0.19 at equal n).
- **Question:** (a) word Figure 3's perturbed group as "ventilated donors" rather than "brain-dead organ donors"?
  (b) register a post hoc test (labelled post hoc) of the offset-correlation difference, organ donors vs ventilated
  post-mortem, as support for your SCN reading, or leave it descriptive?
- **On hold:** Figure 3 narrative wording; any A132 follow-up. A110 continues (H110.4 already uses organ donors).
- Evidence: results_preliminary.md, section A132.
- **Answer (PI, 2026-10-24):** describe both groups (the 0 / 12 h pull in both, phase relations kept only without
  brain death); a post hoc p for the difference is fine, labelled post hoc -> A132b (rho difference 0.517, p 0.008).
  PI reading: the pull is likely drug-induced; the kept offsets are hard to explain.

### Q1 (2026-10-23): Figure 3 framing after A109
H109.1 passes (organs share a clock-shaped state in every death group); H109.2 fails in ventilated donors (phase
relations collapse to 0 / 12 h; about half of the shared covariance is a displacement) but passes in Hardy 4
(secondary). Registered reading: "organs share a donor state", WP2 out of the main text. Proposed alternative: coupled
running clocks without cues, shown in Hardy 4 / Hardy 3, with ventilated donors as a shared displacement overriding
the clock relations. PI comment (2026-10-23): check first that Hardy 4 clocks are decorrelated from wall time, and
whether the 0 / 12 h push is seen in COVID ventilated patients -> registered as A109c (dab55ab), being scored.
- **A109c H109c.1-2 (2026-10-23):** no detectable wall-time alignment in Hardy 4 (R 0.097, p 0.56; MDE 0.145-0.208);
  offsets preserved after removing wall time (rho 0.750, p 0.0002), also without ventilator cases and in inpatients
  only. COVID (D109c.3-4, descriptive): Delorey organs uninformative (too few shared donors); Melms lung shows no
  push in COVID vs controls. The circularity ratio is confounded by host-phase concentration (sampling window), so the
  offset test (H109.2) is the decisive readout; the ventilated collapse to 0 / 12 h rests on it (p_closer 0.001).
- **Question:** after A109c, which framing (registered default, or the Hardy 4 framing as a registered amendment
  labelled post hoc relative to A109)?
- **On hold:** WP2 (A110); edits of the Figure 3 table and narrative.
- Evidence: results_preliminary.md, sections A109, A109b (and A109c when scored).
- **Answer (PI, 2026-10-23):** reframed claim accepted: when daily cues are lost through illness (Hardy 4), organs
  remain coupled clocks with their phase relations; in ventilated (mostly brain-dead) organ donors a shared perturbation
  overrides them. PI note: brain-dead donors have probably lost the SCN, Hardy 4 donors keep it with perturbed
  zeitgebers; caveat: drugs in both groups, not only cues. Follow-up: A132 (SCN loss vs ventilation), A110 rewritten
  for Hardy 4 / Hardy 3.


## Q15 (2026-10-29): A153 has no primaries as registered — which eligibility rule?
**Finding (A159, results_preliminary "Phase 3 ... A159"):**
- No Figure 4 control arm has a genome-wide table on disk.
- The registered per-compartment flag (c above both pseudo q95) is underpowered: a single compartment's pseudo q95 is
  0.6-0.9, so only 3 of 31 healthy compartments pass, all in breast. Heart fibroblast at c 0.72 fails.

So A153 as registered has an empty primary family. Changing the rule after seeing D159.1 is a design change (rule d).
**Options:**
- (a) **Family-level eligibility from healthy atlases** [recommended]. A family is eligible if its coupling, averaged over
  compartments and atlases, exceeds the matched and responsive q95 of the same average under pseudo sets (set j pooled,
  as H150.1 / H159.2). Disease primaries are then the compartments of eligible families with MDE on d* ≤ 0.2.
  - Uses A150 + A159 atlases (5).
  - Needs no disease genome tables.
  - Likely admits stroma and endothelium; immune would be decided by the test.
- (b) Build genome-wide tables for the non-brain disease datasets (about 18 contrasts; Census/h5ad sources), then apply (a)
  per control arm, or a per-arm family version of it. This is 1-2 days of builder work and disk, and also makes H159
  scorable on new atlases. It can follow (a) as a confirmation.
- (c) Drop eligibility: every estimable compartment with MDE ≤ 0.2 is a primary, and Δc is reported descriptively.
  Larger Holm family, so less power.

**Answer (PI, 2026-10-30):** start with (a); (b) is also good and can follow later if it seems important, after the
other tests. Confirmed in session on 2026-10-30: register the A153 note with rule (a) and run A153.

## Q16 (2026-10-29): fresh tumour-normal stromal pairing (roadmap 7.2) — worth registering on what is available?
**Scout** (`exploratory/scouts_2026-10-29/a151_fresh_scout.md`, metadata only). Few fresh paired single-cell cohorts with
stroma exist:
- Kyushu gastric (16 untreated pairs, verified; Seurat .rds of 18-22 GB, parsed in Python only);
- Li 2022 ccRCC (up to 12, unverified);
- OPSCC (≤ 7, unverified);
- oesophageal GSE173950 (conditional).

That is 16 verified patients and about 45 if the rest hold. Colorectal, liver and breast stay not estimable. The bounded
statistic flags small cohorts (n 4-11), so a fresh A151 would likely be an underpowered descriptive.
**Larger alternative (not fresh):** marker-based fibroblast re-annotation of cohorts already in A149/A151 that have 0
fibroblast pairs: GSE150290 23, GSE199654 13, GSE144240 10, GSE178481 9, GSE242299 8. These have been seen at the
malignant level but not for fibroblasts.
**Options:**
- (a) Register a fresh A151 on Kyushu + ccRCC + OPSCC, as the main-text candidate.
- (b) Register the re-annotation of the used cohorts, labelled "not fresh".
- (c) **[default, recommended]** Keep the pairing descriptive in Figure 5 (lead with A97/A127 coherence vs
  proliferation), and spend no more compute on 7.2.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (c): single-cell stromal pairing stays descriptive. The bulk host-referenced tumour–normal tests (HANDOFF_2026-10-31, C1) take over.

## Q17 (2026-10-29): coupling intervals undercover (A154 V154.3 fails) — how to report c?
**Finding:**
- With every compartment at the same nominal c, ĉ is biased low (c 0.3 → median 0.09-0.26), and 95 % bootstrap
  intervals cover in 22 of 24 cells only at 0.00-0.89.
- A noise-free diagnostic shows the cause is a host built from weakly coupled anchors: with precise anchors ĉ is 0.288
  for 0.3.
- Contrasts at a common host stay calibrated (V154.2 passes). Real anchors are the stromal-vascular compartments, with c
  about 0.6-0.7.

**Options:**
- (a) **[recommended]** A registered follow-up simulation (A162) at realistic anchor coupling (anchors c 0.6-0.8,
  targets 0.1-0.9), measuring coverage and bias. If coverage is adequate there, report c with the caveat "conditional on
  the anchor host".
- (b) Develop an attenuation correction (host reliability) in clocksync and validate it by simulation (several days).
- (c) Report c only as ranks and contrasts (Δc, d*), not absolute values with intervals.

**Default if no answer:** (a). I will register and run it.

**Outcome (2026-10-30):** A162 fails (V162.1 19/32 cells, V162.2 16/32; results_preliminary A162). Coverage holds at
c ≤ 0.3 but not at 0.6-0.9; the estimate is shrunk by a roughly constant factor (0.87-0.92 with anchors at 0.8), as an
errors-in-variables attenuation by the host's reliability. Per the registered reading, absolute c is reported without
intervals (ranks, Δc, d*), and a correction is developed (option b): registered as A163.

**Outcome of A163 (2026-10-31):** the frozen correction fails validation (V163.1 30/32, V163.2 27/32, V163.3 58/64
cells; V163.4 passes). It removes most of the attenuation: coverage at c_target 0.9 goes from 0.02-0.19 to 0.90-0.94
(median), but it misses at c_target ≥ 0.6 with n = 30, K = 4, and with an interval-width shortfall at c_target 0.9.
Per the registered reading, c stays reported without intervals (ranks, Δc, d*). **Decision for you:**
- (1) accept this, which is my default;
- (2) register a restricted claim (c_anchor ≥ 0.8 with n ≥ 100 or K = 8) on fresh seeds;
- (3) develop a better interval (e.g. a bootstrap that re-estimates R per replicate, or BCa) before a new validation.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** option (1): absolute c is reported without intervals (ranks, Δc, d*); no further interval development now.

## Q18 (2026-10-30): A153 under rule (a) still has no primaries — where should the disease positive come from?

- Rule (a) makes fibroblast_stromal the eligible family (endothelial too in the 3-atlas sensitivity).
- Every eligible compartment of the 13 estimable disease targets has an MDE on d* > 0.2. The smallest are
  dcm_healthy heart fibroblast 0.232 and hlca_fibrotic EC_arterial 0.214. So H153.1-3 are not scorable, and the A17
  islets beta effect does not survive (results_preliminary, A153).
- **Nothing depending on this is paused.** A164 (GTEx T2D/obesity, matched) runs as registered, and roadmap 6.6
  (within-person designs) is next in the handoff order.
- **Options for the PI:**
  - (i) Relax the MDE bound to 0.25 in a new registration. That would admit only dcm_healthy heart fibroblast; it is
    post hoc in spirit, so I do not recommend it.
  - (ii) Option (b) of Q15. It does not change the MDEs.
  - (iii) Larger single-cell disease cohorts (an outward action: data requests).
  - (iv) Rely on within-person and GTEx designs for the disease positive.
- My recommendation is (iv) now, with (iii) as a longer-term ask.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (iv): the disease positive is sought in larger bulk within-person designs from routine clinical sampling (HANDOFF_2026-10-31, Priority 2, H1), with cancer C1 as the host-referenced tumour test. Scouts list only estimable bulk designs: at least 2 compartments (tissues, sites or sorted cell types) per person, since one compartment per person has no host and gives L1 coherence only (PI, 2026-10-31).

## Q19 (2026-10-30): within-person designs (roadmap 6.6) — the on-disk paired cohorts are already analysed; acquire Smillie 2019?

- Scout (metadata only, no clock statistic): `exploratory/scouts_2026-10-30/within_person/SCOUT.md`.
- **On disk, paired, but not blind:**
  - TAURUS (UC + CD): 29 patients with inflamed and non-inflamed biopsies; 23-24 have both eligible families
    (fibroblast_stromal and endothelial, >= 50 cells per state). Archive A45 and A53 already read it within patient:
    stromal coupling 0.80-0.91 in both states, phases unchanged, no amplitude loss. Inflamed vs non-inflamed is also
    confounded with pre/post treatment (pre 67/37, post 37/63 samples).
  - Psoriasis single cell (GSE173706, GSE228421): 12 patients with both families in lesional and non-lesional
    skin; used in A54/A58.
  - Kong 2025 (Crohn vs normal, same donor): 3-4 donors, too few.
- **New, independent:** Smillie 2019 UC atlas (18 UC patients with inflamed and adjacent non-inflamed biopsies,
  12 healthy; 366,650 cells; not on disk, pairing to be verified). Hosted on the Broad Single Cell Portal, which needs
  a login: an access step for the PI.
- **Nothing is paused.** I have not registered a TAURUS or psoriasis re-read: any new readout there would be chosen
  after A45/A53/A54 were seen, so it could only be a secondary.
- **Options:** (i) acquire Smillie 2019 and register a blind within-patient contrast (stroma and endothelium:
  d*, phase, z2, psi), with TAURUS as the not-blind secondary; (ii) register the TAURUS readouts as an explicitly
  non-blind secondary only; (iii) drop 6.6, since gut stroma looks resilient to inflammation (A45).
- My recommendation is (i).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** bulk first: closed for now. H1 scouts paired bulk IBD series; Smillie 2019 is revisited only if no bulk design reaches MDE ≤ 0.15 on d* (the portal login stays a PI step). Scouts list only estimable bulk designs: at least 2 compartments (tissues, sites or sorted cell types) per person, since one compartment per person has no host and gives L1 coherence only (PI, 2026-10-31).

## Q20 (2026-10-31): proliferation coherence loss in gastric malignant cells — worth our own malignant-cell calling?
- **Context:** the proliferation clock-coherence loss is now replicated in breast (A97 discovery, A127 replication; A166
  H166.3 failed, underpowered), lung adenocarcinoma
  (A171) and gastric cancer (A173, two of three cohorts, clock-specific). It is not pan-cancer (A172 fails, power OK).
  In breast single cells, ER− malignant cells are less coherent (A170, p 0.044) but not clock-specifically, so
  the registered cell-intrinsic reading is not supported.
- **Gastric single cell (scout `scouts_2026-10-31/sc_gastric/`):** no public set has malignant-cell labels together
  with raw counts. The candidates are GSE270680 (25 tumour donors), GSE183904 (26) and GSE150290 (23). All three are
  public 10x / csv matrices, about 2.8 GB in total, and none has cell annotations. A test would need our own epithelial
  annotation and CNV-based malignant calling, a methodological choice made before scoring. By analogy with Wu 2021,
  perhaps 10-20 patients per set would have ≥ 50 malignant cells.
- **Nothing is paused.**
- **Options:**
  - (i) register a pooled three-set gastric single-cell test with a pre-specified CNV caller;
  - (ii) stop the gastric line at bulk;
  - (iii) instead, test whether "adenocarcinoma vs squamous" holds where both are large. This would mean the TCGA
    ESCA / HNSC / CESC rows, which are seen, so only as a secondary.
- **My recommendation is (ii) for now.** The bulk replication is the paper-relevant result. Single-cell n would be
  marginal, and the CNV step adds a degree of freedom.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (ii): gastric stops at bulk. Single-cell cancer work (C4) uses sets with author malignant labels (Wu + Pal breast, a labelled LUAD atlas after a scout).

## Q23 (2026-11-04): ARNTL was zero in the A18-era per-cell files — rebuild all and rerun the per-cell readouts?

- **Finding** (results_preliminary, integrity 2026-11-04): seven per-cell files (hlca_full, islets_t2d, natri2024_ild,
  parkinson_regions, psychad_mssm, psychad_radc, seaad_dlpfc) have ARNTL = 0 (Census symbol BMAL1). Pseudobulk results,
  the A152 files and the cancer per-cell files are not affected.
- **Affected readouts:** A18 part 2 per-cell amplitude / synchrony, A144 part 2 (descriptive, inconclusive), the
  Extended Data per-cell calibration. None is in the claims ledger.
- **Done:** builder fixed; HLCA rebuilt for A203.
- **Nothing else is paused.**
- **Options:** (i) rebuild the six other files (about 3 h of Census streaming, about 9 GB) and rerun the three readouts,
  with the old and new numbers side by side; (ii) drop the per-cell readouts from the paper and keep the erratum.
- **My recommendation is (i)** if any per-cell number stays in the Extended Data, else (ii).
- **Answer (user, 2026-11-04):** (ii) drop the per-cell readouts (A18 part 2 per-cell amplitude and synchrony, A144 part 2,
  the Extended Data per-cell calibration) and keep the erratum. Amended the same day: the malignant per-cell D of WP0 item 2
  stays (cancer files have ARNTL; malignant depth 4.3 clock UMI per cell, where D resolves gene-phase coherence). The
  six other per-cell files are not rebuilt. Basis: the dropped readouts rest on the ARNTL-zero files and were inconclusive; none is in the claims ledger. A203
  bears on donor-level precision, not on within-sample readouts such as D.
