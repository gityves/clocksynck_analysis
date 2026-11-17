# Results section: slide by slide

Order (agreed 2026-11-11): R1 model validation, R2 organisation in health (with ageing), R3 development, R4 loss of
daily cues, R5 chronic disease, R6 where clocks break: locally (psoriasis, cancer), R7 misalignment and limits,
R8 synthesis. Every number below was re-read from the table named in its slide's **Source** line. GTEx figures show
per-tissue aggregates only, never per-donor times.

Status labels: **registered** (primary or secondary, scored as written), **confirmation** (archive result known,
re-run on current code), **seen data** (registered after the data were looked at), **post hoc**.

---

# R1. Model validation

**Chapter message:** where the true time is known, the untimed method recovers it. Organ phase map from
working-hours deaths, patient order in living surgery, organ offsets in a full-day primate atlas, and the two
perturbations that should move it (a liver without its clock, a feeding shift) are detected.

---

## R1.1 GTEx: the organ map from deaths between 08 and 18 h

**Message:** from donors who all died in working hours, the method recovers the organ phase map and the organ
amplitudes that whole-day timing gives.

**Design:** GTEx v10 fast deaths (Hardy 1-2). Test set: the 111 donors who died 08-18 h, 40 tissues with ≥ 20 such
donors. Truth (i): cosinor on the recorded time of death of all fast deaths. Truth (iii): the same on the 133 donors who
died outside 08-18 h only (disjoint from the test set).

**Bullets**
1. Offsets against truth (i): circular r 0.945 [0.907, 0.972], median error 0.35 h. Amplitudes: Pearson (log2) 0.976.
2. Against the disjoint truth (iii): amplitudes 0.965 (pass). Offsets 0.679 [0.467, 0.839], below the registered bar
   of 0.7 (**fail as registered**).
3. Why (post hoc): the estimate reproduces across disjoint donor halves (0.873), the night-only truth does not (its
   halves agree at 0.544). The shortfall comes from the truth's noise, not from the estimate.

**Figure:** `figures/gtex_organ_map/gtex_organ_map.png` (`scripts/gtex_organ_map.py`), the CHIRAL-paper organ map
(angle: tissue offset, radius: relative amplitude, upper half of the clock face, -6 to +6 h) for A. the timed truth (i),
B. clocksync and C. CHIRAL, both untimed on the 08-18 h donors; in B and C each tissue is joined to its truth (grey).
To add if wanted: the split-half agreement (estimate halves vs truth (iii) halves).

**Speaker notes**
- Why the best direction and rotation: an orderer's cycle has an arbitrary direction and origin; scored as run, a
  reversed draw negates every tissue offset (CYCLOPS was reversed in all full-day draws), which made the orderers look
  worse with the full day than in hospital hours. Aligned for every method alike, the full day helps everyone, and the
  advantage of clocksync is specific to hospital hours (post hoc scoring, 2026-11-17).
- The orderers' working-hours failure is not a reflection or stretch artefact: allowing the better direction and a
  least-squares stretch, CHIRAL stays at |r| 0.165 (post hoc check).
- CoFE genome-wide and CYCLOPS also fail with the whole day (0.139, 0.115): they order donors along non-clock
  variation.
- The distinctive property on real data is partial-day coverage, not coupling estimation as such: on idealised
  whole-day simulations two-step orderer pipelines recover coupling nearly as well (A145).
- DIP vs internal time (post hoc diagnostic, `scripts/gtex_dip_diag.py`, full day): the two disagree at DIP 12-18 h,
  where CHIRAL places only 35 donors (62 died in that interval). These bins hold the weakest shared clock states (0.75-0.80
  of the median donor) and their recorded times are scattered (R 0.32-0.44); at DIP 15-18 h the DIP is > 6 h off the
  recorded time for 32 % of donors (clocksync 21 %). Weak clocks are harder for both (Spearman of strength with error
  −0.27 clocksync, −0.40 DIP); strongest third 1.41 vs 1.44 h, weakest third 2.87 vs 4.18 h.
- Why not circular r: over the full day the DIP's |circular r| with time of death is higher (0.677 vs 0.600), but the
  difference is not resolved (paired bootstrap 95 % CI [−0.46, +0.40]): with deaths spread over the day the
  Jammalamadaka-SenGupta mean direction is ill-defined. The median error difference is (DIP − clocksync +0.53 h
  [0.16, 0.86]).
- Status: H123.4 registered, pass.

**Source:** `exploratory/A123_working_hours/output/a123b/{tests,h123_4,posthoc_reflect}.csv`.

---

## R1.3 Living patients: internal time follows the surgery hour

**Message:** in living humans, each patient's untimed internal time follows the recorded hour, and this is carried by
the clock genes.

**Bullets**
1. Breast surgery (42 patients, tumour + adjacent normal, resection 08-18 h): Spearman of internal time with the
   recorded hour 0.558; 200 expression-matched non-clock gene sets: median −0.053, none reaches it (p 0.005).
   Median error with no recorded time used 1.42 h (shuffled hours 2.28 h).
2. Heart surgery (16 patients, 8 operated 06-10 h, 8 at 14-18 h; four cell types): the afternoon host sits away from
   the morning host (p 0.0004), in the direction forward time predicts (19.7 h against 18 h, within the registered
   ±3 h).
3. Subcutaneous adipose (12 men, 48 biopsies over the day): circular r with clock time 0.822 (p 0.001), median error
   1.11 h.

**Figure:** `figures/living_patients/living_patients.png` (`scripts/living_breast_compute.py`, `scripts/living_patients.py`).
A. breast: internal time (13:00 rule, no recorded time used) vs resection hour, one dot per patient, same range on both
axes. B. heart: each patient's internal time by surgery window (06-10 h vs 14-18 h), circular means. C. adipose:
single-biopsy internal time vs zeitgeber time (common offset removed), circular mean and SD per time point. D. per-sample
absolute error after removing one common offset: breast (08-18 h) 1.12 h, adipose (full day) 1.11 h.

**Speaker notes**
- Breast: the ordering does not rest on glucocorticoid-induced PER1 (0.560 without it); it leans on NR1D1 (0.233
  without PER1 and NR1D1). This cohort was seen before the registration (V25.5): it is a robustness test, not a
  replication. A fresh cohort with recorded resection times is the open step.
- Breast, why most patients look early (post hoc, `output/living/breast_calibration.csv`): the 13:00 rule assumes the
  cohort's mean hour; the recorded mean is 13.8 h, so every estimate starts 0.86 h early (64 % underestimated). The
  internal time is also compressed (slope 0.59 on the recorded hour): late patients come out about 2 h early. The order is
  right (Spearman 0.558); with the centre fixed, the per-patient error (1.12 h) equals the full-day adipose biopsies (1.11 h).
- Compartments: breast = two bulk compartments per patient (tumour, adjacent normal), internal time = phase of their shared
  host. Heart = four snRNA-seq pseudobulks (cardiomyocyte, fibroblast, endothelial, neuron; anchors fibroblast,
  cardiomyocyte, neuron). Adipose = one bulk biopsy, no host model: the biopsy's phase on the GTEx clock plane.
- Heart: the truth is a 4-h window per patient, so the error is on the group displacement, not per patient. Holds in
  each surgical arm (p 0.032, 0.033), with fixed anchors and with a left-ventricle waveform. Archive result known:
  confirmation on current code.
- Adipose: muscle (GSE108539) fails on the registered statistic; it passes only with an estimator chosen after the
  failure (Fisher-Lee 0.576), so the muscle pass is post hoc and not shown.
- Status: A200 H200.1-2 registered, pass (seen data); A174 V174.1-3 pass (confirmation); A113 H113.1 adipose
  registered, pass.

**Source:** `exploratory/A200_breast_surgery_rerun/output/{tests,secondary}.csv`;
`exploratory/A174_v3_heart/output/v3_scores.csv`; `exploratory/A113_living_time/output/a113_summary.csv`.

---

## R1.4 A primate sampled around the clock: 63 tissue offsets

**Message:** in a baboon atlas sampled every 2 h over a day, the untimed tissue offsets match the timed ones within
half an hour.

**Bullets**
1. 63 tissues, untimed offset against the timed offset: circular r (Fisher-Lee) 0.912, median error 0.45 h
   (peripheral tissues 0.38 h, neural and ocular 0.50 h).
2. The untimed host phase orders the 12 sampling times: R 0.987, median error 0.39 h.

**Figure:** `figures/baboon_offsets/baboon_offsets.png` (`scripts/baboon_compute.py`, `scripts/baboon_offsets.py`; wide).
A. organ map as the GTEx one: untimed offset (angle) and relative amplitude (radius) per tissue, coloured by organ system,
grey rings: the timed values. B. 63 tissues, untimed vs timed offset after one common shift; dashed: identity. C. the 12
animals' internal time (host phase) vs zeitgeber time after one common shift. Tissue codes: Mure 2018
(`figures/baboon_offsets/tissue_mapping.jpeg`).

**Speaker notes**
- Archive amendment (A36), re-run on current code with leave-one-out anchors (median error 0.455 → 0.452 h). The
  registered statistic (Jammalamadaka-SenGupta on the host phase) is degenerate for an even 12-point grid: the numbers
  here are descriptive.
- One animal per time point: the 12 host phases are 12 animals.
- In 8-h windows of the same atlas (12 windows of 5 time points): order |Spearman| ≥ 0.9 in every window, direction
  right in 11 of 12; offsets degrade to median r 0.55, median error 1.2 h (A40).

**Source:** `results/figure_1/output/species/a36_readouts.csv`; `results/figure_1/output/windows/baboon_windows_summary.csv`.

---

## R1.5 Positive controls: the readout moves when the clock is broken or shifted

**Message:** removing the liver's clock decouples the liver, and only the liver; shifting the liver with daytime
feeding moves its untimed offset as the timed one moves.

**Bullets**
1. Hepatocyte Bmal1 knockout vs control (24 vs 24 mice, ad libitum): liver coupling 0.967 → 0.676, d* −0.287
   (p 0.001). Same permutations: white fat −0.009, lung +0.019; the liver drop exceeds both (p 0.001).
2. Night to day restricted feeding (Xin, 4 tissues): untimed offset changes track the timed ones, Fisher-Lee 0.988
   (exact p 0.042, the minimum attainable with 4 tissues). Liver moves +8.8 h untimed, +9.2 h timed.

**Figure:** `figures/positive_controls/positive_controls.png` (`scripts/positive_controls.py`; half slide). A. change in
coupling at the control's host spread (d*) per tissue, knockout minus control, 95 % interval: liver −0.287 [−0.411, −0.118],
white fat −0.009 [−0.089, 0.055], lung +0.019 [−0.045, 0.087]. B. inverted feeding (food in the day instead of the night): each
tissue's shift relative to the other three, untimed vs timed; dashed: identity.

**Speaker notes**
- The knockout liver keeps a clock-shaped shared state (its gate still passes) and sits −1.31 h from the fat + lung
  host (timed −0.98 h).
- Not shown but stated: daytime feeding in the Manella cohort (3 tissues) cannot reach p < 0.05 by design (minimum
  exact p 0.167; r 0.986); and the registered prediction that daytime feeding lowers liver coupling failed (white fat
  and lung drop as much). The readout detected the shift, not a liver-specific decoupling.
- Mouse controls: the whole day is sampled, so these test the readout, not the hospital-hours design.
- Status: A111/A129 H111.3 registered, pass; H111.1 Xin registered, pass; Manella fail by design; H111.2 fail.

**Source:** `exploratory/A111_feeding_controls/output/{a111_h111_3,a111_h111_1,a111_h111_2,a111_clockloss_liver_offset}.csv`.

---

## R1 what is not claimed

- The disjoint-truth bar for GTEx offsets was missed (0.679 < 0.7); the noise explanation is post hoc.
- Living-patient validation rests on cohorts seen before registration (breast) or known from the archive (heart);
  no fresh timed resection cohort yet.
- Gate passes alone are not claims: the gate makes 13 to 26 false claims in 150 no-clock simulations, the gate plus
  the matched pseudo-clock companion none (A209, A209b; method slide M6).

---

# R2. Organisation in health

**Chapter message:** a person's organs share one clock state, and the stronger an organ's clock, the more it follows
it. Inside tissues the shared clock is carried by stromal-vascular cells; lymphocytes carry none of their own. With
age, internal time advances while the organisation stays.

---

## R2.1 Organs share one clock state; weak clocks are the loosely coupled ones, the brain included

**Message:** every organ follows the body's clock state, and how closely it follows is set by the strength of its own
clock: strong clocks are the coupled ones, the opposite of what oscillator theory predicts. The brain's lower
coupling is this rule, not a separate brain clock. Within a person, a tissue stays within 1 to 3 h of the body's time.

**Bullets**
1. Across the 27 GTEx peripheral tissues, stronger clocks are more coupled (Spearman rho +0.734; registered prediction,
   from oscillator theory, rho < 0: not supported).
2. Specific to the clock genes: 200 expression-matched and 200 stress-responsive gene sets give rho near 0 (medians
   −0.160 and 0.109; q95 0.258 and 0.453). Calibrated against the stress-responsive sets: A290, running.
3. The brain follows the rule. In the baboon, where each animal's time is known, coupling rises with each tissue's
   rhythm amplitude against known time (Spearman 0.82 over 57 tissues; periphery 0.77, brain 0.74), and the weakly
   rhythmic brain regions fall on the peripheral line (registered on the untimed amplitude, A289: median residual
   −0.001). Human brain couples at 0.406 against 0.679 in the periphery.
4. A tissue's clock departs from its person's body time with an SD of 1.0 h (thyroid) to 3.2 h (lung), 11 body
   tissues; an upper bound (not above the companion gene sets).

**Figure:** `figures/organs_health/organs_health.png` (`scripts/organs_compute.py`, `scripts/organs_health.py`; wide).
A. human (GTEx fast deaths): clock-gene variation (A104 "own amplitude": all clock-plane variation between donors,
shared and own; log scale) against coupling to the body's clock state, one dot per tissue, coloured by organ system
(brain pink); grey: least-squares line over the peripheral tissues. B. baboon: rhythm amplitude against known time
(A36 timed gain, estimated without the host; log scale) against coupling, 35 peripheral tissues and 22 brain regions;
grey: the peripheral line. C. variation-coupling correlation of the clock genes (blue) against 200
expression-matched and 200 stress-responsive gene sets (grey dots). D. SD of a tissue's departure from its person's
body time, with 90 % interval.

**Speaker notes**
- Status: A104 registered (H104.1, rho < 0, not supported; the positive relation is a post hoc reading made
  registered by A144 and A271 Part B); A144 part 1 (registered): not noise shared within donors, it lies in the shared
  part; A271 Part B H271.2 pass (q95 rule; the calibrated form, A290, is running); A289 H289.3 (registered, fail:
  brain on the line, the outcome shown); A259 registered (H259.1 fail, H259.2 no tissue above its companions).
- Why panel B uses the timed amplitude (post hoc choice of display): the untimed "own amplitude" counts every
  clock-gene variation between animals, not only the rhythm, and shares a term with coupling. Baboon tissues at
  coupling 0 with sizeable untimed amplitude (lateral hypothalamus, thalamus, testis, pineal; smooth muscle 0.06) have
  weak rhythms against known time (0.29 to 0.76, smooth muscle 1.67; periphery median 2.46): their untimed amplitude
  is non-rhythmic variation. The timed amplitude is estimated without the host, so the relation in B cannot come from
  a shared term; it also shows that a coupling of 0 marks a tissue without a detectable shared rhythm.
- Why strong clocks couple more is open. A144's registered candidate: a tissue's amplitude is a population amplitude;
  body-wide signals keep its cells in phase, raising amplitude and coupling together (part 2 inconclusive).
- Panel A plots the A263 fit (peripheral host, pituitary in neither group: 26 peripheral tissues, rho 0.745); the
  registered A104 statistic counts the pituitary as peripheral (27 tissues, 0.734). Human brain regions sit above the
  peripheral line, but only by extrapolation (brain amplitudes 0.96-1.33 at the periphery's lower edge, post hoc); the
  baboon, with overlapping amplitudes, decides (A289).
- The brain's lower coupling is not specific to the clock genes (A263): every gene set drops in the brain, and the
  clock is the program the brain shares most with the body (brain coupling 0.406 against stress-responsive and
  matched set medians 0.064 and 0.029; post hoc). Brain vs periphery as blocks (A142, registered): r* 0.639 [0.530,
  0.728] against random 13-tissue peripheral blocks (median 0.978). Wording: "the brain is less coupled to the body's
  shared state", never "the brain's clock is the exception".
- Not a measurement artefact (post hoc): the brain's own unshared variation is smaller than the periphery's (psi
  median 0.38 vs 0.94); it shares far less of the body's state (z2 0.16 vs 1.51); coupling does not follow clock-gene
  counts (A271 post hoc, rho 0.09); the fresh, timed baboon brain shows the same.
- Not the blood-brain barrier (A289, registered): pineal and pituitary, outside it, couple less than brain regions in
  baboon (0.12 vs 0.42) and sheep (0.08 vs 0.20), the pineal not at all (post hoc reading: its clock is set by the
  SCN's sympathetic input). They are glands, so the barrier is not isolated.
- Compared with CHIRAL (Talamanca 2023): CHIRAL gives each tissue a mean amplitude and phase relative to one donor
  phase, assuming every tissue follows it fully. Coupling estimates how fully; panel D the spread between people
  around a tissue's mean offset (a variance across donors, not a per-person estimate). With donors over the whole day,
  a two-step CHIRAL pipeline could approach coupling (A145); hospital hours and the gene-set comparison are what it
  lacks.
- Panel D bounds how far a tissue can run on its own time; it does not show that it does (H259.2).

**Source:** `presentation/output/organs/couplings.csv`; `results/figure_2/output/organs/r2_couplings.csv`;
`results/figure_1/output/species/a36_tissues.csv`;
`exploratory/A271_responsive_companions/output/{a104_sets,a104_tests}.csv`; `exploratory/A289_brain_barrier/output/tests.csv`;
`exploratory/A263_companions/output/sets_human.csv`; `exploratory/A142_organ_blocks/output/tests.csv`;
`exploratory/A259_tissue_deviation/output/score/{targets,tests}.csv`.
